"""Portable speed / solver benchmark for the 4D aligner staging model.

Run on any machine to answer two questions:
  1. how long would the 70-step staging run take here?
  2. is a faster linear solver available on this machine than AMG-preconditioned CG?

The hot path is ~30 nonlinear iterations per solve, 140 solves per group.  Each
iteration rebuilds the PDL stiffness, assembles the contact operator, builds an
AMG hierarchy and runs CG.  On the reference machine (4-core laptop) the AMG
setup alone is ~75% of the wall time, so that is where any win has to come from.

    python benchmark.py              # full benchmark
    python benchmark.py --quick      # skip the alternative-solver trials
    python benchmark.py --probe-only # just report the environment
"""
import os, sys, time, platform, argparse, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, 'ofj4d'))
sys.path.insert(0, os.path.join(HERE, 'ofj'))
os.chdir(HERE)
from solve4d import solve_config, MAX_GAP          # noqa: E402

BAR = '=' * 78


def human(n):
    return f"{n/1e6:.2f}M" if n >= 1e6 else f"{n/1e3:.1f}k"


# --------------------------------------------------------------- environment --
def probe_env():
    print(BAR); print("ENVIRONMENT"); print(BAR)
    print(f"  platform      {platform.platform()}")
    print(f"  python        {platform.python_version()}  ({platform.machine()})")
    try:
        import psutil
        print(f"  cpu           {psutil.cpu_count(logical=False)} physical / "
              f"{psutil.cpu_count(logical=True)} logical")
        vm = psutil.virtual_memory()
        print(f"  memory        {vm.total/2**30:.1f} GiB total, "
              f"{vm.available/2**30:.1f} GiB available")
    except ImportError:
        print(f"  cpu           os.cpu_count() = {os.cpu_count()} (pip install psutil for detail)")
    print(f"  processor     {platform.processor()[:60]}")
    for v in ('OMP_NUM_THREADS', 'MKL_NUM_THREADS', 'OPENBLAS_NUM_THREADS'):
        if os.environ.get(v):
            print(f"  {v:<13} {os.environ[v]}")
    import scipy
    print(f"  numpy {np.__version__}   scipy {scipy.__version__}")
    try:
        cfg = np.__config__.show(mode='dicts')
        blas = cfg.get('Build Dependencies', {}).get('blas', {})
        print(f"  BLAS          {blas.get('name','?')} {blas.get('version','')}")
    except Exception:
        pass

    print("\n  candidate accelerators / solvers")
    found = {}
    checks = [
        ('pyamg',       'AMG preconditioner (required baseline)'),
        ('pypardiso',   'Intel MKL PARDISO direct solver  <-- usually the big win'),
        ('sksparse',    'SuiteSparse CHOLMOD sparse Cholesky'),
        ('scikits.umfpack', 'UMFPACK LU'),
        ('petsc4py',    'PETSc (many solvers, MPI)'),
        ('cupy',        'CUDA arrays + cuSPARSE CG'),
        ('cupyx',       'CuPy scipy-compatible sparse'),
        ('torch',       'PyTorch (GPU SpMV)'),
        ('pyamgx',      'NVIDIA AMGX: GPU AMG *including setup*'),
        ('numba',       'JIT for the assembly loops'),
    ]
    for mod, what in checks:
        try:
            m = __import__(mod)
            v = getattr(m, '__version__', '')
            print(f"    [available] {mod:<18} {v:<10} {what}")
            found[mod] = v or True
        except Exception:
            print(f"    [  --     ] {mod:<18} {'':<10} {what}")
    if 'torch' in found:
        try:
            import torch
            if torch.cuda.is_available():
                p = torch.cuda.get_device_properties(0)
                print(f"    GPU: {p.name}, {p.total_memory/2**30:.1f} GiB, "
                      f"compute {p.major}.{p.minor}")
                found['gpu'] = p.name
                bw, fl = gpu_balance()
                print(f"      fp64 bandwidth {bw[0]:6.1f} GB/s vs host {bw[1]:5.1f} GB/s"
                      f"   ({bw[0]/bw[1]:5.1f}x)  <- assembly is bandwidth-bound")
                print(f"      fp64 flops     {fl[0]:6.2f} TF/s  vs host {fl[1]:5.2f} TF/s"
                      f"   ({fl[0]/fl[1]:5.1f}x)  <- factorisation is flop-bound")
                if fl[0] < 2*fl[1]:
                    print("      => assemble on the device, factorise on the host "
                          "(what ofj4d/gpu4d.py does)")
            else:
                print("    torch present but CUDA not available")
        except Exception as e:
            print(f"    torch CUDA probe failed: {e}")
    return found


def gpu_balance():
    """The two numbers that decide what belongs on the device: fp64 memory
    bandwidth and fp64 arithmetic, each against the host."""
    import torch
    d = torch.device('cuda')
    n = 20_000_000
    a = torch.ones(n, dtype=torch.float64, device=d); b = torch.ones_like(a)
    for _ in range(2): c = a + b
    torch.cuda.synchronize(); t0 = time.time()
    for _ in range(5): c = a + b
    torch.cuda.synchronize(); g_bw = 5*3*n*8/(time.time()-t0)/1e9
    an = np.ones(n); bn = np.ones(n)
    t0 = time.time()
    for _ in range(3): cn = an + bn
    h_bw = 3*3*n*8/(time.time()-t0)/1e9
    m = 3072
    A = torch.randn(m, m, dtype=torch.float64, device=d); B = torch.randn_like(A)
    for _ in range(2): C = A @ B
    torch.cuda.synchronize(); t0 = time.time()
    for _ in range(3): C = A @ B
    torch.cuda.synchronize(); g_fl = 3*2*m**3/(time.time()-t0)/1e12
    An = np.random.randn(m, m); Bn = np.random.randn(m, m)
    t0 = time.time(); Cn = An @ Bn
    h_fl = 2*m**3/(time.time()-t0)/1e12
    return (g_bw, h_bw), (g_fl, h_fl)


def bench_loop(use_gpu, n_short=3, n_long=9):
    """Seconds per nonlinear iteration, measured on the real loop rather than
    modelled: the slope of wall time against iteration count, which cancels the
    once-per-solve setup."""
    import solve4d
    from model4d import Static, State
    from arch4d import ArchForm
    from morph4d import Morpher
    from loads4d import Elastic
    from contact import build_pairs
    solve4d.USE_GPU = bool(use_gpu)
    S = Static(); st = State(S); A = ArchForm(S); M = Morpher(S, A)
    for u in S.keep:
        S.tooth[u]['centroid_cur'] = st.centroid(u)
    load = Elastic(S, st, group='control')
    slave = np.unique(S.aligner_inner_faces)
    shape = M.shape(1)
    pre = S.prescribed_nodes(st, A, 1, shape)
    Cp = build_pairs(pre, slave, S.master_tris, k=16, max_gap=MAX_GAP)
    gap_pre, _ = S.gap_for_pairing(pre, slave, Cp['tri'], Cp['bary'])
    cref = dict(slave=slave, tri=Cp['tri'], bary=Cp['bary'], gap_pre=gap_pre)
    st.aligner = M.seat(shape, st, 1)
    ts, ns = [], []
    for n in (n_short, n_long):
        t0 = time.time()
        r = solve_config(S, st, load, cache={}, contact_ref=cref, maxit=n)
        ts.append(time.time()-t0)
        ns.append(len(r['hist']))       # the loop can converge before the cap
    if ns[1] > ns[0]:
        per_iter = (ts[1]-ts[0])/(ns[1]-ns[0])
        setup = ts[0] - ns[0]*per_iter
    else:                               # both runs did the same work: no slope
        per_iter, setup = ts[1]/max(ns[1], 1), 0.0
    return per_iter, max(setup, 0.0)


# ------------------------------------------------------------------ hot path --
def build_system():
    """Assemble one representative system (step 1, control group)."""
    from model4d import Static, State
    from arch4d import ArchForm
    from morph4d import Morpher
    from loads4d import Elastic
    from solve4d import Operators, contact_ops, K_N, K_T, MAX_GAP
    from contact import build_pairs
    import spec
    from scipy.sparse import coo_matrix, diags

    t = {}
    t0 = time.time()
    S = Static(); st = State(S); A = ArchForm(S); M = Morpher(S, A)
    for u in S.keep:
        S.tooth[u]['centroid_cur'] = st.centroid(u)
    L = Elastic(S, st, group='control')
    t['model setup (once)'] = time.time()-t0

    slave = np.unique(S.aligner_inner_faces)
    shape = M.shape(1)
    pre = S.prescribed_nodes(st, A, 1, shape)
    Cp = build_pairs(pre, slave, S.master_tris, k=16, max_gap=MAX_GAP)
    gap_pre, _ = S.gap_for_pairing(pre, slave, Cp['tri'], Cp['bary'])
    st.aligner = M.seat(shape, st, 1)
    nodes = st.nodes(); ndof = 3*S.n_global

    t0 = time.time(); P = S.build_P(nodes); t['build_P'] = time.time()-t0
    t0 = time.time(); Brm = S.rigid_modes(nodes); t['rigid modes'] = time.time()-t0
    t0 = time.time()
    opA = Operators(nodes, S.aligner_tets, spec.NU_ALIGNER)
    opP = Operators(nodes, S.pdl_tets, spec.PDL_NU, reuse=True)
    t['element operators (once/config)'] = time.time()-t0
    t0 = time.time()
    KAr = (P.T @ opA.K(np.full(len(S.aligner_tets), spec.E_ALIGNER), ndof) @ P).tocsr()
    t['aligner K + triple product'] = time.time()-t0

    E = np.full(len(S.pdl_tets), spec.PDL_E1)
    # first call also builds the fixed CSR pattern; the per-iteration cost is the
    # second call, which is what the nonlinear loop actually pays
    t0 = time.time(); KP = opP.K(E, ndof); t['PDL K + pattern (once/config)'] = time.time()-t0
    t0 = time.time(); KP = opP.K(E, ndof); t['PDL K (per iteration)'] = time.time()-t0
    t0 = time.time(); Kr = KAr + (P.T @ KP @ P).tocsr(); t['PDL triple product'] = time.time()-t0

    gap_now, nrm = S.gap_for_pairing(nodes, slave, Cp['tri'], Cp['bary'])
    C = dict(slave=slave, tri=Cp['tri'], bary=Cp['bary'], normal=nrm, gap0=gap_now)
    dofs, R, a = contact_ops(C)
    gref = gap_now - gap_pre
    act = np.where(gref <= 0)[0]
    Pt = np.eye(3)[None] - np.einsum('ni,nj->nij', nrm, nrm)
    t0 = time.time()
    Ke = (K_N*np.einsum('ni,nj->nij', a[act], a[act])
          + K_T*np.einsum('nki,nkl,nlj->nij', R[act], Pt[act], R[act]))
    d = dofs[act]
    Kc = coo_matrix((Ke.ravel(), (np.repeat(d, 12, axis=1).ravel(),
                                  np.tile(d, (1, 12)).ravel())),
                    shape=(ndof, ndof)).tocsr()
    t['contact operator'] = time.time()-t0
    fc = np.zeros(ndof)
    np.add.at(fc, d.ravel(), ((-K_N*gref[act])[:, None]*a[act]).ravel())
    t0 = time.time(); Am = (Kr + (P.T @ Kc @ P)).tocsr(); t['assemble A'] = time.time()-t0
    b = P.T @ (L.vector(S, st, ndof) + fc)
    dg = np.sqrt(Am.diagonal()); dg[dg <= 0] = 1.0
    Di = diags(1.0/dg)
    t0 = time.time(); As = (Di @ Am @ Di).tocsr(); t['diagonal scaling'] = time.time()-t0
    return As, b/dg, Brm*dg[:, None], t, len(act), len(gref)


def bench_baseline(As, bs, B):
    from solve4d import make_amg
    from scipy.sparse.linalg import cg
    print(f"\n{BAR}"); print("BASELINE: AMG-preconditioned CG (what the model uses)")
    print(BAR)
    t0 = time.time(); ml = make_amg(As, B); ts = time.time()-t0
    n = [0]; t0 = time.time()
    x, _ = cg(As, bs, rtol=1e-6, maxiter=600, M=ml.aspreconditioner(cycle='V'),
              callback=lambda z: n.__setitem__(0, n[0]+1))
    tc = time.time()-t0
    res = np.linalg.norm(bs - As@x)/np.linalg.norm(bs)
    print(f"  AMG setup      {ts:8.2f} s     <-- dominates; rebuilt EVERY iteration")
    print(f"  CG             {tc:8.2f} s  ({n[0]} iterations, residual {res:.1e})")
    print(f"  total          {ts+tc:8.2f} s per nonlinear iteration")
    return ts+tc, x


def bench_alternatives(As, bs, xref, found):
    print(f"\n{BAR}"); print("ALTERNATIVE LINEAR SOLVERS on the same matrix")
    print(BAR)
    nrm = max(np.linalg.norm(xref), 1e-30)
    rows = []

    def report(name, fn):
        try:
            t0 = time.time(); x = fn(); dt = time.time()-t0
            err = np.linalg.norm(x - xref)/nrm
            res = np.linalg.norm(bs - As@x)/np.linalg.norm(bs)
            print(f"  {name:<34}{dt:8.2f} s   residual {res:.1e}   "
                  f"rel.diff vs baseline {err:.1e}")
            rows.append((name, dt, res))
        except Exception as e:
            print(f"  {name:<34}   failed: {type(e).__name__}: {str(e)[:44]}")

    if 'pypardiso' in found:
        import pypardiso
        report("MKL PARDISO (direct)", lambda: pypardiso.spsolve(As.tocsr(), bs))
    if 'sksparse' in found:
        from sksparse.cholmod import cholesky
        def chol():
            f = cholesky(As.tocsc())
            return f(bs)
        report("CHOLMOD Cholesky (direct)", chol)
    if 'scikits.umfpack' in found:
        from scipy.sparse.linalg import spsolve
        report("UMFPACK LU (direct)",
               lambda: spsolve(As.tocsc(), bs, use_umfpack=True))
    from scipy.sparse.linalg import splu
    report("SciPy SuperLU (direct)", lambda: splu(As.tocsc()).solve(bs))
    if 'cupy' in found:
        import cupy, cupyx.scipy.sparse as csp
        import cupyx.scipy.sparse.linalg as csl
        def gpu_cg():
            Ag = csp.csr_matrix(As.astype(np.float64))
            bg = cupy.asarray(bs)
            xg, _ = csl.cg(Ag, bg, tol=1e-6, maxiter=2000)
            return cupy.asnumpy(xg)
        report("CuPy CG (GPU, no preconditioner)", gpu_cg)
    if 'pyamgx' in found:
        print("  pyamgx present: NVIDIA AMGX can do the AMG SETUP on the GPU, which is")
        print("    the 75% that CuPy cannot touch. Worth wiring in if this machine has it.")
    if not rows:
        print("  (no alternative direct solvers installed -- see requirements-optional.txt)")
    return rows


def extrapolate(rows, n_iter=26, solves_per_step=2, steps=70):
    """`rows` is a list of (name, seconds per nonlinear iteration, setup per solve)."""
    print(f"\n{BAR}"); print("PROJECTED RUNTIME ON THIS MACHINE"); print(BAR)
    print(f"  assumed {n_iter} nonlinear iterations per solve, "
          f"{solves_per_step} solves per staging step")
    for name, per_it, setup in rows:
        per_solve = setup + n_iter*per_it
        print(f"\n  {name} ({per_it:.2f} s per nonlinear iteration)")
        for s in (10, 20, 30, 70):
            tot = s*solves_per_step*per_solve
            mark = '   <- full protocol' if s == steps else ''
            print(f"    {s:3d} steps  {tot/3600:7.2f} h per group "
                  f"({2*tot/3600:6.2f} h for buccal TAD + control){mark}")
    print(f"\n  reference (4-core i7 laptop, 2 threads, AMG+CG): ~0.15 h/step "
          f"=> 30 steps ~ 5-6 h per group")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--quick', action='store_true')
    ap.add_argument('--probe-only', action='store_true')
    ap.add_argument('--json', default=None)
    A = ap.parse_args()

    found = probe_env()
    if A.probe_only:
        return

    print(f"\n{BAR}"); print("ASSEMBLY (per nonlinear iteration unless noted)"); print(BAR)
    As, bs, B, t, nact, npair = build_system()
    for k, v in t.items():
        print(f"  {k:<34}{v:8.2f} s")
    setup_once = t['model setup (once)'] + t['element operators (once/config)'] \
        + t['aligner K + triple product'] + t['build_P'] + t['rigid modes'] \
        + t['PDL K + pattern (once/config)']
    per_assembly = (t['PDL K (per iteration)'] + t['PDL triple product']
                    + t['contact operator'] + t['assemble A'] + t['diagonal scaling'])
    print(f"\n  system: {As.shape[0]} dof, {human(As.nnz)} nnz, "
          f"{nact}/{npair} contact pairs active")
    print(f"  per-iteration assembly total: {per_assembly:.2f} s")

    from solve4d import Pardiso
    per_lin, xref = bench_baseline(As, bs, B)
    rows = [] if A.quick else bench_alternatives(As, bs, xref, found)

    best = min(rows, key=lambda r: r[1]) if rows else None
    proj = [('modelled, AMG+CG (the original path)',
             per_assembly + per_lin, setup_once)]
    if best is not None and best[1] < per_lin and best[2] < 1e-8:
        proj.append((f'modelled, host assembly + {best[0]}',
                     per_assembly + best[1], setup_once))

    # the configured path, measured on the real loop rather than modelled
    if not A.quick:
        import gpu4d
        print(f"\n{BAR}"); print("MEASURED ON THE REAL NONLINEAR LOOP"); print(BAR)
        for label, use_gpu in (('host assembly', False), ('GPU assembly', True)):
            if use_gpu and not (gpu4d.available() and Pardiso.available()):
                continue
            per_it, setup = bench_loop(use_gpu)
            print(f"  {label:<16}{per_it:6.2f} s per nonlinear iteration "
                  f"(+{setup:5.2f} s once per solve)")
            proj.append((f'measured, {label}', per_it, setup))

    extrapolate(proj)

    if best is not None and best[1] < per_lin and best[2] < 1e-8:
        print(f"\n  *** {best[0]} is {per_lin/best[1]:.1f}x faster than AMG+CG here.")
        if 'PARDISO' in best[0] and Pardiso.available():
            print(f"      ofj4d/solve4d.py already uses it (OFJ_LINSOLVE=auto); "
                  f"set OFJ_LINSOLVE=amg to force the AMG path back.")
        else:
            print(f"      Wire it into ofj4d/solve4d.py (see the Pardiso class there).")
    if A.json:
        json.dump(dict(env=platform.platform(), gpu=found.get('gpu'),
                       dof=As.shape[0], nnz=int(As.nnz),
                       assembly_s=per_assembly, linear_s=per_lin,
                       alternatives=[(r[0], r[1]) for r in rows],
                       per_iteration=[(n, it, su) for n, it, su in proj]),
                  open(A.json, 'w'), indent=1)
        print(f"\n  wrote {A.json}")


if __name__ == '__main__':
    main()
