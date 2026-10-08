"""Validate and report the three 70-step staging runs."""
import sys, os, pickle
import numpy as np
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(HERE,'ofj4d')); sys.path.insert(0, os.path.join(HERE,'ofj')); os.chdir(HERE)
import spec

RUNS = {'control':     'build4d/run_control.pkl',
        'buccal_tad':  'build4d/run_buccal_tad.pkl',
        'palatal_tad': 'build4d/run_palatal_tad.pkl'}
# the paper's reported values, for the comparison at the end
PAPER = {
    'step10_17':   {'control': (0.97, 5.92), 'buccal_tad': (1.11, 6.75),
                    'palatal_tad': (1.11, 6.75)},          # 1.09-1.13 mm, 6.63-6.87 deg
    'step70_17':   {'control': 0.73, 'palatal_tad': 1.31},
    'step70_16':   {'control': 0.62, 'palatal_tad': 1.00},
    'step70_inc':  {'control': 0.87, 'buccal_tad': 1.71, 'palatal_tad': 1.88},
    'drift_20_40': {'control': 0.89},                       # others 0.61-0.76
}
BAR = '='*78
problems = []


def load(path):
    return pickle.load(open(path, 'rb'))


def check_complete(name, D, n_steps=70):
    steps = [h['step'] for h in D['hist']]
    ok = steps == list(range(1, n_steps+1))
    print(f"  {name:<12} steps {min(steps)}-{max(steps)} ({len(steps)} of {n_steps})"
          f"   {'complete' if ok else 'INCOMPLETE'}")
    if not ok:
        problems.append(f"{name}: step list is not 1..{n_steps} ({len(steps)} entries)")
    teeth = set(D['hist'][-1]['meas'].keys())
    want = set(spec.HALF_ARCH_UNN)
    if teeth != want:
        problems.append(f"{name}: teeth {sorted(teeth)} != {sorted(want)}")
    # every measured quantity finite, at every step
    bad = 0
    for h in D['hist']:
        for u, m in h['meas'].items():
            for k, v in m.items():
                a = np.asarray(v, dtype=float)
                if not np.isfinite(a).all():
                    bad += 1
        for i in h['info']:
            if not all(np.isfinite([i['Fc'], i['du'], i['emax'], i['emax_abs']])):
                bad += 1
    if bad:
        problems.append(f"{name}: {bad} non-finite measured values")
    print(f"  {'':<12} non-finite values: {bad}")
    return ok


def convergence(name, D):
    it = np.array([[i['nit'] for i in h['info']] for h in D['hist']])       # (70,2)
    du = np.array([[i['du'] for i in h['info']] for h in D['hist']])
    cap = np.array([[not i['converged'] for i in h['info']] for h in D['hist']])
    Fc = np.array([[i['Fc'] for i in h['info']] for h in D['hist']])
    e99 = np.array([[i['emax'] for i in h['info']] for h in D['hist']])
    print(f"\n  {name}")
    print(f"    iterations/solve   min {it.min():3d}  median {int(np.median(it)):3d}  "
          f"max {it.max():3d}   (cap is 35)")
    print(f"    solves at the cap  {cap.sum()} of {cap.size} "
          f"({100*cap.sum()/cap.size:.0f}%)")
    print(f"    final du/solve     median {np.median(du):.2e}  max {du.max():.2e}"
          f"   (tolerance 1e-3)")
    print(f"    contact force      median {np.median(Fc):6.3f} N  max {Fc.max():6.3f} N")
    print(f"    PDL p99 strain     median {100*np.median(e99):5.2f}%  max {100*e99.max():5.2f}%")
    conv = du[~cap]
    if len(conv):
        print(f"    of solves that converged, final du: max {conv.max():.2e}")
    if du.max() > 1e-2:
        problems.append(f"{name}: a solve ended with du={du.max():.1e}, far above tolerance")
    if e99.max() > 0.5:
        problems.append(f"{name}: p99 PDL strain reached {100*e99.max():.0f}%")
    return dict(it=it, du=du, cap=cap, Fc=Fc, e99=e99)


def movement(D):
    """Cumulative crown/root movement per tooth at the final step."""
    return D['hist'][-1]['meas']


def at_step(D, step):
    return D['hist'][step-1]['meas']


def expression(name, D):
    print(f"\n  {name}: prescribed vs achieved (crown, end of run)")
    print(f"    {'tooth':<7}{'steps':>10}{'prescribed':>12}{'achieved':>11}"
          f"{'rate':>8}{'  direction'}")
    m = movement(D)
    rows = []
    for u in spec.HALF_ARCH_UNN:
        a, b = spec.V_PATTERN[u]
        want = spec.total_prescribed(u)
        if u in (7, 8):                       # incisors move palatally
            got, lab = -m[u]['crown_buccal'], 'palatal'
        else:
            got, lab = m[u]['crown_distal'], 'distal'
        rate = got/want if want else float('nan')
        rows.append((spec.FDI[u], want, got, rate))
        print(f"    FDI {spec.FDI[u]:<3}{f'{a}-{b}':>10}{want:>10.2f} mm{got:>8.3f} mm"
              f"{100*rate:>7.0f}%  {lab}")
    return rows


def anchorage(runs, steps=(10, 20, 30, 70)):
    """The paper's comparison: molar distalisation and anterior anchorage loss."""
    print(f"\n{BAR}\nANCHORAGE COMPARISON (the paper's endpoint)\n{BAR}")
    print("  crown movement, mm; distal +, palatal + means palatal for incisors\n")
    hdr = f"  {'step':<6}{'tooth':<8}"
    for g in runs: hdr += f"{g:>14}"
    hdr += f"{'best - control':>16}"
    print(hdr)
    out = {}
    ctrl = 'control' if 'control' in runs else list(runs)[0]
    for st in steps:
        for u in spec.HALF_ARCH_UNN:
            vals = {}
            for g, D in runs.items():
                m = at_step(D, st)[u]
                vals[g] = -m['crown_buccal'] if u in (7, 8) else m['crown_distal']
            gain = max(vals.values()) - vals[ctrl]
            out[(st, u)] = vals
            print(f"  {st:<6}{'FDI '+str(spec.FDI[u]):<8}"
                  + ''.join(f"{vals[g]:>14.4f}" for g in runs)
                  + f"{gain:>16.4f}")
        print()
    return out


def tipping(name, D):
    """Crown vs root movement: bodily movement or tipping?"""
    print(f"\n  {name}: movement type at the end of the run")
    print(f"    {'tooth':<8}{'crown':>10}{'root':>10}{'root/crown':>12}{'tip (deg)':>11}")
    m = movement(D)
    for u in spec.HALF_ARCH_UNN:
        if u in (7, 8):
            c, r = -m[u]['crown_buccal'], -m[u]['root_buccal']
        else:
            c, r = m[u]['crown_distal'], m[u]['root_distal']
        ratio = r/c if abs(c) > 1e-6 else float('nan')
        print(f"    FDI {spec.FDI[u]:<4}{c:>10.3f}{r:>10.3f}{ratio:>12.2f}"
              f"{m[u]['rot_coronal']:>11.2f}")


def compare_paper(runs):
    """Side by side with the values Mao et al. report."""
    print(f"\n{BAR}\nAGAINST THE PAPER\n{BAR}")

    def row(label, paper, ours, unit='mm'):
        p = f"{paper:7.2f}" if paper is not None else "      -"
        print(f"  {label:<44}{p} {unit:<4}{ours:8.3f} {unit}")

    print("\n  second molar (FDI 17) at step 10 -- distal movement")
    for g in runs:
        pv = PAPER['step10_17'].get(g)
        row(f"{g}", pv[0] if pv else None, runs[g]['hist'][9]['meas'][2]['crown_distal'])
    print("\n  second molar (FDI 17) at step 10 -- tipping")
    for g in runs:
        pv = PAPER['step10_17'].get(g)
        row(f"{g}", pv[1] if pv else None,
            runs[g]['hist'][9]['meas'][2]['rot_coronal'], 'deg')
    print("\n  mesial drift of FDI 17 over steps 20-40 (paper: 0.89 control, 0.61-0.76 others)")
    for g in runs:
        a = runs[g]['hist'][19]['meas'][2]['crown_distal']
        b = runs[g]['hist'][39]['meas'][2]['crown_distal']
        row(f"{g}", PAPER['drift_20_40'].get(g), a-b)
    print("\n  at step 70 -- FDI 17 / FDI 16 distal, and incisor retraction")
    for g in runs:
        m = runs[g]['hist'][-1]['meas']
        row(f"{g}  FDI 17", PAPER['step70_17'].get(g), m[2]['crown_distal'])
        row(f"{g}  FDI 16", PAPER['step70_16'].get(g), m[3]['crown_distal'])
        row(f"{g}  FDI 11 retraction", PAPER['step70_inc'].get(g), -m[8]['crown_buccal'])
    print("\n  group ordering on incisor retraction at step 70")
    order = sorted(((g, -runs[g]['hist'][-1]['meas'][8]['crown_buccal']) for g in runs),
                   key=lambda x: -x[1])
    print(f"    paper:  palatal_tad > buccal_tad > control")
    print(f"    ours:   {' > '.join(g for g, _ in order)}")


def main():
    # optional prefix, e.g. build4d/variants/tip_run_  ->  build4d/variants/tip_run_<group>.pkl
    paths = ({g: f"{sys.argv[1]}{g}.pkl" for g in RUNS} if len(sys.argv) > 1 else RUNS)
    print(BAR); print("COMPLETENESS"); print(BAR)
    runs = {}
    for g, p in paths.items():
        if not os.path.exists(p):
            problems.append(f"{g}: {p} missing"); print(f"  {g:<12} MISSING {p}"); continue
        D = load(p)
        if D.get('group') != g:
            problems.append(f"{g}: checkpoint says group={D.get('group')}")
        check_complete(g, D)
        runs[g] = D

    print(f"\n{BAR}"); print("CONVERGENCE"); print(BAR)
    diag = {g: convergence(g, D) for g, D in runs.items()}

    print(f"\n{BAR}"); print("EXPRESSION"); print(BAR)
    for g, D in runs.items():
        expression(g, D)

    print(f"\n{BAR}"); print("MOVEMENT TYPE"); print(BAR)
    for g, D in runs.items():
        tipping(g, D)

    if len(runs) >= 2:
        anchorage(runs)
    compare_paper(runs)

    print(f"{BAR}"); print("PROBLEMS FOUND"); print(BAR)
    if problems:
        for p in problems: print(f"  ! {p}")
    else:
        print("  none: both runs complete, all values finite, all solves within tolerance")
    return runs, diag


if __name__ == '__main__':
    main()
