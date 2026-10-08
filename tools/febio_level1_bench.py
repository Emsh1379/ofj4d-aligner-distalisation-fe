"""Head-to-head wall-clock benchmark on the Level-1 problem: this project's solver
(tools/febio_level1_custom.py) against FEBio 4 + MKL PARDISO on the identical mesh and
load (tools/febio_level1_export.py).  Run on an IDLE machine.

For each thread count, each code is run REPEAT times as a fresh process, sequentially.
Reported: process wall time (includes interpreter / input parsing), and the solver's own
in-solve breakdown.

usage: python tools/febio_level1_bench.py FEB_FILE [--threads 1 8] [--repeat 3] [--json out.json]
"""
import os, sys, subprocess, time, json, re, argparse, statistics
HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FEBIO = os.environ.get('FEBIO', os.path.expanduser('~/febio/febio.sh'))

ap = argparse.ArgumentParser()
ap.add_argument('feb'); ap.add_argument('--threads', type=int, nargs='+', default=[1, 8])
ap.add_argument('--repeat', type=int, default=3); ap.add_argument('--json')
a = ap.parse_args()
feb = os.path.abspath(a.feb); wd = os.path.dirname(feb)
sys.path.insert(0, os.path.join(HERE, 'tools'))


def env_for(n):
    e = dict(os.environ, OMP_NUM_THREADS=str(n), MKL_NUM_THREADS=str(n), OFJ_GPU='off')
    return e


def febio_times(log):
    txt = open(log, errors='ignore').read()
    g = lambda label: float(re.search(r'\s' + label + r'\s*\.*\s*:\s*[\d:]+\s*\(([\d.eE+-]+) sec\)', txt).group(1))
    it = int(re.search(r'Total number of equilibrium iterations\s*\.*\s*:\s*(\d+)', txt).group(1))
    rf = int(re.search(r'Total number of stiffness reformations\s*\.*\s*:\s*(\d+)', txt).group(1))
    return dict(total=g('Total elapsed time'), input=g('Input time'), init=g('Initialization time'),
                solve=g('Solve time'), linear_solver=g('time in linear solver'),
                stiffness=g('evaluating stiffness'), residual=g('evaluating residual'),
                update=g('model update'), io=g(r'IO-time \(plot, dmp, data\)'), iters=it, reforms=rf)


out = {}
for n in a.threads:
    rec = dict(custom=[], febio=[])
    for r in range(a.repeat):
        t0 = time.perf_counter()
        p = subprocess.run([sys.executable, os.path.join(HERE, 'tools/febio_level1_custom.py'), '--repeat', '1',
                            '--out', os.path.join(wd, f'bench_custom_t{n}.pkl')],
                           env=env_for(n), capture_output=True, text=True, cwd=HERE)
        wall = time.perf_counter() - t0
        tm = json.loads([l for l in p.stdout.splitlines() if l.startswith('TIMING ')][-1][7:])
        tm['process_wall'] = wall
        rec['custom'].append(tm)
        t0 = time.perf_counter()
        subprocess.run([FEBIO, '-i', feb, '-silent'], env=env_for(n), capture_output=True, cwd=wd)
        wall = time.perf_counter() - t0
        ft = febio_times(feb[:-4] + '.log'); ft['process_wall'] = wall
        rec['febio'].append(ft)
        print(f"threads {n} rep {r}: custom to-first-solution {tm['to_first_solution']:.2f} s "
              f"(assemble {tm['assemble']:.2f}, factor+solve {tm['solve_first']:.2f}), process {tm['process_wall']:.1f} s | "
              f"FEBio total {ft['total']:.2f} s (stiffness {ft['stiffness']:.2f}, linear {ft['linear_solver']:.2f}, "
              f"{ft['iters']} it / {ft['reforms']} reforms), process {ft['process_wall']:.1f} s", flush=True)
    med = lambda L, k: statistics.median(x[k] for x in L)
    rec['median'] = dict(custom={k: med(rec['custom'], k) for k in rec['custom'][0]},
                         febio={k: med(rec['febio'], k) for k in rec['febio'][0]})
    out[n] = rec
    print(f"threads {n} MEDIAN custom {json.dumps({k: round(v, 3) for k, v in rec['median']['custom'].items()})}")
    print(f"threads {n} MEDIAN febio  {json.dumps({k: round(v, 3) for k, v in rec['median']['febio'].items()})}")
if a.json:
    json.dump(out, open(a.json, 'w'), indent=1); print('->', a.json)
