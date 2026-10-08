"""Two staging steps of the final model, compared with the shipped checkpoint.

Run once on a new machine, after installing, before committing hours of compute:

    python pipeline/smoke_test.py

Passes when every tooth's crown/root movement matches the reference run
(build4d/final/run_control.pkl) within 1 um after steps 1 and 2.
"""
import os, sys, subprocess, pickle, tempfile, time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(HERE)
REF = 'build4d/final/run_control.pkl'
TOL_MM = 1e-3
KEYS = ('crown_distal', 'crown_buccal', 'crown_occlusal', 'root_distal', 'root_buccal')


def main():
    out = os.path.join(tempfile.mkdtemp(), 'smoke_control.pkl')
    env = dict(os.environ, OFJ_PRESCRIBE='tip', OFJ_TIP_BLEND='0.40')
    t0 = time.time()
    r = subprocess.run([sys.executable, '-u', 'ofj4d/run4d.py', '--group', 'control',
                        '--steps', '2', '--out', out, '--no-resume'],
                       env=env, capture_output=True, text=True)
    for line in (r.stdout + r.stderr).splitlines():
        if 'Resetting the LHS' in line:
            continue
        if (line.startswith(('bone model', 'global nodes', 'reduced dofs', '[control]', '    ['))
                or 'Traceback' in line or 'Error' in line):
            print(line)
    if r.returncode != 0 or not os.path.exists(out):
        print('\n'.join((r.stdout + r.stderr).splitlines()[-25:]))
        print('\nSMOKE TEST FAILED: the run did not complete (output above)')
        return 1
    new = pickle.load(open(out, 'rb'))['hist']
    ref = pickle.load(open(REF, 'rb'))['hist']
    worst, where = 0.0, None
    for s in (0, 1):
        for u, m in ref[s]['meas'].items():
            for k in KEYS:
                d = abs(new[s]['meas'][u][k] - m[k])
                if d > worst:
                    worst, where = d, (s + 1, u, k)
    print(f"\nlargest difference from the reference: {worst*1e3:.4f} um "
          f"(step {where[0]}, tooth index {where[1]}, {where[2]})   [{time.time()-t0:.0f} s]")
    if worst < TOL_MM:
        print('SMOKE TEST PASSED')
        return 0
    print(f'SMOKE TEST FAILED: differs from the reference by more than {TOL_MM*1e3:.0f} um')
    return 1


if __name__ == '__main__':
    sys.exit(main())
