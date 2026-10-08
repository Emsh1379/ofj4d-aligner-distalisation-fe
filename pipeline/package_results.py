"""Package the pipeline and its results as ofj4d_final_pipeline/ in a zip (same layout as the input package).

usage: python pipeline/package_results.py [--out ../ofj4d_final_pipeline_results.zip]

Included: all code, docs, figures, run results and logs, repaired meshes, bone mesh, miniscrew positions, job files.
Excluded (regenerable, machine-specific or very large): FEBio working models and node dumps, compiled plug-ins,
Python caches, refined mesh-convergence meshes (rebuild with build4d/meshconv/rebuild.sh), raw-data links and caches.
"""
import os, sys, zipfile, argparse, fnmatch, time

HERE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ap = argparse.ArgumentParser()
ap.add_argument('--out', default=os.path.join(os.path.dirname(HERE), 'ofj4d_final_pipeline_results.zip'))
a = ap.parse_args()

EXCLUDE_DIRS = ['__pycache__', '.venv', 'build4d/febio_run', 'pipeline/febio_cudss/build', 'pipeline/febio_pdl/build',
                'pipeline/febio_contact/build', 'pipeline/febio_pairs/build', 'data/P8', 'data/cache', '.ipynb_checkpoints']
EXCLUDE_FILES = ['build4d/meshconv/*.pkl', 'build4d/febio_val/*.feb', 'build4d/febio_val/*_nodes.txt', 'build4d/febio_val/*_rb.txt',
                 'build4d/febio_val/*.log', 'build4d/febio_val/*.stdout', 'build4d/febio_val/*.time', 'build4d/febio_val/*.xplt',
                 'build4d/febio_val/custom.pkl', 'build4d/febio_val/bench_custom_*.pkl', '*.pyc', '*.so', '*.o', 'build4d/static_info.pkl']

n, total = 0, 0
with zipfile.ZipFile(a.out, 'w', zipfile.ZIP_DEFLATED, compresslevel=6) as z:
    for root, dirs, files in os.walk(HERE):
        rel_root = os.path.relpath(root, HERE)
        dirs[:] = [d for d in dirs if not any(os.path.normpath(os.path.join(rel_root, d)) == os.path.normpath(e) or d == e for e in EXCLUDE_DIRS)]
        for f in files:
            rel = os.path.normpath(os.path.join(rel_root, f))
            if any(fnmatch.fnmatch(rel, pat) for pat in EXCLUDE_FILES):
                continue
            full = os.path.join(root, f)
            if os.path.islink(full):
                continue
            z.write(full, os.path.join('ofj4d_final_pipeline', rel))
            n += 1; total += os.path.getsize(full)
size = os.path.getsize(a.out)
with zipfile.ZipFile(a.out) as z:
    bad = z.testzip()
print(f"{a.out}: {n} files, {total/1e6:.1f} MB uncompressed, {size/1e6:.1f} MB zipped; integrity {'OK' if bad is None else 'FAILED at ' + bad}")
