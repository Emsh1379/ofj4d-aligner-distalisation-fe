#!/usr/bin/env bash
# Rebuild the refined meshes of the mesh-convergence study (not shipped: ~120 MB). A few minutes on a CPU.
cd "$(dirname "$0")/../.."
python tools/smooth_pdl.py --iterations 20 --collapse 0.35 --layers 6 --subdivide 0 --out build4d/meshconv/geom_fe_t6.pkl
python tools/smooth_pdl.py --iterations 20 --collapse 0.35 --layers 3 --subdivide 1 --out build4d/meshconv/geom_fe_h2.pkl
python tools/smooth_pdl.py --iterations 20 --collapse 0.35 --layers 6 --subdivide 1 --out build4d/meshconv/geom_fe_h2t6.pkl
OFJ_GPU=off python tools/volmesh_aligner.py --switches 'pq1.5/15YQ' --step 2 --out build4d/meshconv/aligner_fe_s2.pkl
