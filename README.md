# Clear-aligner maxillary distalisation: patient-specific 4D finite element study

Code, run scripts, result tables and figures supporting the manuscript "Miniscrew position, elastic force,
precision-cut location and staging protocol in clear-aligner maxillary whole-arch distalisation:
a patient-specific 4D finite element study".

- `ofj/`, `ofj4d/`: model construction and 4D staging solver
- `pipeline/`: FEBio 4 plug-ins, build and run scripts
- `build4d/`: job queue and supervisor scripts (`build4d/v2/queue.sh` sets the run environment)
- `tools/`, `results/scripts/`: analysis, figure and table scripts
- `results/`: tables, per-run data, figures
- `docs/`: setup, specification, findings log, replication notes

Production-run environment (`build4d/v2/queue.sh`): `OFJ_K_N=1000`, `OFJ_PRESCRIBE=tip`,
`OFJ_TIP_BLEND=0.40`, `OFJ_BONE=rigid`.

Not included: patient-derived meshes, geometry pickles (*.pkl), stress fields (*.npz), run logs.
