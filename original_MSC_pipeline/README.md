# Project TARRAGON — original MSC processing pipeline

Original reference scripts for applying the volumetric, unnormalised Vogt–Bailey (VB) Index to the Midnight Scan Club (MSC, OpenNeuro ds000224) fMRIPrep derivatives. These are the **three original source scripts** supplied with the handover, preserved without editing; this README was added for GitHub.

## Workflow

1. `01_preflight_inventory.py` — discover/check fMRIPrep derivatives, masks, confounds and acquisition structure; produce `Scripts/reports/preflight_manifest.csv`.
2. `02_denoise_36P.py` — 36-parameter nuisance regression and 128-second discrete-cosine high-pass filtering in a joint regression model; no added spatial smoothing; retain provenance and residual quality-control output.
3. `03_run_vb_simple.py` — invoke `vb_tool volumetric` on denoised BOLD to produce unnormalised VB7 and VB27 maps; `--little-neighbourhood` selects VB7.

## Methodological details

The 36P model includes the six motion parameters, white matter, CSF, and global signals, plus first derivatives, squares, and squared first derivatives. Initial non-steady-state volumes are removed where flagged. Voxelwise VB calculation is performed on the denoised images without additional spatial smoothing.

## Important reuse notes

The scripts contain original MSC-specific default paths and participant/session defaults (`sub-MSC01`, `ses-func01`), and expect the corresponding original fMRIPrep directory structure. Update the defaults or use the available command-line options when working in another environment. The original scripts have **not** been generalised, and the files alone are not a self-contained reproduction of the source dataset, fMRIPrep outputs, VB toolbox, or Python environment. The original research `config.yaml` was provided separately and is not counted among these three source scripts.

No datasets, participant images, passwords, SSH keys, or local machine secrets are included.
