# Pilot Test Scripts

## Destination layout

Every Test 2 generated item is under `~/TARRAGON/Pilot/Pilot_Results_fMRIPrep/`:

```text
Pilot_Results_fMRIPrep/
├── BIDS/HCA8889115_V1/                 # derived BIDS-format inputs (hardlinks)
└── HCA8889115/V1/
    ├── fMRIPrep_output/                # fMRIPrep derivatives, FreeSurfer
    ├── fMRIPrep_work/                  # fMRIPrep intermediates
    ├── fMRIPrep_logs/                  # fMRIPrep console logs
    └── fMRIPrep_VB/                    # 36P + VB analysis
        ├── Denoising_Output/
        ├── VB_Output/
        └── Scripts/{reports,logs,tmp}/
```

The scripts themselves remain in `Pilot/Scripts/Test2_fMRIPrep/`. The raw AABC downloads and prior `Pilot_Results/` (Test 1) remain untouched. The `config.yaml` documents the same paths; all executable scripts are set to these locations.

## Run order

**Do not run fMRIPrep until BIDS assembly and space/quota/license checks pass.**

1. `module purge; module load python/3.12.9; unset PYTHONHOME PYTHONPATH`
2. `python 01_prepare_bids.py` (inventory only)
3. Inspect output, then `python 01_prepare_bids.py --build` (hard-links NIfTI, writes copied JSON metadata)
4. `command -v bids-validator && bids-validator "$HOME/TARRAGON/Pilot/Pilot_Results_fMRIPrep/BIDS/HCA8889115_V1"` if validator installed. Do not skip errors; fMRIPrep validates BIDS itself as well.
5. `export FS_LICENSE=/path/to/your/valid/FreeSurfer/license.txt` if not already defined.
6. `bash 02_run_fmriprep.sh` (loads fMRIPrep 22.1.1 on Pavie, quota prompt)
7. After successful fMRIPrep: open a new Python-3.12 terminal and run `python 03_preflight_derivatives.py`.
8. Verify all eight `READY` rows. Run `python 04_denoise_36P.py --all --space both`. Do not run if the manifest needs attention.
9. Run `python 05_run_vb_fmriprep.py --run REST1_PA --space MNI152NLin2009cAsym --neighbourhood 7 --dry-run`.
10. Run `python 05_run_vb_fmriprep.py --run REST1_PA --space MNI152NLin2009cAsym --neighbourhood 7`.
11. Run `python 05_run_vb_fmriprep.py --all-runs --space both --neighbourhood both --heartbeat-seconds 10`.

## Scientific choices

- REST1/REST2 and AP/PA are four independent BOLD runs. Matching fieldmaps are separate acquisitions, grouped as `SpinEchoFieldMap1` for REST1 and `SpinEchoFieldMap3` for REST2, as indicated in the source inventory.
- The BIDS conversion uses source JSONs with *no guessed timing/geometry*. It fails if known AP/PA polarity and BIDS phase-encoding values disagree. Images are hardlinked on the same filesystem (no duplicated payload).
- BIDS names `sub-HCA8889115_ses-V1_task-rest_dir-AP_run-1_bold.nii.gz` etc distinguish all four acquisitions.
- fMRIPrep 22.1.1 is run directly using the Pavie module, not the MSC `fmriprep-docker` wrapper. It uses 20 processes, 8 OMP threads, 100000 MB limit, `func`, `T1w`, and 2mm `MNI152NLin2009cAsym` outputs.
- The original 36P method is adapted from the user-provided MSC script: motion (6), WM/CSF/global (3), derivatives/power2/derivative-power2 = 36 regressors, joint 128s DCT high-pass. Initial NSS frames are removed if present; no spatial smoothing. This is not the same cleaning pipeline as HCP's FIX-based data.
- The `func` output and MNI output are not always identical spatial grids. Both are processed independently so the future native-vs-MNI VB comparison can be made without inverse warping already-normalized data.
- **HCP MNI and fMRIPrep MNI152NLin2009cAsym are not necessarily the same reference grid or registration. Do not compute voxelwise HCP–fMRIPrep differences directly without spatial harmonisation.**

## Source preservation

Scripts only read `AABC/HCA8889115/Resting State fMRI Unprocessed` and `AABC/HCA8889115/Structural Unprocessed`. Original downloads and Test 1 results remain unchanged. The BIDS folder only introduces hardlinks to source images and new metadata sidecars. Be aware that hardlinks share bytes with originals; do not edit image files through either path.

## Resource safety

Home disk quota was previously reported as 200 GB, distinct from server free space. Raising the fMRIPrep memory limit to 100000 MB does not increase available disk quota. fMRIPrep intermediate results can occupy many GB. `02_run_fmriprep.sh` prints `quota -s` and `du -sh` and asks for explicit confirmation. The CPU/memory values are fMRIPrep limits, not exclusive reservations; do not approve without enough headroom. No /data/extra usage is configured, per your storage preference.

## Limitations

The actual Pavie files and fMRIPrep module have not been tested by the assistant. Preflight and BIDS validation should be reviewed before triggering resource-intensive preprocessing. `03_preflight_derivatives.py` deliberately fails if expected outputs are missing or if the 36P columns are unusable. Final VB existence checks are not full statistical quality control.


## GitHub handover note (25 November 2026)

These eight source files are the **original uploaded AABC HCA8889115 V1 pilot files**, retained unchanged, apart from this appended README note. They were configured for **the independently run V1 fMRIPrep 22.1.1 processing**, not the later V2–V4 visit processing from fMRIPrep 25.2.5 derivatives.

To process other participants or visits, check hardcoded participant/visit IDs, acquisition labels, BIDS entity patterns, expected manifests, FMRIPrep version and flags, and spatial transforms before running. `06_register_native_vb_to_MNI.sh` is a V1-specific comparison step and should be verified for each new spatial reference.

The preflight, 36P denoising and VB calculation scripts implement the same *downstream* scientific method as the original MSC TARRAGON scripts, with dataset-specific input/output adaptation. This does not imply identical upstream fMRIPrep processing across datasets/visits.

Only source scripts and a parameter configuration are uploaded; no actual fMRI data or outputs are included.
