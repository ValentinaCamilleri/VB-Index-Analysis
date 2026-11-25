# Project TARRAGON — Processing Pipelines

Research code for computing unnormalised volumetric Vogt–Bailey Index (VB7 and VB27) maps from resting-state fMRI.

## Repository structure

- `original_MSC_pipeline/`: three original Midnight Scan Club (MSC) processing scripts, plus their README.
- `AABC_V1_pilot/`: eight files from the HCA8889115 Visit 1 AABC adaptation, plus its README.

## Processing overview

Input fMRIPrep derivatives → preflight/manifest validation → joint 36-parameter nuisance regression and 128-second DCT high-pass filtering, with no additional spatial smoothing → unnormalised volumetric VB7/VB27 calculation.

The AABC V1 folder additionally contains data preparation, independent fMRIPrep execution and native-VB-to-MNI registration scripts. Consult each folder's README before execution. The MSC and AABC scripts are **not** interchangeable by changing only the input paths: participant/session/acquisition discovery and configuration must also be reviewed.

## Reproducibility and limitations

- `original_MSC_pipeline` represents the original MSC implementation.
- `AABC_V1_pilot` represents the initial HCA8889115 V1 pilot and is not a generalised V1–V4 pipeline.
- Later V2–V4 analyses used fMRIPrep derivatives provided separately and adapted scripts not included here.
- This archive contains code and documentation only, not neuroimaging data or processed maps.
- Absolute local file paths reflect the original computational environments and require adaptation.
- Confirm software versions, data permissions and registration transforms before re-running.

## Data

Midnight Scan Club: OpenNeuro `ds000224`. AABC/HCP-Aging data access and usage are governed by the dataset provider. Do not commit raw MRI or large derivatives to this repository.
