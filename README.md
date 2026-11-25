# Project Processing Pipelines

Research code for computing unnormalised volumetric Vogt–Bailey Index (VB7 and VB27) maps from resting-state fMRI.

## Repository structure

- `original_MSC_pipeline/`: three original Midnight Scan Club (MSC) processing scripts, plus their README.
- `AABC_V1_pilot/`: eight files from the HCA8889115 Visit 1 AABC adaptation, plus its README.

## Processing overview

Input fMRIPrep derivatives → preflight/manifest validation → joint 36-parameter nuisance regression and 128-second DCT high-pass filtering, with no additional spatial smoothing → unnormalised volumetric VB7/VB27 calculation.

The AABC V1 folder additionally contains data preparation, independent fMRIPrep execution and native-VB-to-MNI registration scripts. Consult each folder's README before execution. The MSC and AABC scripts are **not** interchangeable by changing only the input paths: participant/session/acquisition discovery and configuration must also be reviewed.

