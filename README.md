# VB Index Analysis

This repository contains research scripts for preprocessing resting-state functional MRI (fMRI) data and calculating the unnormalised volumetric Vogt–Bailey (VB) Index using VB7 and VB27 neighbourhood configurations.

The scripts were originally developed for the Midnight Scan Club (MSC) dataset and subsequently adapted for a pilot investigation using the Aging Adult Brain Connectome (AABC) dataset.

## Repository Structure

- **`original_MSC_pipeline/`** — The three original processing scripts developed for the Midnight Scan Club dataset, together with their README.
- **`AABC_V1_pilot/`** — Eight files documenting the adaptation of the processing workflow for participant HCA8889115, Visit 1 (V1), including preprocessing, VB calculation, configuration and accompanying documentation.

## Processing Overview

The principal workflow comprises three stages:

1. **Preflight validation:** Verification of fMRIPrep derivatives, brain masks, acquisition metadata and nuisance regressors.
2. **Nuisance regression:** Application of a 36-parameter nuisance regression model combined with 128-second discrete cosine transform (DCT) high-pass filtering, without additional spatial smoothing.
3. **VB Index calculation:** Computation of unnormalised volumetric VB7 and VB27 maps from the denoised functional time series.

The AABC V1 implementation additionally includes scripts for preparing the original imaging data in BIDS format, running fMRIPrep and transforming completed native-space VB maps into MNI standard space.

Although the AABC scripts were adapted from the original MSC implementation, the underlying denoising and VB calculation methodology remains consistent. Dataset-specific modifications include participant and session identification, acquisition handling, file paths and processing configurations.

Please consult the README within each directory for further information regarding processing parameters, execution instructions and adaptation to other datasets.

## Scope

This repository provides the preprocessing and VB calculation components used in the research investigations. It does not constitute a complete implementation of the wider Project TARRAGON research framework, which includes additional analyses beyond the processing stages presented here.

## Citation

If these scripts or the Vogt–Bailey Index are used in research leading to a publication, please cite the original scientific publication describing the VB Index and acknowledge the relevant software and methods.

**Original VB Index publication:** https://doi.org/10.1016/j.neuroimage.2020.117140
