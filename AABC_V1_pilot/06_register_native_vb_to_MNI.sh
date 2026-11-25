#!/usr/bin/env bash
set -euo pipefail

BASE="$HOME/TARRAGON/Pilot/Pilot_Results_fMRIPrep/HCA8889115/V1"

FUNC="$BASE/fMRIPrep_output/sub-HCA8889115/ses-V1/func"
ANAT="$BASE/fMRIPrep_output/sub-HCA8889115/ses-V1/anat"
VB="$BASE/fMRIPrep_VB/VB_Output"

RUN="REST1_PA"
PREFIX="sub-HCA8889115_ses-V1_task-rest_dir-PA_run-1"

DEST="$BASE/fMRIPrep_VB/Comparison/$RUN"
mkdir -p "$DEST"

# Existing fMRIPrep transforms
BOLD_TO_T1="$FUNC/${PREFIX}_from-scanner_to-T1w_mode-image_xfm.txt"

T1_TO_MNI="$ANAT/sub-HCA8889115_ses-V1_from-T1w_to-MNI152NLin2009cAsym_mode-image_xfm.h5"

# MNI reference grid
REFERENCE="$FUNC/${PREFIX}_space-MNI152NLin2009cAsym_desc-brain_mask.nii.gz"

# Native functional brain mask
NATIVE_MASK="$FUNC/${PREFIX}_desc-brain_mask.nii.gz"

echo "=== Checking inputs ==="

for f in "$BOLD_TO_T1" "$T1_TO_MNI" "$REFERENCE" "$NATIVE_MASK"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Missing $f"
        exit 1
    fi
done

command -v antsApplyTransforms >/dev/null || {
    echo "ERROR: ANTs not loaded"
    exit 1
}

echo "=== Transforming native brain mask ==="

antsApplyTransforms \
    -d 3 \
    -i "$NATIVE_MASK" \
    -r "$REFERENCE" \
    -o "$DEST/REST1_PA_nativeMask_to_MNI.nii.gz" \
    -n NearestNeighbor \
    -t "$T1_TO_MNI" \
    -t "$BOLD_TO_T1"

for N in 7 27; do

    INPUT="$VB/$RUN/native/${PREFIX}_desc-vb${N}_volumetric_output.unnorm.vbi-vol.nii.gz"

    OUTPUT="$DEST/REST1_PA_VB${N}_nativeVB_to_MNI.nii.gz"

    if [[ ! -f "$INPUT" ]]; then
        echo "ERROR: Missing VB${N} input: $INPUT"
        exit 1
    fi

    if [[ -e "$OUTPUT" ]]; then
        echo "Existing output preserved: $OUTPUT"
        continue
    fi

    echo "=== Transforming native VB${N} ==="

    antsApplyTransforms \
        -d 3 \
        -i "$INPUT" \
        -r "$REFERENCE" \
        -o "$OUTPUT" \
        -n Linear \
        -t "$T1_TO_MNI" \
        -t "$BOLD_TO_T1"

    echo "Completed VB${N}"
done

echo "=== Registration finished ==="
echo "Output directory: $DEST"
