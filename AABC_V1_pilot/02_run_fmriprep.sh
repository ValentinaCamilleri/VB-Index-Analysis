#!/usr/bin/env bash
# Test 2: fMRIPrep 22.1.1 on HCA8889115 V1 BIDS (source images never edited)
set -euo pipefail
if ! type module >/dev/null 2>&1; then
  if [[ -r /etc/profile.d/modules.sh ]]; then source /etc/profile.d/modules.sh; fi
fi
if ! type module >/dev/null 2>&1; then
  echo 'Environment Modules command not found; launch this script from the Pavie login/VNC shell.' >&2
  exit 2
fi
module purge
unset PYTHONHOME PYTHONPATH
module load fmriprep/22.1.1
BASE="$HOME/TARRAGON/Pilot/Pilot_Results_fMRIPrep"
BIDS="$BASE/BIDS/HCA8889115_V1"
OUT="$BASE/HCA8889115/V1/fMRIPrep_output"
WORK="$BASE/HCA8889115/V1/fMRIPrep_work"
LOGS="$BASE/HCA8889115/V1/fMRIPrep_logs"
for path in "$BIDS/dataset_description.json" "$BIDS/sub-HCA8889115/ses-V1/anat/sub-HCA8889115_ses-V1_T1w.nii.gz"; do
  if [[ ! -f "$path" ]]; then echo "MISSING BIDS input: $path; run 01_prepare_bids.py --build first" >&2; exit 2; fi
done
for run in 1 2; do
  for dir in AP PA; do
    bold="$BIDS/sub-HCA8889115/ses-V1/func/sub-HCA8889115_ses-V1_task-rest_dir-${dir}_run-${run}_bold.nii.gz"
    [[ -f "$bold" ]] || { echo "MISSING BOLD: $bold" >&2; exit 2; }
  done
done
LICENSE="${FS_LICENSE:-}"
if [[ -z "$LICENSE" && -n "${FREESURFER_HOME:-}" && -f "$FREESURFER_HOME/license.txt" ]]; then
  LICENSE="$FREESURFER_HOME/license.txt"
fi
if [[ -z "$LICENSE" || ! -f "$LICENSE" ]]; then
  echo 'No FreeSurfer license found. Set FS_LICENSE to your valid accessible license.txt path and rerun.' >&2
  echo 'For example: export FS_LICENSE="$HOME/.licenses/freesurfer/license.txt" (if you have it there)' >&2
  exit 2
fi
# Important: home has a per-user quota; df -h alone does not reveal that quota.
echo '=== Available filesystem space (not quota) ==='
df -h "$BASE"
echo '=== User quota, if supported ==='
quota -s 2>/dev/null || true
echo '=== Current TARRAGON disk usage ==='
du -sh "$HOME/TARRAGON" 2>/dev/null || true
read -r -p 'Confirm sufficient HOME quota and permission to use 20 processes/8 OpenMP threads/100GB on this shared workstation? Type yes: ' response
[[ "$response" == yes ]] || { echo 'Aborted before running fMRIPrep.'; exit 1; }
mkdir -p "$OUT" "$WORK" "$LOGS"
LOG="$LOGS/fmriprep_HCA8889115_V1_$(date -u +%Y%m%dT%H%M%SZ).log"
echo "fMRIPrep version: $(fmriprep --version)"
echo "Logging to: $LOG"
# Keep raw BIDS source and work/output separately on the same HOME filesystem.
# Do not enable ICA-AROMA, filtering, or smoothing here; 36P comes afterward.
fmriprep "$BIDS" "$OUT" participant \
  --participant-label HCA8889115 \
  --fs-license-file "$LICENSE" \
  -w "$WORK" \
  --output-spaces func T1w MNI152NLin2009cAsym:res-2 \
  --nprocs 20 \
  --omp-nthreads 8 \
  --mem_mb 100000 \
  --stop-on-first-crash \
  --resource-monitor \
  --notrack \
  -v 2>&1 | tee "$LOG"
echo "Completed. fMRIPrep derivatives: $OUT"
