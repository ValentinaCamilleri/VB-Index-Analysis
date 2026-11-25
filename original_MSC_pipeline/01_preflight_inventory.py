#!/usr/bin/env python3
"""MSC rest preflight v2: create layout, inspect fMRIPrep files,
check against the acquisition inventory, and produce read-only reports.

No source data, denoised images, or VB maps are modified or generated.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import re
import sys
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

DEFAULT_FMRIPREP = Path(
    '/mnt/SlowData/sgal0013/Midnight_Scan_Club/'
    'ds000224_fMRIprep_git/OpenNeuroDerivatives/ds000224-fmriprep'
)
DEFAULT_EXP = Path(
    '/mnt/SlowData/sgal0013/Midnight_Scan_Club/'
    'VB_Work_Output/Sephora/20261026_VB_MSC'
)
MNI = 'MNI152NLin2009cAsym'
SPACES = ('native', MNI)
BASE = ('trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z',
        'white_matter', 'csf', 'global_signal')
COLS36 = tuple(f'{base}{suffix}' for base in BASE for suffix in
               ('', '_derivative1', '_power2', '_derivative1_power2'))
ENTITIES = ('sub', 'ses', 'task', 'acq', 'ce', 'dir', 'rec', 'run',
            'echo', 'part', 'recording', 'proc')
FILE_ENTITIES = ENTITIES + ('space', 'res', 'den')
MANIFEST_COLUMNS = [
    'subject', 'session', 'task', 'run', 'space', 'status', 'issues', 'notes',
    'bold', 'mask', 'confounds_tsv', 'confounds_json', 'bold_json',
    'n_volumes', 'tr_seconds', 'voxel_sizes_mm', 'n_mask_voxels',
    'confounds_rows', 'n_nonsteady_state', 'n_retained',
    'n_cosine_regressors', 'n_planned_dct_regressors',
    'n_first_row_derivative_nan_fills', 'n_36p_regressors',
]


def parse_entities(filename: str) -> dict[str, str]:
    out = {}
    for part in filename.split('_'):
        if '-' in part:
            k, v = part.split('-', 1)
            if k in FILE_ENTITIES and v:
                out[k] = v
    return out


def acq_key(ent: dict[str, str]) -> tuple[str, ...]:
    return tuple(ent.get(k, '') for k in ENTITIES)


def src_prefix(ent: dict[str, str]) -> str:
    order = ('sub', 'ses', 'task', 'acq', 'ce', 'dir', 'rec', 'run',
             'echo', 'part', 'recording', 'proc')
    return '_'.join(f'{k}-{ent[k]}' for k in order if k in ent)


def report_csv(path: Path, rows: list[dict], fieldnames: list[str]):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', encoding='utf-8', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(rows)


def existing_json(path: Path):
    if not path.is_file():
        return None
    try:
        with path.open(encoding='utf-8') as f:
            return json.load(f)
    except (OSError, ValueError):
        return None


def init_layout(root: Path):
    for directory in ('Scripts', 'Scripts/reports', 'Scripts/logs', 'Scripts/qc',
                      'Denoising_Output', 'VB_Output'):
        (root / directory).mkdir(parents=True, exist_ok=True)
    for name, title, generator in (
        ('Denoising_Output', 'MSC 36P + 128s DCT high-pass BOLD derivatives',
         'MSC36PDCT'),
        ('VB_Output', 'MSC volumetric unnormalised VB7 / VB27 derivatives',
         'MSCVolumetricVB'),
    ):
        descriptor = root / name / 'dataset_description.json'
        if not descriptor.exists():
            value = {
                'Name': title,
                'BIDSVersion': '1.11.2',
                'DatasetType': 'derivative',
                'GeneratedBy': [{'Name': generator}],
            }
            descriptor.write_text(json.dumps(value, indent=2) + '\n', encoding='utf-8')
        readme = root / name / 'README.md'
        if not readme.exists():
            readme.write_text(
                f'# {title}\n\n'
                'This dataset will be generated from fMRIPrep outputs. '
                'Do not interpret directory creation as successful processing.\n\n'
                'VB map filenames retain the original `.vbi-vol.nii.gz` and '
                '`.vbi-neigh.nii.gz` toolbox extensions; these are custom '
                'nonstandard BIDS suffixes documented for this project.\n',
                encoding='utf-8'
            )


def select_one(folder: Path, prefix: str, pattern_suffix: str) -> list[Path]:
    return sorted(folder.glob(f'{prefix}*{pattern_suffix}'))


def validate_image(bold_path: Path, mask_path: Path, nib) -> tuple[dict, list[str]]:
    """Reads headers and the 3D mask; does NOT decompress entire 4D image."""
    issues = []
    result = {}
    try:
        img = nib.load(str(bold_path))
        if len(img.shape) != 4 or img.shape[3] < 50:
            issues.append('BOLD_SHAPE_INVALID')
        else:
            result['n_volumes'] = int(img.shape[3])
            result['tr_seconds'] = float(img.header.get_zooms()[3])
            if not math.isfinite(result['tr_seconds']) or result['tr_seconds'] <= 0:
                issues.append('BOLD_TR_INVALID')
        result['voxel_sizes_mm'] = ';'.join(f'{v:.5g}' for v in img.header.get_zooms()[:3])
        if not np.isfinite(img.affine).all():
            issues.append('BOLD_AFFINE_INVALID')
    except Exception as e:
        issues.append(f'BOLD_UNREADABLE:{type(e).__name__}')
        return result, issues

    try:
        mask_img = nib.load(str(mask_path))
        if len(mask_img.shape) != 3:
            issues.append('MASK_NOT_3D')
        if mask_img.shape != img.shape[:3]:
            issues.append('MASK_SHAPE_MISMATCH')
        if not np.allclose(mask_img.affine, img.affine, rtol=1e-5, atol=1e-3):
            issues.append('MASK_AFFINE_MISMATCH')
        if len(mask_img.shape) == 3:
            mdata = np.asanyarray(mask_img.dataobj)
            if not np.all(np.isfinite(mdata)):
                issues.append('MASK_NONFINITE')
            if not np.all(np.isin(np.unique(mdata), (0, 1))):
                issues.append('MASK_NONBINARY')
            result['n_mask_voxels'] = int(np.count_nonzero(mdata))
            if result['n_mask_voxels'] < 100:
                issues.append('MASK_TOO_SMALL')
    except Exception as e:
        issues.append(f'MASK_UNREADABLE:{type(e).__name__}')
    return result, issues


def validate_confounds(path: Path) -> tuple[dict, list[str], list[str]]:
    """Check 36P columns, NSS outliers, numerical validity, and DCT availability.

    NOTES are expected conditions to handle during denoising, not blockers.
    We generate an exact 128-second DCT basis *after* initial-volume removal.
    """
    issues = []
    notes = []
    result = {'n_first_row_derivative_nan_fills': 0}
    try:
        with path.open(encoding='utf-8', newline='') as f:
            reader = csv.DictReader(f, delimiter='\t')
            names = reader.fieldnames or []
            rows = list(reader)
    except Exception as e:
        return result, [f'CONFOUNDS_UNREADABLE:{type(e).__name__}'], notes

    result['confounds_rows'] = len(rows)
    absent = [c for c in COLS36 if c not in names]
    result['n_36p_regressors'] = 36 - len(absent)
    if absent:
        issues.append('36P_COLUMNS_MISSING:' + ';'.join(absent))

    cosine_cols = sorted(c for c in names if re.fullmatch(r'cosine_\d+', c))
    result['n_cosine_regressors'] = len(cosine_cols)
    if not cosine_cols:
        notes.append('DCT128_TO_BE_GENERATED_AFTER_TRIM')
    else:
        notes.append('SUPPLIED_DCT_NOT_USED_GENERATE_EXACT_128S_AFTER_TRIM')

    nss_cols = sorted(c for c in names if re.fullmatch(r'non_steady_state_outlier\d+', c))
    flagged = set()
    for col in nss_cols:
        positions = []
        for i, row in enumerate(rows):
            try:
                if float(row[col]) > 0.5:
                    positions.append(i)
            except (ValueError, TypeError):
                issues.append('NSS_COLUMN_NONNUMERIC:' + col)
                break
        flagged.update(positions)
        if len(positions) != 1:
            issues.append('NSS_NOT_ONE_HOT:' + col)

    nskip = len(flagged)
    result['n_nonsteady_state'] = nskip
    result['n_retained'] = len(rows) - nskip
    if flagged != set(range(nskip)):
        issues.append('NSS_NOT_INITIAL_CONTIGUOUS')
    if result['n_retained'] < 50:
        issues.append('INSUFFICIENT_TIMEPOINTS')
    if nskip:
        notes.append(f'INITIAL_NONSTEADY_VOLUMES_TO_REMOVE:{nskip}')

    # Check the *retained* rows only. In runs with no NSS volumes, the first
    # derivative and squared-derivative entries are naturally undefined.
    # Denoising will replace those exact first-row values with 0.0 and record
    # the action. We do NOT silently fill any other missing or invalid cell.
    for col in list(COLS36):
        if col not in names:
            continue
        bad = []
        for i in range(nskip, len(rows)):
            cell = rows[i][col]
            try:
                value = float(cell)
            except (ValueError, TypeError):
                value = float('nan')
            if not math.isfinite(value):
                bad.append(i)
        if bad:
            derivative = col.endswith(('_derivative1', '_derivative1_power2'))
            if nskip == 0 and bad == [0] and derivative:
                result['n_first_row_derivative_nan_fills'] += 1
            else:
                issues.append('CONFOUNDS_NONFINITE:' + col)

    if result['n_first_row_derivative_nan_fills']:
        notes.append('FIRST_ROW_DERIVATIVE_NANS_ZERO_FILL_REQUIRED:' +
                     str(result['n_first_row_derivative_nan_fills']))
    return result, issues, notes


def load_reference(path: Path):
    if not path.is_file():
        return None
    with path.open(encoding='utf-8', newline='') as f:
        ref = list(csv.DictReader(f))
    return {(r['subject'], r['session'], r.get('run', '')): r for r in ref}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fmriprep', type=Path, default=DEFAULT_FMRIPREP)
    parser.add_argument('--experiment', type=Path, default=DEFAULT_EXP)
    parser.add_argument('--expected', type=Path,
                        help='Reference CSV (subject/session/run/n_volumes/tr_seconds); '
                             'defaults to Scripts/expected_rest_acquisitions.csv')
    args = parser.parse_args()
    root = args.fmriprep.expanduser().resolve()
    experiment = args.experiment.expanduser().resolve()
    if not root.is_dir():
        parser.error('fMRIPrep directory does not exist: ' + str(root))
    try:
        import nibabel as nib
    except ImportError:
        parser.error('nibabel is required: use the Python environment running fMRIPrep/VB')
    init_layout(experiment)
    reference_path = args.expected or (experiment / 'Scripts' / 'expected_rest_acquisitions.csv')
    expected = load_reference(reference_path)
    if expected is None:
        print(f'NOTE: reference CSV not found ({reference_path}); checking directly from files.')

    discovered = defaultdict(dict)
    for bold in sorted(root.glob('sub-*/ses-*/func/*_task-rest*_desc-preproc_bold.nii.gz')):
        ent = parse_entities(bold.name)
        space = ent.get('space', 'native')
        if space not in SPACES:
            continue
        key = acq_key(ent)
        discovered[key].setdefault(space, []).append(bold)

    rows = []
    conf_cache = {}
    observed = set()
    for key in sorted(discovered):
        f = dict(zip(ENTITIES, key))
        subj, ses, run = f['sub'], f['ses'], f['run']
        ref_key = ('sub-' + subj, 'ses-' + ses, run)
        observed.add(ref_key)
        prefix = src_prefix({k:v for k,v in f.items() if v})
        folder = root / ('sub-' + subj) / ('ses-' + ses) / 'func'
        confound_matches = select_one(folder, prefix, '_desc-confounds_timeseries.tsv')
        if not confound_matches:
            confound_matches = select_one(folder, prefix, '_desc-confounds_regressors.tsv')
        for space in SPACES:
            matched = discovered[key].get(space, [])
            row = {
                'subject': 'sub-' + subj, 'session': 'ses-' + ses,
                'task': f['task'], 'run': run,
                'space': space, 'status': 'CHECK', 'issues': '', 'notes': '',
            }
            issues = []
            notes = []
            if len(matched) != 1:
                issues.append(f'BOLD_MATCH_COUNT:{len(matched)}')
            if len(confound_matches) != 1:
                issues.append(f'CONFOUNDS_MATCH_COUNT:{len(confound_matches)}')
            if matched:
                bold = matched[0]
                row['bold'] = str(bold)
                bold_prefix = bold.name.rsplit('_desc-preproc_bold.nii.gz', 1)[0]
                candidate_masks = [folder / (bold_prefix + '_desc-brain_mask.nii.gz')]
                # Rebuild suffix from exactly the same entities, excluding a
                # possible BOLD-only 'desc-preproc' part.
                masks = [p for p in candidate_masks if p.is_file()]
                if len(masks) != 1:
                    issues.append(f'MASK_MATCH_COUNT:{len(masks)}')
                else:
                    row['mask'] = str(masks[0])
                    metrics, subissues = validate_image(bold, masks[0], nib)
                    row.update(metrics)
                    issues.extend(subissues)
                metadata_file = folder / (bold_prefix + '_desc-preproc_bold.json')
                if metadata_file.is_file():
                    row['bold_json'] = str(metadata_file)
                    metadata = existing_json(metadata_file)
                    if metadata is None:
                        issues.append('BOLD_JSON_UNREADABLE')
                    elif 'RepetitionTime' in metadata and row.get('tr_seconds'):
                        if abs(float(metadata['RepetitionTime']) - row['tr_seconds']) > 0.001:
                            issues.append('TR_HEADER_JSON_MISMATCH')
                else:
                    # A BIDS inheritance rule may supply this information.
                    row['bold_json'] = ''
            if confound_matches:
                conf = confound_matches[0]
                row['confounds_tsv'] = str(conf)
                sidecar = conf.with_suffix('.json')
                row['confounds_json'] = str(sidecar) if sidecar.is_file() else ''
                if str(conf) not in conf_cache:
                    conf_cache[str(conf)] = validate_confounds(conf)
                metrics, subissues, subnotes = conf_cache[str(conf)]
                row.update(metrics)
                issues.extend(subissues)
                notes.extend(subnotes)
                if row.get('n_volumes') is not None and row.get('confounds_rows') is not None:
                    if row['n_volumes'] != row['confounds_rows']:
                        issues.append('CONFOUNDS_BOLD_LENGTH_MISMATCH')
            if expected is not None:
                ref = expected.get(ref_key)
                if ref is None:
                    issues.append('UNEXPECTED_REST_ACQUISITION')
                elif row.get('n_volumes') and int(ref['n_volumes']) != row['n_volumes']:
                    issues.append('REFERENCE_VOLUME_MISMATCH')
            # DCT regressors will be built from the retained time series (post
            # NSS removal) using this BOLD image's TR. 128 seconds is a period.
            # Formula agrees with the commonly used cosine drift order.
            n_keep = row.get('n_retained')
            tr = row.get('tr_seconds')
            if n_keep is not None and tr is not None and tr > 0:
                row['n_planned_dct_regressors'] = int(math.floor(
                    2.0 * n_keep * tr / 128.0
                ))
            row['issues'] = '|'.join(sorted(set(issues)))
            row['notes'] = '|'.join(sorted(set(notes)))
            row['status'] = 'READY' if not issues else 'CHECK'
            rows.append(row)

    if expected is not None:
        for ref_key, ref in sorted(expected.items()):
            if ref_key not in observed:
                for space in SPACES:
                    rows.append({
                        'subject': ref_key[0], 'session': ref_key[1],
                        'task': 'rest', 'run': ref_key[2], 'space': space,
                        'status': 'CHECK', 'issues': 'EXPECTED_REST_ACQUISITION_MISSING',
                    })

    reportdir = experiment / 'Scripts' / 'reports'
    report_csv(reportdir / 'preflight_manifest.csv', rows, MANIFEST_COLUMNS)
    counts = Counter(r['status'] for r in rows)
    summary = {
        'timestamp_utc': datetime.now(timezone.utc).isoformat(),
        'source_fmriprep': str(root),
        'experiment': str(experiment),
        'reference_csv': str(reference_path) if expected else None,
        'expected_reference_acquisitions': len(expected) if expected else None,
        'discovered_acquisitions': len(discovered),
        'discovered_subjects': len({k[0] for k in discovered}),
        'space_level_rows': len(rows),
        'counts': dict(counts),
        'all_ready': len(rows) == 200 and counts.get('READY', 0) == 200,
        'notes_count': dict(Counter(note for row in rows for note in
                                    str(row.get('notes', '')).split('|') if note)),
        'comments': [
            'This is a header/mask/confound-data preflight; full 4D payload '
            'and residual-QC checks occur at denoising.',
            'Missing cosine columns are NOT a blocker: generate an exact '
            '128-second DCT basis after removing initial NSS volumes.',
            'Expected first-row derivative NaNs are NOT a blocker: denoising '
            'must set only these values to zero and document the adjustment.',
        ],
    }
    (reportdir / 'preflight_summary.json').write_text(json.dumps(summary, indent=2) + '\n',
                                                       encoding='utf-8')
    print(json.dumps(summary, indent=2))
    print('Manifest:', reportdir / 'preflight_manifest.csv')
    if not summary['all_ready']:
        print('Preflight needs inspection. See issues in the manifest.', file=sys.stderr)
        return 2
    print('Preflight complete: 100 acquisitions × 2 spaces = 200 READY rows.')
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
