#!/usr/bin/env python3
"""AABC V1: 36P + 128-second DCT high-pass denoising of fMRIPrep outputs.

Consumes the READY rows of Scripts/reports/preflight_manifest.csv produced by
03_preflight_derivatives.py. Defaults to HCA8889115 session V1 (both native and MNI spaces).

This script never writes to the source fMRIPrep dataset. A run is written to
Denoising_Output only if both data and residual QC checks pass. Uses on-disk
memory maps for large MNI images and validates provenance on resume.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import shutil
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

DEFAULT_EXP = Path.home() / 'TARRAGON/Pilot/Pilot_Results_fMRIPrep/HCA8889115/V1/fMRIPrep_VB'
MNI = 'MNI152NLin2009cAsym'
BASE = ('trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z',
        'white_matter', 'csf', 'global_signal')
COLS36 = tuple(f'{base}{suffix}' for base in BASE for suffix in
               ('', '_derivative1', '_power2', '_derivative1_power2'))
METHOD = 'aabc-fmriprep36p-dct128-nss-trim-v1'
DESC = '36P128s'
CUT_SEC = 128.0


def utcnow():
    return datetime.now(timezone.utc).isoformat()


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path: Path):
    with path.open(encoding='utf-8') as f:
        return json.load(f)


def atomic_json(path: Path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.partial')
    tmp.write_text(json.dumps(obj, indent=2, allow_nan=False) + '\n', encoding='utf-8')
    os.replace(tmp, path)


def load_manifest(path: Path):
    with path.open(newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError('Empty preflight manifest: ' + str(path))
    if not all(r.get('status') == 'READY' for r in rows):
        raise RuntimeError('Preflight manifest contains non-READY rows. Resolve before processing.')
    seen = set()
    for row in rows:
        k = (row['subject'], row['session'], row.get('run_id', ''), row['space'])
        if k in seen:
            raise RuntimeError('Duplicate manifest entry: ' + str(k))
        seen.add(k)
    return rows


def intended_rows(rows, args):
    if args.all:
        selected = rows
    else:
        selected = [r for r in rows if r['subject'] == args.subject and
                    (args.all_sessions or r['session'] == args.session)]
    if args.space != 'both':
        selected = [r for r in selected if r['space'] == args.space]
    return sorted(selected, key=lambda r: (r['subject'], r['session'], r.get('run_id',''), r['space']))


def read_design(confounds_path: Path, n_volumes: int, tr: float,
                expected_skip: int, expected_retained: int):
    with confounds_path.open(newline='', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')
        names = reader.fieldnames or []
        rows = list(reader)
    if len(rows) != n_volumes:
        raise ValueError(f'Confounds rows ({len(rows)}) != image volumes ({n_volumes})')
    missing = sorted(set(COLS36) - set(names))
    if missing:
        raise ValueError('Missing 36P columns: ' + ','.join(missing))
    nss_cols = [c for c in names if re.fullmatch(r'non_steady_state_outlier\d+', c)]
    detected = set()
    for c in nss_cols:
        flagged = []
        for t, row in enumerate(rows):
            try:
                v = float(row[c])
            except (ValueError, TypeError):
                raise ValueError(f'NSS field non-numeric: {c} volume {t}') from None
            if not np.isfinite(v):
                raise ValueError(f'NSS field nonfinite: {c} volume {t}')
            if v > 0.5:
                flagged.append(t)
        if len(flagged) != 1:
            raise ValueError(f'NSS not one-hot: {c}')
        detected.add(flagged[0])
    nskip = len(detected)
    if detected != set(range(nskip)):
        raise ValueError('NSS volumes are not all initial and contiguous')
    if nskip != expected_skip or n_volumes - nskip != expected_retained:
        raise ValueError('NSS/retained count differs from the validated manifest')
    nt = n_volumes - nskip
    if nt < 50:
        raise ValueError('Insufficient retained volumes')

    x = np.empty((nt, len(COLS36)), dtype=np.float64)
    fills = []
    for j, col in enumerate(COLS36):
        for ti, source_row in enumerate(range(nskip, n_volumes)):
            try:
                v = float(rows[source_row][col])
            except (ValueError, TypeError):
                v = float('nan')
            if not np.isfinite(v):
                allowed = (nskip == 0 and source_row == 0 and
                           col.endswith(('_derivative1', '_derivative1_power2')))
                if not allowed:
                    raise ValueError(f'Unexpected nonfinite confound {col} at index {source_row}')
                v = 0.0
                fills.append(col)
            x[ti, j] = v
    # Scale nuisance variables to avoid ill-conditioning caused by differences
    # in units and quadratic expansions. Demeaning removes their intercept.
    mu = x.mean(axis=0)
    sd = x.std(axis=0)
    constant = sd < 1.e-12
    z = (x - mu) / np.where(constant, 1.0, sd)
    z[:, constant] = 0.0

    # DCT-II basis: frequencies k/(2*T) up to 1/128 Hz. Identical convention
    # to nilearn.glm.first_level.design_matrix.make_first_level_design_matrix.
    ndct = int(np.floor(2 * nt * tr / CUT_SEC))
    ndct = min(ndct, nt - 1)
    t = np.arange(nt, dtype=np.float64)
    cos = (np.sqrt(2.0 / nt) * np.cos(
        np.pi * (t[:, None] + 0.5) *
        np.arange(1, ndct + 1, dtype=np.float64)[None, :] / nt))
    dct_cols = [f'cosine_{i:02d}' for i in range(ndct)]
    columns = list(COLS36) + dct_cols + ['intercept']
    design = np.column_stack((z, cos, np.ones(nt)))
    u, s, _ = np.linalg.svd(design, full_matrices=False)
    rank_tolerance = max(design.shape) * np.finfo(float).eps * s[0]
    rank = int(np.count_nonzero(s > rank_tolerance))
    q = np.ascontiguousarray(u[:, :rank])
    if rank < 5 or nt <= rank + 10:
        raise ValueError(f'Invalid design rank {rank}, retained volumes {nt}')
    stats = dict(removed_initial_volumes=nskip, retained_volumes=nt,
                 first_row_derivative_zero_fills=fills,
                 dct_regressors=ndct, regressors_36p=list(COLS36),
                 constant_confounds=[COLS36[i] for i in np.flatnonzero(constant)],
                 rank=rank, total_design_columns=design.shape[1],
                 effective_residual_degrees_of_freedom=nt-rank,
                 column_names=columns, confound_means=mu.tolist(),
                 confound_std=sd.tolist())
    return design, q, stats


def expected_output(row, exp):
    source = Path(row['bold'])
    if not source.name.endswith('_desc-preproc_bold.nii.gz'):
        raise ValueError('Unexpected source BOLD naming: ' + source.name)
    name = source.name.replace('_desc-preproc_bold.nii.gz',
                               f'_desc-{DESC}_bold.nii.gz')
    folder = exp / 'Denoising_Output' / row['subject'] / row['session'] / 'func'
    bold_path = folder / name
    json_path = folder / name.replace('.nii.gz', '.json')
    return bold_path, json_path


def series_prefix(row):
    # Preserve both run and direction entities: AP and PA of REST1 have
    # distinct confounds and must NEVER share a design-regressor file.
    filename = Path(row['confounds_tsv']).name
    for ending in ('_desc-confounds_timeseries.tsv', '_desc-confounds_regressors.tsv'):
        if filename.endswith(ending):
            return filename[:-len(ending)]
    raise ValueError('Unexpected confounds filename: ' + filename)


def shared_design_path(row, exp):
    folder = exp / 'Denoising_Output' / row['subject'] / row['session'] / 'func'
    return folder / f'{series_prefix(row)}_desc-{DESC}_regressors.tsv'


def get_mask(path, image):
    import nibabel as nib
    mask_img = nib.load(str(path))
    if len(mask_img.shape) != 3 or mask_img.shape != image.shape[:3]:
        raise ValueError('Mask/image shape mismatch')
    if not np.allclose(mask_img.affine, image.affine, atol=1.e-3, rtol=1.e-5):
        raise ValueError('Mask/image affine mismatch')
    data = np.asarray(mask_img.dataobj)
    if not np.all(np.isfinite(data)) or not np.all(np.isin(np.unique(data), [0, 1])):
        raise ValueError('Mask contains nonbinary or nonfinite values')
    mask_ix = np.flatnonzero(data.reshape(-1, order='F') == 1)
    if len(mask_ix) < 100:
        raise ValueError('Mask contains too few valid voxels')
    return mask_ix


def candidate_is_complete(bold_path, json_path, expected_shape, fingerprint):
    if not bold_path.is_file() or not json_path.is_file():
        return False
    try:
        import nibabel as nib
        meta = read_json(json_path)
        img = nib.load(str(bold_path))
        if (meta.get('ProcessingFingerprint') != fingerprint or
                tuple(img.shape) != tuple(expected_shape) or
                int(meta.get('RetainedVolumes', -1)) != expected_shape[-1] or
                not meta.get('QC', {}).get('passed', False)):
            return False
        # Read a few slices; full payload validation only occurs at generation.
        for t in {0, expected_shape[-1]//2, expected_shape[-1]-1}:
            if not np.isfinite(np.asarray(img.dataobj[..., t])).all():
                return False
        return True
    except (OSError, ValueError, KeyError, TypeError):
        return False


def decompress_if_needed(source_path: Path, tempdir: Path):
    path = tempdir / 'source_bold.nii'
    if source_path.name.endswith('.nii.gz'):
        with gzip.open(source_path, 'rb') as inp, path.open('wb') as dest:
            shutil.copyfileobj(inp, dest, length=16 * 1024 * 1024)
    elif source_path.suffix == '.nii':
        # The source preflight currently supplies .nii.gz; future-proofing.
        shutil.copyfile(source_path, path)
    else:
        raise ValueError('Unsupported source image format')
    return path


def write_design(path, design, names, fingerprint):
    # The same acquisition-level model is used in both spaces. If it exists,
    # ensure it is identical rather than silently overwriting it.
    import io
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter='\t', lineterminator='\n')
    writer.writerow(names)
    for values in design:
        writer.writerow([f'{float(v):.17g}' for v in values])
    payload = buf.getvalue()
    if path.exists():
        existing = path.read_text(encoding='utf-8')
        if existing != payload:
            raise RuntimeError('Existing design TSV differs from current model: ' + str(path))
    else:
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_name(path.name + '.partial')
        tmp.write_text(payload, encoding='utf-8')
        os.replace(tmp, path)
    return digest(payload.encode('utf-8'))


def ensure_scratch_space(scratch, raw_nbytes, mask_voxels, nret):
    # Approximate scratch need = decompressed input + denoised full 4D
    # + time x masked-voxel signal matrix, plus conservative overhead.
    minimum = int(2 * raw_nbytes + 4 * mask_voxels * nret + 1_000_000_000)
    free = shutil.disk_usage(scratch).free
    if free < minimum:
        raise RuntimeError(
            f'Insufficient scratch space: {free/1e9:.1f} GB free; '
            f'approximately {minimum/1e9:.1f} GB required. '
            'Pass --scratch-dir to a filesystem with sufficient space.')
    return minimum


def perform_one(row, args, exp, scratch):
    import nibabel as nib
    src_path = Path(row['bold'])
    mask_path = Path(row['mask'])
    conf_path = Path(row['confounds_tsv'])
    input_img = nib.load(str(src_path))
    if len(input_img.shape) != 4:
        raise ValueError('Input BOLD is not 4D')
    original_vol = int(input_img.shape[-1])
    nt_expected = int(row['n_retained'])
    tr = float(input_img.header.get_zooms()[3])
    if (original_vol != int(row['n_volumes']) or
            not np.isclose(tr, float(row['tr_seconds']), atol=1.e-4)):
        raise ValueError('Input header differs from manifest')
    if not np.isfinite(tr) or tr <= 0:
        raise ValueError('TR invalid')
    mask_ix = get_mask(mask_path, input_img)
    design, q, stats = read_design(conf_path, original_vol, tr,
                                   int(row['n_nonsteady_state']), nt_expected)
    if stats['dct_regressors'] != int(row['n_planned_dct_regressors']):
        raise ValueError('Planned DCT count differs from manifest')

    bold_path, json_path = expected_output(row, exp)
    input_paths = [str(src_path), str(mask_path), str(conf_path)]
    fingerprint_payload = dict(method=METHOD, input_files=[
        {'path': str(p), 'size': p.stat().st_size,
         'mtime_ns': p.stat().st_mtime_ns} for p in (src_path, mask_path, conf_path)],
        tr=tr, nskip=stats['removed_initial_volumes'],
        original_volumes=original_vol, dct_period=128,
        desc=DESC, mask_voxels=len(mask_ix))
    fingerprint = digest(json.dumps(fingerprint_payload, sort_keys=True).encode())
    expected_shape = input_img.shape[:3] + (nt_expected,)
    if candidate_is_complete(bold_path, json_path, expected_shape, fingerprint):
        print(f'SKIP validated: {bold_path}')
        return dict(status='SKIP_VALID', output=str(bold_path),
                    subject=row['subject'], session=row['session'], space=row['space'])
    if (bold_path.exists() or json_path.exists()) and not args.overwrite:
        raise RuntimeError('Existing incomplete or mismatched output; refusing overwrite: '
                           + str(bold_path) + '. Inspect or use --overwrite.')
    raw_nbytes = int(np.prod(input_img.shape) * 4)
    min_scratch = ensure_scratch_space(scratch, raw_nbytes, len(mask_ix), nt_expected)
    print(f'PROCESS {row["subject"]} {row["session"]} {row["space"]}: '
          f'{original_vol} -> {nt_expected} volumes, '
          f'{len(mask_ix):,} in-mask voxels, '
          f'DCT={stats["dct_regressors"]}, rank={stats["rank"]}; '
          f'scratch allowance ~{min_scratch/1e9:.1f} GB', flush=True)
    if args.dry_run:
        return dict(status='DRY_RUN', output=str(bold_path),
                    subject=row['subject'], session=row['session'], space=row['space'])

    out_dir = bold_path.parent
    out_dir.mkdir(parents=True, exist_ok=True)
    design_path = shared_design_path(row, exp)
    design_hash = write_design(design_path, design, stats['column_names'], fingerprint)
    design_json_path = design_path.with_suffix('.json')
    if not design_json_path.exists():
        atomic_json(design_json_path, {
            'Model': METHOD, 'Description': 'Joint 36P nuisance and 128-second DCT high-pass model',
            'Columns': stats['column_names'], 'ColumnScaling': '36P demeaned and divided by population SD; constants set to zero',
            'TRSeconds': tr, 'RemovedInitialVolumes': stats['removed_initial_volumes'],
            'RetainedVolumes': nt_expected, 'DCTPeriodSeconds': CUT_SEC,
            'DCTColumns': stats['dct_regressors'], 'DesignRank': stats['rank'],
            'FirstRowDerivativeNaNZeroFillColumns': stats['first_row_derivative_zero_fills'],
            'DesignSHA256': design_hash, 'SourceConfounds': str(conf_path)
        })
    elif read_json(design_json_path).get('DesignSHA256') != design_hash:
        raise RuntimeError('Existing design metadata hash does not match design TSV')
    nxyz = int(np.prod(input_img.shape[:3]))
    chunk = args.chunk_voxels
    with tempfile.TemporaryDirectory(prefix='aabc36p_', dir=str(scratch)) as temp_name:
        temp = Path(temp_name)
        local_nii = decompress_if_needed(src_path, temp)
        in_img = nib.load(str(local_nii), mmap='r')
        nmask = len(mask_ix)
        y = np.memmap(str(temp / 'Y_timeseries.dat'), dtype=np.float32,
                      mode='w+', shape=(nt_expected, nmask))
        for t in range(nt_expected):
            frame = np.asarray(in_img.dataobj[..., t+stats['removed_initial_volumes']],
                               dtype=np.float32)
            selected = frame.reshape(-1, order='F')[mask_ix]
            if not np.isfinite(selected).all():
                raise ValueError(f'Nonfinite masked BOLD data at retained frame {t}')
            y[t, :] = selected
        y.flush()
        o = np.memmap(str(temp / 'output_4d.dat'), dtype=np.float32,
                      mode='w+', shape=(nxyz, nt_expected), order='F')
        max_relative_projection = 0.0
        count_constant = 0
        std_before_sum, std_after_sum = 0.0, 0.0
        valid_samples = 0
        for start in range(0, nmask, chunk):
            stop = min(start+chunk, nmask)
            original = np.asarray(y[:, start:stop], dtype=np.float64)
            res = original - q @ (q.T @ original)
            if not np.isfinite(res).all():
                raise ValueError(f'Nonfinite residuals in mask chunk {start}:{stop}')
            original_sd = np.std(original, axis=0)
            residual_sd = np.std(res, axis=0)
            count_constant += int(np.count_nonzero(residual_sd < 1e-8))
            std_before_sum += float(original_sd.sum())
            std_after_sum += float(residual_sd.sum())
            valid_samples += stop-start
            proj = q.T @ res
            proj_rms = np.sqrt(np.mean(proj**2, axis=0))
            # Normalized projection of residuals onto modeled nuisance space.
            rel = proj_rms / np.maximum(original_sd, 1e-8)
            max_relative_projection = max(max_relative_projection, float(rel.max()))
            if float(rel.max()) > 1e-5:
                raise ValueError(f'Nuisance projection QC failed: relative={rel.max():.6g}')
            spatial_idx = mask_ix[start:stop]
            for t in range(nt_expected):
                o[spatial_idx, t] = res[t, :].astype(np.float32)
        o.flush()
        header = in_img.header.copy()
        if 'toffset' in header:
            header['toffset'] = float(header['toffset']) + stats['removed_initial_volumes'] * tr
        out_data = o.reshape(expected_shape, order='F')
        final_img = nib.Nifti1Image(out_data, in_img.affine, header=header)
        final_img.header.set_slope_inter(1.0, 0.0)
        # Save first to a .nii.gz file in the destination directory, then
        # perform post-write checks and atomically promote to final filename.
        pending = bold_path.with_name(bold_path.name.replace('.nii.gz', '.partial.nii.gz'))
        try:
            nib.save(final_img, str(pending))
            test_img = nib.load(str(pending))
            if tuple(test_img.shape) != expected_shape:
                raise RuntimeError('Output saved with incorrect NIfTI shape')
            if not np.allclose(test_img.affine, input_img.affine, atol=1e-3, rtol=1e-5):
                raise RuntimeError('Output spatial affine changed')
            for ti in (0, nt_expected//2, nt_expected-1):
                if not np.isfinite(np.asarray(test_img.dataobj[..., ti])).all():
                    raise RuntimeError('Saved NIfTI contains nonfinite samples')
            qc = {
                'passed': True,
                'MeanVoxelTemporalSDPre': std_before_sum/valid_samples,
                'MeanVoxelTemporalSDPost': std_after_sum/valid_samples,
                'NearConstantResidualVoxels': count_constant,
                'MaxRelativeNuisanceProjection': max_relative_projection,
                'MaskedVoxelCount': nmask,
                'OutputGeometryMatchesSource': True,
                'SampledSavedFramesFinite': True,
            }
            source_meta = read_json(Path(row['bold_json'])) if Path(row['bold_json']).is_file() else {}
            meta = dict(source_meta)
            meta.update({
                'Description': 'AABC V1 fMRIPrep 36-parameter nuisance regression with joint 128-second DCT high-pass',
                'Sources': input_paths,
                'SpatialReference': ('native BOLD' if row['space']=='native' else row['space']),
                'RepetitionTime': tr,
                'DenoisingModel': METHOD,
                'GlobalSignalRegression': True,
                'DCTHighPassPeriodSeconds': CUT_SEC,
                'HighPassRegressors': stats['dct_regressors'],
                'RemovedInitialVolumes': stats['removed_initial_volumes'],
                'RemovedInitialVolumeIndices': list(range(stats['removed_initial_volumes'])),
                'OriginalVolumes': original_vol,
                'RetainedVolumes': nt_expected,
                'NonSteadyStateDerivativeNaNFills': stats['first_row_derivative_zero_fills'],
                'MotionScrubbing': False,
                'SpatialSmoothingFWHM': 0.0,
                'SpatialProcessing': 'No additional resampling was performed during denoising',
                'MaskSource': str(mask_path),
                'ConfoundsSource': str(conf_path),
                'ModelRegressors': str(design_path),
                'ModelRegressorsSHA256': design_hash,
                'DesignRank': stats['rank'],
                'ResidualDegreesOfFreedom': stats['effective_residual_degrees_of_freedom'],
                'ProcessingFingerprint': fingerprint,
                'ProcessingUTC': utcnow(),
                'QC': qc,
            })
            # Commit image then JSON. A partial pair is treated as incomplete
            # and never skipped on resume.
            os.replace(pending, bold_path)
            atomic_json(json_path, meta)
        finally:
            if pending.exists():
                pending.unlink()
        del out_data, final_img, o, y
        print(f'COMPLETE: {bold_path}', flush=True)
        print(f'  QC: residual SD {qc["MeanVoxelTemporalSDPost"]:.4g}; '
              f'max nuisance projection {max_relative_projection:.2g}', flush=True)
        return dict(status='COMPLETE', output=str(bold_path),
                    subject=row['subject'], session=row['session'], space=row['space'],
                    masked_voxels=nmask, retained_volumes=nt_expected,
                    dct_regressors=stats['dct_regressors'], design_rank=stats['rank'],
                    max_relative_projection=max_relative_projection)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--experiment', type=Path, default=DEFAULT_EXP)
    ap.add_argument('--manifest', type=Path, help='Override Scripts/reports/preflight_manifest.csv')
    ap.add_argument('--subject', default='sub-HCA8889115')
    ap.add_argument('--session', default='ses-V1')
    ap.add_argument('--all-sessions', action='store_true',
                    help='Process all sessions of --subject (after successful pilot).')
    ap.add_argument('--all', action='store_true', help='Process whole dataset (after validation).')
    ap.add_argument('--space', choices=['both','native',MNI], default='both')
    ap.add_argument('--scratch-dir', type=Path, default=None,
                    help='Filesystem for large temporary files (default: <experiment>/Scripts/tmp).')
    ap.add_argument('--chunk-voxels', type=int, default=512)
    ap.add_argument('--dry-run', action='store_true', help='Validate inputs/design, do not write outputs')
    ap.add_argument('--overwrite', action='store_true', help='Allow replacement of existing output pairs')
    args = ap.parse_args()
    if args.chunk_voxels < 16 or args.chunk_voxels > 10000:
        ap.error('--chunk-voxels must be between 16 and 10000')
    exp = args.experiment.expanduser().resolve()
    manifest = args.manifest or (exp/'Scripts/reports/preflight_manifest.csv')
    scratch = args.scratch_dir or (exp/'Scripts/tmp')
    scratch.mkdir(parents=True, exist_ok=True)
    rows = intended_rows(load_manifest(manifest), args)
    if not rows:
        ap.error('No matching READY rows found in the preflight manifest')
    print(f'Selected {len(rows)} space-specific rows; outputs -> {exp/"Denoising_Output"}')
    print(f'Temporary work directory: {scratch}')
    results = []
    for row in rows:
        try:
            results.append(perform_one(row, args, exp, scratch))
        except Exception as exc:
            print(f'FAILED {row["subject"]} {row["session"]} {row["space"]}: {exc}',
                  file=sys.stderr, flush=True)
            results.append(dict(status='FAILED', subject=row['subject'],
                                session=row['session'], space=row['space'], error=str(exc)))
            # For the initial pilot stop early to make debugging straightforward.
            break
    report = {
        'timestamp_utc': utcnow(), 'method': METHOD,
        'dry_run': args.dry_run, 'requested_rows': len(rows),
        'processed_rows': len(results), 'results': results,
    }
    report_path = exp/'Scripts/reports/denoising_latest_run.json'
    atomic_json(report_path, report)
    print(json.dumps(report, indent=2))
    print('Run report:', report_path)
    if len(results) != len(rows) or any(r['status'] == 'FAILED' for r in results):
        sys.exit(2)


if __name__ == '__main__':
    main()
