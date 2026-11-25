#!/usr/bin/env python3
"""Run the original VB Toolbox on MSC 36P-denoised rest data.

This wrapper intentionally does NOT validate the denoised image or VB results.
It selects files from the existing preflight manifest, calls `vb_tool volumetric`,
and records only command execution status and expected-file presence.

Defaults: sub-MSC01 / ses-func01 / native / VB7 / one VB worker.
Use --space both --neighbourhood both for four VB conditions per acquisition.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

DEFAULT_EXP = Path(
    '/mnt/SlowData/sgal0013/Midnight_Scan_Club/'
    'VB_Work_Output/Sephora/20261026_VB_MSC'
)
MNI = 'MNI152NLin2009cAsym'
SOURCE_SUFFIX = '_desc-preproc_bold.nii.gz'
DENOISED_SUFFIX = '_desc-36P128s_bold.nii.gz'


def load_selected_rows(path: Path, args):
    if not path.is_file():
        raise FileNotFoundError('Run 01_preflight_inventory.py first: ' + str(path))
    with path.open(newline='', encoding='utf-8') as f:
        rows = list(csv.DictReader(f))
    if not rows:
        raise RuntimeError('Empty preflight manifest: ' + str(path))
    selected = []
    for row in rows:
        if row.get('task') != 'rest':
            continue
        if not args.all and row.get('subject') != args.subject:
            continue
        if not args.all and not args.all_sessions and row.get('session') != args.session:
            continue
        if args.space != 'both' and row.get('space') != args.space:
            continue
        if row.get('status') != 'READY':
            raise RuntimeError('Non-READY preflight row selected: ' +
                               str((row.get('subject'), row.get('session'), row.get('space'))))
        selected.append(row)
    if not selected:
        raise RuntimeError('No matching rest acquisitions found')
    # Avoid accidentally executing duplicate conditions.
    keys = [(r['subject'], r['session'], r.get('run', ''), r['space']) for r in selected]
    if len(keys) != len(set(keys)):
        raise RuntimeError('Duplicate acquisition/space in preflight manifest')
    return sorted(selected, key=lambda r: (r['subject'], r['session'], r.get('run', ''), 0 if r['space'] == 'native' else 1))


def plan_condition(exp: Path, row: dict, neighbourhood: int):
    source_name = Path(row['bold']).name
    if not source_name.endswith(SOURCE_SUFFIX):
        raise RuntimeError('Unexpected fMRIPrep filename: ' + source_name)
    entities = source_name[:-len(SOURCE_SUFFIX)]
    denoised = (exp / 'Denoising_Output' / row['subject'] /
                row['session'] / 'func' / (entities + DENOISED_SUFFIX))
    mask = Path(row['mask'])
    # Let the original toolbox append .unnorm.vbi-vol.nii.gz and
    # .unnorm.vbi-neigh.nii.gz. No renaming or rewriting of its NIfTIs.
    out_dir = exp / 'VB_Output' / row['subject'] / row['session'] / 'func'
    prefix = out_dir / (entities + f'_desc-vb{neighbourhood}_volumetric_output')
    vol = Path(str(prefix) + '.unnorm.vbi-vol.nii.gz')
    neigh = Path(str(prefix) + '.unnorm.vbi-neigh.nii.gz')
    return denoised, mask, prefix, vol, neigh


def execute_one(exp: Path, row: dict, neighbourhood: int, args, vb_exec: str):
    denoised, mask, prefix, vb, neigh = plan_condition(exp, row, neighbourhood)
    label = f'{row["subject"]} {row["session"]} {row["space"]} VB{neighbourhood}'
    cmd = [vb_exec, 'volumetric', '--data', str(denoised),
           '--volmask', str(mask), '--jobs', str(args.jobs),
           '--output', str(prefix)]
    if neighbourhood == 7:
        cmd.append('--little-neighbourhood')
    print('\n' + label)
    print('COMMAND:', shlex.join(cmd), flush=True)
    print('EXPECTED:', vb, '\n         ', neigh, flush=True)

    if args.dry_run:
        return {'status': 'DRY_RUN', 'label': label, 'command': cmd,
                'vb_map': str(vb), 'neighbourhood_map': str(neigh)}

    # File existence checks are basic execution safeguards, NOT image QC.
    if not denoised.is_file():
        raise FileNotFoundError('Denoised BOLD not found; run 02_denoise_36P.py: ' + str(denoised))
    if not mask.is_file():
        raise FileNotFoundError('Matching fMRIPrep mask not found: ' + str(mask))
    existing = (vb.is_file(), neigh.is_file())
    if all(existing) and not args.overwrite:
        print('SKIP_EXISTING (outputs not quality-checked)', flush=True)
        return {'status': 'SKIP_EXISTING_UNCHECKED', 'label': label,
                'vb_map': str(vb), 'neighbourhood_map': str(neigh)}
    if any(existing) and not args.overwrite:
        raise FileExistsError('Partial VB outputs exist. Review them, then use --overwrite if desired.')
    if args.overwrite:
        # Remove these exact two outputs before rerunning so stale files cannot
        # be mistaken for new outputs. Other processing conditions are untouched.
        vb.unlink(missing_ok=True)
        neigh.unlink(missing_ok=True)

    prefix.parent.mkdir(parents=True, exist_ok=True)
    logdir = exp / 'Scripts' / 'logs'
    logdir.mkdir(parents=True, exist_ok=True)
    logfile = logdir / (prefix.name + '.log')
    env = os.environ.copy()
    env['OPENBLAS_NUM_THREADS'] = '1'
    env['MKL_NUM_THREADS'] = '1'
    env['OMP_NUM_THREADS'] = '1'
    start = time.monotonic()
    with logfile.open('w', encoding='utf-8') as log:
        log.write('COMMAND: ' + shlex.join(cmd) + '\n\n')
        log.flush()
        proc = subprocess.run(cmd, stdout=log, stderr=subprocess.STDOUT, env=env)
    duration = time.monotonic() - start
    if proc.returncode != 0:
        raise RuntimeError(f'vb_tool returned exit code {proc.returncode}; see {logfile}')
    if not (vb.is_file() and neigh.is_file()):
        raise RuntimeError('vb_tool returned success but one or both expected files are absent; see ' + str(logfile))
    print(f'COMPLETE: {label}, {duration:.1f}s', flush=True)
    return {'status': 'COMPLETE', 'label': label, 'seconds': round(duration, 3),
            'vb_map': str(vb), 'neighbourhood_map': str(neigh),
            'log': str(logfile), 'command': cmd}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--experiment', type=Path, default=DEFAULT_EXP)
    p.add_argument('--manifest', type=Path, default=None)
    p.add_argument('--subject', default='sub-MSC01')
    p.add_argument('--session', default='ses-func01')
    p.add_argument('--all-sessions', action='store_true')
    p.add_argument('--all', action='store_true', help='All subjects and rest sessions')
    p.add_argument('--space', choices=('native', MNI, 'both'), default='native')
    p.add_argument('--neighbourhood', choices=('7', '27', 'both'), default='7')
    p.add_argument('--jobs', type=int, default=1)
    p.add_argument('--vb-executable', default='vb_tool')
    p.add_argument('--dry-run', action='store_true')
    p.add_argument('--overwrite', action='store_true',
                   help='Replace only the selected VB condition(s) if output files exist')
    args = p.parse_args()
    if args.jobs < 1:
        p.error('--jobs must be >=1')
    exp = args.experiment.expanduser().resolve()
    manifest = args.manifest or exp / 'Scripts' / 'reports' / 'preflight_manifest.csv'
    try:
        rows = load_selected_rows(manifest, args)
    except Exception as error:
        p.error(str(error))
    neighbourhoods = [7, 27] if args.neighbourhood == 'both' else [int(args.neighbourhood)]
    executable = shutil.which(args.vb_executable)
    if executable is None:
        p.error(f'Cannot find VB executable: {args.vb_executable}')
    report = {'timestamp_utc': datetime.now(timezone.utc).isoformat(),
              'experiment': str(exp), 'vb_executable': executable,
              'vb_qc_performed': False, 'requested_conditions': len(rows) * len(neighbourhoods),
              'results': []}
    error_found = False
    for row in rows:
        for neighbourhood in neighbourhoods:
            try:
                result = execute_one(exp, row, neighbourhood, args, executable)
            except Exception as error:
                error_found = True
                result = {'status': 'FAILED',
                          'label': f'{row["subject"]} {row["session"]} {row["space"]} VB{neighbourhood}',
                          'error': str(error)}
                print('FAILED:', result['label'], error, file=sys.stderr, flush=True)
            report['results'].append(result)
            if error_found:
                break
        if error_found:
            break
    reportdir = exp / 'Scripts' / 'reports'
    reportdir.mkdir(parents=True, exist_ok=True)
    report_path = reportdir / 'vb_simple_latest_run.json'
    report_path.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print('REPORT:', report_path)
    if error_found:
        sys.exit(2)


if __name__ == '__main__':
    main()
