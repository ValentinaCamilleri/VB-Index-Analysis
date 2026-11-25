#!/usr/bin/env python3
"""Inventory fMRIPrep 22.1.1 BOLD/mask/confounds for original MSC 36P step.

Read-only on fMRIPrep derivatives. Writes a READY/CHECK manifest in
Pilot_Results_fMRIPrep/HCA8889115/V1/fMRIPrep_VB/Scripts/reports.
"""
from __future__ import annotations
import csv, json, math, re, sys
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import nibabel as nib

BASE=Path.home()/'TARRAGON/Pilot/Pilot_Results_fMRIPrep/HCA8889115/V1'
FMRIPREP=BASE/'fMRIPrep_output'
EXP=BASE/'fMRIPrep_VB'
MNI='MNI152NLin2009cAsym'
RUNS={(r,pe):f'REST{r}_{pe}' for r in (1,2) for pe in ('AP','PA')}
BASE_COLS=('trans_x','trans_y','trans_z','rot_x','rot_y','rot_z','white_matter','csf','global_signal')
COLS36=tuple(f'{a}{b}' for a in BASE_COLS for b in ('','_derivative1','_power2','_derivative1_power2'))
FIELDS=['subject','session','task','run','run_id','dir','space','status','issues','notes','bold','mask','confounds_tsv','confounds_json','bold_json','n_volumes','tr_seconds','voxel_sizes_mm','n_mask_voxels','confounds_rows','n_nonsteady_state','n_retained','n_cosine_regressors','n_planned_dct_regressors','n_first_row_derivative_nan_fills','n_36p_regressors']

def get_info(bold,mask,confounds,row):
    issues=[]; notes=[]
    if not bold.is_file(): issues.append('BOLD_MISSING')
    if not mask.is_file(): issues.append('MASK_MISSING')
    if not confounds.is_file(): issues.append('CONFOUNDS_MISSING')
    if issues:return issues,notes
    b=nib.load(str(bold)); m=nib.load(str(mask))
    if len(b.shape)!=4 or b.shape[3]<50: issues.append('BOLD_INVALID_4D')
    if m.shape!=b.shape[:3]: issues.append('MASK_SHAPE_MISMATCH')
    if not np.allclose(m.affine,b.affine,rtol=1e-5,atol=1e-3): issues.append('MASK_AFFINE_MISMATCH')
    try:
        v=np.asarray(m.dataobj)
        if not np.isfinite(v).all() or not np.all(np.isin(np.unique(v),[0,1])): issues.append('MASK_NOT_BINARY')
        row['n_mask_voxels']=int(np.count_nonzero(v))
        if row['n_mask_voxels']<100:issues.append('MASK_TOO_SMALL')
    except Exception as exc:issues.append('MASK_READ_ERROR:'+str(exc))
    row['n_volumes']=int(b.shape[-1]); row['tr_seconds']=float(b.header.get_zooms()[3])
    row['voxel_sizes_mm']=';'.join(map(str,b.header.get_zooms()[:3]))
    if row['tr_seconds']<=0 or not math.isfinite(row['tr_seconds']): issues.append('INVALID_TR')
    with confounds.open(encoding='utf-8',newline='') as fh:
        reader=csv.DictReader(fh,delimiter='\t'); columns=reader.fieldnames or []; lines=list(reader)
    row['confounds_rows']=len(lines)
    if row['n_volumes']!=len(lines):issues.append('CONFOUNDS_ROWS_MISMATCH')
    missing=sorted(set(COLS36)-set(columns))
    row['n_36p_regressors']=36-len(missing)
    if missing:issues.append('MISSING_36P_COLUMNS:'+','.join(missing))
    nss_cols=[c for c in columns if re.fullmatch(r'non_steady_state_outlier\d+',c)]
    flagged=set()
    for col in nss_cols:
        hits=[]
        for idx,data in enumerate(lines):
            try:v=float(data[col])
            except (ValueError,TypeError):v=float('nan')
            if not np.isfinite(v):issues.append('INVALID_NSS_COLUMN:'+col);break
            if v>0.5:hits.append(idx)
        if len(hits)!=1:issues.append('NSS_NOT_ONE_HOT:'+col)
        flagged.update(hits)
    nskip=len(flagged)
    if flagged!=set(range(nskip)):issues.append('NSS_NOT_INITIAL_CONTIGUOUS')
    row['n_nonsteady_state']=nskip
    row['n_retained']=len(lines)-nskip
    if row['n_retained']<50: issues.append('INSUFFICIENT_RETAINED_VOLUMES')
    row['n_cosine_regressors']=sum(bool(re.fullmatch(r'cosine_\d+',c)) for c in columns)
    row['n_planned_dct_regressors']=int(math.floor(2*row['n_retained']*row['tr_seconds']/128))
    if row['n_cosine_regressors']: notes.append('Existing fMRIPrep cosines not reused; regenerate DCT after NSS removal')
    fills=0
    for col in COLS36:
        if col not in columns:continue
        bad=[]
        for idx in range(nskip,len(lines)):
            try:x=float(lines[idx][col])
            except (ValueError,TypeError):x=float('nan')
            if not math.isfinite(x):bad.append(idx)
        if bad:
            if nskip==0 and bad==[0] and col.endswith(('_derivative1','_derivative1_power2')):
                fills+=1
            else:issues.append('INVALID_36P_VALUES:'+col)
    row['n_first_row_derivative_nan_fills']=fills
    if fills:notes.append(f'EXPECTED_FIRST_ROW_DERIVATIVE_NANS:{fills}')
    return issues,notes

def has_entities(name, run, pe):
    parts = name.split('_')
    return (f'run-{run}' in parts and f'dir-{pe}' in parts and 'task-rest' in parts)


def confound_for(folder, run, pe):
    matches = [p for p in folder.glob('*_desc-confounds_timeseries.tsv')
               if has_entities(p.name, run, pe)]
    if not matches:
        matches = [p for p in folder.glob('*_desc-confounds_regressors.tsv')
                   if has_entities(p.name, run, pe)]
    if len(matches)!=1:
        raise ValueError(f'Expected one confounds TSV for REST{run}_{pe}; found {len(matches)}: {matches}')
    return matches[0]


def main():
    rows=[]
    # Identify images by their BIDS entities rather than assuming a fixed ordering.
    for r in (1,2):
        for pe in ('AP','PA'):
            stem=f'sub-HCA8889115_ses-V1_task-rest_dir-{pe}_run-{r}'
            folder=FMRIPREP/'sub-HCA8889115/ses-V1/func'
            if not folder.is_dir():
                print('MISSING fMRIPrep output directory:',folder,file=sys.stderr)
                return 2
            try:
                conf=confound_for(folder,r,pe)
            except ValueError as exc:
                print('CONFOUNDS ERROR:',exc,file=sys.stderr)
                return 2
            for space in ('native',MNI):
                available=[]
                for path in folder.glob('*_desc-preproc_bold.nii.gz'):
                    if not has_entities(path.name,r,pe):
                        continue
                    fname=path.name
                    match=re.search(r'_space-([^_]+)',fname)
                    found=match.group(1) if match else 'native'
                    if space=='native' and found in ('native','func','bold','run'):
                        available.append(path)
                    elif space==MNI and found==MNI:
                        available.append(path)
                row={'subject':'sub-HCA8889115','session':'ses-V1',
                     'task':'rest','run':str(r),'run_id':RUNS[(r,pe)],'dir':pe,
                     'space':space,'status':'CHECK','issues':'','notes':'',
                     'confounds_tsv':str(conf),'confounds_json':str(conf.with_suffix('.json')) if conf.with_suffix('.json').exists() else ''}
                issues=[];notes=[]
                if len(available)!=1:
                    issues.append(f'BOLD_MATCH_COUNT:{len(available)}')
                else:
                    bold=available[0]
                    mask=bold.with_name(bold.name.replace('_desc-preproc_bold.nii.gz','_desc-brain_mask.nii.gz'))
                    row.update(bold=str(bold),mask=str(mask),
                               bold_json=str(bold.with_name(bold.name.replace('.nii.gz','.json'))))
                    try:
                        a,b=get_info(bold,mask,conf,row); issues.extend(a); notes.extend(b)
                    except Exception as exc:
                        issues.append(f'VALIDATION_ERROR:{type(exc).__name__}:{exc}')
                row['issues']='|'.join(issues);row['notes']='|'.join(notes)
                row['status']='READY' if not issues else 'CHECK'
                rows.append(row)
                print(f"{row['run_id']} {space}: {row['status']}"+(f" | {row['issues']}" if issues else ''))
    reportdir=EXP/'Scripts/reports';reportdir.mkdir(parents=True,exist_ok=True)
    manifest=reportdir/'preflight_manifest.csv'
    with manifest.open('w',newline='',encoding='utf-8') as fh:
        wr=csv.DictWriter(fh,fieldnames=FIELDS);wr.writeheader();wr.writerows(rows)
    summary={'timestamp_utc':datetime.now(timezone.utc).isoformat(),
             'derivatives':str(FMRIPREP),'manifest':str(manifest),
             'expected_conditions':8,'ready_conditions':sum(x['status']=='READY' for x in rows)}
    (reportdir/'preflight_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(json.dumps(summary,indent=2))
    return 0 if summary['ready_conditions']==8 else 2
if __name__=='__main__':sys.exit(main())
