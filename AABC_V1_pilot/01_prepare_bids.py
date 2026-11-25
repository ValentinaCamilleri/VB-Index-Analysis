#!/usr/bin/env python3
"""HCA8889115 V1: inspect and optionally assemble BIDS from original HCP packages.

The source packages are NOT modified. NIfTI files are hard-linked to avoid
copying large data and to ensure links remain usable inside a container bind.
All JSON values are copied from acquisition metadata (never guessed).

Usage: python 01_prepare_bids.py             # dry inventory
       python 01_prepare_bids.py --build     # create only after reviewing plan
"""
from __future__ import annotations
import argparse, json, os, sys
from pathlib import Path
from datetime import datetime, timezone

HOME = Path.home()
SOURCE = HOME/'TARRAGON/AABC/HCA8889115'
BIDS = HOME/'TARRAGON/Pilot/Pilot_Results_fMRIPrep/BIDS/HCA8889115_V1'
SUB = 'sub-HCA8889115'
SES = 'ses-V1'
RUNS = ((1,'AP'),(1,'PA'),(2,'AP'),(2,'PA'))

def unique(root, needle, parent_name=None):
    matches = [p for p in root.rglob('*') if p.is_file() and p.name.endswith(needle)
               and (parent_name is None or p.parent.name == parent_name)]
    if len(matches) != 1:
        raise ValueError(f'Expected exactly one source for {needle} in {root}, found {len(matches)}:\n'+'\n'.join(map(str,matches)))
    return matches[0]

def sidecar(path):
    if not path.name.endswith('.nii.gz'):
        raise ValueError(f'Expected a NIfTI file: {path}')
    candidate = path.with_name(path.name[:-7]+'.json')
    if not candidate.is_file():
        raise FileNotFoundError(f'Corresponding JSON sidecar missing: {candidate}')
    d = json.loads(candidate.read_text(encoding='utf-8-sig'))
    if not isinstance(d,dict): raise ValueError(f'Not JSON dictionary: {candidate}')
    return d

def positive_number(metadata, key, source):
    v=metadata.get(key)
    if not isinstance(v,(int,float)) or isinstance(v,bool) or v <=0:
        raise ValueError(f'{source}: missing or invalid {key} ({v!r})')

def validate_epi(metadata, source, pe):
    positive_number(metadata,'TotalReadoutTime',source)
    positive_number(metadata,'RepetitionTime',source)
    direction=metadata.get('PhaseEncodingDirection')
    if direction not in ('i','i-','j','j-','k','k-'):
        raise ValueError(f'{source}: invalid PhaseEncodingDirection {direction!r}')
    # This polarity-to-BIDS assignment is established by the source inventory.
    expected={'AP':'j-','PA':'j'}[pe]
    if direction != expected:
        raise ValueError(f'{source}: expected {expected!r} for {pe}; got {direction!r}. Inspect orientation before continuing.')

def prepare_plan(func_raw, anat_raw):
    plan=[]
    subj = Path(SUB)/SES
    funcs={}
    for r,pe in RUNS:
        label=f'REST{r}_{pe}'
        stem=f'{SUB}_{SES}_task-rest_dir-{pe}_run-{r}'
        src=unique(func_raw, f'_rfMRI_{label}.nii.gz', parent_name=f'rfMRI_{label}')
        meta=sidecar(src)
        validate_epi(meta,src,pe)
        if not isinstance(meta.get('SliceTiming'),list) or len(meta['SliceTiming'])==0:
            raise ValueError(f'{src}: SliceTiming missing (not inferred)')
        bold_rel=subj/'func'/f'{stem}_bold.nii.gz'
        plan.append((src,bold_rel,meta))
        funcs[label]=bold_rel.relative_to(Path(SUB)).as_posix()
        sbref=unique(func_raw, f'_rfMRI_{label}_SBRef.nii.gz', parent_name=f'rfMRI_{label}')
        sbref_meta=sidecar(sbref)
        validate_epi(sbref_meta,sbref,pe)
        plan.append((sbref,subj/'func'/f'{stem}_sbref.nii.gz',sbref_meta))
    for label,dest in [('T1w_MPR_vNav_4e_e1e2_mean',f'{SUB}_{SES}_T1w.nii.gz'),
                       ('T2w_SPC_vNav',f'{SUB}_{SES}_T2w.nii.gz')]:
        src=unique(anat_raw, f'_{label}.nii.gz', parent_name=label)
        plan.append((src,subj/'anat'/dest,sidecar(src)))
    # Independently acquired spin-echo pairs, per original REST1 / REST2 group.
    # Source package may include duplicates in several folders. Select a copy
    # from REST1_AP / REST2_AP to unambiguously bind to the matching session.
    for r in (1,2):
        group=1 if r==1 else 3
        intended=[funcs[f'REST{r}_{pe}'] for pe in ('AP','PA')]
        for pe in ('AP','PA'):
            src=unique(func_raw, f'_SpinEchoFieldMap{group}_{pe}.nii.gz',parent_name=f'rfMRI_REST{r}_AP')
            meta=sidecar(src)
            validate_epi(meta,src,pe)
            meta=dict(meta)
            meta['IntendedFor']=intended
            dest=subj/'fmap'/f'{SUB}_{SES}_acq-rest{r}_dir-{pe}_epi.nii.gz'
            plan.append((src,dest,meta))
    # Check correspondence of BOLD fieldmap readout times.
    for r in (1,2):
        row=[(path,meta) for path,rel,meta in plan if f'_run-{r}_bold' in rel.name or f'_acq-rest{r}_' in rel.name]
        times={round(float(meta['TotalReadoutTime']),8) for path,meta in row}
        if len(times)!=1:
            raise ValueError(f'REST{r}: fieldmap/BOLD TotalReadoutTime differs: {row}')
    return plan

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-root',type=Path,default=SOURCE)
    p.add_argument('--bids-root',type=Path,default=BIDS)
    p.add_argument('--build',action='store_true',help='Create BIDS using hardlinks. Default only inventories and prints plan.')
    args=p.parse_args()
    root=args.source_root.expanduser().resolve()
    dest=args.bids_root.expanduser().resolve()
    func=root/'Resting State fMRI Unprocessed'
    anat=root/'Structural Unprocessed'
    if not func.is_dir() or not anat.is_dir():
        p.error('Missing raw folder(s). Check '+str(func)+' AND '+str(anat))
    try:
        plan=prepare_plan(func,anat)
    except Exception as exc:
        p.error(str(exc))
    print(f'SOURCE {root}\nBIDS TARGET {dest}\nPLAN: {len(plan)} NIfTI + metadata pairs')
    for src,rel,meta in plan:
        print(f'  {src}\n      -> {dest/rel}')
    if not args.build:
        print('\nDRY RUN ONLY. If paths and metadata look correct, rerun with --build.')
        return 0
    dest.mkdir(parents=True,exist_ok=True)
    # Fail early if source and BIDS root are different file systems.
    if os.stat(func).st_dev != os.stat(dest).st_dev or os.stat(anat).st_dev != os.stat(dest).st_dev:
        p.error('Source and BIDS target not on same filesystem; hardlinking would fail. No data copied.')
    desc={'Name':'HCA8889115 V1 locally assembled from AABC/HCP unprocessed packages',
          'BIDSVersion':'1.8.0','DatasetType':'raw'}
    dd=dest/'dataset_description.json'
    if dd.exists() and json.loads(dd.read_text()) != desc:
        p.error(f'Existing dataset_description differs: {dd}')
    else:
        dd.write_text(json.dumps(desc,indent=2)+'\n',encoding='utf-8')
    counts={'linked':0,'existing':0,'metadata':0}
    for src,rel,meta in plan:
        target=dest/rel
        target.parent.mkdir(parents=True,exist_ok=True)
        if target.exists():
            if not os.path.samefile(src,target):
                raise FileExistsError(f'Existing target not hardlinked to expected source: {target}')
            counts['existing']+=1
        else:
            os.link(src,target)
            counts['linked']+=1
        json_target=target.with_name(target.name[:-7]+'.json')
        content=json.dumps(meta,indent=2,ensure_ascii=False)+'\n'
        if json_target.exists():
            if json_target.read_text(encoding='utf-8') != content:
                raise FileExistsError(f'Existing metadata differs: {json_target}')
        else:
            json_target.write_text(content,encoding='utf-8')
            counts['metadata']+=1
    print('BIDS CREATED (zero MRI payload copies):',counts)
    print('Now check: find',str(dest),'-type f -o -type l')
    print('Run BIDS Validator if installed; otherwise fMRIPrep performs validation.')
    return 0
if __name__=='__main__': sys.exit(main())
