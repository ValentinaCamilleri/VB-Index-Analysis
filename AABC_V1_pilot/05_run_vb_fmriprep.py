#!/usr/bin/env python3
"""Run official unnormalised VB7/VB27 on MSC 36P-denoised AABC V1 fMRIPrep BOLD.

No additional denoising, smoothing, filtering, trimming or registration here.
Requires successful 03_preflight_derivatives.py and 04_denoise_36P.py.
"""
from __future__ import annotations
import argparse,csv,json,os,shlex,shutil,subprocess,sys,time
from datetime import datetime, timezone
from pathlib import Path
EXP=Path.home()/'TARRAGON/Pilot/Pilot_Results_fMRIPrep/HCA8889115/V1/fMRIPrep_VB'
RUN_IDS=['REST1_AP','REST1_PA','REST2_AP','REST2_PA']
MNI='MNI152NLin2009cAsym'

def elapsed(sec):
    h, rem=divmod(int(sec),3600);m,s=divmod(rem,60)
    return f'{h:02}:{m:02}:{s:02}'

def progress(done,total,detail=''):
    width=24; n=round(width*done/total) if total else 0
    print(f'OVERALL [{"#"*n}{"-"*(width-n)}] {done}/{total} ({100*done//total if total else 100}%) {detail}',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    group=p.add_mutually_exclusive_group()
    group.add_argument('--run',choices=RUN_IDS,default=None)
    group.add_argument('--all-runs',action='store_true')
    p.add_argument('--space',choices=['native',MNI,'both'],default=MNI)
    p.add_argument('--neighbourhood',choices=['7','27','both'],default='7')
    p.add_argument('--jobs',type=int,default=1)
    p.add_argument('--heartbeat-seconds',type=int,default=10)
    p.add_argument('--vb-executable',default=None)
    p.add_argument('--dry-run',action='store_true')
    p.add_argument('--overwrite',action='store_true')
    args=p.parse_args()
    if args.jobs<1 or args.heartbeat_seconds<1:p.error('jobs and heartbeat-seconds must be positive')
    manifest=EXP/'Scripts/reports/preflight_manifest.csv'
    if not manifest.exists():p.error('Missing manifest: run 03_preflight_derivatives.py first')
    with manifest.open(newline='',encoding='utf-8') as fh:rows=list(csv.DictReader(fh))
    selected_runs=RUN_IDS if args.all_runs else [args.run or 'REST1_PA']
    selected_spaces=['native',MNI] if args.space=='both' else [args.space]
    sizes=[7,27] if args.neighbourhood=='both' else [int(args.neighbourhood)]
    keys={}
    for row in rows:
        key=(row.get('run_id'),row.get('space'))
        if key in keys:p.error(f'Duplicate preflight row for {key}')
        keys[key]=row
    conditions=[]
    for run in selected_runs:
        for space in selected_spaces:
            row=keys.get((run,space))
            if not row or row['status']!='READY':p.error(f'Missing READY preflight row for {run} / {space}')
            for size in sizes:conditions.append((row,size))
    exe=args.vb_executable or shutil.which('vb_tool') or str(Path.home()/'.local/bin/vb_tool')
    if not Path(exe).is_file():p.error('Cannot locate vb_tool. Pass --vb-executable /path/to/vb_tool')
    results=[]; began=time.monotonic(); failed=False
    print(f'RUNNING fMRIPrep–36P–VB: {len(conditions)} condition(s)',flush=True)
    progress(0,len(conditions))
    for index,(row,size) in enumerate(conditions,1):
        run,space=row['run_id'],row['space']
        source=Path(row['bold'])
        denoised=EXP/'Denoising_Output'/row['subject']/row['session']/'func'/source.name.replace('_desc-preproc_bold.nii.gz','_desc-36P128s_bold.nii.gz')
        json_sidecar=denoised.with_name(denoised.name.replace('.nii.gz','.json'))
        outdir=EXP/'VB_Output'/run/space
        prefix=outdir/(source.name.replace('_desc-preproc_bold.nii.gz',f'_desc-vb{size}_volumetric_output'))
        vb=Path(str(prefix)+'.unnorm.vbi-vol.nii.gz')
        neigh=Path(str(prefix)+'.unnorm.vbi-neigh.nii.gz')
        mask=Path(row['mask'])
        cmd=[str(exe),'volumetric','--data',str(denoised),'--volmask',str(mask),'--jobs',str(args.jobs),'--output',str(prefix)]
        if size==7:cmd.append('--little-neighbourhood')
        label=f'{run} {space} VB{size}'
        print(f'\n[{index}/{len(conditions)}] {label}',flush=True)
        print('COMMAND:',shlex.join(cmd),flush=True)
        try:
            if not args.dry_run:
                if not denoised.is_file() or not json_sidecar.is_file():
                    raise FileNotFoundError(f'Missing 36P BOLD/JSON; run 04_denoise_36P.py: {denoised}')
                meta=json.loads(json_sidecar.read_text(encoding='utf-8'))
                if meta.get('DenoisingModel') != 'aabc-fmriprep36p-dct128-nss-trim-v1' or not meta.get('QC',{}).get('passed'):
                    raise ValueError('36P metadata does not verify the expected model and passed QC')
                if str(source) not in meta.get('Sources',[]):
                    raise ValueError('36P metadata does not identify the manifest BOLD input')
                if not mask.is_file():raise FileNotFoundError(f'Missing original fMRIPrep mask: {mask}')
                old=[vb.exists(),neigh.exists()]
                if all(old) and not args.overwrite:
                    print('SKIPPED: both output files already exist (not full QC)',flush=True)
                    results.append({'status':'SKIP_EXISTING_UNCHECKED','condition':label})
                    progress(index,len(conditions),label)
                    continue
                if any(old) and not args.overwrite:raise FileExistsError('Partial output pair already exists')
                if args.overwrite:
                    vb.unlink(missing_ok=True);neigh.unlink(missing_ok=True)
                outdir.mkdir(parents=True,exist_ok=True)
                logs=EXP/'Scripts/logs';logs.mkdir(parents=True,exist_ok=True)
                log=logs/(prefix.name+'.log')
                env=os.environ.copy()
                for k in ('OPENBLAS_NUM_THREADS','MKL_NUM_THREADS','OMP_NUM_THREADS'):env[k]='1'
                start=time.monotonic()
                with log.open('w',encoding='utf-8') as fh:
                    fh.write('COMMAND: '+shlex.join(cmd)+'\n');fh.flush()
                    with subprocess.Popen(cmd,stdout=fh,stderr=subprocess.STDOUT,env=env) as proc:
                        try:
                            while True:
                                try:ret=proc.wait(timeout=args.heartbeat_seconds);break
                                except subprocess.TimeoutExpired:
                                    print(f'STAGE: RUNNING {label} | elapsed {elapsed(time.monotonic()-start)} | condition {index}/{len(conditions)}',flush=True)
                        except KeyboardInterrupt:
                            proc.terminate();proc.wait(timeout=20);raise
                if ret!=0:raise RuntimeError(f'VB exited {ret}; see {log}')
                if not vb.is_file() or not neigh.is_file():raise RuntimeError('One or both expected VB output files missing')
                status='COMPLETE'; seconds=round(time.monotonic()-start,3)
                print(f'COMPLETE in {seconds:.1f}s: {vb}',flush=True)
                results.append({'condition':label,'status':status,'duration_seconds':seconds,'vb':str(vb),'neighbours':str(neigh),'log':str(log)})
            else:
                print('DRY RUN: no processing',flush=True)
                results.append({'condition':label,'status':'DRY_RUN','command':cmd})
            progress(index,len(conditions),label)
        except Exception as exc:
            failed=True;results.append({'condition':label,'status':'FAILED','error':str(exc)})
            print(f'FAILED {label}: {exc}',file=sys.stderr,flush=True);break
    reportdir=EXP/'Scripts/reports';reportdir.mkdir(parents=True,exist_ok=True)
    report={'timestamp_utc':datetime.now(timezone.utc).isoformat(),'requested_conditions':len(conditions),
            'elapsed_seconds':round(time.monotonic()-began,3),'results':results}
    rpath=reportdir/'vb_latest_run.json';rpath.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('REPORT:',rpath)
    return 2 if failed else 0
if __name__=='__main__':sys.exit(main())
