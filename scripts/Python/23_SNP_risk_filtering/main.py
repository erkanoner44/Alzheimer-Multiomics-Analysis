#!/usr/bin/env python3
"""Module 23: pre-specified SNP risk/cross-hybridisation exclusion and sensitivity.
No cross-reactive blacklist is invented. BH recalculated across retained tested CpGs.
"""
import argparse, csv, gzip, json, sys
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

PRIORITY=['FKBP5','C1QA','GAD1','SST','VGF','BDNF','C4A','C4B','TNF','GAD2','PVALB','S100A8','S100A9','CRH','SLC32A1']

def read_ids(path):
    if not path: return set()
    p=Path(path)
    if not p.is_file():raise FileNotFoundError(str(p))
    opener=gzip.open if p.suffix=='.gz' else open
    ids=set()
    with opener(p,'rt',encoding='utf-8-sig') as f:
        for line in f:
            token=line.strip().split(',')[0].split('\t')[0].strip().strip('"')
            if token.lower() in ('cpg','name','probe_id','illumina_id','probe') or not token or token.startswith('#'):continue
            if not token.startswith('cg'):raise ValueError(f'{p.name}: unexpected ID {token!r}; input must contain one cg ID per line')
            ids.add(token)
    return ids

def main():
    pa=argparse.ArgumentParser()
    pa.add_argument('--adjusted',required=True,help='Module18 genome-wide all result CSV.GZ')
    pa.add_argument('--annotation',required=True,help='450K CpG-to-gene CSV')
    pa.add_argument('--snp-risk',help='R generated CpG/SBE SNP risk list (one cg ID per line)')
    pa.add_argument('--cross-reactive',help='VALIDATED external list of cross-reactive probes, one cg ID per line')
    pa.add_argument('--out',default='Module23_results')
    pa.add_argument('--min-effect',type=float,default=.02)
    args=pa.parse_args();out=Path(args.out);out.mkdir(parents=True,exist_ok=True)
    audit={'status':'STARTED','started_utc':datetime.now(timezone.utc).isoformat(),'inputs':{'adjusted':args.adjusted,'annotation':args.annotation,'snp_risk':args.snp_risk,'cross_reactive':args.cross_reactive},'notes':[]}
    try:
        df=pd.read_csv(args.adjusted,low_memory=False)
        required={'CpG','adjusted_p','adjusted_AD_beta','BH_FDR_genomewide'}
        if not required.issubset(df.columns):raise ValueError('Missing Module18 fields: '+str(required-set(df.columns)))
        if df['CpG'].duplicated().any(): raise ValueError('Duplicate probe IDs in Module18')
        if not df['CpG'].astype(str).str.match(r'^cg\d{8}$').all():raise ValueError('Non-standard cg probe IDs in Module18')
        tested=df['adjusted_p'].notna() & df['adjusted_p'].between(0,1)
        if tested.sum()<1000:raise ValueError('Genome-wide results appear incomplete')
        expected=multipletests(df.loc[tested,'adjusted_p'].to_numpy(),method='fdr_bh')[1]
        old=df.loc[tested,'BH_FDR_genomewide'].to_numpy(dtype=float)
        max_error=float(np.max(np.abs(expected-old)))
        if max_error>1e-6:raise ValueError(f'Input genomewide FDR mismatch, max difference {max_error:g}')
        snp=read_ids(args.snp_risk);cross=read_ids(args.cross_reactive)
        if not args.snp_risk and not args.cross_reactive:
            audit['status']='NO_FILTER_APPLIED'
            audit['notes'].append('No validated SNP-risk or cross-reactive probe list supplied; sensitivity analysis not run.')
            audit['n_tested_original']=int(tested.sum());audit['n_original_screen_hits']=int((tested & (df['BH_FDR_genomewide']<.05)&(df['adjusted_AD_beta'].abs()>=args.min_effect)).sum())
            (out/'module23_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
            print('NO_FILTER_APPLIED: provide --snp-risk and/or --cross-reactive list.');return
        if args.cross_reactive: audit['notes'].append('External cross-reactive exclusion list provenance MUST be documented by researcher; code cannot certify it.')
        known=set(df.CpG);snp_in=snp&known;cross_in=cross&known;excluded=snp_in|cross_in
        retained=tested&~df.CpG.isin(excluded)
        if retained.sum()<1000:raise ValueError('Too few probes left after filtering; check inputs')
        df['snp_risk_flag']=df.CpG.isin(snp_in)
        df['cross_reactive_flag']=df.CpG.isin(cross_in)
        df['retained_after_filter']=retained
        df['BH_FDR_filtered_genomewide']=np.nan
        df.loc[retained,'BH_FDR_filtered_genomewide']=multipletests(df.loc[retained,'adjusted_p'].to_numpy(),method='fdr_bh')[1]
        df['old_hit']=(tested)&(df.BH_FDR_genomewide<.05)&(df.adjusted_AD_beta.abs()>=args.min_effect)
        df['new_hit']=retained&(df.BH_FDR_filtered_genomewide<.05)&(df.adjusted_AD_beta.abs()>=args.min_effect)
        cols=['CpG','adjusted_AD_beta','adjusted_p','BH_FDR_genomewide','BH_FDR_filtered_genomewide','snp_risk_flag','cross_reactive_flag','retained_after_filter','old_hit','new_hit']
        df[cols].to_csv(out/'MTG_probe_QC_sensitivity_all.csv.gz',index=False,compression='gzip')
        hits=df.loc[df.new_hit,cols].sort_values('BH_FDR_filtered_genomewide')
        hits.to_csv(out/'MTG_filtered_FDR005_absBeta002.csv',index=False)
        ann=pd.read_csv(args.annotation,usecols=['CpG','UCSC_RefGene_Name','UCSC_RefGene_Group'],dtype=str).drop_duplicates('CpG')
        q=df.loc[df.new_hit|df.old_hit,cols]
        q=q.merge(ann,on='CpG',how='left')
        priority=set(PRIORITY)
        pri=q[q.UCSC_RefGene_Name.fillna('').map(lambda v:any(t.strip() in priority for t in v.split(';')))]
        pri.to_csv(out/'priority_genes_filtered_comparison.csv',index=False)
        with open(out/'filtered_tested_CpGs.txt','w') as f:f.write('\n'.join(df.loc[retained,'CpG'])+'\n')
        with open(out/'filtered_candidate_CpGs.txt','w') as f:f.write('\n'.join(df.loc[df.new_hit,'CpG'])+'\n')
        with open(out/'filtered_hypo_CpGs.txt','w') as f:f.write('\n'.join(df.loc[df.new_hit&(df.adjusted_AD_beta<0),'CpG'])+'\n')
        with open(out/'filtered_hyper_CpGs.txt','w') as f:f.write('\n'.join(df.loc[df.new_hit&(df.adjusted_AD_beta>0),'CpG'])+'\n')
        audit.update({'status':'FILTERED_SENSITIVITY_COMPLETED','n_tested_original':int(tested.sum()),'n_excluded_in_tested':int((tested&df.CpG.isin(excluded)).sum()),'n_excluded_SNP_in_tested':int((tested&df.snp_risk_flag).sum()),'n_excluded_cross_reactive_in_tested':int((tested&df.cross_reactive_flag).sum()),'n_excluded_both_in_tested':int((tested&df.snp_risk_flag&df.cross_reactive_flag).sum()),'n_retained':int(retained.sum()),'n_original_screen_hits':int(df.old_hit.sum()),'n_filtered_screen_hits':int(df.new_hit.sum()),'n_original_hits_excluded':int((df.old_hit&~retained).sum()),'n_new_hits_after_BH':int((df.new_hit&~df.old_hit).sum()),'FDR':'BH recomputed over retained successfully tested CpG probes only','p_value_model':'Original Module18 age+sex+plate OLS p-values reused, not refitted','min_effect':args.min_effect,'input_BH_max_diff':max_error})
        audit['notes']+=['Filtered BH adjusts the testing family but does not fix cell mixture, unavailable detection p QC, absent raw-array QC, or model misspecification.','SNP-risk list uses array annotation, not observed cohort genotype.','Gene symbols are annotation links, not gene-level FDR.']
        (out/'module23_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
        print('FILTERED_SENSITIVITY_COMPLETED: retained=',int(retained.sum()),' hits=',int(df.new_hit.sum()))
    except Exception as e:
        audit['status']='ERROR';audit['error']=str(e)
        (out/'module23_audit.json').write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
        raise
if __name__=='__main__': main()
