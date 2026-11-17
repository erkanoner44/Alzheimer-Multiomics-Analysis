#!/usr/bin/env python3
"""Module 17: real MTG age/sex/plate adjusted *priority CpG* sensitivity analysis.
Does not claim genome-wide adjusted FDR, methylation detection QC, or cell adjustment.
"""
import argparse, csv, gzip, json, math, re, traceback
from pathlib import Path
from collections import Counter
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

DATA_NAMES = {'AD':'GSE134379_processedSamples_mtgAD.txt.gz','ND':'GSE134379_processedSamples_mtgND.txt.gz'}
PROBE_RE = re.compile(r'^Beta value for\s+MTG\s+(.+?)\s*$',re.I)
TITLE_RE = re.compile(r'^(AD|ND),\s*(MTG|CBL) sample\s+(.+)$',re.I)

def load_series(path):
    fields={}; characteristic=[]
    with gzip.open(path,'rt',encoding='utf-8-sig',errors='replace') as f:
        for line in f:
            if not line.startswith('!Sample_'):continue
            rec=next(csv.reader([line],delimiter='\t'))
            if rec[0]=='!Sample_characteristics_ch1':characteristic.append(rec[1:])
            elif rec[0] in ('!Sample_title','!Sample_geo_accession'):fields[rec[0]]=rec[1:]
    if not {'!Sample_title','!Sample_geo_accession'}.issubset(fields):raise ValueError('GEO serisi: sample title / GSM alanlari yok')
    n=len(fields['!Sample_title']); data={'title':fields['!Sample_title'],'gsm':fields['!Sample_geo_accession']}
    for vals in characteristic:
        if len(vals)!=n:raise ValueError('GEO characteristics column count mismatch')
        key=vals[0].split(':',1)[0].strip().lower().replace(' ','_')
        if key in data:raise ValueError('GEO characteristics duplicate: '+key)
        data[key]=[v.split(':',1)[1].strip() if ':' in v else '' for v in vals]
    df=pd.DataFrame(data)
    required=('brain_region','diagnosis','age','gender','plate','bead_chip')
    missing=[x for x in required if x not in df]
    if missing:raise ValueError('GEO missing covariates: '+str(missing))
    df=df[df.brain_region.str.upper()=='MTG'].copy()
    df['subject_id']=df.title.map(lambda s: TITLE_RE.match(s).group(3).strip() if TITLE_RE.match(s) else '')
    if df.subject_id.eq('').any():raise ValueError('Some sample titles could not be matched to subject IDs')
    if df.subject_id.duplicated().any():raise ValueError('Duplicate MTG subject IDs')
    df['age_numeric']=pd.to_numeric(df.age,errors='coerce')
    if df.age_numeric.isna().any() or df.gender.isna().any() or df.plate.isna().any():raise ValueError('Covariate missing')
    if not df.gender.isin(['M','F']).all():raise ValueError('Unexpected gender coding')
    if not df.diagnosis.isin(['AD','ND']).all():raise ValueError('Unexpected diagnosis coding')
    if Counter(df.diagnosis)!=Counter({'AD':225,'ND':179}):raise ValueError('Unexpected MTG group sizes')
    return df

def read_header(path, group, md):
    with gzip.open(path,'rt',encoding='utf-8-sig',errors='replace') as f:
        header=f.readline().rstrip('\r\n').split('\t')
    if len(header)%2:raise ValueError('Odd header width: '+path.name)
    ids=[]
    for i in range(0,len(header),2):
        m=PROBE_RE.match(header[i]);
        if not m or header[i+1].strip().lower()!='detection p value':raise ValueError('Unexpected header pair: '+str(header[i:i+2]))
        ids.append(m.group(1).strip())
    if len(ids)!=len(set(ids)):raise ValueError('Duplicated matrix subject IDs')
    expected=set(md[md.diagnosis==group].subject_id)
    if set(ids)!=expected:raise ValueError(group+' IDs mismatch GEO metadata: '+str(list(set(ids)^expected)[:5]))
    return ids, len(header)+1

def audit_design(md, ids_ad, ids_nd, with_plate=True):
    indexed=md.set_index('subject_id'); ordered=pd.concat([indexed.loc[ids_ad],indexed.loc[ids_nd]])
    X=pd.DataFrame({'AD':(ordered.diagnosis=='AD').astype(float).to_numpy(),
                    'age_per_10_years':((ordered.age_numeric-ordered.age_numeric.mean())/10).astype(float).to_numpy(),
                    'male':(ordered.gender=='M').astype(float).to_numpy()})
    if with_plate:
        dummies=pd.get_dummies(ordered.plate.astype(str),prefix='plate',drop_first=True,dtype=float).reset_index(drop=True)
        X=pd.concat([X,dummies],axis=1)
    X.insert(0,'intercept',1.)
    matrix=X.to_numpy(dtype=float)
    rank=int(np.linalg.matrix_rank(matrix)); cond=float(np.linalg.cond(matrix))
    if rank != matrix.shape[1]:raise ValueError('Design rank-deficient: '+str((rank,matrix.shape[1])))
    if not np.isfinite(cond) or cond>1e8:raise ValueError('Design numerically unstable, condition='+str(cond))
    tab=pd.crosstab(md.plate,md.diagnosis)
    return matrix,list(X.columns),ordered,{'rank':rank,'n_columns':matrix.shape[1],'condition_number':cond,
      'plate_diagnosis_counts':tab.to_dict(),'covariates':list(X.columns),'age_summary':{
      g:{'n':int(len(z)),'mean':float(z.age_numeric.mean()),'sd':float(z.age_numeric.std()),'min':float(z.age_numeric.min()),'max':float(z.age_numeric.max())} for g,z in md.groupby('diagnosis')},
      'sex_by_group':pd.crosstab(md.gender,md.diagnosis).to_dict()}

def test_one(y,x):
    valid=np.isfinite(y)&(y>=0)&(y<=1)
    total=int(valid.sum()); ad=int(valid[:225].sum()); nd=int(valid[225:].sum())
    row={'n_AD':ad,'n_ND':nd,'n_nonmissing':total,'mean_beta_AD':float(np.mean(y[:225][valid[:225]])) if ad else None,
         'mean_beta_ND':float(np.mean(y[225:][valid[225:]])) if nd else None}
    if ad<180 or nd<144:
        row.update({'status':'INSUFFICIENT_COVERAGE','adjusted_AD_beta':np.nan,'adjusted_p':np.nan,'CI95_low':np.nan,'CI95_high':np.nan})
        return row
    xv=x[valid]; yv=y[valid]; rank=np.linalg.matrix_rank(xv)
    if rank != xv.shape[1] or len(yv)<=xv.shape[1]:
        row.update({'status':'SUBSET_DESIGN_RANK_FAIL','adjusted_AD_beta':np.nan,'adjusted_p':np.nan,'CI95_low':np.nan,'CI95_high':np.nan})
        return row
    coef=np.linalg.lstsq(xv,yv,rcond=None)[0]
    resid=yv-xv@coef; dof=len(yv)-xv.shape[1]
    sse=float(resid@resid); mse=sse/dof
    var=float(np.linalg.pinv(xv.T@xv)[1,1]*mse)
    if var<=0 or not np.isfinite(var):
        row.update({'status':'NONPOSITIVE_VARIANCE','adjusted_AD_beta':float(coef[1]),'adjusted_p':np.nan,'CI95_low':np.nan,'CI95_high':np.nan})
        return row
    se=math.sqrt(var); tval=float(coef[1]/se); p=float(2*student_t.sf(abs(tval),dof)); cutoff=float(student_t.ppf(.975,dof))
    row.update({'status':'OK','adjusted_AD_beta':float(coef[1]),'adjusted_p':p,
                'CI95_low':float(coef[1]-cutoff*se),'CI95_high':float(coef[1]+cutoff*se),'df':dof})
    return row

def run(data_dir,series,priority,out):
    out.mkdir(parents=True,exist_ok=True)
    audit={'status':'STARTED','scope':'PRIORITY_CPG_ONLY','methods':'OLS beta ~ diagnosis_AD + age + sex + plate; per-probe complete cases',
           'limitations':['Selection based on prior unadjusted results; subset BH is not genome-wide FDR.',
             'No validated detection p values; no cell composition adjustment.',
             'Postmortem samples; associational model, not causality.']}
    report=out/'module17_audit.json'
    try:
        md=load_series(series)
        md.to_csv(out/'MTG_covariates_verified.csv',index=False)
        pa=data_dir/DATA_NAMES['AD'];pn=data_dir/DATA_NAMES['ND']
        ad_ids,nca=read_header(pa,'AD',md);nd_ids,ncn=read_header(pn,'ND',md)
        x,cols,ordered,design=audit_design(md,ad_ids,nd_ids)
        audit.update({'metadata':{'total':len(md),'AD':225,'ND':179},'design':design})
        pri=pd.read_csv(priority,low_memory=False)
        if 'CpG' not in pri:raise ValueError('Priority CSV requires CpG column')
        wanted={str(v) for v in pri.CpG.dropna().unique() if re.fullmatch(r'cg\d+',str(v),re.I)}
        if not wanted:raise ValueError('No cg CpG IDs in priority CSV')
        records=[]; nrows=0; inspected_non_cg=0
        with gzip.open(pa,'rt',encoding='utf-8-sig',errors='replace') as fa,gzip.open(pn,'rt',encoding='utf-8-sig',errors='replace') as fn:
            next(fa);next(fn)
            for lineno,(la,ln) in enumerate(zip(fa,fn),start=2):
                aa=la.rstrip('\r\n').split('\t');nn=ln.rstrip('\r\n').split('\t')
                if len(aa)!=nca or len(nn)!=ncn:raise ValueError('Wrong row width at line '+str(lineno))
                if aa[0]!=nn[0]:raise ValueError('AD ND CpG order disagreement at line '+str(lineno))
                probe=aa[0];nrows+=1
                if not probe.lower().startswith('cg'):inspected_non_cg+=1
                if probe not in wanted:continue
                values=[float(z) if z not in ('','NA','NaN','nan') else float('nan') for z in aa[1::2]+nn[1::2]]
                y=np.array(values,dtype=float)
                res=test_one(y,x);res['CpG']=probe;records.append(res)
            # zip stops at shorter file, ensure same length
            if fa.readline() or fn.readline():raise ValueError('AD ND unequal matrix length')
        out_df=pd.DataFrame(records)
        from statsmodels.stats.multitest import multipletests
        out_df['BH_within_priority_ONLY']=np.nan
        eligible=out_df.adjusted_p.notna()
        if eligible.any():out_df.loc[eligible,'BH_within_priority_ONLY']=multipletests(out_df.loc[eligible,'adjusted_p'].to_numpy(),method='fdr_bh')[1]
        out_df.to_csv(out/'MTG_priority_CpG_age_sex_plate_adjusted.csv',index=False)
        columns=['CpG']
        extra=[k for k in ['gene_symbol','annotation_region','delta_beta_AD_minus_ND','BH_FDR'] if k in pri.columns]
        detail=pri[columns+extra].merge(out_df,on='CpG',how='left',validate='many_to_one')
        detail.to_csv(out/'MTG_priority_gene_CpG_adjusted_join.csv',index=False)
        audit.update({'status':'COMPLETED_PRIORITY_ONLY','n_all_probe_rows_seen':nrows,'n_non_cg_rows_seen':inspected_non_cg,
                      'n_priority_ids_requested':len(wanted),'n_priority_ids_found':int(len(out_df)),
                      'n_priority_ids_missing':len(wanted)-len(out_df),
                      'n_adjusted_p_calculated':int(eligible.sum()),'n_within_priority_BH_lt_005':int((out_df.BH_within_priority_ONLY<.05).sum())})
    except Exception as e:
        audit['status']='FAILED';audit['error']=str(e);audit['traceback']=traceback.format_exc(limit=5)
        raise
    finally:
        report.write_text(json.dumps(audit,ensure_ascii=False,indent=2),encoding='utf-8')
    print('SUCCESS: priority models computed',len(out_df),'CpGs. Audit:',report)

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--data-dir',type=Path,required=True)
    ap.add_argument('--series',type=Path,required=True)
    ap.add_argument('--priority',type=Path,required=True)
    ap.add_argument('--out',type=Path,default=Path('Module17_results'))
    args=ap.parse_args();run(args.data_dir,args.series,args.priority,args.out)
