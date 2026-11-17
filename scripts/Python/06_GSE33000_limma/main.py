"""GSE33000 extraction and limma wrapper. Requires R + limma for differential expression."""
import argparse
import csv
import gzip
import json
import math
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

BASE = Path(__file__).resolve().parent

def locate(arg, names):
    if arg:
        x=Path(arg).expanduser().resolve()
        if not x.is_file(): raise FileNotFoundError(x)
        return x
    roots=[Path.cwd(), BASE, BASE.parent, Path.cwd()/"GSE33000_output", BASE.parent/"GSE33000_output"]
    for root in roots:
        for name in names:
            x=root/name
            if x.is_file(): return x.resolve()
    raise FileNotFoundError(f"Bulunamadı: {names}. --soft ve --metadata argümanları ile yolları verin.")

def metadata_load(path):
    df=pd.read_csv(path,dtype=str).fillna("")
    required={'gsm','group','age_years','sex'}
    if not required.issubset(df.columns): raise ValueError(f"Metadata eksik sütunlar: {required-set(df.columns)}")
    df=df.loc[df.group.isin(['AD','CN'])].copy()
    if df.gsm.duplicated().any(): raise ValueError("Tekrarlanan GSM örnekleri var")
    df['age']=pd.to_numeric(df.age_years.str.extract(r'(\d+(?:\.\d+)?)')[0],errors='coerce')
    df['sex']=df.sex.str.strip().str.lower()
    if df.age.isna().any() or not df.sex.isin(['male','female']).all():
        raise ValueError('Yaş/cinsiyet değeri eksik veya geçersiz; otomatik doldurulmadı.')
    groups=df.group.value_counts().to_dict()
    if not groups.get('AD') or not groups.get('CN'): raise ValueError('AD ve CN her ikisi de gerekli')
    return df

def read_sample_header(handle):
    gsm=None
    for line in handle:
        if line.startswith('^SAMPLE = '):
            gsm=line.strip().split('=',1)[1].strip()
        elif line.strip().lower()=='!sample_table_begin':
            col=next(handle).rstrip('\r\n').split('\t')
            if col[:2]!=['ID_REF','VALUE']:
                raise ValueError(f'{gsm} beklenmeyen sütunlar: {col}')
            return gsm
    return None

def extract(soft,meta,output,max_samples=None,max_probes=None):
    wanted=set(meta.gsm)
    if max_samples is not None:
        wanted=set(meta.gsm.iloc[:max_samples])
    arrays={}; first_ids=None
    opener=gzip.open if str(soft).lower().endswith('.gz') else open
    with opener(soft,'rt',encoding='utf-8',errors='replace',newline='') as f:
        while True:
            gsm=read_sample_header(f)
            if gsm is None: break
            if gsm in wanted:
                ids=[]; vals=[]
                for line in f:
                    if line.strip().lower()=='!sample_table_end': break
                    if not line.strip(): continue
                    fields=line.rstrip('\r\n').split('\t')
                    if len(fields)<2: raise ValueError(f'{gsm} eksik ölçüm satırı')
                    ids.append(fields[0]);
                    try: vals.append(float(fields[1]) * math.log2(10))
                    except ValueError: vals.append(np.nan)
                if max_probes: ids=ids[:max_probes]; vals=vals[:max_probes]
                if first_ids is None:
                    first_ids=ids
                    if len(first_ids)!=len(set(first_ids)): raise ValueError('Tekrarlayan prob ID')
                elif ids!=first_ids: raise ValueError(f'{gsm}: prob sırası/kimliği uyuşmuyor')
                arrays[gsm]=np.asarray(vals,dtype=np.float32)
                print(f'{len(arrays):3d}/{len(wanted)} {gsm}, {len(vals)} prob',flush=True)
            else:
                for line in f:
                    if line.strip().lower()=='!sample_table_end': break
            if len(arrays)==len(wanted): break
    missing=wanted-set(arrays)
    if missing: raise ValueError(f'SOFT dosyasında GSM eksik: {sorted(missing)[:8]}')
    selected=meta[meta.gsm.isin(wanted)].copy()
    selected.to_csv(output/'sample_metadata_analysis.csv',index=False)
    ordered=selected.gsm.tolist()
    # Matrix order: probes x biological samples; log2 sample/reference ratio.
    df=pd.DataFrame(np.stack([arrays[g] for g in ordered],axis=1),index=pd.Index(first_ids,name='ID_REF'),columns=ordered)
    if not np.isfinite(df.to_numpy()).all():
        nbad=int((~np.isfinite(df.to_numpy())).sum())
        print(f'UYARI: {nbad} geçersiz değer var; limma kodunda eksik değer denetlenecek')
    df.to_csv(output/'expression_log2ratio_probes.csv.gz',compression={'method':'gzip','compresslevel':1},float_format='%.7g')
    summary={'input_soft':str(soft),'n_samples':len(ordered),'n_probes':len(df),
             'AD':int((selected.group=='AD').sum()),'CN':int((selected.group=='CN').sum()),
             'scale':'log2(Cy5/Cy3), converted from original normalized log10 ratio',
             'n_nonfinite':int((~np.isfinite(df.to_numpy())).sum()),
             'extraction_is_subset':bool(max_samples or max_probes),
             'r_limma_done':False}
    (output/'preparation_audit.json').write_text(json.dumps(summary,indent=2,ensure_ascii=False),encoding='utf-8')
    return summary

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--soft',default=None,help='GSE33000_family.soft.gz yolu')
    p.add_argument('--metadata',default=None,help='sample_metadata_AD_CN_corrected.csv yolu')
    p.add_argument('--out',default='GSE33000_Modul06_output')
    p.add_argument('--max-samples',type=int,default=None,help='SADECE smoke test')
    p.add_argument('--max-probes',type=int,default=None,help='SADECE smoke test')
    p.add_argument('--prepare-only',action='store_true')
    args=p.parse_args()
    soft=locate(args.soft,['GSE33000_family.soft.gz','GSE33000_family.soft(1).gz'])
    meta_path=locate(args.metadata,['sample_metadata_AD_CN_corrected.csv','sample_metadata_AD_CN_corrected(1).csv'])
    out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True)
    meta=metadata_load(meta_path)
    print('Metadata sayıları:',meta.group.value_counts().to_dict())
    if not args.max_samples and len(meta)!=467: raise ValueError(f'Tam analiz için 467 bekleniyor, {len(meta)} bulundu')
    audit=extract(soft,meta,out,args.max_samples,args.max_probes)
    print('Hazırlama sonucu:',audit)
    if args.max_samples or args.max_probes or args.prepare_only:
        print('Hazırlama tamamlandı. DEG ÇALIŞTIRILMADI.');return
    r=shutil.which('Rscript')
    if not r:
        print('Rscript bulunamadı. R kurulumu sonrası çalıştırın: Rscript limma_analysis.R '+str(out))
        return
    subprocess.run([r,str(BASE/'limma_analysis.R'),str(out)],check=True)
    audit['r_limma_done']=True
    (out/'preparation_audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')

if __name__=='__main__':main()
