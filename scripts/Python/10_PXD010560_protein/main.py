"""Modul 10: PXD010560 supporting protein tables discovery and cautious mRNA-protein matching."""
from __future__ import annotations
import argparse, json, re, sys
from pathlib import Path
import numpy as np
import pandas as pd
import requests

HERE = Path(__file__).resolve().parent
SUPP_URLS = {
    'S1':'https://pmc.ncbi.nlm.nih.gov/articles/instance/7775912/bin/NIHMS1620738-supplement-supp_fileS1.xlsx',
    'S2':'https://pmc.ncbi.nlm.nih.gov/articles/instance/7775912/bin/NIHMS1620738-supplement-supp_fileS2.xlsx',
    'S3':'https://pmc.ncbi.nlm.nih.gov/articles/instance/7775912/bin/NIHMS1620738-supplement-supp_fileS3.xlsx',
    'S4':'https://pmc.ncbi.nlm.nih.gov/articles/instance/7775912/bin/NIHMS1620738-supplement-supp_fileS4.xlsx',
}
PRIDE_API = 'https://www.ebi.ac.uk/pride/ws/archive/v2/projects/PXD010560/files'

def fetch(url, dest):
    if dest.exists() and dest.stat().st_size>1000:
        return 'existing'
    tmp=dest.with_suffix(dest.suffix+'.part')
    try:
        with requests.get(url, timeout=(25,90), stream=True, headers={'User-Agent':'Academic-research-proteomics/1.0'}) as r:
            r.raise_for_status()
            with tmp.open('wb') as f:
                for c in r.iter_content(chunk_size=1024*256):
                    if c:f.write(c)
        if dest.suffix=='.xlsx':
            with tmp.open('rb') as f:
                if f.read(2)!=b'PK':raise ValueError('Not an XLSX / server returned HTML')
        tmp.replace(dest)
        return 'downloaded'
    except Exception as e:
        tmp.unlink(missing_ok=True)
        return 'FAILED: '+str(e)[:300]

def norm(v):
    s=str(v).strip().upper()
    return re.sub(r'\s+', '', s)

def tokenize(s):
    return [x for x in re.split(r'[;,|/\s]+',str(s).upper()) if x]

def detect_symbol_columns(df):
    cols=[]
    for c in df.columns:
        lc=str(c).strip().lower()
        if re.search(r'(gene.?symbol|gene.?name|gene.?id|^symbol$|^gene$|gene.?names)',lc):
            if not any(x in lc for x in ['ontology','description','annotation score','gene set']):cols.append(c)
    return cols

def inspect_table(df,name,sheet,strong):
    info={'file':name,'sheet':str(sheet),'rows':int(len(df)),'columns':int(len(df.columns)),'fields':' | '.join(map(str,df.columns))[:7000]}
    symbol_cols=detect_symbol_columns(df)
    info['symbol_columns']=' | '.join(map(str,symbol_cols))
    output=[]
    genes=set(str(x).upper() for x in strong['gene_symbol'].dropna())
    for col in symbol_cols:
        values=df[col].astype(str)
        for idx,v in values.items():
            hits=genes.intersection(tokenize(v))
            for hit in hits:
                output.append({'file':name,'sheet':str(sheet),'row_number_excel_approx':int(idx)+2,'matched_gene':hit,
                               'source_gene_column':str(col),'source_gene_value':v[:300],
                               'evidence_class':'TABLE_GENE_NAME_MATCH_ONLY',
                               'interpretation':'Detection/match only; NOT AD-vs-control differential protein evidence'})
    return info,output

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--download',action='store_true',help='Download supplementary spreadsheets and PRIDE manifest')
    ap.add_argument('--input',default=str(HERE/'input'))
    ap.add_argument('--out',default=str(HERE/'outputs'))
    args=ap.parse_args()
    inp=Path(args.input).resolve();out=Path(args.out).resolve();out.mkdir(parents=True,exist_ok=True);inp.mkdir(parents=True,exist_ok=True)
    strong_file=inp/'strong_candidate_genes.csv'
    if not strong_file.exists():raise FileNotFoundError(strong_file)
    strong=pd.read_csv(strong_file)
    if 'gene_symbol' not in strong:raise ValueError('strong candidate file lacks gene_symbol')
    audit={'study':'PXD010560','phase':'processed-table discovery and evidence audit',
           'proteomics_statistical_validation':'NOT_ESTABLISHED',
           'source_article':'https://doi.org/10.1111/jnc.15116',
           'critical_warning':'Transcriptomic GSE33000 cohort is independent from proteomic Sathe et al. cohort. Protein matching alone is NOT validation of differential abundance.'}
    if args.download:
        audit['supplement_downloads']={key:fetch(url,inp/('PXD010560_'+key+'.xlsx')) for key,url in SUPP_URLS.items()}
        try:
            response=requests.get(PRIDE_API,timeout=40);response.raise_for_status();data=response.json()
            (out/'PRIDE_files_manifest.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
            audit['pride_manifest']='retrieved'
            objs=data if isinstance(data,list) else data.get('_embedded',{}).get('files',[]) if isinstance(data,dict) else []
            rows=[]
            for x in objs:
                if isinstance(x,dict):rows.append({'fileName':x.get('fileName',''),'fileCategory':str(x.get('fileCategory','')),'fileSize':x.get('fileSizeBytes',x.get('fileSize','')),'downloadLink':str(x.get('publicFileLocations',x.get('downloadLink','')))[:1000]})
            pd.DataFrame(rows).to_csv(out/'PRIDE_files_manifest.csv',index=False)
        except Exception as e:audit['pride_manifest']='FAILED: '+str(e)[:250]
    sheets=[];matches=[];problems=[]
    for p in sorted(list(inp.glob('*.xlsx')) + list(inp.glob('*.xls'))):
        try:
            wb=pd.ExcelFile(p)
            for sheet in wb.sheet_names:
                try:
                    df=pd.read_excel(wb,sheet_name=sheet)
                    info,found=inspect_table(df,p.name,sheet,strong);sheets.append(info);matches.extend(found)
                except Exception as e:problems.append({'file':p.name,'sheet':sheet,'error':str(e)[:350]})
        except Exception as e: problems.append({'file':p.name,'error':str(e)[:350]})
    pd.DataFrame(sheets,columns=['file','sheet','rows','columns','fields','symbol_columns']).to_csv(out/'supplementary_table_inventory.csv',index=False,encoding='utf-8-sig')
    pd.DataFrame(matches,columns=['file','sheet','row_number_excel_approx','matched_gene','source_gene_column','source_gene_value','evidence_class','interpretation']).to_csv(out/'candidate_gene_presence_UNVALIDATED.csv',index=False,encoding='utf-8-sig')
    priority=['SST','BDNF','VGF','GAD1','GAD2','C1QA','C4A','C4B','S100A8','S100A9','FKBP5','PVALB','SYT1']
    matched_genes=set(x['matched_gene'] for x in matches)
    preview=strong[strong['gene_symbol'].astype(str).str.upper().isin(priority)].copy()
    preview['matched_in_supplementary_gene_column']=preview['gene_symbol'].astype(str).str.upper().isin(matched_genes)
    preview['protein_differential_abundance_confirmed']=False
    preview['protein_AD_CN_log2FC']=np.nan
    preview['protein_FDR']=np.nan
    preview.to_csv(out/'priority_gene_protein_evidence_status.csv',index=False,encoding='utf-8-sig')
    audit.update({'supplementary_workbooks_scanned':len(set(x['file'] for x in sheets)), 'sheets_scanned':len(sheets),
        'gene_mentions_matches_not_validations':len(matches),'unique_candidate_genes_in_sheets':len(matched_genes),
        'issues':problems,'protein_differential_abundance_confirmed':0,
        'next_step':'Inspect column definitions and supplement legends; determine which file contains sample-level TMT reporter intensities or verified protein comparisons. Only then calculate a valid AD/CN protein effect and FDR.'})
    (out/'module10_audit.json').write_text(json.dumps(audit,indent=2,ensure_ascii=False),encoding='utf-8')
    print('Inventory saved:',out/'supplementary_table_inventory.csv')
    print('Candidate gene matches (NOT protein verification):',len(matches))
    print('Audit:',out/'module10_audit.json')
    if not sheets:print('No readable supplementary XLSX found yet. Run: python main.py --download; or copy S1-S4 xlsx files into input/')

if __name__=='__main__':main()
