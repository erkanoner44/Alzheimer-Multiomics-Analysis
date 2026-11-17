#!/usr/bin/env python3
"""GSE33000 GPL4372 annotation, conservative gene-level stats and optional g:Profiler ORA."""
import argparse
import gzip
import json
import re
from pathlib import Path
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests


def extract_platform(soft):
    in_table = False
    header = None
    rows = []
    with gzip.open(soft, 'rt', encoding='utf-8', errors='replace') as f:
        for line in f:
            line = line.rstrip('\r\n')
            if line == '!platform_table_begin':
                in_table = True
                continue
            if not in_table:
                continue
            if line == '!platform_table_end':
                break
            if header is None:
                header = line.split('\t')
                continue
            vals = line.split('\t')
            if len(vals) == len(header):
                rows.append(vals)
    if not header or not rows:
        raise ValueError('SOFT dosyasında GPL4372 platform tablosu bulunamadı.')
    plat = pd.DataFrame(rows, columns=header)
    for col in ['ID','ORF','EntrezGeneID']:
        if col not in plat:
            raise ValueError(f'Platformda {col} yok.')
    return plat


def prepare_probe(deg, platform):
    deg = pd.read_csv(deg, dtype={'ID_REF':str})
    required = {'ID_REF','logFC','P.Value','adj.P.Val','AveExpr'}
    if not required.issubset(deg):
        raise ValueError('DEG sütunları eksik: '+str(required-set(deg)))
    if deg.ID_REF.duplicated().any():
        raise ValueError('Aynı prob birden fazla kayıt içeriyor.')
    platform = platform.drop_duplicates('ID').rename(columns={'ID':'ID_REF','ORF':'gene_symbol','EntrezGeneID':'entrez_id'})
    merged = deg.merge(platform, on='ID_REF',how='left',validate='one_to_one',indicator=True)
    merged['gene_symbol'] = merged['gene_symbol'].fillna('').astype(str).str.strip()
    merged['entrez_id'] = merged['entrez_id'].fillna('').astype(str).str.strip()
    # Retain only single Entrez IDs with a platform-provided symbol, avoid ambiguous IDs.
    valid = merged.gene_symbol.str.match(r'^[A-Za-z][A-Za-z0-9_.@-]*$') & merged.entrez_id.str.fullmatch(r'\d+')
    merged['usable_gene_mapping'] = valid
    return merged


def summarize_genes(probes):
    x=probes[probes.usable_gene_mapping].copy()
    x['P.Value'] = pd.to_numeric(x['P.Value'],errors='coerce')
    x['logFC'] = pd.to_numeric(x['logFC'],errors='coerce')
    x=x.dropna(subset=['P.Value','logFC'])
    x=x[(x['P.Value']>=0)&(x['P.Value']<=1)]
    # Use NCBI Entrez ID, rather than potentially outdated ambiguous symbols, as grouping key.
    results=[]
    for eid,g in x.groupby('entrez_id',sort=False):
        g=g.sort_values(['P.Value','ID_REF'],kind='stable')
        best=g.iloc[0]
        m=len(g)
        # Valid under arbitrary dependence within gene. Conservative Bonferroni minP.
        p_gene=min(1.0,float(best['P.Value'])*m)
        signs=np.sign(g.logFC.to_numpy())
        sign_consistent = bool((signs>=0).all() or (signs<=0).all())
        results.append({'entrez_id':eid,'gene_symbol':best.gene_symbol,'n_probes':m,
            'representative_probe':best.ID_REF,'logFC_representative':best.logFC,
            'logFC_median':float(g.logFC.median()),'p_gene_bonferroni_minP':p_gene,
            'probes_direction_consistent':sign_consistent,
            'min_probe_FDR':float(g['adj.P.Val'].min())})
    res=pd.DataFrame(results)
    if res.empty:
        raise ValueError('Platformda hiçbir güvenilir gen eşleştirmesi bulunamadı.')
    res['FDR_gene_BH'] = multipletests(res.p_gene_bonferroni_minP.to_numpy(),method='fdr_bh')[1]
    res['significant_FDR_005'] = res.FDR_gene_BH < .05
    res['strong_FDR_005_FC_05'] = res.significant_FDR_005 & (res.logFC_representative.abs()>=.5)
    res['strong_FDR_005_FC_1'] = res.significant_FDR_005 & (res.logFC_representative.abs()>=1)
    return res.sort_values(['FDR_gene_BH','p_gene_bonferroni_minP'])


def gprofiler_enrichment(genetable, out, threshold=0.05):
    try:
        import requests
        subset=genetable[(genetable.FDR_gene_BH<threshold)&(genetable.logFC_representative.abs()>=0.5)&genetable.probes_direction_consistent]
        background=genetable.entrez_id.astype(str).tolist()
        query=subset.entrez_id.astype(str).tolist()
        if len(query)<5:
            return {'status':'SKIPPED','reason':'Anlamlı gen sayısı <5.'}
        payload={'organism':'hsapiens','query':query,'background':background,
            'sources':['GO:BP','GO:MF','GO:CC','KEGG','REAC'],
            'user_threshold':0.05,'significance_threshold_method':'g_SCS',
            'domain_scope':'custom','no_evidences':True}
        r=requests.post('https://biit.cs.ut.ee/gprofiler/api/gost/profile/',json=payload,timeout=90)
        r.raise_for_status()
        obj=r.json();data=obj.get('result',[])
        cols=['source','native','name','p_value','term_size','query_size','intersection_size','effective_domain_size','significant']
        pd.DataFrame(data).reindex(columns=cols).to_csv(out/'GO_KEGG_Reactome_gprofiler.csv',index=False)
        (out/'gprofiler_response.json').write_text(json.dumps(obj,ensure_ascii=False,indent=2),encoding='utf-8')
        return {'status':'OK','terms':len(data),'query_gene_count':len(query),'background_gene_count':len(background)}
    except Exception as exc:
        return {'status':'NOT_RUN','reason':str(exc),'note':'No enrichment results invented. Check internet/g:Profiler access.'}


def main():
    a=argparse.ArgumentParser()
    a.add_argument('--soft',required=True,type=Path)
    a.add_argument('--adjusted',required=True,type=Path)
    a.add_argument('--unadjusted',type=Path,default=None)
    a.add_argument('--out',type=Path,default=Path('Module07_results'))
    a.add_argument('--enrich',action='store_true',help='Requires internet; g:Profiler GO/KEGG/Reactome with measured gene background')
    args=a.parse_args();out=args.out;out.mkdir(parents=True,exist_ok=True)
    plat=extract_platform(args.soft)
    plat.to_csv(out/'GPL4372_platform_from_SOFT.csv',index=False)
    ap=prepare_probe(args.adjusted,plat)
    ap.to_csv(out/'annotated_adjusted_probes.csv.gz',index=False,compression='gzip')
    genes=summarize_genes(ap)
    genes.to_csv(out/'adjusted_gene_level_results.csv',index=False)
    for cut in ['strong_FDR_005_FC_05','strong_FDR_005_FC_1']:
        genes.loc[genes[cut]].to_csv(out/(cut+'.csv'),index=False)
    summary={'analysis':'GSE33000 AD vs CN (age + sex adjusted)','platform':'GPL4372',
        'n_platform_rows':len(plat),'n_DEG_probes':len(ap),
        'n_probes_mapped_to_unique_entrez_and_symbol':int(ap.usable_gene_mapping.sum()),
        'n_genes_tested':len(genes),
        'n_genes_FDR_lt_005':int(genes.significant_FDR_005.sum()),
        'n_genes_FDR_lt_005_abs_logFC_ge_05':int(genes.strong_FDR_005_FC_05.sum()),
        'n_genes_FDR_lt_005_abs_logFC_ge_1':int(genes.strong_FDR_005_FC_1.sum()),
        'note':'Gene p = min probe raw p * n probes (Bonferroni); BH across genes; representative logFC is exploratory, not fitted gene-level model.'}
    if args.unadjusted:
        up=prepare_probe(args.unadjusted,plat)
        ug=summarize_genes(up)
        ug.to_csv(out/'unadjusted_gene_level_results.csv',index=False)
        compare=genes[['entrez_id','gene_symbol','FDR_gene_BH','logFC_representative']].merge(
            ug[['entrez_id','FDR_gene_BH','logFC_representative']],on='entrez_id',suffixes=('_adjusted','_unadjusted'))
        compare.to_csv(out/'adjusted_vs_unadjusted_gene_comparison.csv',index=False)
        summary['unadjusted_n_genes_FDR_lt_005']=int(ug.significant_FDR_005.sum())
    if args.enrich:
        summary['enrichment']=gprofiler_enrichment(genes,out)
    else:
        summary['enrichment']={'status':'NOT_RUN','reason':'Pass --enrich to run online ORA.'}
    (out/'module07_audit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    print(json.dumps(summary,ensure_ascii=False,indent=2))

if __name__=='__main__':main()
