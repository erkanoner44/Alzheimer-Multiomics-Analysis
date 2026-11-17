#!/usr/bin/env python3
"""Module 25: descriptive gene-level transcript/protein/CpG integration.
Requires only pandas numpy matplotlib. No pooled p-values, no causal assertions.
"""
from pathlib import Path
import argparse, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

PRIORITY = ['VGF','NEUROD6','BDNF','CRH','FKBP5','C4A','SST','C4B','GAD1','GAD2','C1QA','PVALB','S100A9','S100A8','SYT1']

def present(v):
    return not pd.isna(v) and str(v).strip() not in ('','nan','None')

def run(inputs,results,figures):
    results.mkdir(parents=True,exist_ok=True);figures.mkdir(parents=True,exist_ok=True)
    tr=pd.read_csv(inputs/'priority_genes_transcript_evidence.csv')
    pr=pd.read_csv(inputs/'candidate_transcript_protein_comparison.csv')
    met=pd.read_csv(inputs/'MTG_filtered_FDR005_absBeta002.csv')
    annot=pd.read_csv(inputs/'filtered_CpG_gene_annotation.csv')
    targeted=pd.read_csv(inputs/'priority_genes_filtered_comparison.csv')
    audit=json.loads((inputs/'module23_audit.json').read_text(encoding='utf-8'))
    assert len(met)==audit['n_filtered_screen_hits'],(len(met),audit['n_filtered_screen_hits'])
    assert met.CpG.is_unique
    assert set(annot.CpG).issubset(set(met.CpG))
    assert met['new_hit'].fillna(False).astype(bool).all()
    sig=set(met.CpG)
    mapping=[]
    for z in annot.itertuples(index=False):
        genes=sorted({g.strip() for g in str(z.gene_symbol).split(';') if g.strip() and g.strip().lower()!='nan'})
        for gene in genes:
            mapping.append({'gene_symbol':gene,'CpG':z.CpG,'annotation_region':z.annotation_region})
    edges=pd.DataFrame(mapping).drop_duplicates(['gene_symbol','CpG'])
    edges=edges.merge(met[['CpG','adjusted_AD_beta','BH_FDR_filtered_genomewide']],on='CpG',how='left',validate='many_to_one')
    assert set(edges.CpG).issubset(sig)
    edges['direction']=np.where(edges.adjusted_AD_beta<0,'hypo','hyper')
    edges.sort_values(['gene_symbol','BH_FDR_filtered_genomewide','CpG']).to_csv(results/'Module25_CpG_gene_edges_filtered.csv',index=False)
    gp=edges.groupby('gene_symbol').agg(cpg_count=('CpG','nunique'),hypo_cpg=('direction',lambda x:int((x=='hypo').sum())),hyper_cpg=('direction',lambda x:int((x=='hyper').sum())),min_cpg_FDR=('BH_FDR_filtered_genomewide','min'),largest_abs_adj_beta=('adjusted_AD_beta',lambda x:float(abs(x).max())),cpg_ids=('CpG',lambda x:';'.join(sorted(set(x))))).reset_index()
    trsub=tr[['gene_symbol','logFC_representative','FDR_gene_BH','n_probes']].rename(columns={'logFC_representative':'transcript_log2FC','FDR_gene_BH':'transcript_gene_FDR','n_probes':'transcript_n_probes'})
    psub=pr[['gene_symbol','protein_detected_file3','protein_log2FC_best_p_entry','protein_BH_FDR_exploratory','protein_unique_peptides','transcript_protein_same_direction']].copy()
    psub=psub.rename(columns={'protein_log2FC_best_p_entry':'protein_log2FC','protein_BH_FDR_exploratory':'protein_exploratory_FDR'})
    psub=psub.drop_duplicates('gene_symbol')
    genes=sorted(set(trsub.gene_symbol)|set(psub.gene_symbol)|set(gp.gene_symbol))
    joined=pd.DataFrame({'gene_symbol':genes}).merge(trsub,on='gene_symbol',how='left').merge(psub,on='gene_symbol',how='left').merge(gp,on='gene_symbol',how='left')
    for col in ['cpg_count','hypo_cpg','hyper_cpg']:
        joined[col]=joined[col].fillna(0).astype(int)
    joined['priority_gene']=joined.gene_symbol.isin(PRIORITY)
    joined['transcript_significant']=joined.transcript_gene_FDR.lt(0.05).fillna(False)
    joined['protein_exploratory_significant']=joined.protein_exploratory_FDR.lt(0.05).fillna(False)
    joined['has_filtered_candidate_CpG']=joined.cpg_count.gt(0)
    joined['has_three_evidence_layers']=joined.transcript_gene_FDR.notna() & joined.protein_log2FC.notna() & joined.has_filtered_candidate_CpG
    joined['transcript_protein_concordant']=np.where(joined.transcript_log2FC.notna() & joined.protein_log2FC.notna(),np.sign(joined.transcript_log2FC)==np.sign(joined.protein_log2FC),pd.NA)
    joined['interpretation_note']='Transcript cohort, proteome cohort and methylome cohort are independent; CpG-gene annotation is not regulatory validation.'
    joined=joined.sort_values(['priority_gene','has_three_evidence_layers','cpg_count','gene_symbol'],ascending=[False,False,False,True])
    joined.to_csv(results/'Module25_all_gene_descriptive_evidence.csv',index=False)
    pri=joined[joined.priority_gene].copy()
    pri.to_csv(results/'Module25_priority_15_gene_evidence.csv',index=False)
    tgt=targeted.copy()
    tgt['gene_symbol']=tgt.UCSC_RefGene_Name.fillna('').astype(str).str.split(';')
    tgt=tgt.explode('gene_symbol').drop_duplicates(['gene_symbol','CpG'])
    tgt=tgt[tgt.gene_symbol.isin(PRIORITY)]
    tgt.to_csv(results/'Module25_priority_four_CpGs_audit.csv',index=False)
    summary={'status':'COMPLETED_DESCRIPTIVE_INTEGRATION','n_gene_records':len(joined),'n_priority_genes':len(pri),'n_all_filtered_candidate_cpg':len(met),'n_annotated_candidate_cpg':annot.CpG.nunique(),'n_candidate_cpg_gene_edges':len(edges),'n_genes_with_candidate_cpg':int(joined.has_filtered_candidate_CpG.sum()),'n_priority_genes_with_candidate_cpg':int(pri.has_filtered_candidate_CpG.sum()),'n_priority_three_evidence_layers':int(pri.has_three_evidence_layers.sum()),'priority_genes_with_cpg':pri.loc[pri.has_filtered_candidate_CpG,'gene_symbol'].tolist(),'SNP_filter':'Module23 SNP CpG/SBE MAF>=0.05 only; no cross-reactive-probe exclusion','important':'CpG FDR is CpG-level, proteomic FDR exploratory, no cohort-matched multiomic validation; gene FDR only from transcriptome.'}
    (results/'module25_audit.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
    plot=pri.copy().sort_values('gene_symbol')
    fig,ax=plt.subplots(figsize=(10.2,6.0))
    y=np.arange(len(plot))
    ax.barh(y-0.18,plot.hypo_cpg,height=.34,label='Hypomethylated CpGs')
    ax.barh(y+0.18,plot.hyper_cpg,height=.34,label='Hypermethylated CpGs')
    ax.set_yticks(y,plot.gene_symbol);ax.invert_yaxis();ax.set_xlabel('Unique SNP-filtered significant CpGs (annotation-linked)');ax.set_title('Module 25 | Priority genes: filtered methylation candidates');ax.legend(frameon=False);ax.grid(axis='x',alpha=.16);fig.tight_layout()
    for ext in ['png','pdf']: fig.savefig(figures/f'Figure25A_priority_gene_CpGs.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    p=pri[pri.transcript_log2FC.notna()].sort_values('transcript_log2FC').copy()
    fig,ax=plt.subplots(figsize=(10,5.9))
    idx=np.arange(len(p)); ax.barh(idx,p.transcript_log2FC)
    ax.set_yticks(idx,p.gene_symbol);ax.axvline(0,color='gray',lw=.8);ax.set_xlabel('Microarray transcript log2 fold-change (AD vs CN)');ax.set_title('Module 25 | Priority gene transcriptional effects');ax.grid(axis='x',alpha=.14);fig.tight_layout()
    for ext in ['png','pdf']:fig.savefig(figures/f'Figure25B_transcript_priority.{ext}',dpi=220,bbox_inches='tight')
    plt.close(fig)
    return summary

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--inputs',default='Module25_inputs');p.add_argument('--results',default='Module25_results');p.add_argument('--figures',default='Module25_figures');a=p.parse_args();print(json.dumps(run(Path(a.inputs),Path(a.results),Path(a.figures)),indent=2,ensure_ascii=False))
