#!/usr/bin/env python3
"""Independently verify published counts from supplied supplementary files.

Python 3.9+, standard library only. This is NOT a raw-data analysis pipeline.
Usage: python Code/verify_supplementary.py --package revised_package.zip
       python Code/verify_supplementary.py --folder path/to/extracted/package
"""
import argparse, csv, gzip, io, json, zipfile
from pathlib import Path

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    g=ap.add_mutually_exclusive_group(required=True)
    g.add_argument('--package',type=Path)
    g.add_argument('--folder',type=Path)
    ap.add_argument('--write-gometh-inputs',type=Path,help='Create filtered tested/significant/hypo/hyper CpG files')
    a=ap.parse_args()
    arc=zipfile.ZipFile(a.package) if a.package else None
    def get(path):
        if arc: return arc.read(path)
        return (a.folder/path).read_bytes()
    def rows(n):
        raw=get('Supplementary/'+n)
        if n.endswith('.gz'): raw=gzip.decompress(raw)
        return csv.DictReader(io.StringIO(raw.decode('utf-8-sig')))
    results={}; tests=[]
    def check(name,value,expected):
        results[name]={'observed':value,'expected':expected,'pass':value==expected}
        tests.append(value==expected)
    sig=list(rows('ESM_13_S5_filtered_CpG_hits.csv'))
    check('candidate_CpGs',len(sig),11471)
    check('hypomethylated',sum(float(r['adjusted_AD_beta'])<0 for r in sig),7382)
    check('hypermethylated',sum(float(r['adjusted_AD_beta'])>0 for r in sig),4089)
    source={
      'pooled_GO':('ESM_15_S7_pooled_GO.csv',174),
      'pooled_KEGG':('ESM_16_S8_pooled_KEGG.csv',18),
      'hypo_GO':('ESM_17_S9_hypomethylated_GO.csv',25),
      'hypo_KEGG':('ESM_2_S10_hypomethylated_KEGG.csv',1),
      'hyper_GO':('ESM_3_S11_hypermethylated_GO.csv',167),
      'hyper_KEGG':('ESM_4_S12_hypermethylated_KEGG.csv',10),
    }
    for name,(file,count) in source.items():
        check(name,sum(float(x['FDR'])<0.05 for x in rows(file) if x['FDR']),count)
    protein=[x for x in rows('ESM_6_S14_proteome_transcript_comparison.csv') if x['transcript_protein_same_direction'].strip().lower() in ('same','opposite')]
    check('protein_transcript_pairs',len(protein),74)
    check('same_direction',sum(r['transcript_protein_same_direction'].lower()=='same' for r in protein),60)
    check('opposite_direction',sum(r['transcript_protein_same_direction'].lower()=='opposite' for r in protein),14)
    if a.write_gometh_inputs:
        d=a.write_gometh_inputs; d.mkdir(parents=True,exist_ok=True)
        all_ids=[]
        for r in rows('ESM_14_S6_full_CpG_QC_and_statistics.csv.gz'):
            if r['retained_after_filter'].lower()=='true':all_ids.append(r['CpG'])
        check('filtered_tested_CpGs',len(all_ids),384360)
        groups={'all_tested_cpg.txt':all_ids,'significant_cpg.txt':[r['CpG'] for r in sig],
                'hypomethylated_cpg.txt':[r['CpG'] for r in sig if float(r['adjusted_AD_beta'])<0],
                'hypermethylated_cpg.txt':[r['CpG'] for r in sig if float(r['adjusted_AD_beta'])>0]}
        for name,ids in groups.items(): (d/name).write_text('\n'.join(ids)+'\n',encoding='utf-8')
    print(json.dumps(results,indent=2,ensure_ascii=False))
    if not all(tests): raise SystemExit('FAIL: one or more published counts disagree with supplementary tables')
    print('PASS: All checked supplementary-data counts match the manuscript.')
    if arc:arc.close()
if __name__=='__main__':main()
