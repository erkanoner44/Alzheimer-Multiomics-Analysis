#!/usr/bin/env python3
"""Module 24: compare existing Module21 and newly run Module24 g:Profiler-independent missMethyl output."""
import argparse,csv,json,math
from pathlib import Path

def load(path):
    if not path.is_file(): raise FileNotFoundError(f'Missing: {path}')
    with path.open(encoding='utf-8-sig',newline='') as f:
        data=list(csv.DictReader(f))
    if not data or 'FDR' not in data[0] or not ('Description' in data[0] or 'TERM' in data[0]):
        raise ValueError(f'Invalid GO/KEGG result: {path}')
    by_name={}
    for r in data:
        name=(r.get('Description') or r.get('TERM') or '').strip()
        if not name: continue
        if name in by_name: raise ValueError(f'Duplicate pathway in {path}: {name}')
        raw=r.get('FDR','')
        fdr=float(raw) if raw.strip() else math.nan
        by_name[name]={'description':name,'FDR':fdr,'row':r}
    return by_name

def fmt(x):return '' if not math.isfinite(x) else f'{x:.12g}'

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--before',default='Module21_reference')
    p.add_argument('--after',default='Module24_results')
    p.add_argument('--out',default='Module24_comparison')
    a=p.parse_args();before=Path(a.before);after=Path(a.after);out=Path(a.out)
    audit_file=after/'module24_audit.json'
    if not audit_file.is_file(): raise FileNotFoundError(f'Run R first: {audit_file}')
    audit=json.loads(audit_file.read_text(encoding='utf-8'))
    if audit.get('status')!='COMPLETED': raise RuntimeError('Filtered R enrichment not COMPLETED')
    if (audit.get('n_tested'),audit.get('n_sig'))!=(384360,11471): raise RuntimeError('Unexpected input counts')
    out.mkdir(exist_ok=True,parents=True)
    records=[];counts=[];summary=[]
    for group in ('pooled','hypo','hyper'):
        for db in ('GO','KEGG'):
            b=load(before/f'gometh_{group}_{db}.csv')
            n=load(after/f'gometh_filtered_{group}_{db}.csv')
            pre_sig={x for x,r in b.items() if math.isfinite(r['FDR']) and r['FDR']<.05}
            post_sig={x for x,r in n.items() if math.isfinite(r['FDR']) and r['FDR']<.05}
            for name in sorted(set(b)|set(n)):
                old=b.get(name,{}).get('FDR',math.nan)
                new=n.get(name,{}).get('FDR',math.nan)
                records.append({'group':group,'collection':db,'pathway':name,
                    'FDR_before_SNP_filter':fmt(old),'FDR_after_SNP_filter':fmt(new),
                    'significant_before':name in pre_sig,'significant_after':name in post_sig,
                    'category':('retained' if name in pre_sig&post_sig else 'lost' if name in pre_sig else 'gained' if name in post_sig else 'neither')})
            counts.append({'group':group,'collection':db,'tested_terms_before':len(b),'tested_terms_after':len(n),
                'sig_before':len(pre_sig),'sig_after':len(post_sig),'retained':len(pre_sig&post_sig),
                'lost':len(pre_sig-post_sig),'gained':len(post_sig-pre_sig)})
            summary.append(f'{group} {db}: {len(pre_sig)} -> {len(post_sig)}; retained={len(pre_sig&post_sig)}, lost={len(pre_sig-post_sig)}, gained={len(post_sig-pre_sig)}')
    def writecsv(path,rows):
        with path.open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    writecsv(out/'pathway_before_after_ALL.csv',records)
    writecsv(out/'pathway_significance_comparison.csv',counts)
    salient=('Glutamatergic synapse','Dopaminergic synapse','Axon guidance','Cholinergic synapse','GABAergic synapse','Long-term potentiation','Leukocyte transendothelial migration','Synapse','postsynapse')
    select=[r for r in records if any(k.casefold() == r['pathway'].casefold() for k in salient)]
    writecsv(out/'priority_pathways.csv',select)
    result={'status':'COMPARISON_COMPLETED','module24_r_status':'COMPLETED','counts':counts,
       'method_note':'Separate missMethyl::gometh with platform CpG-number bias adjustment; BEFORE and AFTER cohorts identical but SNP-risk universe differs.',
       'limitations':['SNP filter only; no cross-reactive blacklist','Unavailable detection-p QC','No cell composition correction','GO labels need biological specificity review']}
    (out/'module24_comparison_audit.json').write_text(json.dumps(result,ensure_ascii=False,indent=2),encoding='utf-8')
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        labels=[f"{r['group']} {r['collection']}" for r in counts]
        fig,ax=plt.subplots(figsize=(10,5.5));x=list(range(len(counts)));wd=.38
        ax.bar([v-wd/2 for v in x],[r['sig_before'] for r in counts],wd,label='Before SNP QC')
        ax.bar([v+wd/2 for v in x],[r['sig_after'] for r in counts],wd,label='After SNP QC')
        ax.set_xticks(x,labels,rotation=25,ha='right');ax.set_ylabel('Significant enriched terms (FDR < 0.05)')
        ax.set_title('Module 24: missMethyl pathway sensitivity after SNP-risk filtering')
        ax.legend(frameon=False);fig.tight_layout()
        for ext in ('png','pdf'): fig.savefig(out/f'Figure24A_pathway_sensitivity.{ext}',dpi=250)
        plt.close(fig)
        selected=[r for r in select if r['collection']=='KEGG' and r['group']=='pooled']
        if selected:
            fig,ax=plt.subplots(figsize=(10,5.5));names=[r['pathway'] for r in selected]
            x=list(range(len(names)))
            for shift,key,label in ((-.2,'FDR_before_SNP_filter','Before'),(.2,'FDR_after_SNP_filter','After')):
                values=[-math.log10(max(float(r[key]),1e-300)) if r[key] else 0 for r in selected]
                ax.barh([v+shift for v in x],values,height=.38,label=label)
            ax.set_yticks(x,names);ax.invert_yaxis();ax.axvline(-math.log10(.05),linestyle='--',linewidth=1)
            ax.set_xlabel('-log10(FDR)');ax.legend(frameon=False)
            ax.set_title('Selected KEGG pathways: sensitivity to SNP-risk filtering')
            fig.tight_layout()
            for ext in ('png','pdf'): fig.savefig(out/f'Figure24B_selected_KEGG.{ext}',dpi=250)
            plt.close(fig)
    except ImportError:
        print('matplotlib not installed; CSV and audit available, figures skipped')
    print('\n'.join(summary))
    print('Results:',out.resolve())
if __name__=='__main__':main()
