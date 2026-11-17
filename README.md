# Alzheimer Multi-Omics Analysis — Research Code

**Author:** Erkan Öner · Adıyaman University, Faculty of Pharmacy, Department of Biochemistry, Türkiye  
**ORCID:** https://orcid.org/0000-0002-6332-6484

Research code supporting a cross-cohort Alzheimer’s disease study involving publicly available transcriptomic (GSE33000), proteomic (PXD010560), and methylation (GSE134379) datasets.

## What is included

- `scripts/R/`: archived R scripts for limma, Illumina 450K SNP-risk annotation, and `missMethyl` enrichment.
- `scripts/Python/`: archived modules for data preparation, CpG-level regression and filtering, protein integration, and supplementary-count verification.
- `data/CpG_lists/`: filtered CpG ID lists used as inputs to the enrichment step, **not patient-level raw data**.
- `docs/`: original code audit notes and warnings about external dependencies.

## Reproducibility status — read before use

**These are recovered original analysis scripts, not a fully independently validated end-to-end reproducibility pipeline.** A Python supplementary-table count-check was previously executed successfully; some other scripts have not been re-executed against all raw data in this environment. Some scripts require external public data, specific Python/R packages, upstream output files, and local file-path configuration. No validated claim of full re-execution is made.

For an independent run, obtain public data at:

- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33000
- https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE134379
- https://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD010560

### Typical analysis sequence

1. GSE33000 microarray processing and limma (`06_GSE33000_limma`)
2. Probe-to-gene summarization (`07_GSE33000_gene_aggregation`)
3. Proteomics matching (`10_PXD010560_protein`)
4. Methylation regression (`17_methyl_priority_covariates`, `18_methyl_genomewide_OLS`)
5. Common SNP-risk filtering (`23_SNP_risk_filtering`)
6. Filtered GO/KEGG enrichment (`24_filtered_missMethyl`, `run_gometh.R`)
7. Cross-omics integration (`25_three_omics_integration`)
8. Supplementary table audit (`verify_supplementary.py`)

These stages are **not** a turnkey command sequence; review each script's `--help`/arguments and comments.

## Example package setup

Python (version and compatible libraries depend on individual module):

```bash
python -m venv .venv
# activate the environment for your operating system
pip install -r requirements.txt
```

R packages commonly used by these scripts include `limma`, `missMethyl`, `IlluminaHumanMethylation450kanno.ilmn12.hg19`, `ggplot2`, and `jsonlite`; install the Bioconductor packages through `BiocManager`. Run R modules only after reviewing input data and version compatibility.

## Repository notes

- Do not commit raw identifiable participant data, secrets, access tokens, or intermediary download files.
- The manuscript presents independent *cross-cohort* comparison, not matched patient-level multi-omics integration.
- For a stable manuscript citation, create a tagged GitHub release and archive it with Zenodo after code inspection.
- See `docs/KOD_DENETIM_README_TR.md` for prior audit limitations.

## Citation

Associated manuscript: *Synaptic and neuroimmune signatures in Alzheimer's disease: cross-cohort integration of transcriptomics, proteomics and SNP-risk-filtered DNA methylation*. Manuscript submission status and DOI should be updated only when officially available.
