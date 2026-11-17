#!/usr/bin/env Rscript
# Conservative annotation-only SNP risk screen. This is not a genotype-based test.
suppressPackageStartupMessages({
  library(minfi)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})
args <- commandArgs(trailingOnly=TRUE)
outdir <- if(length(args)>=1) args[1] else 'Module23_inputs'
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)
anno <- minfi::getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19::IlluminaHumanMethylation450kanno.ilmn12.hg19)
cols <- colnames(anno)
cat('Available annotation columns:\n',paste(cols,collapse=', '),'\n')
# Only CpG and single-base extension site variants. Do NOT remove every probe with any SNP somewhere in its body.
required <- c('CpG_maf','SBE_maf')
if(!all(required %in% cols)) stop(paste('Required SNP frequency columns missing:',paste(setdiff(required,cols),collapse=', '),'\nCheck annotation version; no SNP risk list generated.'))
ids <- rownames(anno)
if(anyDuplicated(ids)) stop('Duplicate annotation identifiers')
cpg <- suppressWarnings(as.numeric(anno[['CpG_maf']]))
sbe <- suppressWarnings(as.numeric(anno[['SBE_maf']]))
if(length(cpg)!=length(ids) || length(sbe)!=length(ids)) stop('Malformed SNP frequency columns')
# MAF threshold is an operational sensitivity analysis threshold, not a clinical genetic threshold.
maf_cutoff <- 0.05
risk <- (is.finite(cpg)&cpg>=maf_cutoff) | (is.finite(sbe)&sbe>=maf_cutoff)
risk_ids <- ids[risk & grepl('^cg[0-9]{8}$',ids)]
writeLines(risk_ids,file.path(outdir,'450K_SNP_CpG_or_SBE_MAF_ge005.txt'))
meta <- c(
  'Source: IlluminaHumanMethylation450kanno.ilmn12.hg19 via minfi::getAnnotation',
  paste('R date:',Sys.Date()),
  paste('Annotation package version:',as.character(packageVersion('IlluminaHumanMethylation450kanno.ilmn12.hg19'))),
  paste('minfi version:',as.character(packageVersion('minfi'))),
  'Rules: exclude CpGs with CpG_maf >= 0.05 or SBE_maf >= 0.05',
  paste('Probes listed:',length(risk_ids)),
  paste('Annotation rows:',length(ids)),
  paste('CpG_maf nonmissing:',sum(is.finite(cpg))),
  paste('SBE_maf nonmissing:',sum(is.finite(sbe))),
  'This list is SNP-risk based on reference annotation, NOT cross-reactive detection or individual genotyping',
  'Missing MAF does NOT mean SNP-free; probe body SNP risk is not screened here'
)
writeLines(meta,file.path(outdir,'SNP_list_provenance.txt'))
cat(paste(meta,collapse='\n'),'\n')
