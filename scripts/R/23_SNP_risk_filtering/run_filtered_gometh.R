#!/usr/bin/env Rscript
# Run ONLY after Module23 Python returns FILTERED_SENSITIVITY_COMPLETED.
suppressPackageStartupMessages({
  library(missMethyl)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})
args <- commandArgs(trailingOnly=TRUE)
inputdir <- if(length(args)>=1) args[[1]] else 'Module23_results'
outdir <- if(length(args)>=2) args[[2]] else 'Module23_gometh_results'
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
if(!file.exists(file.path(inputdir,'module23_audit.json'))) stop('Module23 audit file missing')
if(!requireNamespace('jsonlite',quietly=TRUE)) stop('Install jsonlite: install.packages("jsonlite")')
audit<-jsonlite::fromJSON(file.path(inputdir,'module23_audit.json'))
if(!identical(audit$status,'FILTERED_SENSITIVITY_COMPLETED')) stop('Prob filtering has not been confirmed; stop')
load_ids <- function(filename) {
  p <- file.path(inputdir,filename)
  if(!file.exists(p)) stop(paste('Missing:',p))
  unique(trimws(readLines(p,warn=FALSE)))
}
background <- load_ids('filtered_tested_CpGs.txt')
if(length(background)<1000) stop('Background too small')
sets <- list(pooled=load_ids('filtered_candidate_CpGs.txt'),
             hypo=load_ids('filtered_hypo_CpGs.txt'),
             hyper=load_ids('filtered_hyper_CpGs.txt'))
for(group in names(sets)){
  ids <- sets[[group]]
  if(length(ids)<20) {warning(paste('Too few CpGs for',group,'; skipped'));next}
  if(!all(ids %in% background)) stop('Candidate IDs outside filtered background')
  for(db in c('GO','KEGG')) {
    message('Running ',group,' ',db,' with ',length(ids),' candidates')
    z <- missMethyl::gometh(sig.cpg=ids,all.cpg=background,
          collection=db,array.type='450K',prior.prob=TRUE,
          equiv.cpg=TRUE,fract.counts=TRUE)
    write.csv(z,file.path(outdir,paste0('gometh_filtered_',group,'_',db,'.csv')),row.names=FALSE)
  }
}
writeLines(c(paste('Analysis completed:',Sys.time()),
             paste('Filtered background:',length(background)),
             paste('Pooled:',length(sets$pooled)),
             paste('Hypo:',length(sets$hypo)),
             paste('Hyper:',length(sets$hyper))),file.path(outdir,'filtered_gometh_run_summary.txt'))
