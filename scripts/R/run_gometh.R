#!/usr/bin/env Rscript
# Module 21 - missMethyl with platform-specific CpG probe-number bias correction.
# GO / KEGG only; Reactome was queried separately via g:Profiler in Module 20.
args <- commandArgs(trailingOnly=TRUE)
argval <- function(flag, default=NULL) {
  w <- which(args == flag)
  if (length(w)==0L) return(default)
  if (w[length(w)]==length(args)) stop(paste('Missing value for',flag))
  args[w[length(w)]+1L]
}
in_dir <- argval('--inputs','Module21_inputs')
out_dir <- argval('--out','Module21_results')
run_directions <- '--directions' %in% args
if (!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)
log_file <- file.path(out_dir,'module21_run.log')
log_msg <- function(...) {
  msg <- paste(format(Sys.time(),'%Y-%m-%d %H:%M:%S'),paste(...,collapse=' '))
  cat(msg,'\n');cat(msg,'\n',file=log_file,append=TRUE)
}
if (!requireNamespace('missMethyl',quietly=TRUE)) stop('Install missMethyl: Rscript install_R_packages.R')
if (!requireNamespace('jsonlite',quietly=TRUE)) stop('Install jsonlite: Rscript install_R_packages.R')
if (!requireNamespace('IlluminaHumanMethylation450kanno.ilmn12.hg19',quietly=TRUE)) {
  stop('Install 450K annotation package: Rscript install_R_packages.R')
}
# missMethyl 1.46.0 uses a platform annotation object from the attached package
# search path; requireNamespace() alone only loads the namespace.
suppressPackageStartupMessages({
  library(missMethyl)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})
if (!("package:IlluminaHumanMethylation450kanno.ilmn12.hg19" %in% search())) {
  stop('450K annotation package was not attached to the R search path')
}
log_msg('450K annotation attached successfully')
readids <- function(name) {
  file <- file.path(in_dir,name)
  if (!file.exists(file)) stop(paste('Input file missing:',file))
  ids <- trimws(readLines(file, warn=FALSE, encoding='UTF-8'))
  ids <- ids[nzchar(ids)]
  if(anyDuplicated(ids)) stop(paste('Duplicate CpGs:',file))
  if(any(!grepl('^cg[0-9]{8}$',ids))) stop(paste('Invalid cg CpG ID:',file))
  ids
}
all_cpg <- readids('all_tested_cpg.txt')
sig_cpg <- readids('significant_cpg.txt')
hypo <- readids('hypomethylated_cpg.txt')
hyper <- readids('hypermethylated_cpg.txt')
if(length(all_cpg)<1000L) stop('Background too small; need genome-wide tested probes')
if(!all(sig_cpg %in% all_cpg) || !all(hypo %in% sig_cpg) || !all(hyper %in% sig_cpg)) stop('Significant sets must be subsets of tested universe')
if(length(intersect(hypo,hyper))>0) stop('Opposite direction lists overlap')
if(!setequal(sig_cpg,c(hypo,hyper))) stop('Direction sets do not partition significant set')
# We deliberately pass all.cpg; never default to full array for a filtered analysis.
log_msg('Universe:',length(all_cpg),'significant:',length(sig_cpg),'hypo:',length(hypo),'hyper:',length(hyper))
log_msg('missMethyl version:',as.character(packageVersion('missMethyl')))

collect <- list()
run_one <- function(group_name, ids, collection) {
  log_msg('Starting',group_name,collection)
  if(length(ids)==0L) stop('Cannot run enrichment with zero significant probes')
  ans <- missMethyl::gometh(sig.cpg=ids, all.cpg=all_cpg,
            collection=collection, array.type='450K',
            prior.prob=TRUE, plot.bias=FALSE)
  ans <- as.data.frame(ans)
  if(!('FDR' %in% names(ans))) stop('gometh output missing FDR column')
  ans$collection <- collection
  ans$input_group <- group_name
  ans <- ans[order(ans$FDR,na.last=TRUE),,drop=FALSE]
  dest <- file.path(out_dir,paste0('gometh_',group_name,'_',collection,'.csv'))
  write.csv(ans,dest,row.names=FALSE,na='')
  n_sig <- sum(!is.na(ans$FDR) & ans$FDR<0.05)
  log_msg('Finished',group_name,collection,'terms:',nrow(ans),'FDR<0.05:',n_sig)
  collect[[paste(group_name,collection,sep='_')]] <<- list(terms=nrow(ans),fdr_lt_005=n_sig,output=basename(dest))
  ans
}

status <- list(status='STARTED',method='missMethyl::gometh',array_type='450K',
               all_cpg=length(all_cpg),sig_cpg=length(sig_cpg),hypo_cpg=length(hypo),hyper_cpg=length(hyper),
               genomewide_FDR_source='Module 18',bias_adjustment='prior.prob=TRUE',
               caveats=c('No detection-p QC; SNP-risk exclusion is supplied through filtered input lists, not performed inside this script; no cross-reactivity blacklist applied',
                         'No cell composition correction',
                         'CpG-level bias-corrected GO/KEGG, not DMR or causal analysis'))
write_status <- function() jsonlite::write_json(status,file.path(out_dir,'module21_audit.json'),pretty=TRUE,auto_unbox=TRUE)
write_status()
tryCatch({
  for (coll in c('GO','KEGG')) run_one('pooled',sig_cpg,coll)
  if(run_directions) {
    for (g in c('hypo','hyper')) {
      subset_ids <- if(g=='hypo') hypo else hyper
      for (coll in c('GO','KEGG')) run_one(g,subset_ids,coll)
    }
  }
  status$status <- 'COMPLETED'
  status$results <- collect
  write_status()
  log_msg('ALL DONE. See:',out_dir)
},error=function(e) {
  status$status <- 'FAILED'
  status$error <- conditionMessage(e)
  status$results <- collect
  write_status()
  log_msg('ERROR:',conditionMessage(e))
  stop(e)
})
