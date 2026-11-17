#!/usr/bin/env Rscript
# Module 24. Recomputed missMethyl GO/KEGG on Mod23 SNP-risk-filtered CpG universe.
args <- commandArgs(trailingOnly=TRUE)
value <- function(flag, default) {
  pos <- which(args == flag)
  if(length(pos)==0) return(default)
  if(pos[length(pos)] >= length(args)) stop(paste('Missing value for',flag))
  args[pos[length(pos)]+1]
}
in_dir <- value('--inputs', 'Module24_inputs')
out_dir <- value('--out', 'Module24_results')
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
log_path <- file.path(out_dir,'module24_run.log')
logit <- function(...) {
  msg <- paste(format(Sys.time(), '%Y-%m-%d %H:%M:%S'), paste(...,collapse=' '))
  cat(msg,'\n');cat(msg,'\n',file=log_path,append=TRUE)
}
readids <- function(file) {
  path <- file.path(in_dir,file)
  if(!file.exists(path)) stop(paste('Missing:',path))
  x <- trimws(readLines(path, warn=FALSE))
  x <- x[nzchar(x)]
  if(anyDuplicated(x)) stop(paste('Duplicated CpG:',file))
  if(any(!grepl('^cg[0-9]{8}$',x))) stop(paste('Invalid CpG:',file))
  x
}
status <- list(status='STARTED',module=24,method='missMethyl::gometh',array_type='450K',
  filter='SNP-risk CpG/SBE MAF>=0.05 only; NO cross-reactive blacklist applied',
  model='Module 18 age+sex+plate original p values; Mod23 BH recomputed over retained probes',
  biases='missMethyl prior.prob=TRUE; annotated CpG multiplicity correction',
  caveats=c('No raw intensity/detection-p QC','No cell-composition correction','Not a DMR analysis'))
write_audit <- function() jsonlite::write_json(status,file.path(out_dir,'module24_audit.json'),pretty=TRUE,auto_unbox=TRUE)
tryCatch({
  for(pkg in c('missMethyl','jsonlite','IlluminaHumanMethylation450kanno.ilmn12.hg19')) {
    if(!requireNamespace(pkg,quietly=TRUE)) stop(paste('Package missing:',pkg))
  }
  suppressPackageStartupMessages({
    library(missMethyl)
    library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
  })
  if(!('package:IlluminaHumanMethylation450kanno.ilmn12.hg19' %in% search())) stop('450K annotation package not attached')
  all_cpg <- readids('all_tested_cpg.txt')
  sig <- readids('significant_cpg.txt')
  hypo <- readids('hypomethylated_cpg.txt')
  hyper <- readids('hypermethylated_cpg.txt')
  if(length(all_cpg)!=384360L || length(sig)!=11471L || length(hypo)!=7382L || length(hyper)!=4089L) stop('Counts disagree with Module23 audit; aborting')
  if(!all(sig %in% all_cpg) || !all(hypo %in% sig) || !all(hyper %in% sig)) stop('CpG subset mismatch')
  if(length(intersect(hypo,hyper)) || !setequal(sig,c(hypo,hyper))) stop('Direction partition mismatch')
  status$n_tested <- length(all_cpg); status$n_sig <- length(sig)
  status$n_hypo <- length(hypo);status$n_hyper <- length(hyper)
  status$missMethyl_version <- as.character(packageVersion('missMethyl'))
  status$results <- list();write_audit()
  for(group in c('pooled','hypo','hyper')) {
    current <- switch(group,pooled=sig,hypo=hypo,hyper=hyper)
    for(db in c('GO','KEGG')) {
      logit('Starting',group,db,'with filtered universe')
      ans <- as.data.frame(missMethyl::gometh(sig.cpg=current,all.cpg=all_cpg,
        collection=db,array.type='450K',prior.prob=TRUE,plot.bias=FALSE))
      if(!('FDR' %in% names(ans))) stop('Expected FDR column absent')
      ans$collection <- db;ans$input_group <- group
      ans <- ans[order(ans$FDR,na.last=TRUE),,drop=FALSE]
      filename <- paste0('gometh_filtered_',group,'_',db,'.csv')
      write.csv(ans,file.path(out_dir,filename),row.names=FALSE,na='')
      n <- sum(!is.na(ans$FDR) & ans$FDR < .05)
      status$results[[paste(group,db,sep='_')]] <- list(terms=nrow(ans),fdr_lt_005=n,filename=filename)
      write_audit()
      logit('Finished',group,db,'terms',nrow(ans),'FDR<0.05',n)
    }
  }
  status$status <- 'COMPLETED';write_audit()
  logit('ALL DONE; next: python compare_pathways.py')
},error=function(e){
  status$status <- 'FAILED';status$error <- conditionMessage(e)
  if(requireNamespace('jsonlite',quietly=TRUE)) write_audit()
  logit('ERROR:',conditionMessage(e))
  stop(e)
})
