suppressPackageStartupMessages({library(limma);library(ggplot2)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=1) stop('Rscript limma_analysis.R OUTPUT_DIR')
out <- normalizePath(args[1],winslash='/',mustWork=TRUE)
meta <- read.csv(file.path(out,'sample_metadata_analysis.csv'),stringsAsFactors=FALSE,check.names=FALSE)
expr <- read.csv(gzfile(file.path(out,'expression_log2ratio_probes.csv.gz')),row.names=1,check.names=FALSE)
stopifnot(all(meta$gsm %in% colnames(expr)),!anyDuplicated(meta$gsm))
expr <- as.matrix(expr[,meta$gsm,drop=FALSE]);storage.mode(expr)<-'numeric'
if(any(!is.finite(expr))) {
  good <- rowSums(is.finite(expr)) == ncol(expr)
  message('Complete-case probe filter (removed): ',sum(!good))
  expr <- expr[good,,drop=FALSE]
}
if(nrow(expr)<1000 || ncol(expr)<30) stop('Too few probes or samples for complete analysis')
meta$group <- factor(meta$group,levels=c('CN','AD'))
meta$sex <- factor(tolower(meta$sex),levels=c('female','male'))
if(anyNA(meta$sex) || anyNA(meta$age)) stop('Age/sex missing')
# Use a continuous covariate; explicit sex adjustment; disease contrast AD - CN.
design <- model.matrix(~ group + age + sex,data=meta)
if(qr(design)$rank < ncol(design)) stop('Design matrix rank deficient')
fit <- eBayes(lmFit(expr,design),trend=TRUE,robust=FALSE)
res <- topTable(fit,coef='groupAD',number=Inf,sort.by='P',adjust.method='BH')
res$ID_REF <- rownames(res)
write.csv(res,file.path(out,'DEG_AD_vs_CN_age_sex_adjusted_probes.csv'),row.names=FALSE)
# Crude model for confounding diagnostics.
base_design <- model.matrix(~ group,data=meta)
base_fit <- eBayes(lmFit(expr,base_design),trend=TRUE)
base_res <- topTable(base_fit,coef='groupAD',number=Inf,adjust.method='BH')
base_res$ID_REF <- rownames(base_res)
write.csv(base_res,file.path(out,'DEG_AD_vs_CN_unadjusted_probes.csv'),row.names=FALSE)
# Note: sample-wise group differences in age/sex can limit identifiability even after adjustment.
qc <- data.frame(model=c('adjusted','unadjusted'),
  n_FDR_005=c(sum(res$adj.P.Val<0.05,na.rm=TRUE),sum(base_res$adj.P.Val<0.05,na.rm=TRUE)),
  n_FDR_005_abs_log2FC_1=c(sum(res$adj.P.Val<0.05 & abs(res$logFC)>=1,na.rm=TRUE),
                           sum(base_res$adj.P.Val<0.05 & abs(base_res$logFC)>=1,na.rm=TRUE)),
  probes_tested=nrow(expr),samples=ncol(expr))
write.csv(qc,file.path(out,'limma_summary.csv'),row.names=FALSE)
vol <- data.frame(logFC=res$logFC,nlog10p=-log10(pmax(res$adj.P.Val,.Machine$double.xmin)),
                  sig=res$adj.P.Val<0.05)
p <- ggplot(vol,aes(logFC,nlog10p,color=sig))+geom_point(size=.7,alpha=.6)+
  geom_vline(xintercept=c(-1,1),linetype='dashed',color='gray50')+
  geom_hline(yintercept=-log10(.05),linetype='dashed',color='gray50')+
  scale_color_manual(values=c('FALSE'='gray65','TRUE'='#B22222'))+
  labs(title='GSE33000 AD vs CN | adjusted for age and sex',
       x='log2 fold change (AD - CN)',y='-log10(BH FDR)',color='FDR < 0.05')+
  theme_minimal()
ggsave(file.path(out,'volcano_adjusted.png'),plot=p,width=9,height=6,dpi=200)
# Unsupervised PCA on top variable probes, not a disease model.
variance <- apply(expr,1,var)
keep <- order(variance,decreasing=TRUE)[seq_len(min(2000,nrow(expr)))]
pca <- prcomp(t(expr[keep,,drop=FALSE]),center=TRUE,scale.=FALSE)
coords <- data.frame(PC1=pca$x[,1],PC2=pca$x[,2],group=meta$group)
variance_explained <- 100 * pca$sdev^2 / sum(pca$sdev^2)
pp <- ggplot(coords,aes(PC1,PC2,color=group))+geom_point(alpha=.75,size=1.8)+theme_minimal()+
  labs(title='GSE33000 PCA (top 2000 variable probes)',
       x=sprintf('PC1 (%.1f%%)',variance_explained[1]),y=sprintf('PC2 (%.1f%%)',variance_explained[2]))
ggsave(file.path(out,'PCA_top_variable_probes.png'),plot=pp,width=8,height=6,dpi=200)
write.csv(coords,file.path(out,'PCA_coordinates.csv'),row.names=FALSE)
message('limma tamamlandı. Output: ',out)
