args <- commandArgs(TRUE)
archive <- args[1]
out <- args[2]
source(file.path(out,"ql_wald_ci.R"))
library(edgeR)
direct <- data.frame(
  contrast=c("Chaffin_DCM_vs_HCM_LV","Reichart_DCM_vs_ARVC_LV",
             "Reichart_DCM_vs_ARVC_RV","DoxCM_vs_NICM_LV"),
  cohort=c(1,7,7,2),region=c("LV","LV","RV","LV"),
  case=c("DCM",rep("dilated cardiomyopathy",2),"DoxCM"),
  control=c("HCM",rep("arrhythmogenic right ventricular cardiomyopathy",2),"NICM")
)
read_tsv <- function(f) read.delim(f,check.names=FALSE,colClasses="character")
rows <- list()
for (i in seq_len(nrow(direct))) {
  c <- direct[i,]
  root <- file.path(archive,paste0("Cohort_",c$cohort))
  obs <- read_tsv(file.path(root,"Pseudobulk_Donors.tsv"))
  ft <- read_tsv(file.path(root,"Pseudobulk_Features.tsv"))
  mat <- as.matrix(read.delim(gzfile(file.path(root,"Pseudobulk_Counts.tsv.gz")),
                             row.names=1,check.names=FALSE))
  stopifnot(identical(colnames(mat),obs$sample),identical(rownames(mat),ft$feature_id))
  for (lin in unique(obs$lineage[obs$cardiac_region==c$region])) {
    use <- obs$lineage==lin & obs$cardiac_region==c$region & obs$disease %in% c(c$case,c$control) &
      as.integer(obs$nuclei)>=20L
    group <- factor(ifelse(obs$disease[use]==c$case,"Case","Control"),levels=c("Control","Case"))
    if (any(table(group)<3L)) next
    design <- model.matrix(~group)
    d <- DGEList(counts=mat[,use,drop=FALSE])
    keep <- filterByExpr(d,design=design)
    if (!any(keep)) next
    d <- calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method="TMM")
    fit <- glmQLFit(estimateDisp(d,design,robust=TRUE),design,robust=TRUE)
    z <- add_ql_wald_ci(fit,glmQLFTest(fit,coef="groupCase"),coef="groupCase")
    gene <- ft$cyp_gene[match(rownames(z),ft$feature_id)]
    ok <- !is.na(gene) & nzchar(gene)
    if (any(ok)) rows[[length(rows)+1L]] <- cbind(
      c[rep(1,sum(ok)),],lineage=lin,cyp_gene=gene[ok],z[ok,]
    )
  }
}
d <- do.call(rbind,rows)
stopifnot(length(unique(d$contrast))==4L)
prior <- subset(read.delim(file.path(out,"..","sex_etiology","CYP_Effects.tsv")),mode=="all_donor_direct")
key <- function(x) with(x,paste(contrast,lineage,cyp_gene,sep="/"))
stopifnot(!anyDuplicated(key(d)),!anyDuplicated(key(prior)),setequal(key(d),key(prior)))
prior <- prior[match(key(d),key(prior)),]
stopifnot(isTRUE(all.equal(d$logFC,prior$logFC,tolerance=1e-7)),
          isTRUE(all.equal(d$PValue,prior$PValue,tolerance=1e-7)))
d$FDR_direct <- prior$FDR
write.table(d,file.path(out,"Direct_Comparison_With_CI.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
focal <- subset(d,(lineage=="cardiomyocyte" & cyp_gene %in% c("CYP2J2","CYP11A1")) |
                  (lineage=="fibroblast" & cyp_gene %in% c("CYP19A1","CYP4X1")))
lims <- range(0,focal$CI_95_low,focal$CI_95_high,na.rm=TRUE)
pdf(file.path(out,"Direct_Comparison_Forest.pdf"),width=7,height=6.1,useDingbats=FALSE)
par(mar=c(4,15,2,1),family="sans",las=1)
nr <- nrow(focal)
y <- rev(seq_len(nr))
col <- ifelse(focal$FDR_direct<.05,"#0072B2","#646D75")
plot(focal$logFC_unshrunk,y,xlim=lims,ylim=c(0,nr+1),type="n",yaxt="n",bty="n",
     xlab=expression("Unshrunk log"[2]*" fold change"),ylab="")
abline(v=0,lty=2,col="grey65")
segments(focal$CI_95_low,y,focal$CI_95_high,y,col=col,lwd=1.4)
points(focal$logFC_unshrunk,y,col=col,pch=19,cex=.8)
labels <- gsub("_"," ",paste(focal$contrast,focal$cyp_gene))
axis(2,at=y,labels=labels,tick=FALSE,cex.axis=.64)
dev.off()
