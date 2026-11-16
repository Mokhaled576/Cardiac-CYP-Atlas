args <- commandArgs(TRUE)
old <- read.delim(args[1])
new <- read.delim(args[2])
out <- args[3]
direct <- subset(new, mode=="all_donor_direct")
stopifnot(length(unique(direct$contrast))==4L)
parent <- list(
  Chaffin_DCM_vs_HCM_LV=c("Chaffin_DCM_LV","Chaffin_HCM_LV"),
  Reichart_DCM_vs_ARVC_LV=c("Reichart_DCM_LV","Reichart_ARVC_LV"),
  Reichart_DCM_vs_ARVC_RV=c("Reichart_DCM_RV","Reichart_ARVC_RV"),
  DoxCM_vs_NICM_LV=c("DoxCM_LV","NICM_LV")
)
rows <- list()
pdf(file.path(out,"Direct_Etiology_Volcanoes.pdf"), width=8,height=6)
for (k in names(parent)) {
  d <- subset(direct,contrast==k)
  y <- -log10(pmax(d$FDR,.Machine$double.xmin))
  plot(d$logFC,y,pch=19,col=ifelse(d$FDR<.05,"#AC4633","#678797"),
       main=k,xlab="Direct disease comparison log2 fold change",ylab="-log10(FDR)")
  abline(v=0,h=-log10(.05),lty=2,col="grey55")
  for (arm in seq_along(parent[[k]])) {
    ref <- subset(old,population=="all_donors" & contrast==parent[[k]][arm])
    joined <- merge(d[c("lineage","cyp_gene","logFC","FDR")],
                    ref[c("lineage","cyp_gene","logFC","CYP_BH_across_tier")],
                    by=c("lineage","cyp_gene"),suffixes=c("_direct","_vs_control"))
    rows[[length(rows)+1L]] <- transform(joined,direct_contrast=k,control_contrast=parent[[k]][arm])
  }
}
dev.off()
joined <- do.call(rbind,rows)
write.table(joined,file.path(out,"Direct_vs_Control_Effects.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
counts <- aggregate(cbind(eligible=rep(1,nrow(joined)),
                          agree=as.integer(sign(joined$logFC_direct)==sign(joined$logFC_vs_control))),
                    joined[c("direct_contrast","control_contrast")],sum)
write.table(counts,file.path(out,"Direct_vs_Control_Concordance.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
