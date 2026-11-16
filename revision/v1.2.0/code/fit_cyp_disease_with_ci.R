args <- commandArgs(TRUE)
if (length(args) != 2L) stop("Usage: Rscript fit_cyp_disease_with_ci.R INPUT_DIR NEW_OUTPUT_DIR")
path <- args[[1]]
destination <- args[[2]]
if (!dir.create(destination, recursive=FALSE)) stop("Output directory must be new")
script_arg <- grep("^--file=", commandArgs(FALSE), value=TRUE)
source(file.path(dirname(sub("^--file=", "", script_arg)), "ql_wald_ci.R"))
suppressPackageStartupMessages(library(edgeR))
read <- function(file) read.delim(file.path(path,file),check.names=FALSE,colClasses="character")
obs <- read("Pseudobulk_Donors.tsv")
obs$nuclei <- as.integer(obs$nuclei)
features <- read("Pseudobulk_Features.tsv")
contrasts <- read("Contrasts.tsv")
counts <- as.matrix(read.delim(gzfile(file.path(path,"Pseudobulk_Counts.tsv.gz")),row.names=1,check.names=FALSE))
stopifnot(identical(rownames(counts),features$feature_id),identical(colnames(counts),obs$sample),
          !anyDuplicated(obs[c("cardiac_region","lineage","donor_id")]),all(obs$nuclei>=20))
result <- list()
audit <- list()
for (i in seq_len(nrow(contrasts))) {
  c <- contrasts[i,]
  for (lin in sort(unique(obs$lineage[obs$cardiac_region==c$region]))) {
    within <- obs$cardiac_region==c$region & obs$lineage==lin & obs$disease %in% c(c$case,c$control)
    for (population in c("all_donors","exclude_flagged_overlap")) {
      selected <- within
      if (population=="exclude_flagged_overlap") selected <- selected & obs$overlap_flag=="no"
      group <- factor(ifelse(obs$disease[selected]==c$case,"Case","Control"),levels=c("Control","Case"))
      n <- table(group)
      reason <- if (any(n<3L)) "fewer_than_3_per_group" else "tested"
      if (population=="exclude_flagged_overlap" && identical(selected,within)) reason <- "identical_to_all_donors"
      key <- paste(c$contrast,lin,population,sep="/")
      audit[[key]] <- data.frame(c, lineage=lin,population=population,control_donors=as.integer(n[1]),
                                  case_donors=as.integer(n[2]),status=reason)
      if (reason!="tested") next
      design <- model.matrix(~group)
      d <- DGEList(counts=counts[,selected,drop=FALSE],group=group)
      keep <- filterByExpr(d,design=design)
      if (!any(keep)) {
        audit[[key]]$status <- "no_features_pass_filter"
        next
      }
      d <- calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method="TMM")
      d <- estimateDisp(d,design,robust=TRUE)
      fit <- glmQLFit(d,design,robust=TRUE)
      test <- glmQLFTest(fit,coef="groupCase")
      z <- add_ql_wald_ci(fit, test, coef="groupCase")
      j <- match(rownames(z),features$feature_id)
      stopifnot(!anyNA(j))
      meta <- audit[[key]][rep(1,nrow(z)),setdiff(names(audit[[key]]),"status"),drop=FALSE]
      result[[key]] <- cbind(meta,data.frame(feature_id=features$feature_id[j],
                             symbol=features$symbol[j],cyp_gene=features$cyp_gene[j],z))
    }
  }
}
write.table(do.call(rbind,audit),file.path(destination,"Lineage_Eligibility.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
stopifnot(length(result)>0)
z <- do.call(rbind,result)
con <- gzfile(file.path(destination,"Genomewide_edgeR.tsv.gz"),"wt")
write.table(z,con,sep="\t",row.names=FALSE,quote=FALSE)
close(con)
cyp <- z[nzchar(z$cyp_gene),]
stopifnot(!anyDuplicated(cyp[c("contrast","lineage","population","cyp_gene")]))
write.table(cyp,file.path(destination,"CYP_edgeR.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines(paste("Finished",nrow(contrasts),"contrasts;",nrow(cyp),"CYP tests including sensitivity"),file.path(destination,"COMPLETED.txt"))

capture.output(sessionInfo(), file=file.path(destination, 'sessionInfo.txt'))
