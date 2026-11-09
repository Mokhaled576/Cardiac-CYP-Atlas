args <- commandArgs(TRUE)
path <- args[1]
out <- args[2]
suppressPackageStartupMessages(library(edgeR))
obs <- read.delim(file.path(out,"Source_Donor_Covariates.tsv"),check.names=FALSE,colClasses="character")
features <- read.delim(file.path(path,"Pseudobulk_Features.tsv"),check.names=FALSE,colClasses="character")
contrasts <- read.delim(file.path(path,"Contrasts.tsv"),check.names=FALSE)
counts <- as.matrix(read.delim(gzfile(file.path(path,"Pseudobulk_Counts.tsv.gz")),row.names=1,check.names=FALSE))
stopifnot(identical(rownames(counts),features$feature_id),identical(colnames(counts),obs$sample))
obs$age <- as.numeric(obs$age)
results <- list()
audit <- list()
for(i in seq_len(nrow(contrasts))) {
  c <- contrasts[i,]
  formula <- if(c$study %in% c("Chaffin","Simonson_ICM")) ~group+age10+sex else ~group+sex
  fields <- if(c$study %in% c("Chaffin","Simonson_ICM")) c("age","sex") else "sex"
  for(lin in sort(unique(obs$lineage[obs$cardiac_region==c$region]))) {
    base <- obs$cardiac_region==c$region & obs$lineage==lin & obs$disease %in% c(c$case,c$control)
    select <- base & complete.cases(obs[,fields,drop=FALSE])
    m <- obs[select,]
    m$group <- factor(ifelse(m$disease==c$case,"Case","Control"),levels=c("Control","Case"))
    m$age10 <- (m$age-50)/10
    m$sex <- factor(m$sex)
    n <- table(m$group)
    rec <- data.frame(c,lineage=lin,case_donors=as.integer(n[2]),control_donors=as.integer(n[1]),
        covariates=paste(fields,collapse="+"),complete_case_donors=nrow(m),excluded_missing=sum(base & !select),
        status="tested")
    if(any(n<3)) rec$status <- "fewer_than_3_complete_donors_per_group"
    if(rec$status=="tested" && nlevels(m$sex)<2) rec$status <- "sex_invariant"
    if(rec$status=="tested") {
      adjusted <- model.matrix(formula,m)
      if(qr(adjusted)$rank<ncol(adjusted) || nrow(adjusted)-ncol(adjusted)<3)
        rec$status <- "rank_deficient_or_fewer_than_3_residual_df"
    }
    audit[[paste(c$contrast,lin)]] <- rec
    if(rec$status!="tested") next
    d <- DGEList(counts=counts[,select,drop=FALSE])
    keep <- filterByExpr(d,design=model.matrix(~group,m))
    d <- calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method="TMM")
    for(mode in c("complete_case_unadjusted","complete_case_adjusted")) {
      design <- if(mode=="complete_case_adjusted") adjusted else model.matrix(~group,m)
      test <- glmQLFTest(glmQLFit(estimateDisp(d,design,robust=TRUE),design,robust=TRUE),coef="groupCase")
      z <- topTags(test,n=Inf,sort.by="none")$table
      gene <- features$cyp_gene[match(rownames(z),features$feature_id)]
      valid <- !is.na(gene) & nzchar(gene)
      if(any(valid)) {
        t <- z[valid,,drop=FALSE]
        results[[paste(c$contrast,lin,mode)]] <- cbind(rec[rep(1,nrow(t)),],mode=mode,cyp_gene=gene[valid],t)
      }
    }
  }
}
write.table(do.call(rbind,audit),file.path(out,"Covariate_Model_Eligibility.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
if(length(results)) {
  write.table(do.call(rbind,results),file.path(out,"Covariate_CYP_Effects.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
}
writeLines(paste(length(results),"fitted model/lineage combinations"),file.path(out,"COMPLETED.txt"))
