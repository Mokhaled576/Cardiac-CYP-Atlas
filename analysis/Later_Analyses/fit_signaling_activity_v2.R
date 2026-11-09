args <- commandArgs(TRUE)
path <- args[1]
out <- args[2]
dir.create(out, recursive=TRUE)
suppressPackageStartupMessages({library(edgeR); library(limma); library(statmod)})
read_meta <- function(file) read.delim(file.path(path,file),check.names=FALSE,colClasses="character")
obs <- read_meta("Pseudobulk_Donors.tsv")
feature <- read_meta("Pseudobulk_Features.tsv")
contrast <- read_meta("Contrasts.tsv")
counts <- as.matrix(read.delim(gzfile(file.path(path,"Pseudobulk_Counts.tsv.gz")),row.names=1,check.names=FALSE))
stopifnot(identical(rownames(counts),feature$feature_id),identical(colnames(counts),obs$sample))
load("source_metadata/Regulatory_Resources/PROGENy_model_human_full.rda")
model <- model_human_full
model <- model[order(model$pathway,model$p.value,model$gene),]
model <- model[!duplicated(model[c("pathway","gene")]),]
model <- do.call(rbind,lapply(split(model,model$pathway),head,100))
stopifnot(length(unique(model$pathway))==14,all(table(model$pathway)==100),
          all(is.finite(model$weight)))
net <- data.frame(key=paste0("footprint:",model$pathway),target=model$gene,weight=model$weight)
write.table(model,file.path(out,"Footprint_Weights.tsv"),row.names=FALSE,sep="\t",quote=FALSE)
audit <- list()
res <- list()
scores_out <- list()
for(i in seq_len(nrow(contrast))) {
  c <- contrast[i,]
  for(lineage in sort(unique(obs$lineage[obs$cardiac_region==c$region]))) {
    within <- obs$cardiac_region==c$region & obs$lineage==lineage & obs$disease %in% c(c$case,c$control)
    for (population in c("all_donors","exclude_flagged_overlap")) {
      select <- if(population=="all_donors") within else within & obs$overlap_flag=="no"
      m <- obs[select,]
      group <- factor(ifelse(m$disease==c$case,"Case","Control"),levels=c("Control","Case"))
      n <- table(group)
      status <- if (any(n<3)) "fewer_than_3_per_group" else "tested"
      if(population!="all_donors" && identical(within,select)) status <- "identical_to_all_donors"
      k <- paste(c$contrast,lineage,population,sep="/")
      rec <- data.frame(c,lineage=lineage,population=population,case_donors=as.integer(n[2]),
                        control_donors=as.integer(n[1]),status=status)
      audit[[k]] <- rec
      if(status!="tested") next
      stopifnot(!anyDuplicated(m$donor_id))
      design <- model.matrix(~group)
      d <- DGEList(counts=counts[,select,drop=FALSE])
      keep <- filterByExpr(d,design=design)
      if(sum(keep)<100) {audit[[k]]$status<-"fewer_than_100_expressed_features"; next}
      d <- calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method="TMM")
      e <- cpm(d,log=TRUE,prior.count=2)
      symbols <- feature$symbol[match(rownames(e),feature$feature_id)]
      good <- !is.na(symbols) & nzchar(symbols) & !duplicated(symbols) & !duplicated(symbols,fromLast=TRUE)
      e <- e[good,,drop=FALSE]
      idx <- match(net$target,symbols[good])
      hit <- split(which(!is.na(idx)),net$key[!is.na(idx)])
      hit <- hit[lengths(hit)>=20]
      if(!length(hit)) {audit[[k]]$status<-"no_eligible_pathways"; next}
      scores <- t(vapply(hit,function(j) colSums(e[idx[j],,drop=FALSE]*net$weight[j])/sum(abs(net$weight[j])),numeric(nrow(m))))
      scale <- apply(scores,1,stats::sd)
      keep_scores <- is.finite(scale)&scale>0
      if(!any(keep_scores)) {audit[[k]]$status<-"no_variable_footprints";next}
      scores <- (scores[keep_scores,,drop=FALSE]-rowMeans(scores[keep_scores,,drop=FALSE]))/scale[keep_scores]
      fit <- eBayes(lmFit(scores,design),robust=TRUE)
      tab <- topTable(fit,coef="groupCase",number=Inf,sort.by="none",confint=TRUE,genelist=data.frame(program_key=rownames(scores)))
      set <- tab$program_key
      tab$program_key <- NULL
      stopifnot(all(set %in% names(hit)))
      meta <- rec[rep(1,nrow(tab)),setdiff(names(rec),"status"),drop=FALSE]
      res[[k]] <- cbind(meta,domain=sub(":.*","",set),program=sub("^[^:]+:","",set),
                         n_targets=lengths(hit[set]),tab)
      if(population=="all_donors") {
        scores_out[[k]] <- data.frame(contrast=c$contrast,lineage=lineage,
          domain=rep(sub(":.*","",rownames(scores)),times=nrow(m)),
          program=rep(sub("^[^:]+:","",rownames(scores)),times=nrow(m)),
          donor_id=rep(m$donor_id,each=nrow(scores)),disease=rep(m$disease,each=nrow(scores)),
          score=as.vector(scores))
      }
    }
  }
}
write.table(do.call(rbind,audit),file.path(out,"Model_Eligibility.tsv"),row.names=FALSE,sep="\t",quote=FALSE)
stopifnot(length(res)>0)
x <- do.call(rbind,res)
x$logFC <- NULL
x$score_delta <- x$CI.L + (x$CI.R-x$CI.L)/2
write.table(x,file.path(out,"Signaling_Activity_Effects.tsv"),row.names=FALSE,sep="\t",quote=FALSE)
write.table(do.call(rbind,scores_out),gzfile(file.path(out,"Donor_Signaling_Scores.tsv.gz")),row.names=FALSE,sep="\t",quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines(paste(nrow(x),"modeled signaling footprint tests"),file.path(out,"COMPLETED.txt"))
