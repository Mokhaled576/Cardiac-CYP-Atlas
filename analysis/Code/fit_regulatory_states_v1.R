args <- commandArgs(TRUE)
path <- args[1]
out <- args[2]
dir.create(out, recursive=TRUE)
suppressPackageStartupMessages({library(edgeR); library(limma); library(org.Hs.eg.db); library(GO.db)})
read_meta <- function(file) read.delim(file.path(path,file),check.names=FALSE,colClasses="character")
obs <- read_meta("Pseudobulk_Donors.tsv")
feature <- read_meta("Pseudobulk_Features.tsv")
contrast <- read_meta("Contrasts.tsv")
counts <- as.matrix(read.delim(gzfile(file.path(path,"Pseudobulk_Counts.tsv.gz")),row.names=1,check.names=FALSE))
stopifnot(identical(rownames(counts),feature$feature_id),identical(colnames(counts),obs$sample))
tf <- read.csv("source_metadata/Regulatory_Resources/CollecTRI_regulons.csv",check.names=FALSE)
tf <- unique(tf[!is.na(tf$source) & nzchar(tf$source) & !is.na(tf$target) & nzchar(tf$target) &
                tf$weight %in% c(-1,1),c("source","target","weight")])
conflict <- duplicated(tf[1:2]) | duplicated(tf[1:2],fromLast=TRUE)
tf <- tf[!conflict,]
stopifnot(nrow(tf)>10000,length(unique(tf$source))>500)
state_ids <- c("GO:0030199","GO:0034340","GO:0001666","GO:0001525","GO:0042775","GO:0030968",
               "GO:0090398","GO:0006635","GO:0030595","GO:0006936","GO:0006954","GO:0006979","GO:0006915")
terms <- unname(AnnotationDbi::Term(GOTERM[state_ids]))
stopifnot(!anyNA(terms),length(terms)==length(state_ids))
mapping <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,keys=state_ids,keytype="GOALL",columns="SYMBOL"))
mapping <- unique(mapping[!is.na(mapping$SYMBOL),c("GOALL","SYMBOL")])
state <- data.frame(source=mapping$GOALL,target=mapping$SYMBOL,weight=1)
tf$domain <- "TF"
state$domain <- "state_program"
net <- rbind(tf,state)
net$key <- paste(net$domain,net$source,sep=":")
write.table(data.frame(GO_ID=state_ids,term=terms,target_count=as.integer(table(factor(state$source,levels=state_ids)))),
            file.path(out,"State_Program_Provenance.tsv"),row.names=FALSE,sep="\t",quote=FALSE)
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
      sd <- apply(e,1,stats::sd)
      good <- !is.na(symbols) & nzchar(symbols) & !duplicated(symbols) &
              !duplicated(symbols,fromLast=TRUE) & is.finite(sd) & sd>0
      z <- (e[good,,drop=FALSE]-rowMeans(e[good,,drop=FALSE]))/sd[good]
      idx <- match(net$target,symbols[good])
      hit <- split(which(!is.na(idx)),net$key[!is.na(idx)])
      hit <- hit[lengths(hit)>=10]
      if(!length(hit)) {audit[[k]]$status<-"no_eligible_sets"; next}
      scores <- t(vapply(hit,function(j) colMeans(z[idx[j],,drop=FALSE]*net$weight[j]),numeric(nrow(m))))
      fit <- eBayes(lmFit(scores,design),robust=TRUE)
      tab <- topTable(fit,coef="groupCase",number=Inf,sort.by="none",confint=TRUE)
      set <- rownames(tab)
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
write.table(x,gzfile(file.path(out,"Regulatory_State_Effects.tsv.gz")),row.names=FALSE,sep="\t",quote=FALSE)
write.table(do.call(rbind,scores_out),gzfile(file.path(out,"Donor_Program_Scores.tsv.gz")),row.names=FALSE,sep="\t",quote=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines(paste(nrow(x),"modeled TF and program tests"),file.path(out,"COMPLETED.txt"))
