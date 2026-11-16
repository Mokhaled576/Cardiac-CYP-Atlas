args <- commandArgs(TRUE)
stopifnot(length(args) == 4L)

read_tsv <- function(f) read.delim(f, check.names=FALSE)
scores <- read_tsv(gzfile(args[1]))
scores <- subset(scores, domain=="footprint" & pathway %in% c("Estrogen","Androgen") &
                   lineage %in% c("cardiomyocyte","fibroblast"))
rows <- list()
for (path in Sys.glob(file.path(args[2],"Cohort_*","Source_Donor_Covariates.tsv"))) {
  cov <- read_tsv(path)
  cov <- subset(cov, toupper(sex_conflict)=="FALSE" & sex %in% c("female","male"))
  contrasts <- read_tsv(file.path(args[3],basename(dirname(path)),"Contrasts.tsv"))
  for (i in seq_len(nrow(contrasts))) {
    c <- contrasts[i,]
    for (lin in c("cardiomyocyte","fibroblast")) {
      source <- subset(cov, cardiac_region==c$region & lineage==lin)
      stopifnot(!anyDuplicated(source$donor_id))
      for (pathway_name in c("Estrogen","Androgen")) {
        d <- subset(scores, contrast==c$contrast & lineage==lin & pathway==pathway_name)
        m <- match(d$donor_id,source$donor_id)
        stopifnot(all(d$disease[!is.na(m)]==source$disease[m[!is.na(m)]]))
        d$sex <- source$sex[m]
        d <- d[!is.na(d$sex) & is.finite(d$score),]
        for (mode in c("matched_unadjusted","sex_adjusted","female","male","interaction")) {
          dat <- if (mode %in% c("female","male")) subset(d,sex==mode) else d
          dat$group <- factor(ifelse(dat$disease==c$case,"Case","Control"),levels=c("Control","Case"))
          dat$sex <- factor(dat$sex,levels=c("female","male"))
          n <- table(dat$group)
          rec <- data.frame(c,lineage=lin,pathway=pathway_name,mode=mode,case_n=unname(n[2]),
                            control_n=unname(n[1]),estimate=NA_real_,PValue=NA_real_,status="ineligible")
          if (all(n>=3L) && (mode!="interaction" || all(table(dat$group,dat$sex)>=3L))) {
            form <- switch(mode,sex_adjusted=score~group+sex,interaction=score~group*sex,score~group)
            X <- model.matrix(form,dat)
            if (qr(X)$rank==ncol(X) && nrow(X)>ncol(X)) {
              # One response per model, so use ordinary donor-level linear inference.
              fit <- summary(lm(form,dat))$coefficients
              term <- if (mode=="interaction") "groupCase:sexmale" else "groupCase"
              rec$estimate <- fit[term,1]
              rec$PValue <- fit[term,4]
              rec$status <- "tested"
            }
          }
          rows[[length(rows)+1L]] <- rec
        }
      }
    }
  }
}
result <- do.call(rbind,rows)
result$FDR <- ave(result$PValue,result$mode,result$tier,FUN=function(p) p.adjust(p,"BH"))
write.table(result,args[4],sep="\t",row.names=FALSE,quote=FALSE)
