args <- commandArgs(TRUE)
stopifnot(length(args) == 3L)
archive <- args[1]
covariates <- args[2]
out <- args[3]
suppressPackageStartupMessages(library(edgeR))
if (!dir.create(out)) stop("Output directory must be new")
read_tsv <- function(path) read.delim(path, check.names=FALSE, colClasses="character")
rows <- list()
audit <- list()
direct <- data.frame(
  contrast=c("Chaffin_DCM_vs_HCM_LV", "Reichart_DCM_vs_ARVC_LV", "Reichart_DCM_vs_ARVC_RV"),
  study=c("Chaffin", "Reichart", "Reichart"), region=c("LV","LV","RV"),
  case=c("DCM", rep("dilated cardiomyopathy",2)),
  control=c("HCM", rep("arrhythmogenic right ventricular cardiomyopathy",2)),
  tier="direct_etiology"
)
for (path in Sys.glob(file.path(archive,"Cohort_*"))) {
  obs <- read_tsv(file.path(path,"Pseudobulk_Donors.tsv"))
  features <- read_tsv(file.path(path,"Pseudobulk_Features.tsv"))
  contrasts <- read_tsv(file.path(path,"Contrasts.tsv"))
  counts <- as.matrix(read.delim(gzfile(file.path(path,"Pseudobulk_Counts.tsv.gz")),
                                 row.names=1,check.names=FALSE))
  stopifnot(identical(colnames(counts),obs$sample), identical(rownames(counts),features$feature_id))
  contrasts <- rbind(contrasts, direct[direct$study %in% contrasts$study,names(contrasts),drop=FALSE])
  contrasts$tier[contrasts$contrast=="DoxCM_vs_NICM_LV"] <- "direct_etiology"
  sex_path <- file.path(covariates,basename(path),"Source_Donor_Covariates.tsv")
  obs$sex <- NA_character_
  if (file.exists(sex_path)) {
    cov <- read_tsv(sex_path)
    stopifnot(!anyDuplicated(cov$sample))
    m <- match(obs$sample,cov$sample)
    verified <- !is.na(m) & cov$sex_conflict[m] %in% c("False","FALSE") &
      cov$sex[m] %in% c("female","male")
    stopifnot(all(obs$donor_id[verified]==cov$donor_id[m[verified]]))
    obs$sex[verified] <- cov$sex[m[verified]]
  }
  for (i in seq_len(nrow(contrasts))) {
    c <- contrasts[i,]
    for (lin in unique(obs$lineage[obs$cardiac_region==c$region])) {
      within <- obs$cardiac_region==c$region & obs$lineage==lin &
        obs$disease %in% c(c$case,c$control) & as.integer(obs$nuclei)>=20L
      modes <- c("complete_case_unadjusted","sex_adjusted")
      if (c$tier=="direct_etiology") modes <- c("all_donor_direct",modes)
      if (lin %in% c("fibroblast","cardiomyocyte")) modes <- c(modes,"female","male")
      for (mode in modes) {
        use <- within
        if (mode %in% c("complete_case_unadjusted","sex_adjusted")) use <- use & !is.na(obs$sex)
        if (mode %in% c("female","male")) use <- use & obs$sex %in% mode
        group <- factor(ifelse(obs$disease[use]==c$case,"Case","Control"),levels=c("Control","Case"))
        sex <- factor(obs$sex[use],levels=c("female","male"))
        n <- table(group)
        record <- data.frame(c,lineage=lin,mode=mode,case_donors=as.integer(n["Case"]),
                             control_donors=as.integer(n["Control"]),status="tested")
        key <- paste(c$contrast,lin,mode,sep="/")
        audit[[key]] <- record
        if (any(n<3L)) {
          audit[[key]]$status <- "fewer_than_3_per_arm"
          next
        }
        base_design <- model.matrix(~group)
        design <- if (mode=="sex_adjusted") model.matrix(~group+sex) else base_design
        if (qr(design)$rank<ncol(design) || nrow(design)-ncol(design)<3L) {
          audit[[key]]$status <- "insufficient_rank_or_residual_df"
          next
        }
        d <- DGEList(counts=counts[,use,drop=FALSE])
        keep <- filterByExpr(d,design=base_design)
        if (!any(keep)) {
          audit[[key]]$status <- "no_features_pass_filter"
          next
        }
        d <- calcNormFactors(d[keep,,keep.lib.sizes=FALSE],method="TMM")
        fit <- glmQLFit(estimateDisp(d,design,robust=TRUE),design,robust=TRUE)
        z <- topTags(glmQLFTest(fit,coef="groupCase"),n=Inf,sort.by="none")$table
        gene <- features$cyp_gene[match(rownames(z),features$feature_id)]
        valid <- !is.na(gene) & nzchar(gene)
        if (mode %in% c("female","male")) {
          valid <- valid & gene==if (lin=="fibroblast") "CYP19A1" else "CYP2J2"
        }
        if (!any(valid)) {
          audit[[key]]$status <- "no_target_passes_filter"
          next
        }
        rows[[key]] <- cbind(record[rep(1,sum(valid)),],cyp_gene=gene[valid],z[valid,,drop=FALSE])
      }
    }
  }
}
write.table(do.call(rbind,audit),file.path(out,"Eligibility.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
if (!length(rows)) stop("No CYP fits were eligible")
results <- do.call(rbind,rows)
results$FDR_genomewide <- results$FDR
results$FDR <- ave(results$PValue,results$tier,results$mode,FUN=function(p) p.adjust(p,"BH"))
write.table(results,file.path(out,"CYP_Effects.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
capture.output(sessionInfo(),file=file.path(out,"sessionInfo.txt"))
