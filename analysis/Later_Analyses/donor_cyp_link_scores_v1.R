args <- commandArgs(TRUE)
source <- args[[1]]
out <- args[[2]]
dir.create(out,recursive=TRUE,showWarnings=FALSE)
suppressPackageStartupMessages({library(edgeR);library(org.Hs.eg.db);library(GO.db)})
read_source <- function(f) read.delim(file.path(source,f),check.names=FALSE,colClasses="character")
obs <- read_source("Pseudobulk_Donors.tsv")
features <- read_source("Pseudobulk_Features.tsv")
contrasts <- read_source("Contrasts.tsv")
study <- unique(contrasts$study)
obs$nuclei <- as.integer(obs$nuclei)
counts <- as.matrix(read.delim(gzfile(file.path(source,"Pseudobulk_Counts.tsv.gz")),row.names=1,check.names=FALSE))
stopifnot(length(study)==1,identical(rownames(counts),features$feature_id),identical(colnames(counts),obs$sample),
          !anyDuplicated(obs[c("cardiac_region","lineage","donor_id")]))
obs$study <- study
obs$age <- NA_real_
obs$sex <- NA_character_
cov <- file.path("results/Covariate_CYP_2965628",basename(source),"Source_Donor_Covariates.tsv")
if(file.exists(cov)) {
  meta <- read.delim(cov,check.names=FALSE,colClasses="character")
  ix <- match(obs$sample,meta$sample)
  stopifnot(!anyNA(ix),identical(obs$donor_id,meta$donor_id[ix]),identical(obs$disease,meta$disease[ix]))
  obs$age <- as.numeric(meta$age[ix])
  obs$sex <- meta$sex[ix]
}

ids <- c("GO:0006690","GO:0019369","GO:0006954","GO:0006979","GO:0006915","GO:0030199","GO:0034340",
         "GO:0001666","GO:0001525","GO:0042775","GO:0030968","GO:0090398","GO:0006635","GO:0030595","GO:0006936")
terms <- unname(AnnotationDbi::Term(GOTERM[ids]))
stopifnot(!anyNA(terms))
go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db,keys=ids,keytype="GOALL",columns="SYMBOL"))
go <- unique(go[!is.na(go$SYMBOL),c("GOALL","SYMBOL")])
load("source_metadata/Regulatory_Resources/PROGENy_model_human_full.rda")
net <- model_human_full[order(model_human_full$pathway,model_human_full$p.value,model_human_full$gene),]
net <- net[!duplicated(net[c("pathway","gene")]),]
net <- do.call(rbind,lapply(split(net,net$pathway),head,100))
paths <- rbind(data.frame(key=paste0("GO:",go$GOALL),symbol=go$SYMBOL,weight=1),
               data.frame(key=paste0("PROGENy:",net$pathway),symbol=net$gene,weight=net$weight))
reference_cyp <- unique(read.delim("results/CYP_Landscape_2958286/CYP_Source_Coverage.tsv",check.names=FALSE)$approved_symbol)
exclude <- unique(c(reference_cyp,"TBXAS1","PTGIS"))
paths <- paths[!paths$symbol %in% exclude & !grepl("^CYP[0-9]",paths$symbol),]
write.table(paths,file.path(out,"Non_CYP_Program_Targets.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(data.frame(id=ids,term=terms),file.path(out,"GO_Program_Names.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
feature_cyp <- !is.na(features$cyp_gene) & nzchar(features$cyp_gene)
stopifnot(!anyDuplicated(features$cyp_gene[feature_cyp]))
scores <- list(); gene_scores <- list(); diagnostics <- list()
for(k in unique(paste(obs$cardiac_region,obs$lineage,sep="|"))) {
  idx <- which(paste(obs$cardiac_region,obs$lineage,sep="|")==k)
  m <- obs[idx,]
  if(nrow(m)<8) next
  d <- calcNormFactors(DGEList(counts=counts[,idx,drop=FALSE]),method="TMM")
  e <- cpm(d,log=TRUE,prior.count=2)
  linear <- cpm(d)
  present <- rowSums(linear>=1)>=3
  valid <- !is.na(features$symbol) & nzchar(features$symbol) & !duplicated(features$symbol) &
           !duplicated(features$symbol,fromLast=TRUE) & present & !feature_cyp &
           !features$symbol %in% exclude & !grepl("^CYP[0-9]",features$symbol)
  sd <- apply(e,1,stats::sd)
  valid <- valid & is.finite(sd) & sd>0
  standard <- (e[valid,,drop=FALSE]-rowMeans(e[valid,,drop=FALSE]))/sd[valid]
  hit <- match(paths$symbol,features$symbol[valid])
  hit <- split(which(!is.na(hit)),paths$key[!is.na(hit)])
  for(key in names(hit)) {
    ids <- hit[[key]]
    weighted <- startsWith(key,"PROGENy:")
    if(length(ids)<if(weighted) 20 else 10) next
    values <- if(weighted) e[valid,,drop=FALSE] else standard
    target <- match(paths$symbol[ids],features$symbol[valid])
    score <- colSums(values[target,,drop=FALSE]*paths$weight[ids])/sum(abs(paths$weight[ids]))
    if(stats::sd(score)==0) next
    scores[[length(scores)+1L]] <- data.frame(sample=m$sample,program=key,score=score)
    diagnostics[[length(diagnostics)+1L]] <- data.frame(study=study,region=m$cardiac_region[1],lineage=m$lineage[1],
        program=key,n_targets=length(ids))
  }
  for(g in which(feature_cyp)) {
    gene_scores[[length(gene_scores)+1L]] <- data.frame(sample=m$sample,gene=features$cyp_gene[g],
        log_CPM=e[g,],CPM=linear[g,])
  }
}
stopifnot(length(scores)>0,length(gene_scores)>0)
write.table(obs,file.path(out,"Donor_Metadata.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,scores),file.path(out,"CYP_Free_Program_Scores.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,gene_scores),file.path(out,"Donor_CYP_Expression.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,diagnostics),file.path(out,"Program_Coverage.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
writeLines("Verified score inputs ready.",file.path(out,"SCORES_COMPLETE.txt"))
