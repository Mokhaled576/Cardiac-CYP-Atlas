args <- commandArgs(TRUE)
stopifnot(length(args)==4)
suppressPackageStartupMessages({library(Matrix);library(SeuratObject);library(edgeR);library(limma)})
source <- args[1]
donor_column <- args[2]
lineage_column <- args[3]
out <- args[4]
dir.create(out,showWarnings=FALSE)
x <- readRDS(source)
m <- x[[]]
m$published_ident <- as.character(Idents(x)[rownames(m)])
stopifnot(all(c(donor_column,lineage_column)%in%names(m)))
canon <- function(z) toupper(gsub("[^A-Za-z0-9]","",trimws(as.character(z))))
meta <- read.delim("samples.csv",check.names=FALSE,colClasses="character")
meta <- meta[meta$region %in% c("DCM","NF") & meta$ventricle=="LV" & meta$location=="midmyocardium",]
meta$author_donor <- canon(meta$Sample)
meta$library <- canon(meta[["Sample ID"]])
meta$key <- paste(meta$library,trimws(meta$Gender),sep=":")
stopifnot(nrow(meta)>0,!anyDuplicated(meta$author_donor),!anyDuplicated(meta$key),all(nzchar(meta$author_donor)))
use <- as.character(m$region) %in% c("DCM","NF")
samples <- unique(data.frame(sample1=as.character(m[[donor_column]][use]),
    heart_number=as.character(m$heart_number[use]),region=as.character(m$region[use]),
    library=canon(m$library_id[use]),sex=as.character(m$sex[use]),batch=as.character(m$run_number[use])))
stopifnot(nrow(samples)==10,!anyDuplicated(samples$sample1),!anyDuplicated(samples$heart_number),
          sum(samples$region=="DCM")==6,sum(samples$region=="NF")==4)
link <- match(paste(samples$library,samples$sex,sep=":"),meta$key)
stopifnot(!anyNA(link),!anyDuplicated(link),identical(samples$region,meta$region[link]),
          identical(samples$batch,meta$Batch[link]))
meta <- meta[link,]
meta$donor <- samples$sample1
samples$author_donor <- meta$author_donor
stopifnot(length(unique(samples$library[samples$region=="NF"]))==2)
write.table(samples,file.path(out,"Leuven_Sample_Crosswalk.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
label <- read.delim("reviewed_external_lineages.tsv",colClasses="character",check.names=FALSE)
stopifnot(!anyDuplicated(label$published_label),all(label$lineage %in% c("cardiomyocyte","fibroblast")))
m$donor <- as.character(m[[donor_column]])
m$lineage <- label$lineage[match(as.character(m[[lineage_column]]),label$published_label)]
m$condition <- meta$region[match(m$donor,meta$donor)]
candidate <- use & !is.na(m$condition) & !is.na(m$lineage)
stopifnot(any(candidate),length(unique(m$lineage[candidate]))==2,
          all(as.character(m$region[candidate])==m$condition[candidate]))
a <- x[["RNA"]]
layers <- Layers(a)
stopifnot(sum(grepl("^counts($|[.])",layers))==1,"counts"%in%layers)
counts <- LayerData(a,layer="counts")
stopifnot(identical(rownames(m),colnames(counts)),!anyDuplicated(rownames(counts)),
    inherits(counts,"sparseMatrix"),all(is.finite(counts@x)),all(counts@x>=0),all(counts@x==floor(counts@x)))
results <- list(); profile <- list(); audit <- list()
targets <- c(cardiomyocyte="CYP2J2",fibroblast="CYP19A1")
for(lin in names(targets)) {
    gene <- targets[[lin]]
    use <- which(candidate & m$lineage==lin)
    donor <- factor(m$donor[use])
    bins <- Matrix::sparse.model.matrix(~0+donor)
    bulk <- as.matrix(counts[,use,drop=FALSE] %*% bins)
    colnames(bulk) <- levels(donor)
    donors <- meta[match(levels(donor),meta$donor),c("donor","region","Gender","Batch","library")]
    donors$n_nuclei <- as.integer(table(donor))
    donors$lineage <- lin
    donors$pass_nuclei <- donors$n_nuclei>=20
    profile[[lin]] <- donors
    n <- table(factor(donors$region[donors$pass_nuclei],levels=c("NF","DCM")))
    for (mode in c("unadjusted","batch_adjusted","shared_batch_only")) {
        shared <- intersect(donors$Batch[donors$region=="NF"&donors$pass_nuclei],
                            donors$Batch[donors$region=="DCM"&donors$pass_nuclei])
        selected <- donors$pass_nuclei
        if(mode=="shared_batch_only") selected <- selected & donors$Batch%in%shared
        d <- donors[selected,]
        d$group <- factor(d$region,levels=c("NF","DCM"))
        d$batch <- factor(d$Batch)
        n <- table(d$group)
        row <- data.frame(gene=gene,lineage=lin,model=mode,control_n=as.integer(n["NF"]),case_n=as.integer(n["DCM"]),
                          tested=FALSE,reason="",log2FC=NA_real_,CI_low=NA_real_,CI_high=NA_real_,p_value=NA_real_)
        reason <- ""
        if(any(n<3)) reason <- "fewer than three donors per group"
        if(reason=="" && mode=="batch_adjusted" && nlevels(d$batch)<2) reason <- "batch invariant"
        if(reason=="") {
            design <- model.matrix(if(mode=="batch_adjusted") ~group+batch else ~group,d)
            if(qr(design)$rank<ncol(design)||nrow(design)-ncol(design)<3) reason <- "model rank or residual degrees of freedom"
        }
        if(reason=="") {
            y <- DGEList(counts=bulk[,selected,drop=FALSE])
            keep <- filterByExpr(y,design=design)
            if(!gene%in%rownames(y)[keep]) reason <- "target absent or below expression filter"
        }
        if(reason=="") {
            y <- calcNormFactors(y[keep,,keep.lib.sizes=FALSE],method="TMM")
            v <- voom(y,design,plot=FALSE)
            fit <- eBayes(lmFit(v,design),robust=TRUE)
            est <- topTable(fit,coef="groupDCM",number=Inf,sort.by="none",confint=TRUE)[gene,]
            row[1,c("log2FC","CI_low","CI_high","p_value")] <- est[1,c("logFC","CI.L","CI.R","P.Value")]
            row$tested <- TRUE
        }
        row$reason <- reason
        audit[[length(audit)+1L]] <- row
    }
    pb <- data.frame(feature=rownames(bulk),bulk,check.names=FALSE)
    write.table(pb,gzfile(file.path(out,paste0(lin,"_Pseudobulk.tsv.gz"))),sep="\t",row.names=FALSE,quote=FALSE)
}
res <- do.call(rbind,audit)
res$FDR_BH <- ave(res$p_value,res$model,FUN=function(p) p.adjust(p,method="BH",n=2))
write.table(res,file.path(out,"Leuven_DCM_CYP_Validation.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,profile),file.path(out,"Leuven_Donor_Coverage.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
print(res,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(out,"sessionInfo.txt"))
writeLines("Donor-level target tests complete.",file.path(out,"COMPLETED.txt"))
