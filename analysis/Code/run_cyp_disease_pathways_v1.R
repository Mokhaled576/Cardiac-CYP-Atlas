args <- commandArgs(TRUE)
path <- args[1]
study <- args[2]
suppressPackageStartupMessages({
  library(GO.db); library(org.Hs.eg.db); library(fgsea)
})
ids <- c("GO:0006690", "GO:0019369", "GO:0006954", "GO:0006979", "GO:0006915")
families <- unname(AnnotationDbi::Term(GOTERM[ids]))
expected <- c("icosanoid metabolic process", "arachidonate metabolic process",
              "inflammatory response", "response to oxidative stress", "apoptotic process")
if (!identical(sub("^eicosanoid", "icosanoid", families), expected))
  stop("Installed GO.db terms differ from the expected GO IDs")
map <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys=ids,
             keytype="GOALL", columns=c("SYMBOL")))
sets <- lapply(ids,function(id) unique(na.omit(map$SYMBOL[map$GOALL==id])))
names(sets) <- families
write.table(data.frame(GO_ID=ids,term=families,symbol_count=lengths(sets)),
            file.path(path,"Pathway_Gene_Set_Provenance.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
z <- read.delim(gzfile(file.path(path,"Genomewide_edgeR.tsv.gz")),check.names=FALSE)
z <- z[z$population=="all_donors",]
rankings <- split(z,list(z$contrast,z$lineage),drop=TRUE)
out <- list()
set.seed(314159)
for (key in names(rankings)) {
  t <- rankings[[key]]
  t <- t[!is.na(t$symbol) & nzchar(t$symbol) & is.finite(t$F) & is.finite(t$logFC),]
  duplicate <- duplicated(t$symbol) | duplicated(t$symbol,fromLast=TRUE)
  t <- t[!duplicate,]
  rank <- setNames(sign(t$logFC)*sqrt(t$F),t$symbol)
  if (length(rank)<100) next
  f <- fgseaMultilevel(sets,rank,minSize=10,maxSize=2500,eps=1e-10,nproc=1)
  if (!nrow(f)) next
  f$leadingEdge <- vapply(f$leadingEdge,paste,collapse=";",FUN.VALUE=character(1))
  f$GO_ID <- ids[match(f$pathway,families)]
  f$contrast <- t$contrast[1]
  f$study <- t$study[1]
  f$region <- t$region[1]
  f$lineage <- t$lineage[1]
  f$tier <- t$tier[1]
  f$ranked_features <- length(rank)
  out[[key]] <- as.data.frame(f)
}
stopifnot(length(out)>0)
write.table(do.call(rbind,out),file.path(path,"CYP_Disease_Pathways.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
