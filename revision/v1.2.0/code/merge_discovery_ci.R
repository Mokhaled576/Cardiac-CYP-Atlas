# Attach CI columns only after confirming the refit reproduces archived outputs.
args <- commandArgs(TRUE)
if (length(args) < 3L) {
  stop("Usage: Rscript merge_discovery_ci.R ARCHIVED_EFFECTS.tsv NEW_OUTPUT.tsv REFIT_CYP.tsv ...")
}
if (file.exists(args[2])) stop("Output already exists")
read_table <- function(path) read.delim(path, check.names = FALSE, colClasses = "character")
archived <- read_table(args[1])
refit <- do.call(rbind, lapply(args[-c(1, 2)], read_table))
archived <- archived[archived$population == "all_donors", , drop = FALSE]
refit <- refit[refit$population == "all_donors", , drop = FALSE]
key_fields <- c("contrast", "lineage", "population", "feature_id", "cyp_gene")
key <- function(x) do.call(paste, c(x[key_fields], sep = "\t"))
stopifnot(nrow(archived) == 1234L, !anyDuplicated(key(archived)),
          !anyDuplicated(key(refit)), setequal(key(archived), key(refit)))
refit <- refit[match(key(archived), key(refit)), , drop = FALSE]
for (name in c("case_donors", "control_donors", "logFC", "F", "PValue")) {
  a <- as.numeric(archived[[name]])
  b <- as.numeric(refit[[name]])
  tolerance <- 1e-12 + 1e-7 * pmax(abs(a), abs(b))
  if (any(!is.finite(a) | !is.finite(b) | abs(a - b) > tolerance)) {
    stop("Refit differs from the archived results in ", name,
         ". Check the edgeR environment and inputs before combining.")
  }
}
ci_columns <- c("logFC_unshrunk", "SE", "CI_95_low", "CI_95_high", "QL_df_total", "CI_status")
merged <- cbind(archived, refit[ci_columns])
# Preserve the original genome-wide FDR separately from the manuscript's tier-wide CYP FDR.
merged$FDR_edgeR_genomewide <- merged$FDR
merged$FDR <- merged$CYP_BH_across_tier
write.table(merged, args[2], sep = "\t", row.names = FALSE, quote = FALSE)
message("Merged ", nrow(merged), " discovery tests. Unavailable intervals: ",
        sum(merged$CI_status == "unavailable"))
