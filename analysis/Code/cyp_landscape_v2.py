"""Summarize original CYP counts by study, condition, region, lineage, and donor."""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

STUDIES = [
    ("AS", "2941738", "2941738"), ("Chaffin", "2941480", "2941480"),
    ("DoxCM", "2941739", "2941739"), ("Kanemaru", "2941448", "2941448"),
    ("Koenig", "2941737", "2941737"), ("Kuppe", "2941482", "2941482"),
    ("Litvinukova", "2941778", "2941778"), ("Reichart", "2941481", "2941481"),
    ("Simonson_ARVC", "2937872", "2938548"), ("Simonson_ICM", "2937871", "2938547"),
    ("Tucker", "2941483", "2941483"),
]

CYP = """CYP1A1 CYP1A2 CYP1B1 CYP2A6 CYP2A7 CYP2A13 CYP2B6 CYP2C8 CYP2C9 CYP2C18 CYP2C19
CYP2D6 CYP2E1 CYP2F1 CYP2J2 CYP2R1 CYP2S1 CYP2U1 CYP2W1 CYP3A4 CYP3A5 CYP3A7 CYP3A43
CYP4A11 CYP4A22 CYP4B1 CYP4F2 CYP4F3 CYP4F8 CYP4F11 CYP4F12 CYP4F22 CYP4V2 CYP4X1 CYP4Z1
CYP5A1 CYP7A1 CYP7B1 CYP8A1 CYP8B1 CYP11A1 CYP11B1 CYP11B2 CYP17A1 CYP19A1 CYP20A1 CYP21A2
CYP24A1 CYP26A1 CYP26B1 CYP26C1 CYP27A1 CYP27B1 CYP27C1 CYP39A1 CYP46A1 CYP51A1""".split()
assert len(CYP) == len(set(CYP)) == 57

study, prep_job, cluster_job = STUDIES[int(sys.argv[1])]
destination = Path(sys.argv[2])
destination.mkdir(parents=True, exist_ok=False)
qc = Path("qc")
cluster_dir = qc / f"FreshCluster_{study}_{cluster_job}"
assert (cluster_dir / "COMPLETED.txt").is_file()
adata = ad.read_h5ad(qc / f"ClusteringInput_{study}_{prep_job}" / "Source_Counts.h5ad")
assert adata.obs_names.is_unique and adata.var_names.is_unique and adata.X.dtype.kind in "iu"
clusters = pd.read_csv(cluster_dir / "Nucleus_Clusters.tsv.gz", sep="\t",
                       converters={0: str}, dtype=str, keep_default_na=False)
clusters = clusters.set_index(clusters.columns[0], verify_integrity=True)
source_rows = adata.obs_names.get_indexer(clusters.index)
if (source_rows < 0).any():
    raise ValueError(f"{study}: unmatched cluster barcodes")
source_metadata = adata.obs.iloc[source_rows]
fields = ["donor_id", "biosample_id", "cardiac_region", "disease"]
assert source_metadata.provisional_qc_flag.all()
for field in fields:
    assert np.array_equal(clusters[field].to_numpy(), source_metadata[field].astype(str).to_numpy()), field

audit = Path("results/Marker_Donor_Review_2957389")
review = pd.read_csv(audit / "Cluster_Marker_Donor_Review.tsv", sep="\t", dtype={"cluster":str})
review = review.loc[review.study.eq(study)].set_index("cluster", verify_integrity=True)
evidence = pd.read_csv(cluster_dir / "Lineage_Marker_Evidence.tsv", sep="\t", dtype={"cluster":str})
evidence = evidence.sort_values(["cluster", "score"], ascending=[True, False]).groupby("cluster").head(2)
best = evidence.drop_duplicates("cluster").set_index("cluster")
runner_up = evidence.loc[evidence.duplicated("cluster")].set_index("cluster")
review["top_panel_markers"] = best.enriched_markers
review["panel_score_margin"] = best.score - runner_up.score
review["analysis_lineage"] = review.marker_candidate.where(
    review.marker_candidate.ne("unresolved") & review.marker_candidate.eq(best.lineage)
    & review.top_panel_markers.ge(3) & review.panel_score_margin.ge(0.10)
    & review.donors.ge(3) & review.dominant_donor_fraction.le(0.50), "unresolved"
)
review["review_status"] = np.where(review.analysis_lineage.eq("unresolved"),
                                   "held from lineage estimates", "marker and donor support screen passed")
review.to_csv(destination / "Cluster_Landscape_Review.tsv", sep="\t")
clusters["lineage"] = clusters.cluster.map(review.analysis_lineage)
assert clusters.lineage.notna().all()
clusters["library_umi"] = source_metadata.total_counts.to_numpy(dtype=np.int64)

if study in {"Koenig", "DoxCM"}:
    anatomy = pd.read_csv(audit / "Published_Cohort_Anatomy.tsv", sep="\t",
                          dtype={"biosample_id":str,"prep_job":str})
    anatomy = anatomy.loc[anatomy.study.eq(study) & anatomy.prep_job.eq(prep_job)].set_index("biosample_id", verify_integrity=True)
    assert clusters.biosample_id.isin(anatomy.index).all()
    assert clusters.biosample_id.map(anatomy.condition).equals(clusters.disease)
    clusters["cardiac_region"] = clusters.biosample_id.map(anatomy.cardiac_region)
elif study == "Kanemaru":
    clusters["cardiac_region"] = clusters.cardiac_region.replace(
        {"atrioventriculeft cardiac atriumr node": "atrioventricular node"}
    )

symbols = pd.Index(adata.var.get("feature_name", adata.var_names)).astype(str).str.strip().str.upper()
approved = {"CYP5A1": "TBXAS1", "CYP8A1": "PTGIS"}
ensembl = {"CYP5A1": "ENSG00000059377", "CYP8A1": "ENSG00000124212"}
stable_ids = pd.Index(adata.var.get("gene_ids", adata.var_names)).astype(str).str.split(".").str[0]
positions = {}
for gene in CYP:
    matched = symbols.isin([gene, approved.get(gene, gene)])
    if gene in ensembl:
        matched |= stable_ids == ensembl[gene]
    positions[gene] = np.flatnonzero(matched)
catalog = pd.DataFrame({
    "gene": CYP,
    "approved_symbol": [approved.get(gene, gene) for gene in CYP],
    "matched_source_symbols": [";".join(symbols[positions[g]]) for g in CYP],
    "feature_status": ["unique" if len(positions[g]) == 1 else "absent" if len(positions[g]) == 0 else "ambiguous" for g in CYP],
    "source_features": [";".join(adata.var_names[positions[g]]) for g in CYP],
}).assign(study=study)
catalog.to_csv(destination / "CYP_Source_Coverage.tsv", sep="\t", index=False)
measured = catalog.loc[catalog.feature_status.eq("unique"), "gene"].tolist()
matrix = adata.X[:, [positions[g][0] for g in measured]][source_rows].tocsr().astype(np.int64)
if matrix.data.size and matrix.data.min() < 0:
    raise ValueError("Negative raw count")
matrix.eliminate_zeros()
group_cols = ["cardiac_region", "disease", "lineage", "donor_id"]
keys = pd.MultiIndex.from_frame(clusters[group_cols])
codes, groups = pd.factorize(keys, sort=True)
metadata = groups.to_frame(index=False)
metadata.columns = group_cols
metadata["nuclei"] = np.bincount(codes)
operator = csr_matrix((np.ones(len(codes), dtype=np.int64), (codes, np.arange(len(codes)))),
                      shape=(len(groups), len(codes)))
metadata["library_umi"] = np.asarray(operator @ clusters.library_umi.to_numpy()).ravel()
assert metadata.library_umi.gt(0).all()
metadata["study"] = study
metadata["donor_eligible"] = metadata.nuclei.ge(20) & metadata.lineage.ne("unresolved")

totals = (operator @ matrix).toarray()
matrix.data.fill(1)
positive = (operator @ matrix).toarray()
assert (positive <= metadata.nuclei.to_numpy()[:, None]).all()
long = metadata.loc[metadata.index.repeat(len(measured))].reset_index(drop=True)
long["gene"] = np.tile(measured, len(metadata))
long["gene_umi"] = totals.ravel()
long["nuclei_detected"] = positive.ravel()
long["cpm"] = 1e6 * long.gene_umi / long.library_umi
long["nuclei_fraction_detected"] = long.nuclei_detected / long.nuclei
long.to_csv(destination / "CYP_Donor_Landscape.tsv.gz", sep="\t", index=False, compression="gzip")

strata = ["study", "cardiac_region", "disease", "lineage"]
donor_coverage = metadata.groupby(strata, as_index=False).agg(
    donors=("donor_id","size"), donors_with_20_nuclei=("donor_eligible","sum"), nuclei=("nuclei","sum"))
donor_coverage.to_csv(destination / "Landscape_Strata_Coverage.tsv", sep="\t", index=False)
eligible = long.loc[long.donor_eligible].copy()
eligible["detected"] = eligible.gene_umi.gt(0)
summary = eligible.groupby(strata+["gene"], as_index=False).agg(
    donors_evaluated=("donor_id","size"), donors_detected=("detected","sum"),
    median_donor_cpm=("cpm","median"), median_nuclei_fraction_detected=("nuclei_fraction_detected","median"),
    total_gene_umi=("gene_umi","sum"))
summary["donor_fraction_detected"] = summary.donors_detected / summary.donors_evaluated
summary = donor_coverage.loc[donor_coverage.lineage.ne("unresolved")].merge(
    catalog.drop(columns="study"), how="cross").merge(summary, on=strata+["gene"], how="left", validate="one_to_one")
summary.to_csv(destination / "CYP_Study_Landscape.tsv", sep="\t", index=False, na_rep="NA")
(destination / "COMPLETED.txt").write_text(
    f"{study}: {len(clusters)} clustered nuclei, {len(metadata)} donor strata, {len(measured)}/57 unique CYP features.\n")
print((destination / "COMPLETED.txt").read_text(), flush=True)
