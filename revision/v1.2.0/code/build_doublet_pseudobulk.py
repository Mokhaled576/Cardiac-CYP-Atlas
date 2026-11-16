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
review = pd.read_csv(Path("results/CYP_Landscape_2958286/Cluster_Landscape_Review.tsv"), sep="\t", dtype={"cluster":str})
review = review.loc[review.study.eq(study)].set_index("cluster", verify_integrity=True)
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

from doublet_sensitivity import screen_doublets
singlets = screen_doublets(adata, clusters, source_rows, destination)
source_rows = source_rows[singlets]
clusters = clusters.iloc[np.flatnonzero(singlets)].copy()

design_table = pd.read_csv("CYP_Disease_Contrasts_v1.tsv", sep="\t", dtype=str)
design_table = design_table.loc[design_table.study.eq(study)]
assert len(design_table) and design_table.contrast.is_unique
for row in design_table.itertuples(index=False):
    in_region = clusters.loc[clusters.cardiac_region.eq(row.region)]
    assert {row.case, row.control} <= set(in_region.disease), row
mask = np.zeros(len(clusters), dtype=bool)
for row in design_table.itertuples(index=False):
    mask |= (clusters.cardiac_region.eq(row.region) & clusters.disease.isin([row.case,row.control])).to_numpy()
mask &= clusters.lineage.ne("unresolved").to_numpy()
source_rows = source_rows[mask]
clusters = clusters.iloc[np.flatnonzero(mask)].copy()
assert len(clusters)
assert clusters.groupby("donor_id").disease.nunique().le(1).all()
fields = ["cardiac_region","disease","lineage","donor_id"]
keys = pd.MultiIndex.from_frame(clusters[fields])
codes, groups = pd.factorize(keys, sort=True)
samples = groups.to_frame(index=False)
samples.columns = fields
samples["nuclei"] = np.bincount(codes)
retained = samples.nuclei.ge(20).to_numpy()
selected = retained[codes]
operator = csr_matrix((np.ones(selected.sum(), dtype=np.int64),
                       (codes[selected], np.arange(selected.sum()))), shape=(len(groups),selected.sum()))
# Keep the full expression matrix for normalization, dispersion and downstream enrichment.
raw = adata.X[source_rows[selected]]
counts = (operator[retained] @ raw).toarray().T
samples = samples.loc[retained].copy()
assert counts.dtype.kind in "iu" and counts.min() >= 0 and counts.sum(axis=0).min() > 0
expected_total = np.asarray(operator[retained] @ clusters.loc[selected,"library_umi"].to_numpy()).ravel()
assert np.array_equal(counts.sum(axis=0),expected_total)
samples["sample"] = [f"PB{i:04d}" for i in range(len(samples))]
crosswalk = pd.read_csv(audit / "Shared_Identifier_Screen.tsv", sep="\t", dtype=str)
crosswalk = crosswalk.loc[crosswalk.study_a.eq(study) | crosswalk.study_b.eq(study)]
donor_ids = set(crosswalk.loc[crosswalk.field.eq("donor_id"),"shared_identifier"])
sample_ids = set(crosswalk.loc[crosswalk.field.eq("biosample_id"),"shared_identifier"])
donor_ids.update(clusters.loc[clusters.biosample_id.str.removesuffix("_premrna").isin(sample_ids),"donor_id"])
samples["overlap_flag"] = np.where(samples.donor_id.isin(donor_ids), "yes","no")
samples.to_csv(destination/"Pseudobulk_Donors.tsv",sep="\t",index=False)
design_table.to_csv(destination/"Contrasts.tsv",sep="\t",index=False)
symbols = pd.Index(adata.var.get("feature_name",adata.var_names)).astype(str).str.strip().str.upper()
features = pd.DataFrame({"feature_id":adata.var_names,"symbol":symbols})
catalog = pd.read_csv("results/CYP_Landscape_2958286/CYP_Source_Coverage.tsv",sep="\t",dtype=str)
catalog = catalog.loc[catalog.study.eq(study)]
assert len(catalog)==57
catalog.to_csv(destination/"CYP_Source_Coverage.tsv",sep="\t",index=False)
mapping = catalog.loc[catalog.feature_status.eq("unique")].set_index("source_features")["gene"]
assert mapping.index.is_unique and mapping.index.isin(features.feature_id).all()
features["cyp_gene"] = features.feature_id.map(mapping).fillna("")
features.to_csv(destination/"Pseudobulk_Features.tsv",sep="\t",index=False)
pd.DataFrame(counts,index=features.feature_id,columns=samples["sample"]).to_csv(
    destination/"Pseudobulk_Counts.tsv.gz",sep="\t",compression="gzip",index_label="feature_id")
print(f"{study}: {len(samples)} donor/lineage/region profiles and {len(design_table)} matched contrasts",flush=True)
