import importlib.util
import sys
from pathlib import Path
import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import MiniBatchKMeans

out = Path(sys.argv[1])
out.mkdir(parents=True, exist_ok=False)
cohorts = ["Chaffin","Reichart","Kuppe","Tucker","Litvinukova","Kanemaru","AS","DoxCM","Koenig","Simonson_ICM","Simonson_ARVC"]
existing = {"Simonson_ICM":"2937871","Simonson_ARVC":"2937872"}
qc = Path("qc")
rng = np.random.default_rng(314159)
inputs, records = {}, []
audit = Path("results/Marker_Donor_Review_2957389")
review = pd.read_csv("results/CYP_Landscape_2958286/Cluster_Landscape_Review.tsv",sep="\t",dtype=str)
overlap = pd.read_csv(audit/"Shared_Identifier_Screen.tsv",sep="\t",dtype=str)
anatomy = pd.read_csv(audit/"Published_Cohort_Anatomy.tsv",sep="\t",dtype=str)
pair = overlap.study_a.eq("Litvinukova") & overlap.study_b.eq("Reichart")
shared = set(overlap.loc[pair & overlap.field.eq("biosample_id"),"shared_identifier"])
library_links = {}
for study in ("Litvinukova","Reichart"):
    job, = review.loc[review.study.eq(study),"cluster_job"].unique()
    source = qc/f"FreshCluster_{study}_{job}"/"Nucleus_Clusters.tsv.gz"
    sample = pd.read_csv(source,sep="\t",usecols=["donor_id","biosample_id"],dtype=str,keep_default_na=False).drop_duplicates()
    sample["library"] = sample.biosample_id.str.removesuffix("_premrna")
    sample = sample.loc[sample.library.isin(shared),["library","donor_id"]].drop_duplicates()
    assert set(sample.library)==shared and sample.library.is_unique
    library_links[study] = sample.set_index("library").donor_id
links = pd.concat(library_links,axis=1).reset_index().sort_values("library")
donor_links = links[["Litvinukova","Reichart"]].drop_duplicates()
assert donor_links.Litvinukova.is_unique and donor_links.Reichart.is_unique
links.assign(retained_source="Litvinukova",excluded_source="Reichart",
             basis="matching HCA biosample ID after removing Reichart _premrna suffix").to_csv(
                 out/"Verified_Library_Crosswalk.tsv",sep="\t",index=False)
def initial_clusters(data, n_clusters):
    return MiniBatchKMeans(n_clusters=n_clusters,batch_size=2048,n_init=3,max_iter=75,
                           random_state=314159).fit(data).cluster_centers_
for name in cohorts:
    paths = sorted(qc.glob(f"FreshCluster_{name}_*"), key=lambda p: int(p.name.rsplit("_",1)[-1]))
    finished = [path for path in paths if (path / "COMPLETED.txt").is_file()]
    if not finished:
        records.append((name, "missing completed clustering", 0, 0, 0))
        continue
    folder = finished[-1]
    run = folder.name.rsplit("_",1)[-1]
    source = qc / f"ClusteringInput_{name}_{existing.get(name,run)}" / "Source_Counts.h5ad"
    if not source.is_file():
        records.append((name, "missing matching source-count object", 0, 0, 0))
        continue
    meta = pd.read_csv(folder / "Nucleus_Clusters.tsv.gz", sep="\t",
                       converters={0: str}, dtype=str, keep_default_na=False)
    meta = meta.set_index(meta.columns[0], verify_integrity=True)
    if (meta.index == "").any():
        raise ValueError(f"{name}: empty nucleus barcode")
    labels = review.loc[review.study.eq(name) & review.cluster_job.eq(run)].set_index("cluster",verify_integrity=True)
    assert len(labels)
    meta["analysis_lineage"] = meta.cluster.map(labels.analysis_lineage)
    assert meta.analysis_lineage.notna().all()
    initial = len(meta)
    if name == "Reichart":
        meta = meta.loc[~meta.biosample_id.str.removesuffix("_premrna").isin(shared)].copy()
    removed_library = initial-len(meta)
    meta["reference_donor"] = name + "::" + meta.donor_id
    if name=="Reichart":
        remap = donor_links.set_index("Reichart").Litvinukova.map(lambda x:"Litvinukova::"+x)
        meta["reference_donor"] = meta.donor_id.map(remap).fillna(meta.reference_donor)
    row = overlap.loc[overlap.study_a.eq(name)|overlap.study_b.eq(name)]
    meta["candidate_overlap"] = meta.donor_id.isin(row.loc[row.field.eq("donor_id"),"shared_identifier"])
    meta["reviewed_region"] = meta.cardiac_region
    if name in {"DoxCM","Koenig"}:
        mapped = anatomy.loc[anatomy.study.eq(name)].set_index("biosample_id",verify_integrity=True)
        assert meta.biosample_id.map(mapped.condition).equals(meta.disease)
        meta["reviewed_region"] = meta.biosample_id.map(mapped.cardiac_region)
    meta["reviewed_region"] = meta.reviewed_region.replace(
        {"atrioventriculeft cardiac atriumr node":"atrioventricular node"})
    held = int(meta.analysis_lineage.eq("unresolved").sum())
    meta = meta.loc[meta.analysis_lineage.ne("unresolved")].copy()
    keep = []
    for key, donor in meta.groupby("reference_donor", observed=True):
        limit = 60 if name in library_links and key.startswith("Litvinukova::") and key.split("::",1)[1] in set(donor_links.Litvinukova) else 120
        keep.extend(rng.choice(donor.index.to_numpy(),size=min(limit,len(donor)),replace=False))
    b = ad.read_h5ad(source, backed="r")
    missing = meta.index[~meta.index.isin(b.obs_names)]
    if len(missing):
        raise ValueError(f"{name}: {len(missing)} cluster barcodes missing from {source}; examples {missing[:5].tolist()}")
    index = np.sort(b.obs_names.get_indexer(keep))
    raw = b[index, :].to_memory()
    b.file.close()
    symbols = pd.Index(raw.var["feature_name"] if "feature_name" in raw.var else raw.var_names).astype(str).str.upper()
    good = ~symbols.duplicated(keep=False) & symbols.str.fullmatch(r"[A-Z][A-Z0-9.-]*")
    raw = raw[:,good].copy()
    raw.var_names = symbols[good]
    raw.obs["source_barcode"] = raw.obs_names.astype(str)
    for key in ("analysis_lineage","reference_donor","reviewed_region","candidate_overlap"):
        raw.obs[key] = meta.loc[raw.obs_names,key].to_numpy()
    raw.obs_names = pd.Index(name + "::" + raw.obs_names.astype(str))
    raw.obs = raw.obs[["source_barcode","donor_id","biosample_id","cardiac_region","reviewed_region","disease","analysis_lineage","reference_donor","candidate_overlap"]].copy()
    inputs[name] = raw
    records.append((name,"included",raw.n_obs,removed_library,held))
pd.DataFrame(records,columns=["study","status","sampled_nuclei","verified_duplicate_library_nuclei_removed","unresolved_lineage_nuclei_held"]).to_csv(out/"Cohort_Coverage.tsv",sep="\t",index=False)
if len(inputs) != len(cohorts):
    raise ValueError("All 11 cohorts must have matching source counts and completed clustering")
joint = ad.concat(inputs, join="inner", label="study", index_unique=None)
assert joint.var_names.is_unique and joint.n_vars >= 5000
joint.X = joint.X.astype(np.float32)
sc.pp.normalize_total(joint,target_sum=1e4)
sc.pp.log1p(joint)
sc.pp.highly_variable_genes(joint,n_top_genes=2500,flavor="seurat",batch_key="study")
b = joint[:,joint.var.highly_variable.to_numpy()].copy()
pd.Series(b.var_names,name="feature_symbol").to_csv(out/"Reference_HVGs.tsv",sep="\t",index=False)
sc.pp.pca(b,n_comps=35,random_state=314159)
coordinates = "X_pca"
method = "Unadjusted shared PCA"
if importlib.util.find_spec("harmonypy") is not None:
    try:
        sc.external.pp.harmony_integrate(b,key="study",nclust=40,cluster_fn=initial_clusters,random_state=314159)
        assert b.obsm["X_pca_harmony"].shape == b.obsm["X_pca"].shape
        coordinates = "X_pca_harmony"
        method = "Harmony by source study"
    except Exception as error:
        method += f"; Harmony unavailable for this run ({type(error).__name__})"
sc.pp.neighbors(b,use_rep=coordinates,n_neighbors=15,random_state=314159)
sc.tl.umap(b,random_state=314159)
if coordinates != "X_pca_harmony":
    raise RuntimeError(method)
embedding = joint.obs.copy()
assert embedding.groupby("reference_donor",observed=True).size().le(120).all()
embedding[["UMAP1","UMAP2"]] = b.obsm["X_umap"]
embedding.to_csv(out/"Reference_Embedding.tsv.gz",sep="\t",compression="gzip")
(out/"Reference_Status.txt").write_text(f"{method}; {joint.n_obs} nuclei, {joint.n_vars} shared symbols, {len(inputs)} cohorts.\n"
    "Verified duplicate HCA library identifiers excluded. Other donor matches remain flagged. Marker-supported broad lineages only. Reference coordinates are for visualization.\n")
fig, axes = plt.subplots(1,2,figsize=(13,6),layout="constrained")
for ax, key in zip(axes,("study","analysis_lineage")):
    for label, frame in embedding.groupby(key,observed=True):
        ax.scatter(frame.UMAP1,frame.UMAP2,s=1.2,label=str(label),rasterized=True)
    ax.set_title(key.replace("_"," ").capitalize())
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(fontsize=6,markerscale=3,loc="best")
fig.suptitle("Cardiac reference after verified library reconciliation | " + method,fontsize=11)
fig.savefig(out/"Reconciled_Reference_UMAP.pdf",dpi=180)
(out/"COMPLETED.txt").write_text("Integrated 11 cohorts after verified HCA library reconciliation.\n")
