"""Cluster one cohort de novo and save evidence for broad annotation review."""
import importlib.util
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import csr_matrix

source, destination = map(Path, sys.argv[1:])
destination.mkdir(parents=True, exist_ok=False)
seed = 314159
a = ad.read_h5ad(source)
assert {"donor_id", "disease", "provisional_qc_flag"} <= set(a.obs)
assert not set(a.obsm)
assert not {"cell_type", "cluster", "leiden", "louvain"} & set(a.obs)
a = a[a.obs.provisional_qc_flag.to_numpy()].copy()
a.X = a.X.astype(np.float32)
sc.pp.normalize_total(a, target_sum=1e4)
sc.pp.log1p(a)
sc.pp.highly_variable_genes(a, n_top_genes=2500, flavor="seurat", batch_key="donor_id")
b = a[:, a.var.highly_variable.to_numpy()].copy()
sc.pp.pca(b, n_comps=35, random_state=seed)
sc.pp.neighbors(b, n_neighbors=15, n_pcs=35, random_state=seed)
if all(importlib.util.find_spec(pkg) for pkg in ("igraph", "leidenalg")):
    sc.tl.leiden(b, resolution=0.6, random_state=seed)
    cluster = b.obs["leiden"].astype(str).to_numpy()
    algorithm = "Leiden"
else:
    import networkx as nx
    graph = nx.from_scipy_sparse_array(b.obsp["connectivities"])
    communities = nx.community.louvain_communities(graph, resolution=0.6, seed=seed)
    cluster = np.full(a.n_obs, -1)
    for label, members in enumerate(sorted(communities, key=min)):
        cluster[list(members)] = label
    assert (cluster >= 0).all()
    cluster = cluster.astype(str)
    algorithm = "Louvain"
sc.tl.umap(b, random_state=seed)
a.obs["cluster"] = cluster
a.obs[["donor_id", "biosample_id", "disease", "cardiac_region", "cluster"]].assign(
    umap_1=b.obsm["X_umap"][:, 0], umap_2=b.obsm["X_umap"][:, 1]
).to_csv(destination / "Nucleus_Clusters.tsv.gz", sep="\t", compression="gzip")

groups, index = np.unique(cluster, return_inverse=True)
indicator = csr_matrix((np.ones(a.n_obs), (index, np.arange(a.n_obs))), shape=(len(groups), a.n_obs))
sizes = np.bincount(index)
sums = (indicator @ a.X).toarray()
means = sums / sizes[:, None]
rest = (sums.sum(axis=0) - sums) / (a.n_obs - sizes[:, None])
delta = means - rest
symbols = pd.Index(a.var["feature_name"] if "feature_name" in a.var else a.var_names).astype(str).str.upper()
panel = {
    "cardiomyocyte": ["TNNT2", "TNNI3", "MYH7", "RYR2", "ACTC1", "TTN"],
    "fibroblast": ["DCN", "LUM", "COL1A1", "COL3A1", "PDGFRA", "DPT"],
    "endothelial": ["PECAM1", "VWF", "CDH5", "KDR", "ERG", "CLDN5"],
    "lymphatic_endothelial": ["PROX1", "PDPN", "FLT4", "CCL21", "LYVE1"],
    "pericyte": ["RGS5", "PDGFRB", "CSPG4", "ABCC9", "KCNJ8"],
    "smooth_muscle": ["MYH11", "ACTA2", "TAGLN", "CNN1", "MYL9"],
    "myeloid": ["LST1", "TYROBP", "AIF1", "CSF1R", "FCER1G", "CD68"],
    "T_NK": ["CD3D", "CD3E", "TRAC", "NKG7", "KLRD1", "CD247"],
    "B_cell": ["MS4A1", "CD79A", "CD79B", "BANK1", "CD74"],
    "mast_cell": ["KIT", "TPSAB1", "TPSB2", "CPA3", "MS4A2"],
    "adipocyte": ["ADIPOQ", "PLIN1", "CIDEC", "LPL", "PPARG"],
    "neural": ["SOX10", "PLP1", "MPZ", "S100B", "GFRA3"],
    "mesothelial": ["WT1", "MSLN", "UPK3B", "KRT19", "ALDH1A2"],
}
markers = []
scores = []
for row, label in enumerate(groups):
    top = np.argsort(delta[row])[-40:][::-1]
    markers.extend((label, rank, a.var_names[g], symbols[g], delta[row, g], means[row, g])
                   for rank, g in enumerate(top, 1))
    for lineage, genes in panel.items():
        columns = np.flatnonzero(symbols.isin(genes) & ~symbols.duplicated(keep=False))
        if len(columns) >= 3:
            evidence = delta[row, columns]
            scores.append((label, lineage, len(columns), int((evidence > 0.1).sum()), np.mean(evidence)))
pd.DataFrame(markers, columns=["cluster", "rank", "feature_id", "symbol", "logmean_delta", "logmean"]).to_csv(
    destination / "Cluster_Markers.tsv", sep="\t", index=False)
evidence = pd.DataFrame(scores, columns=["cluster", "lineage", "available_markers", "enriched_markers", "score"])
evidence.to_csv(destination / "Lineage_Marker_Evidence.tsv", sep="\t", index=False)
coverage = a.obs.groupby("cluster", observed=True).agg(nuclei=("donor_id", "size"), donors=("donor_id", "nunique"))
candidates = []
for label, data in evidence.groupby("cluster"):
    best = data.sort_values("score", ascending=False).reset_index(drop=True)
    first = best.iloc[0]
    margin = first.score - best.iloc[1].score if len(best) > 1 else np.nan
    candidate = first.lineage if first.enriched_markers >= 2 and margin >= 0.05 else "unresolved"
    candidates.append((label, candidate, first.score, margin))
table = pd.DataFrame(candidates, columns=["cluster", "marker_candidate", "marker_score", "runner_up_margin"]).set_index("cluster")
coverage.join(table).to_csv(destination / "Cluster_Annotation_Candidates.tsv", sep="\t")
a.obs.groupby(["cluster", "donor_id", "disease"], observed=True).size().rename("nuclei").to_csv(
    destination / "Cluster_Donor_Coverage.tsv", sep="\t")
(destination / "COMPLETED.txt").write_text(
    f"{algorithm}; seed {seed}; resolution 0.6; {a.n_obs} nuclei; {len(groups)} clusters. Candidates require review.\n")
print((destination / "COMPLETED.txt").read_text(), flush=True)
