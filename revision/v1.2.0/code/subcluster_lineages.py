"""Subcluster retained lineages inside the archived pseudobulk builder."""
import anndata as ad
import importlib.util
import numpy as np
import pandas as pd
import scanpy as sc


def assign_subclusters(adata, clusters, source_rows, destination):
    symbols = pd.Index(adata.var.get("feature_name", adata.var_names)).astype(str).str.upper()
    panels = {
        "cardiomyocyte": {"stress": ["NPPA", "NPPB"], "contractile": ["MYH6", "ATP2A2", "PLN"]},
        "fibroblast": {"activated": ["POSTN", "FAP", "ACTA2"], "resting": ["DCN", "LUM", "PDGFRA"]},
    }
    result = np.full(len(clusters), "unresolved", dtype=object)
    marker_tables = []
    for lineage, markers in panels.items():
        ix = np.flatnonzero(clusters.lineage.eq(lineage) & clusters.cardiac_region.eq("LV"))
        if not len(ix):
            continue
        profile = ad.AnnData(adata.X[source_rows[ix]].astype(np.float32))
        profile.obs["donor"] = clusters.iloc[ix].donor_id.to_numpy()
        sc.pp.normalize_total(profile, target_sum=10000)
        sc.pp.log1p(profile)
        sc.pp.highly_variable_genes(profile, n_top_genes=2500, batch_key="donor")
        embedding = profile[:, profile.var.highly_variable].copy()
        sc.pp.pca(embedding, n_comps=35, random_state=314159)
        sc.pp.neighbors(embedding, n_neighbors=15, n_pcs=35, random_state=314159)
        if all(importlib.util.find_spec(pkg) for pkg in ("igraph", "leidenalg")):
            sc.tl.leiden(embedding, resolution=0.6, random_state=314159)
            codes = embedding.obs.leiden.astype(str).to_numpy()
        else:
            import networkx as nx
            graph = nx.from_scipy_sparse_array(embedding.obsp["connectivities"])
            communities = nx.community.louvain_communities(graph, resolution=0.6, seed=314159)
            codes = np.empty(len(ix), dtype=int)
            for code, members in enumerate(sorted(communities, key=min)):
                codes[list(members)] = code
            codes = codes.astype(str)
        result[ix] = [f"{lineage}_cluster_{code}" for code in codes]
        for group in np.unique(codes):
            members = codes == group
            record = {"lineage": lineage, "subcluster": group, "nuclei": int(members.sum()),
                      "donors": profile.obs.loc[members, "donor"].nunique()}
            for marker_set in markers.values():
                for gene in marker_set:
                    index = np.flatnonzero(symbols == gene)
                    record[gene] = (float(profile.X[members, int(index[0])].mean())
                                    if len(index) == 1 else np.nan)
            marker_tables.append(record)
    if not marker_tables:
        raise ValueError("No retained focal lineage nuclei in LV")
    pd.DataFrame(marker_tables).to_csv(destination / "Subtype_Marker_Review.tsv", sep="\t", index=False)
    metadata = clusters[["donor_id", "disease", "cardiac_region"]].assign(subcluster=result)
    coverage = metadata.loc[metadata.subcluster.ne("unresolved")].groupby(
        ["cardiac_region", "disease", "subcluster", "donor_id"], observed=True
    ).size().rename("nuclei").reset_index()
    coverage["eligible_profile"] = coverage.nuclei >= 20
    coverage.groupby(["cardiac_region", "disease", "subcluster"], observed=True).agg(
        nuclei=("nuclei", "sum"), eligible_donors=("eligible_profile", "sum")
    ).to_csv(destination / "Subtype_Coverage.tsv", sep="\t")
    return result
