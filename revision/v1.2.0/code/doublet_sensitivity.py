"""Flag predicted doublets per library among archived QC-passing nuclei."""
import numpy as np
import pandas as pd
import scanpy as sc
from anndata import AnnData


def screen_doublets(adata, clusters, source_rows, destination):
    keep = np.ones(len(clusters), dtype=bool)
    audited = np.zeros(len(clusters), dtype=bool)
    summary = []
    for sample, indices in clusters.groupby("biosample_id", sort=True).indices.items():
        indices = np.asarray(indices)
        if len(indices) < 100:
            summary.append((sample, len(indices), 0, "too_few_nuclei"))
            continue
        library = AnnData(adata.X[source_rows[indices]].copy())
        sc.pp.scrublet(library, expected_doublet_rate=0.05, random_state=314159)
        flags = library.obs["predicted_doublet"].to_numpy(dtype=bool)
        threshold = library.uns["scrublet"].get("threshold")
        if threshold is None or not np.isfinite(threshold):
            summary.append((sample, len(indices), 0, "no_automatic_threshold"))
            continue
        keep[indices] = ~flags
        audited[indices] = True
        summary.append((sample, len(indices), int(flags.sum()), "screened"))
    if not audited.any():
        raise ValueError("No library produced a usable Scrublet threshold")
    pd.DataFrame(summary, columns=["biosample", "nuclei", "predicted_doublets", "status"]).to_csv(
        destination / "Doublet_Screen_Audit.tsv", sep="\t", index=False
    )
    return keep
