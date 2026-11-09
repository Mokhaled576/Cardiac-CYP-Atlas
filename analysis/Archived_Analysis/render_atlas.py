"""Render HeartMap figures from frozen result tables, without refitting models."""

import argparse
import html
from pathlib import Path
import re

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import pandas as pd


INK = "#172b40"
TEAL = "#007a79"
BLUE = "#295a94"
ORANGE = "#b84c35"
GREY = "#e6ebef"
PALETTE = list(plt.get_cmap("tab20").colors) + list(plt.get_cmap("Dark2").colors) + list(plt.get_cmap("Set1").colors)
plt.rcParams.update({"font.size": 10, "axes.titlesize": 12, "axes.spines.top": False,
                     "axes.spines.right": False, "pdf.fonttype": 42, "text.color": INK,
                     "axes.labelcolor": INK, "axes.titlecolor": INK})


def read(path, **kwargs):
    return pd.read_csv(path, sep="\t", **kwargs)


def natural(value):
    return [int(token) if token.isdigit() else token.lower()
            for token in re.split(r"(\d+)", str(value))]


def clean(value):
    return str(value).replace("_", " ")


def slug(value):
    return re.sub(r"[^A-Za-z0-9_-]+", "_", str(value)).strip("_")


def wide(data, row, col, value):
    return data.pivot(index=row, columns=col, values=value)


class Atlas:
    def __init__(self, out, name):
        self.out = out
        self.out.mkdir(parents=True, exist_ok=True)
        self.pdf = PdfPages(out / f"{name}.pdf")
        self.records = []
        self.name = name

    def add(self, fig, category, key, title, caption, source):
        number = len(self.records) + 1
        identifier = f"{self.name.split('_')[0]}_{number:03d}_{key}"
        folder = self.out / category
        folder.mkdir(exist_ok=True)
        fig.suptitle(title, x=.06, y=.985, ha="left", fontsize=17, fontweight="bold")
        fig.supxlabel(caption, x=.06, y=.008, ha="left", fontsize=8, color="#40556c", wrap=True)
        fig.tight_layout(rect=(.02, .055, .99, .94))
        fig.savefig(folder / f"{identifier}.png", dpi=145, facecolor="white")
        self.pdf.savefig(fig, dpi=180)
        plt.close(fig)
        self.records.append({"figure": identifier, "section": category, "title": title,
                             "caption": caption, "source": str(source),
                             "image": f"{category}/{identifier}.png"})

    def finish(self):
        self.pdf.close()
        pd.DataFrame(self.records).to_csv(self.out / "Figure_Index.tsv", sep="\t", index=False)
        entries = []
        for row in self.records:
            fields = {key: html.escape(value, quote=True) for key, value in row.items()}
            entries.append('<article><h2><a href="{image}">{title}</a></h2>'
                           '<p>{figure} · {caption}</p><a href="{image}">'
                           '<img loading="lazy" src="{image}" alt="{title}"></a></article>'.format(**fields))
        self.out.joinpath("Open_Figure_Gallery.html").write_text(
            '<!doctype html><meta charset="utf-8"><title>HeartMap Figure Atlas</title>'
            '<style>body{font:16px system-ui;margin:3vw;color:#172b40}'
            'main{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:24px}'
            'article{border:1px solid #ddd;padding:12px}h2{font-size:18px}img{width:100%}</style>'
            f'<h1>{clean(self.name)}</h1><p><a href="{self.name}.pdf">Download figure book (PDF)</a>'
            ' · <a href="Figure_Index.tsv">Figure index and source paths</a></p><main>'
            + "".join(entries) + "</main>", encoding="utf-8")


def heat(ax, matrix, label, *, cmap="viridis", vmin=0, vmax=None, diverging=False, significant=None):
    values = matrix.to_numpy(float)
    palette = plt.get_cmap(cmap).copy()
    palette.set_bad(GREY)
    if diverging:
        limit = vmax or max(float(np.nanquantile(np.abs(values), .98)), .01)
        settings = {"norm": TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit)}
    else:
        settings = {"vmin": vmin, "vmax": vmax}
    image = ax.imshow(np.ma.masked_invalid(values), aspect="auto", cmap=palette,
                      interpolation="none", **settings)
    ax.set_xticks(range(len(matrix.columns)), [clean(x) for x in matrix.columns],
                  rotation=35, ha="right", fontsize=8)
    ax.set_yticks(range(len(matrix.index)), [clean(x) for x in matrix.index], fontsize=8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if significant is not None:
        rows, cols = np.where(significant.reindex_like(matrix).eq(True).to_numpy())
        ax.scatter(cols, rows, color=INK, edgecolor="white", linewidth=.35, s=12)
    ax.figure.colorbar(image, ax=ax, fraction=.035, pad=.025, label=label)


def points(ax, data, column, *, cluster_ids=False):
    groups = sorted(data[column].fillna("Unspecified").astype(str).unique(), key=natural)
    for i, group in enumerate(groups):
        selected = data[data[column].fillna("Unspecified").astype(str).eq(group)]
        color = PALETTE[i % len(PALETTE)]
        ax.scatter(selected.umap_1, selected.umap_2, s=3, lw=0, alpha=.55,
                   color=color, label=clean(group), rasterized=True)
        if cluster_ids:
            ax.text(selected.umap_1.median(), selected.umap_2.median(), group,
                    fontsize=8, ha="center", va="center", fontweight="bold",
                    bbox={"facecolor": "white", "alpha": .78, "edgecolor": "none", "pad": 1})
    if not cluster_ids:
        ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), markerscale=2.2,
                  frameon=False, fontsize=8)
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", aspect="equal")
    ax.set_xticks([])
    ax.set_yticks([])


def clustering(root, out):
    atlas = Atlas(out, "S1_Clustering_and_Annotation")
    paths = sorted((root / "qc").glob("FreshCluster_*"))
    assert len(paths) == 11, f"Expected 11 clustering checkpoints, found {len(paths)}"
    overview = []
    for folder in paths:
        study, job = folder.name.removeprefix("FreshCluster_").rsplit("_", 1)
        status = (folder / "COMPLETED.txt").read_text().strip()
        cat = read(folder / "Cluster_Annotation_Candidates.tsv", dtype={"cluster": str})
        donors = read(folder / "Cluster_Donor_Coverage.tsv", dtype={"cluster": str})
        markers = read(folder / "Cluster_Markers.tsv", dtype={"cluster": str})
        evidence = read(folder / "Lineage_Marker_Evidence.tsv", dtype={"cluster": str})
        cells = read(folder / "Nucleus_Clusters.tsv.gz", dtype={"cluster": str}, usecols=[
            "disease", "cardiac_region", "cluster", "umap_1", "umap_2"])
        assert len(cells) == cat.nuclei.sum(), folder
        assert cells.groupby("cluster").size().sort_index().equals(
            cat.set_index("cluster").nuclei.sort_index()), folder
        valid = cells.dropna(subset=["umap_1", "umap_2"])
        sample = pd.concat([g.sample(min(len(g), 1400), random_state=314159)
                            for _, g in valid.groupby("cluster")], ignore_index=True)
        display = clean(study)
        order = sorted(cat.cluster, key=natural)
        title = f"{display} | job {job}"
        overview.append({"study": study, "cluster_job": job, "nuclei": len(cells),
                         "clusters": len(cat), "plotted_nuclei": len(sample), "checkpoint": status})
        sample["candidate"] = sample.cluster.map(cat.set_index("cluster").marker_candidate)
        source = Path("qc") / folder.name
        fig, ax = plt.subplots(figsize=(10, 7))
        points(ax, sample, "cluster", cluster_ids=True)
        atlas.add(fig, study, f"{study}_clusters", f"{display} · De novo cluster UMAP",
                  f"{title}. {len(cells):,} assigned nuclei, {len(cat)} Louvain clusters. "
                  f"Display sample: {len(sample):,} nuclei, capped per cluster.", source)
        fig, ax = plt.subplots(figsize=(12, 7))
        points(ax, sample, "candidate")
        atlas.add(fig, study, f"{study}_annotation", f"{display} · Provisional annotation UMAP",
                  f"{title}. Colors show marker-based candidates. "
                  "Review marker evidence and donor support before accepting an assignment.", source)
        fig, axes = plt.subplots(1, 2, figsize=(15, 6))
        for ax, column, label in zip(axes, ["disease", "cardiac_region"], ["Source condition", "Source region"]):
            points(ax, sample, column)
            ax.set_title(label)
        atlas.add(fig, study, f"{study}_sample_context", f"{display} · UMAP by condition and source region",
                  f"{title}. Uses deposited metadata from the clustering checkpoint. "
                  "See the regional crosswalk for reconciled disease models.", source)

        top = markers.dropna(subset=["symbol"]).sort_values("rank").groupby("cluster").head(3)
        symbols = list(dict.fromkeys(top.sort_values("cluster", key=lambda s: s.astype(int)).symbol))
        matrix = markers.pivot_table(index="cluster", columns="symbol", values="logmean_delta", aggfunc="max")
        fig, ax = plt.subplots(figsize=(max(12, len(symbols) * .22), 7.5))
        heat(ax, matrix.reindex(index=order, columns=symbols), "Marker log-mean difference", cmap="magma")
        atlas.add(fig, study, f"{study}_marker_genes", f"{display} · Leading genes across de novo clusters",
                  f"{title}. Up to three leading marker symbols per cluster. Gray indicates an unreported "
                  "cluster-gene combination. See the full ranked marker table.", source)
        matrix = evidence.pivot(index="cluster", columns="lineage", values="score").reindex(order)
        fig, ax = plt.subplots(figsize=(12, 7.5))
        heat(ax, matrix, "Marker panel score", cmap="viridis")
        assigned = cat.set_index("cluster").marker_candidate
        ax.set_yticklabels([f"{c}  {clean(assigned.loc[c])}" for c in matrix.index], fontsize=8)
        atlas.add(fig, study, f"{study}_lineage_evidence", f"{display} · Evidence for marker-based annotations",
                  f"{title}. Each row shows the provisional candidate and scores across competing lineage panels.", source)
        by_donor = donors.groupby(["cluster", "donor_id"]).nuclei.sum()
        max_share = by_donor.groupby(level=0).max() / by_donor.groupby(level=0).sum()
        ordered = cat.set_index("cluster").loc[order]
        fig, axes = plt.subplots(1, 3, figsize=(13, 7.5), sharey=True)
        for ax, series, label, color in zip(
            axes, [ordered.nuclei, ordered.donors, max_share.reindex(order) * 100],
            ["Assigned nuclei (log scale)", "Contributing donors", "Largest donor share (%)"],
            [BLUE, TEAL, ORANGE],
        ):
            ax.barh(np.arange(len(order)), series, color=color)
            ax.set_xlabel(label)
        axes[0].set_xscale("log")
        axes[0].set_yticks(np.arange(len(order)), [f"{c}  {clean(assigned.loc[c])}" for c in order], fontsize=8)
        axes[0].invert_yaxis()
        axes[1].axvline(3, color=INK, ls="--", lw=1)
        axes[2].axvline(50, color=INK, ls="--", lw=1)
        atlas.add(fig, study, f"{study}_donor_support", f"{display} · Cluster size and donor representation",
                  f"{title}. Full assignment counts. Dashed lines show 3 donors and 50% donor dominance "
                  "for reference to the operational lineage screen.", source)
        print(f"{study}: six clustering and annotation figures", flush=True)
    pd.DataFrame(overview).to_csv(out / "Cohort_Checkpoints.tsv", sep="\t", index=False)
    atlas.finish()


def expression(root, out):
    atlas = Atlas(out, "S2_CYP_Expression")
    source = Path("results/CYP_Landscape_2958286/CYP_Study_Landscape.tsv")
    data = read(root / source)
    genes = sorted(data.gene.unique(), key=natural)
    assert len(genes) == 57 and data.study.nunique() == 11
    for (study, region, disease), subset in data.groupby(["study", "cardiac_region", "disease"]):
        copy = subset.copy()
        copy["stratum"] = copy.lineage + " (n=" + copy.donors_with_20_nuclei.astype(str) + ")"
        copy["expression"] = np.log10(1 + copy.median_donor_cpm)
        fig, axes = plt.subplots(1, 2, figsize=(max(14, 10 + .48 * copy.lineage.nunique()), 14))
        for ax, value, label, cmap, max_value in [
            (axes[0], "expression", "log10(1 + median donor CPM)", "magma", 3),
            (axes[1], "donor_fraction_detected", "Fraction of donors with CYP RNA", "viridis", 1),
        ]:
            matrix = wide(copy, "gene", "stratum", value).reindex(genes)
            heat(ax, matrix, label, cmap=cmap, vmax=max_value)
            ax.set_yticklabels(["CYP5A1 / TBXAS1" if x == "CYP5A1" else
                                "CYP8A1 / PTGIS" if x == "CYP8A1" else x for x in genes], fontsize=8)
        key = slug(f"{study}_{region}_{disease}")
        atlas.add(fig, study, key, f"{clean(study)} · {clean(region)} · {clean(disease)}",
                  "All 57 CYP genes. n is donors with ≥20 nuclei per lineage. Gray denotes unavailable estimates. "
                  "CPM colors saturate at log10(1 + CPM) = 3. Source condition names are retained.", source)
    atlas.finish()


def model_heat(atlas, subset, contrast, category, row, effect, q, label, title, source,
               *, row_order=None, limit=3, caption=""):
    matrix = wide(subset, row, "lineage", effect)
    if row_order is not None:
        matrix = matrix.reindex(row_order)
    significance = wide(subset.assign(hit=subset[q].lt(.05)), row, "lineage", "hit")
    fig, ax = plt.subplots(figsize=(max(10, 6.5 + .45 * len(matrix.columns)), max(6, 3 + .19 * len(matrix))))
    heat(ax, matrix, label, cmap="RdBu_r", vmax=limit, diverging=True, significant=significance)
    atlas.add(fig, category, f"{slug(contrast)}_{category}", f"{clean(contrast)} · {title}",
              f"{caption} Dots mark tier-wide BH q < 0.05. Gray = unavailable or untested. "
              f"Color saturates beyond ±{limit:g}.", source)


def disease_signaling(root, out):
    atlas = Atlas(out, "S3_Disease_and_Signaling")
    d = Path("results/CYP_Disease_2962292")
    s = Path("results/Signaling_Activity_2967196")
    sources = {"CYP_Disease_Effects": d / "CYP_Disease_Effects.tsv",
               "GO_Enrichment": d / "CYP_Disease_Pathways.tsv",
               "Signaling_Footprints": s / "Signaling_Activity_Effects.tsv"}
    cyp = read(root / sources["CYP_Disease_Effects"]).query("population == 'all_donors'")
    go = read(root / sources["GO_Enrichment"])
    sig = read(root / sources["Signaling_Footprints"]).query("population == 'all_donors'")
    genes = sorted(read(root / "results/CYP_Landscape_2958286/CYP_Study_Landscape.tsv").gene.unique(), key=natural)
    contrasts = sorted(cyp.contrast.unique())
    assert len(contrasts) == 13
    assert set(contrasts) == set(go.contrast) == set(sig.contrast)
    for contrast in contrasts:
        c = cyp[cyp.contrast.eq(contrast)]
        row = c.iloc[0]
        caption = f"{clean(row.study)} | {row.region} | {row.case} vs {row.control} | {row.tier} tier."
        model_heat(atlas, c, contrast, "CYP_Disease_Effects", "cyp_gene", "logFC", "CYP_BH_across_tier",
                   "Donor pseudobulk log2 fold change", "CYP disease effects", sources["CYP_Disease_Effects"],
                   row_order=genes, limit=4, caption=caption)
        model_heat(atlas, go[go.contrast.eq(contrast)], contrast, "GO_Enrichment", "pathway", "NES",
                   "GO_BH_across_tier", "Normalized enrichment score", "Ranked GO enrichment",
                   sources["GO_Enrichment"], limit=2.5, caption=caption)
        model_heat(atlas, sig[sig.contrast.eq(contrast)], contrast, "Signaling_Footprints", "pathway", "score_delta",
                   "FDR_BH", "Case–control difference (score SD)", "Inferred PROGENy signaling",
                   sources["Signaling_Footprints"], limit=2.5, caption=caption)

    for (region, lineage), subset in cyp.groupby(["region", "lineage"]):
        if subset.contrast.nunique() < 2:
            continue
        matrix = wide(subset, "cyp_gene", "contrast", "logFC").reindex(genes)
        significance = wide(subset.assign(hit=subset.CYP_BH_across_tier.lt(.05)), "cyp_gene", "contrast", "hit")
        fig, ax = plt.subplots(figsize=(max(11, 7 + .46 * len(matrix.columns)), 14))
        heat(ax, matrix, "Donor pseudobulk log2 fold change", cmap="RdBu_r", vmax=4, diverging=True,
             significant=significance)
        atlas.add(fig, "Region_Matched_Comparisons", slug(f"{region}_{lineage}"),
                  f"{region} · {clean(lineage)} · CYP effects across matched contrasts",
                  "Columns are within-study effects. Comparison across columns is descriptive. "
                  "Direct disease-vs-disease effects are identified by the contrast name. "
                  "Dots mark tier-wide BH q < 0.05. Gray = untested. Colors saturate at ±4.",
                  sources["CYP_Disease_Effects"])
    atlas.finish()


def effects_forest(atlas, data, category, key, title, caption, source, effect, lower, upper,
                   row_label, x_label, q="FDR_BH"):
    rows = data.dropna(subset=[effect, lower, upper]).sort_values(q).head(18).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(12, max(4.5, len(rows) * .33 + 2.2)))
    y = np.arange(len(rows))
    colors = [TEAL if value < .05 else "#788b9a" for value in rows[q]]
    ax.hlines(y, rows[lower], rows[upper], color=colors, lw=1.5)
    ax.scatter(rows[effect], y, c=colors, s=35, zorder=3)
    ax.axvline(0, color=INK, ls="--", lw=1)
    ax.set_yticks(y, rows[row_label], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel(x_label)
    atlas.add(fig, category, slug(key), title, caption, source)


def programs_validation(root, out):
    atlas = Atlas(out, "S4_Programs_Communication_and_Validation")
    full = Path("results/Full_Analysis_2965831")
    regulatory = full / "Regulatory_State_Effects.tsv.gz"
    reg = read(root / regulatory).query("population == 'all_donors'")
    for (contrast, domain), subset in reg.groupby(["contrast", "domain"]):
        subset = subset.copy()
        category = "TF_Regulons" if domain == "TF" else "Cell_State_Programs"
        if domain == "TF":
            selected = subset.groupby("program").FDR_BH.min().nsmallest(18).index
            subset = subset[subset.program.isin(selected)]
            subset["program_label"] = subset.program
            note = "Up to 18 regulons ranked by the smallest q across lineages. Full results are included."
        else:
            subset["program_label"] = subset.term.fillna(subset.program)
            note = "Prespecified cell-state programs with estimated scores."
        model_heat(atlas, subset, contrast, category, "program_label", "score_delta", "FDR_BH",
                   "Case–control score difference", clean(category), regulatory,
                   limit=1, caption=note)

    occupancy = full / "Cluster_Occupancy_Effects.tsv"
    clusters = read(root / occupancy).query("population == 'all_donors'").copy()
    clusters["label"] = clusters.lineage.map(clean) + " · cluster " + clusters.cluster.astype(str)
    for contrast, subset in clusters.groupby("contrast"):
        effects_forest(atlas, subset, "Cluster_Occupancy", contrast,
                       f"{clean(contrast)} · Relative cluster occupancy",
                       "Up to 18 clusters ranked by q. Teal indicates q < 0.05. "
                       "Effects are relative to the modeled lineage and cluster IDs are study-specific.",
                       occupancy, "log2_odds_delta", "ci_low", "ci_high", "label",
                       "Case–control log2 odds difference (95% CI)")

    communication = full / "Top_Ligand_Receptor.tsv"
    pairs = read(root / communication).query("population == 'all_donors'").copy()
    pairs["label"] = (pairs.sender.map(clean) + " → " + pairs.receiver.map(clean) + " · "
                      + pairs.ligand + " / " + pairs.receptor)
    for contrast, subset in pairs.groupby("contrast"):
        effects_forest(atlas, subset, "Communication", contrast,
                       f"{clean(contrast)} · Leading ligand–receptor scores",
                       "Up to 18 pairs from the project’s ranked export. Teal indicates q < 0.05. "
                       "RNA co-expression supports a predicted interaction.",
                       communication, "score_delta", "ci_low", "ci_high", "label",
                       "Case–control pair score difference (95% CI)")

    focal = Path("results/Donor_CYP_Links_2969410/Permutation_Validated/Prespecified_DCM_Associations.tsv")
    for gene, subset in read(root / focal).groupby("gene"):
        fig, ax = plt.subplots(figsize=(11, 6))
        y = np.arange(len(subset))
        labels = [f"{r.study} · {clean(r.mode)} (n={r.n_donors})" for r in subset.itertuples()]
        ax.scatter(subset.rho, y, s=50, color=TEAL)
        for i, row in enumerate(subset.itertuples()):
            ax.text(1.02, i, f"q={row.Focal_FDR_BH:.2g}", fontsize=8, va="center")
        ax.set(xlim=(-1.05, 1.26), xlabel="Spearman or covariate-adjusted partial rank correlation")
        ax.set_yticks(y, labels, fontsize=9)
        ax.invert_yaxis()
        ax.axvline(0, color=INK, ls="--", lw=1)
        atlas.add(fig, "Donor_Associations", gene, f"{gene} · Inferred estrogen score within LV DCM",
                  "Cardiomyocyte donor-level results. All-donor and available covariate sensitivity "
                  "models are shown separately.", focal)

    leuven = Path("external_validation/GSE298023/Validation_3016861/Leuven_DCM_CYP_Validation.tsv")
    validation = read(root / leuven)
    validation["label"] = validation.gene + " · " + validation.model.map(clean)
    effects_forest(atlas, validation, "External_Validation", "Leuven",
                   "Leuven · Published single-nucleus DCM validation",
                   "External case–control effects across available models. Teal indicates q < 0.05.",
                   leuven, "log2FC", "CI_low", "CI_high", "label", "log2 fold change (95% CI)")
    bulk = Path("external_validation/GSE116250/Bulk_DCM_External_Screen.tsv")
    effects_forest(atlas, read(root / bulk), "External_Validation", "Bulk_LV",
                   "GSE116250 · Published bulk LV DCM screen",
                   "Bulk LV RNA provides tissue-level evidence. Teal indicates q < 0.05.",
                   bulk, "log2_RPKM_plus_0_1_delta", "ci95_low", "ci95_high", "gene",
                   "Case–control difference in log2(RPKM + 0.1) (95% CI)")
    atlas.finish()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("stage", choices=["1", "2", "3", "4"])
    args = parser.parse_args()
    {"1": clustering, "2": expression, "3": disease_signaling, "4": programs_validation}[args.stage](
        args.root, args.out)
