"""Combine complete donor-level CYP landscape outputs and plot descriptive recurrence."""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages

root = Path(sys.argv[1])
cohorts = [root / f"Cohort_{i}" for i in range(11)]
assert all((folder / "COMPLETED.txt").is_file() for folder in cohorts)

def combined(name):
    return pd.concat([pd.read_csv(folder / name, sep="\t", dtype={"donor_id":str}) for folder in cohorts],
                     ignore_index=True)

summary = combined("CYP_Study_Landscape.tsv")
coverage = combined("CYP_Source_Coverage.tsv")
strata = combined("Landscape_Strata_Coverage.tsv")
reviews = combined("Cluster_Landscape_Review.tsv")
assert coverage.study.nunique() == 11 and coverage.groupby("study").gene.nunique().eq(57).all()
assert not summary.duplicated(["study","cardiac_region","disease","lineage","gene"]).any()
assert reviews.shape[0] == 239 and reviews.study.nunique() == 11
summary.to_csv(root / "CYP_Study_Landscape.tsv", sep="\t", index=False, na_rep="NA")
coverage.to_csv(root / "CYP_Source_Coverage.tsv", sep="\t", index=False, na_rep="NA")
strata.to_csv(root / "Landscape_Strata_Coverage.tsv", sep="\t", index=False, na_rep="NA")
reviews.to_csv(root / "Cluster_Landscape_Review.tsv", sep="\t", index=False, na_rep="NA")
combined("CYP_Donor_Landscape.tsv.gz").to_csv(root / "CYP_Donor_Landscape.tsv.gz", sep="\t", index=False, compression="gzip")

controls = {"AS":"Control","Chaffin":"NF","DoxCM":"Donor","Kanemaru":"normal","Koenig":"Donor",
            "Litvinukova":"Control","Reichart":"normal","Simonson_ARVC":"PATO:0000461",
            "Simonson_ICM":"NF","Tucker":"Control"}
regions = {"heart left ventricle":"LV", "heart right ventricle":"RV",
           "left cardiac atrium":"LA", "right cardiac atrium":"RA", "interventricular septum":"IVS",
           "apical region of left ventricle":"AX", "atrioventricular node":"AV node", "sinoatrial node":"SA node"}
baseline = summary.loc[summary.disease.eq(summary.study.map(controls))].copy()
baseline["region"] = baseline.cardiac_region.replace(regions)
baseline = baseline.loc[baseline.donors_evaluated.ge(3) & baseline.feature_status.eq("unique")].copy()
assert not baseline.duplicated(["study","region","lineage","gene"]).any()
baseline["detected_in_majority"] = baseline.donor_fraction_detected.ge(0.5)
baseline["detected_in_any"] = baseline.donors_detected.ge(1)
recurrence = baseline.groupby(["region","lineage","gene"], as_index=False).agg(
    study_strata=("study","nunique"),
    study_strata_detected_in_majority=("detected_in_majority","sum"),
    study_strata_detected_in_any=("detected_in_any","sum"),
    median_of_study_median_cpm=("median_donor_cpm","median"))
recurrence["study_strata_fraction_detected_in_majority"] = (
    recurrence.study_strata_detected_in_majority / recurrence.study_strata)
recurrence.to_csv(root / "CYP_Control_Regional_Recurrence.tsv", sep="\t", index=False)

genes = coverage.loc[coverage.study.eq(coverage.study.iloc[0]), "gene"].tolist()
lineages = sorted(summary.lineage.unique())
cmap = plt.get_cmap("viridis").copy()
cmap.set_bad("#e7e7e7")
with PdfPages(root / "CYP_Control_Regional_Landscape.pdf") as pdf:
    for region in sorted(recurrence.region.unique()):
        table = recurrence.loc[recurrence.region.eq(region)]
        numerator = table.pivot(index="gene", columns="lineage", values="study_strata_detected_in_majority")
        denominator = table.pivot(index="gene", columns="lineage", values="study_strata")
        columns = [lineage for lineage in lineages if lineage in denominator]
        numerator = numerator.reindex(index=genes, columns=columns)
        denominator = denominator.reindex(index=genes, columns=columns)
        fraction = numerator / denominator
        fig, ax = plt.subplots(figsize=(max(8.7, 0.59 * len(columns) + 2.1), 13.1))
        fig.subplots_adjust(left=.18, bottom=.135, right=.86, top=.92)
        image = ax.imshow(fraction.to_numpy(), aspect="auto", cmap=cmap, vmin=0, vmax=1)
        ax.set_xticks(np.arange(len(columns)), [s.replace("_"," ") for s in columns], rotation=48, ha="right", fontsize=8)
        axis_names = {"CYP5A1": "CYP5A1 / TBXAS1", "CYP8A1": "CYP8A1 / PTGIS"}
        ax.set_yticks(np.arange(len(genes)), [axis_names.get(g, g) for g in genes], fontsize=7)
        for i in range(len(genes)):
            for j in range(len(columns)):
                d = denominator.iat[i,j]
                if pd.notna(d):
                    n = numerator.iat[i,j]
                    ax.text(j, i, f"{int(n)}/{int(d)}", ha="center", va="center", fontsize=5.6,
                            color="black" if n/d >= .58 else "white")
        ax.set_title(f"Cardiac CYP detection in control donors | {region}", fontsize=14, pad=15)
        fig.colorbar(image, ax=ax, fraction=.05, pad=.02, label="Fraction of evaluable study strata")
        fig.text(.13, .025, "n/N = study strata with CYP detected in at least half of eligible donors / evaluable strata.\n"
                 "Gray = no evaluable stratum. Minimum 3 donors, each with at least 20 nuclei. Study strata can share donors.",
                 fontsize=8)
        pdf.savefig(fig)
        plt.close(fig)

matched = summary.loc[summary.feature_status.eq("unique") & summary.donors_evaluated.ge(3)].copy()
qc = coverage.assign(unique=coverage.feature_status.eq("unique")).groupby("study", as_index=False).agg(
    unique_CYP_features=("unique","sum"))
stats = reviews.assign(supported=reviews.analysis_lineage.ne("unresolved"))
counts = stats.assign(used_nuclei=stats.nuclei.where(stats.supported, 0)).groupby("study", as_index=False).agg(
    clustered_nuclei=("nuclei","sum"), marker_supported_nuclei=("used_nuclei","sum"),
    clusters=("cluster","size"), supported_clusters=("supported","sum"))
detected = matched.loc[matched.donors_detected.ge(1)].groupby("study").gene.nunique().rename("CYP_detected_in_evaluable_strata")
qc = qc.merge(counts, validate="one_to_one").join(detected, on="study")
qc.to_csv(root / "CYP_Landscape_QC.tsv", sep="\t", index=False)
notes = (
    "Descriptive nuclear RNA landscape using deposited raw integer counts. CPM = CYP UMI divided by all UMI "
    "within each study/region/condition/lineage/donor stratum. Detection = at least one UMI. "
    "The landscape uses broad operational lineages with at least 3 enriched marker genes, score margin >=0.10, "
    "at least 3 contributing donors and no donor exceeding 50% of a cluster. Held clusters are preserved as "
    "unresolved in donor-level data and excluded from lineage summaries. These labels are separate from provisional "
    "and published labels. Descriptive summaries require at least 20 nuclei per donor stratum; recurrence additionally "
    "requires 3 donors per study stratum. The 57-gene roster follows the canonical human CYP literature roster, "
    "https://doi.org/10.3390/ijms17071020 . Uniquely mapped symbols, absent features and ambiguous source symbols are "
    "tracked separately. The CYP5A1 and CYP8A1 rows resolve the approved symbols TBXAS1 and PTGIS and their stable Ensembl IDs. No pooling of donors across studies. Study strata may share donors or source libraries. "
    "Kuppe is excluded from control recurrence because its ischemic, fibrotic and myogenic labels do not establish "
    "nonfailing donor controls in this extract. Region labels are combined only through explicit anatomical synonyms; "
    "AS IVS and LV apex strata remain separate from general LV. This is descriptive RNA expression and does not "
    "estimate disease effects or enzymatic activity.\n"
)
(root / "Landscape_Methods.txt").write_text(notes)
assert (root / "CYP_Control_Regional_Landscape.pdf").stat().st_size > 0
(root / "COMPLETED.txt").write_text(f"11 studies; {int(qc.clusters.sum())} clusters audited; 57 CYP genes.\n")
print(qc.to_string(index=False), flush=True)
