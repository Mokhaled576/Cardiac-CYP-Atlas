"""Assemble verified within-study cardiac CYP contrasts."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

root = Path(sys.argv[1])
folders = [root / f"Cohort_{i}" for i in (0,1,2,4,7,8,9)]
assert all((p / "COMPLETED.txt").is_file() for p in folders)
design = pd.read_csv("CYP_Disease_Contrasts_v1.tsv",sep="\t")
audit = pd.concat([pd.read_csv(p/"Lineage_Eligibility.tsv",sep="\t") for p in folders],ignore_index=True)
results = pd.concat([pd.read_csv(p/"CYP_edgeR.tsv",sep="\t") for p in folders],ignore_index=True)
source = pd.concat([pd.read_csv(p/"CYP_Source_Coverage.tsv",sep="\t") for p in folders],ignore_index=True)
keys = ["contrast","lineage","cyp_gene"]
assert set(audit.contrast)==set(design.contrast) and len(design)==13
assert not results.duplicated(keys+["population"]).any()
assert results.PValue.notna().all() and results.PValue.between(0,1).all()
assert set(results.cyp_gene)<=set(source.gene)
results = results.merge(source[["study","gene","approved_symbol"]],left_on=["study","cyp_gene"],
                        right_on=["study","gene"],validate="many_to_one").drop(columns="gene")
for _, frame in results.groupby(["contrast","population"]):
    results.loc[frame.index,"CYP_BH_within_contrast"] = multipletests(frame.PValue,method="fdr_bh")[1]
for _, frame in results.groupby(["tier","population"]):
    results.loc[frame.index,"CYP_BH_across_tier"] = multipletests(frame.PValue,method="fdr_bh")[1]
results = results.sort_values(["tier","population","CYP_BH_across_tier","PValue","contrast"])
results.to_csv(root/"CYP_Disease_Effects.tsv",sep="\t",index=False)
audit.to_csv(root/"CYP_Disease_Eligibility.tsv",sep="\t",index=False)

base = results.loc[results.population.eq("all_donors")].copy()
retry = results.loc[results.population.eq("exclude_flagged_overlap"),keys+["logFC","PValue","CYP_BH_across_tier"]]
retry = retry.rename(columns={name:f"{name}_overlap_screen" for name in retry if name not in keys})
screen = audit.loc[audit.population.eq("exclude_flagged_overlap"),
                   ["contrast","lineage","status","case_donors","control_donors"]]
screen = screen.rename(columns={c:f"{c}_overlap_screen" for c in ("status","case_donors","control_donors")})
paired = base.merge(screen,on=["contrast","lineage"],validate="many_to_one").merge(retry,on=keys,how="left",validate="one_to_one")
paired["direction_retained"] = pd.NA
present = paired.logFC_overlap_screen.notna()
paired.loc[present,"direction_retained"] = (
    np.sign(paired.loc[present,"logFC"])==np.sign(paired.loc[present,"logFC_overlap_screen"])).to_numpy()
paired.to_csv(root/"CYP_Disease_Overlap_Sensitivity.tsv",sep="\t",index=False,na_rep="NA")

tally = base.groupby(["contrast","tier"],as_index=False).agg(
    CYP_tests=("cyp_gene","size"), CYP_genes_tested=("cyp_gene","nunique"), lineages_tested=("lineage","nunique"),
    tier_FDR_below_005=("CYP_BH_across_tier",lambda x:x.lt(0.05).sum()))
tally = design.merge(tally,on=["contrast","tier"],how="left",validate="one_to_one")
for col in ("CYP_tests","CYP_genes_tested","lineages_tested","tier_FDR_below_005"):
    tally[col]=tally[col].fillna(0).astype(int)
tally.to_csv(root/"CYP_Disease_Contrast_Summary.tsv",sep="\t",index=False)
source.to_csv(root/"CYP_Disease_Feature_Coverage.tsv",sep="\t",index=False)
planned = audit.loc[audit.population.eq("all_donors"),["contrast","study","lineage","status"]].merge(
    source[["study","gene","feature_status"]],on="study",how="left",validate="many_to_many")
planned = planned.rename(columns={"gene":"cyp_gene"})
planned = planned.merge(base[keys+["PValue"]],on=keys,how="left",validate="one_to_one")
planned["test_status"] = np.select(
    [planned.PValue.notna(),planned.status.ne("tested"),planned.feature_status.ne("unique")],
    ["tested",planned.status,planned.feature_status],default="filtered_low_expression")
planned.drop(columns="PValue").to_csv(root/"CYP_Disease_Test_Coverage.tsv",sep="\t",index=False)

(root/"Disease_Methods.txt").write_text(
    "Within-study, region-matched case versus control using genome-wide raw integer pseudobulk counts. "
    "Counts were summed by donor, region, condition and marker-supported lineage. Eligible donor-lineage groups "
    "had at least 20 nuclei and contrasts required three donors per arm. edgeR filterByExpr, TMM and robust "
    "quasi-likelihood models used a case indicator. Estimates are unadjusted for age, sex, medication, "
    "procurement and clinical covariates. Positive logFC indicates higher RNA abundance in the case group. "
    "Primary and secondary CYP tests received separate BH adjustments, and overlap-screen fits were adjusted "
    "separately. Within-contrast BH includes all tested CYP-lineage pairs. The overlap screen excludes flagged "
    "shared donor and sample identifiers within cohorts as sensitivity analysis. Shared donors or controls can "
    "still make studies and contrasts dependent, and study effects were not pooled. "
    "Published cell labels were not used. Lineage assignments remain marker-supported broad calls and held "
    "clusters were excluded. These are discovery analyses, and RNA changes do not establish enzyme activity. "
    "Kanemaru, Litvinukova and Tucker contribute to the descriptive atlas. Kuppe contributes descriptively "
    "because the available tissue-zone groups do not establish a matched nonfailing donor comparison. "
    "HFpEF excluded. Adami screening only.\n"
)
(root/"COMPLETED.txt").write_text(f"{len(folders)} cohorts; {len(design)} planned contrasts; {len(base)} CYP tests.\n")
print(tally.to_string(index=False),flush=True)
print(base.loc[base.tier.eq("primary"),["contrast","lineage","cyp_gene","logFC","PValue","CYP_BH_across_tier"]].head(20).to_string(index=False),flush=True)
