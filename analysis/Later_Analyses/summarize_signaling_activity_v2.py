from pathlib import Path
from zipfile import ZipFile, ZIP_DEFLATED
import sys
import pandas as pd
from statsmodels.stats.multitest import multipletests

out = Path(sys.argv[1])
folders = [out/f"Cohort_{i}" for i in (0,1,2,4,7,8,9)]
assert all((folder/"COMPLETED.txt").is_file() for folder in folders)
def merge(name):
    return pd.concat([pd.read_csv(f/name,sep="\t") for f in folders],ignore_index=True)
results = merge("Signaling_Activity_Effects.tsv").rename(columns={"program":"pathway"})
assert results.contrast.nunique()==13 and results.pathway.nunique()==14
assert results.n_targets.ge(20).all()
assert not results.duplicated(["contrast","lineage","pathway","population"]).any()
assert results["P.Value"].between(0,1).all()
for _, t in results.groupby(["tier","population"]):
    results.loc[t.index,"FDR_BH"] = multipletests(t["P.Value"],method="fdr_bh")[1]
results = results.sort_values(["tier","population","FDR_BH","P.Value"])
results.to_csv(out/"Signaling_Activity_Effects.tsv",sep="\t",index=False)
merge("Model_Eligibility.tsv").to_csv(out/"Model_Eligibility.tsv",sep="\t",index=False)
merge("Donor_Signaling_Scores.tsv.gz").rename(columns={"program":"pathway"}).to_csv(
    out/"Donor_Signaling_Scores.tsv.gz",sep="\t",index=False,compression="gzip")
weights = merge("Footprint_Weights.tsv").drop_duplicates()
assert len(weights)==1400
weights.to_csv(out/"Footprint_Weights.tsv",sep="\t",index=False)
x = results.loc[results.population.eq("all_donors")].copy()
phenotypes={"Chaffin_DCM_LV":"DCM_LV","Koenig_DCM_LV":"DCM_LV","Reichart_DCM_LV":"DCM_LV",
            "Reichart_ARVC_RV":"ARVC_RV","Simonson_ARVC_RV":"ARVC_RV"}
x["phenotype_region"] = x.contrast.map(phenotypes)
rows=[]
for (ph,lin,pt), t in x.dropna(subset=["phenotype_region"]).groupby(["phenotype_region","lineage","pathway"]):
    if len(t)<2: continue
    assert t.study.is_unique
    rows.append(dict(phenotype_region=ph,lineage=lin,pathway=pt,tested_studies=len(t),
        direction_concordant=t.score_delta.gt(0).all() or t.score_delta.lt(0).all(),worst_FDR=t.FDR_BH.max(),
        effects="; ".join(f"{p.study}: {p.score_delta:+.2f} (q={p.FDR_BH:.3g})" for p in t.itertuples(index=False))))
conc=pd.DataFrame(rows)
conc.to_csv(out/"Signaling_Cross_Cohort_Concordance.tsv",sep="\t",index=False)
(out/"Methods.md").write_text(
    "# HeartMap signaling activity\n\n"
    "Donor-level pathway footprints are weighted scores from the human PROGENy perturbation model. "
    "For each of 14 pathways, the 100 targets with the smallest published model p-values were considered. "
    "Scores use TMM log-CPM and signed coefficients divided by the sum of absolute included weights. "
    "At least 20 expressed, uniquely mapped targets are required. Scores are standardized across the modeled donors "
    "within each contrast and lineage. Robust limma compares case and control scores using the same eligible donors, "
    "regions and broad lineages as the original CYP analysis. Effect units are within-model score SDs. "
    "BH correction is separate for primary and secondary tiers and for each population. "
    "The overlap screen excludes candidate shared identifiers. The models are unadjusted for clinical covariates. "
    "The results infer pathway footprints from RNA. They cannot measure pathway flux, CYP activity, phosphorylation "
    "or the direction of cell-cell signaling. PROGENy uses downstream response genes and complements the "
    "existing ranked GO enrichment analysis. The two answer different questions.\n\n"
    "Model source: https://github.com/saezlab/progeny/blob/master/data/model_human_full.rda\n\n"
    "Model SHA-256: 3094f1fa1bb5395c1074280402c0e30ef47ba9ec36f4829ca3e539d019255ee5\n")
name=out/"HeartMap_Signaling_Activity_2967196.zip"
with ZipFile(name,"w",ZIP_DEFLATED) as z:
    for f in ("Signaling_Activity_Effects.tsv","Model_Eligibility.tsv","Donor_Signaling_Scores.tsv.gz",
              "Footprint_Weights.tsv","Signaling_Cross_Cohort_Concordance.tsv","Methods.md"):
        z.write(out/f,f)
    for f in ("HeartMap_Research_Question.md","fit_signaling_activity_v2.R","run_signaling_activity_v2.sh",
              "summarize_signaling_activity_v2.py"):
        z.write(f,f)
    z.write("results/CYP_Disease_2962292/CYP_Disease_Pathways.tsv","Pathway_Enrichment/CYP_Disease_Pathways.tsv")
    z.write("results/CYP_Disease_2962292/Pathway_Methods.txt","Pathway_Enrichment/Pathway_Methods.txt")
with ZipFile(name) as z:
    assert z.testzip() is None
print(results.groupby(["tier","population"]).FDR_BH.agg(tests="size",q_lt_005=lambda y:y.lt(.05).sum()).to_string())
print(conc.loc[conc.direction_concordant & conc.worst_FDR.lt(.05)].sort_values("worst_FDR").head(14).to_string(index=False))
(out/"COMPLETED.txt").write_text(f"{len(results)} signaling tests, 14 pathway models, 13 contrasts.\n")
