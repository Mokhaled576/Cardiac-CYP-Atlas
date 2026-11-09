"""Combine completed HeartMap model families without pooling cohort effects."""
from hashlib import sha256
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile
import sys

import pandas as pd
from statsmodels.stats.multitest import multipletests

root = Path("results")
reference = root/"Reconciled_Reference_2965412"
tf = root/"Regulatory_State_2965168"
lr = root/"Communication_Occupancy_2965363"
cv = root/"Covariate_CYP_2965628"
disease = root/"CYP_Disease_2962292"
landscape = root/"CYP_Landscape_2958286"
out = root/f"Full_Analysis_{sys.argv[1]}"
out.mkdir(exist_ok=False)
idx = (0,1,2,4,7,8,9)
cv_idx = (1,7,8,9)
assert (reference/"COMPLETED.txt").is_file()
assert all((base/f"Cohort_{i}"/"COMPLETED.txt").is_file() for base in (tf,lr) for i in idx)
assert all((cv/f"Cohort_{i}"/"COMPLETED.txt").is_file() for i in cv_idx)

def combined(base, filename, indexes=idx, optional=False):
    paths = [base/f"Cohort_{i}"/filename for i in indexes]
    if optional:
        paths = [p for p in paths if p.is_file()]
    return pd.concat([pd.read_csv(p,sep="\t",low_memory=False) for p in paths],ignore_index=True)

records = []
shortlists = []

def adjust(frame, groups, pcol, keys, name, domain=None):
    assert not frame.duplicated(keys+groups).any(),name
    valid = frame[pcol].notna()
    assert frame.loc[valid,pcol].between(0,1).all(),name
    frame["FDR_BH"] = float("nan")
    for _, m in frame.loc[valid].groupby(groups,observed=True):
        frame.loc[m.index,"FDR_BH"] = multipletests(m[pcol],method="fdr_bh")[1]
    scope = frame
    if domain is not None:
        scope = frame.loc[frame.domain.eq(domain)]
    tally = scope.groupby(groups,observed=True).agg(
        rows=(pcol,"size"),finite_tests=(pcol,"count"),fdr_below_005=("FDR_BH",lambda x:x.lt(.05).sum())).reset_index()
    records.append(tally.assign(analysis=name))
    return frame

def write(frame, name):
    p = out/name
    frame.to_csv(p,sep="\t",index=False,na_rep="NA")
    return p

def highlights(frame, key, name, n=150):
    select = frame.tier.eq("primary")
    if "population" in frame:
        select &= frame.population.eq("all_donors")
    if "mode" in frame:
        select &= frame["mode"].eq("complete_case_adjusted")
    top = frame.loc[select & frame.FDR_BH.notna()].nsmallest(n,["FDR_BH",key])
    shortlists.append(write(top,name))

reg = combined(tf,"Regulatory_State_Effects.tsv.gz")
reg = adjust(reg,["tier","population","domain"],"P.Value",["contrast","lineage","program"],"Regulatory or state program")
state = combined(tf,"State_Program_Provenance.tsv").drop_duplicates()
assert len(state)==13
reg = reg.merge(state[["GO_ID","term"]],left_on="program",right_on="GO_ID",how="left",validate="many_to_one").drop(columns="GO_ID")
reg["term"] = reg.term.fillna(reg.program)
write(reg,"Regulatory_State_Effects.tsv.gz")
for kind in ("TF","state_program"):
    highlights(reg.loc[reg.domain.eq(kind)],"P.Value",f"Top_{kind}.tsv")
write(state,"State_Program_Provenance.tsv")
write(combined(tf,"Model_Eligibility.tsv"),"Regulatory_State_Eligibility.tsv")
del reg

ligand = combined(lr,"Ligand_Receptor_Effects.tsv.gz")
ligand = adjust(ligand,["tier","population"],"p_value",["contrast","sender","receiver","ligand","receptor"],"Ligand receptor")
highlights(ligand,"p_value","Top_Ligand_Receptor.tsv")
write(ligand,"Ligand_Receptor_Effects.tsv.gz")
write(combined(lr,"Ligand_Receptor_Eligibility.tsv"),"Ligand_Receptor_Eligibility.tsv")
del ligand

occ = combined(lr,"Cluster_Occupancy_Effects.tsv")
occ = adjust(occ,["tier","population"],"p_value",["contrast","lineage","cluster"],"Cluster occupancy")
write(occ,"Cluster_Occupancy_Effects.tsv")
write(combined(lr,"Cluster_Occupancy_Eligibility.tsv"),"Cluster_Occupancy_Eligibility.tsv")
highlights(occ,"p_value","Top_Cluster_Occupancy.tsv")

cov = combined(cv,"Covariate_CYP_Effects.tsv",cv_idx,optional=True)
cov = adjust(cov,["tier","mode"],"PValue",["contrast","lineage","cyp_gene"],"Covariate CYP")
write(cov,"Covariate_CYP_Effects.tsv")
write(combined(cv,"Covariate_Model_Eligibility.tsv",cv_idx),"Covariate_CYP_Eligibility.tsv")
highlights(cov,"PValue","Top_Covariate_CYP.tsv")

summary = pd.concat(records,ignore_index=True)
write(summary,"Analysis_Test_Counts.tsv")
cov_status = cov.loc[cov.tier.eq("primary") & cov["mode"].eq("complete_case_adjusted")]
cyp_anchor = cov_status.loc[
    ((cov_status.cyp_gene.eq("CYP2J2") & cov_status.lineage.eq("cardiomyocyte"))
     | (cov_status.cyp_gene.eq("CYP19A1") & cov_status.lineage.eq("fibroblast")))
    & cov_status.contrast.str.contains("DCM_LV")]
print(cyp_anchor[["contrast","lineage","cyp_gene","logFC","covariates","FDR_BH"]].to_string(index=False),flush=True)
print(summary.to_string(index=False),flush=True)

check = pd.read_csv(reference/"Cohort_Coverage.tsv",sep="\t")
assert len(check)==11 and check.status.eq("included").all()
links = pd.read_csv(reference/"Verified_Library_Crosswalk.tsv",sep="\t")
assert len(links)==23
print(f"{check.sampled_nuclei.sum()} integrated nuclei, {len(links)} shared HCA biospecimen IDs, "
      f"{check.verified_duplicate_library_nuclei_removed.sum()} Reichart nuclei excluded before reference sampling.",flush=True)

notes = """# HeartMap full analysis

- The CYP landscape covers 57 CYP genes across 11 cohorts. The disease models use 13 within-study, region-matched contrasts across seven cohorts, with donor-level pseudobulk from original RNA counts. Three control-only cohorts and Kuppe remain descriptive. The source-count layer, reviewed broad marker lineages, HFpEF exclusion and Adami screening rules carry forward.
- The integrated reference excludes Reichart entries sharing exact HCA biospecimen identifiers with the original Litvinukova source, preserves the crosswalk, and applies a 120-nucleus cap to reconciled donor groups. Other donor-ID overlaps remain candidates. The coordinates use shared unique gene symbols, study-corrected Harmony and UMAP. The reference is donor-sampled and is for visualization.
- Regulatory scores average signed, gene-standardized TMM log-CPM across at least 10 expressed CollecTRI targets in each donor profile. State programs use unsigned human GO target sets. Robust limma models compare donor scores within the same region and study.
- The ligand-receptor screen uses the LIANA human consensus resource. Every complex subunit must map uniquely, and expression is the minimum subunit CPM. The score is the average log2(1 + CPM) of ligand and receptor from the same donor's sender and receiver profiles. Pairs need both partners at 1 CPM or more in at least 3 donors and 20% of the modeled donors. Two-sided donor-level Welch tests compare disease and control. This measures RNA co-expression potential. It cannot establish cell proximity, secretion, receptor engagement, or signaling flux.
- Cluster occupancy compares donor-level log2 odds of an existing de novo cluster within its reviewed broad lineage, using a 0.5 pseudocount. Occupancy reflects within-lineage composition and can be affected by sampling and dissection. Clusters and GO programs are retained separately.
- Age plus sex sensitivity was fitted in Chaffin and Simonson ICM. Sex-only sensitivity was fitted where eligible in Reichart and Simonson ARVC. Eligibility uses explicitly recorded covariates, at least three donors per arm, a full rank design and at least three residual degrees of freedom. Adjusted and unadjusted fits use the same complete-case donors. The other cohorts lacked verified covariates in the source-count metadata used here.
- TF, state, ligand-receptor, occupancy and adjusted CYP families receive separate BH corrections by primary/secondary tier and sensitivity population or model, as applicable. Cohort effects were kept separate. RNA associations and overlap screens provide discovery evidence.
- Full results on Fir and laptop include per-cohort raw pseudobulk, complete ligand-receptor tests, donor program scores, integrated coordinates, lineage reviews and fit eligibility. This compact archive contains summary tables, selected ligand-receptor results, figures and the analysis code. The top tables include exploratory rows above FDR 0.05 when fewer than 150 pass.
"""
(out/"README.md").write_text(notes)
resources = Path("source_metadata/Regulatory_Resources")
checks = ["CollecTRI\tprior HeartMap frozen human CollecTRI signed snapshot\thttps://doi.org/10.1093/nar/gkad841",
          "LIANA\thuman consensus rows from LIANA omni_resource.csv\thttps://github.com/saezlab/liana-py/blob/main/src/liana/resource/omni_resource.csv"]
with (out/"Resource_Provenance.tsv").open("w") as f:
    f.write("resource\tdescription\tsource\tsha256\n")
    for row,name in zip(checks,("CollecTRI_regulons.csv","LIANA_omni_resource.csv")):
        f.write(row+"\t"+sha256((resources/name).read_bytes()).hexdigest()+"\n")
sheets = {
    "Landscape": (landscape, ["CYP_Control_Regional_Landscape.pdf","CYP_Study_Landscape.tsv","CYP_Control_Regional_Recurrence.tsv",
                             "CYP_Source_Coverage.tsv","CYP_Ventricular_Recurrence_Sensitivity.tsv","Landscape_Methods.txt"]),
    "Disease": (disease, ["CYP_Disease_Effects.tsv","CYP_Disease_Contrast_Summary.tsv","CYP_Disease_Overlap_Sensitivity.tsv",
                         "CYP_Disease_Test_Coverage.tsv","CYP_Disease_Cross_Cohort_Concordance.tsv","CYP_Disease_Pathways.tsv",
                         "Disease_Methods.txt","Pathway_Methods.txt"]),
    "Reference": (reference, ["Reconciled_Reference_UMAP.pdf","Cohort_Coverage.tsv","Verified_Library_Crosswalk.tsv",
                              "Reference_Status.txt","Reference_HVGs.tsv"]),
    "Downstream": (out, ["Analysis_Test_Counts.tsv","Regulatory_State_Effects.tsv.gz","State_Program_Provenance.tsv",
                         "Regulatory_State_Eligibility.tsv","Ligand_Receptor_Eligibility.tsv","Cluster_Occupancy_Effects.tsv",
                         "Cluster_Occupancy_Eligibility.tsv","Covariate_CYP_Effects.tsv","Covariate_CYP_Eligibility.tsv",
                         "Resource_Provenance.tsv",*[p.name for p in shortlists]])
}
code_files = ["Source_Manifest.tsv","CYP_Disease_Contrasts_v1.tsv","cyp_landscape_v2.py","summarize_cyp_landscape_v2.py",
              "build_cyp_disease_v1.py","fit_cyp_disease_v1.R","summarize_cyp_disease_v1.py",
              "run_cyp_disease_pathways_v1.R","summarize_cyp_disease_pathways_v1.py","make_reconciled_reference_v2.py",
              "fit_regulatory_states_v1.R","fit_communication_occupancy_v1.py","extract_cyp_covariates_v1.py",
              "fit_cyp_covariates_v1.R","summarize_heartmap_full_v1.py"]
package = out/"HeartMap_Full_Analysis_2965412.zip"
with ZipFile(package,"x",compression=ZIP_DEFLATED,compresslevel=5) as z:
    z.write(out/"README.md","README.md")
    for section,(folder,names) in sheets.items():
        for name in names:
            z.write(folder/name,f"{section}/{name}")
    for name in code_files:
        z.write(name,f"Code/{name}")
    z.write("logs/HeartMap_Full_v1.tsv","Code/HeartMap_Full_v1.tsv")
    z.write("logs/HeartMap_CYP_Disease_v1.tsv","Code/HeartMap_CYP_Disease_v1.tsv")
    z.write("logs/HeartMap_CYP_Landscape_v2.tsv","Code/HeartMap_CYP_Landscape_v2.tsv")
with ZipFile(package) as z:
    assert z.testzip() is None
(out/"COMPLETED.txt").write_text("Complete integrated reference and staged donor-level analysis. See README.md for scope and limits.\n")
print(package,package.stat().st_size,flush=True)
