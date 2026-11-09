"""Within-study, within-condition donor associations for CYP expression and pathway scores."""
import sys
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

import numpy as np
import pandas as pd
from scipy.stats import rankdata, t as tdist
from statsmodels.stats.multitest import multipletests

source_root = Path(sys.argv[1])
root = source_root / "Permutation_Validated"
root.mkdir(exist_ok=False)
cohorts = (0,1,2,4,7,8,9)
assert all((source_root/f"Cohort_{i}"/"SCORES_COMPLETE.txt").is_file() for i in cohorts)
review = pd.read_csv("results/CYP_Landscape_2958286/Cluster_Landscape_Review.tsv",sep="\t",dtype=str)
anatomy = pd.read_csv("results/Marker_Donor_Review_2957389/Published_Cohort_Anatomy.tsv",sep="\t",dtype=str)
rng = np.random.default_rng(314159)
screen_rng = np.random.default_rng(271828)
condition_map = {"NF":"NF","Control":"NF","Donor":"NF","normal":"NF","PATO:0000461":"NF",
    "DCM":"DCM","dilated cardiomyopathy":"DCM","arrhythmogenic right ventricular cardiomyopathy":"ARVC",
    "MONDO:0016587":"ARVC","AS":"AS","DoxCM":"DoxCM","NICM":"NICM","HCM":"HCM","ICM":"ICM"}
modes = ["all_donors","exclude_flagged_overlap","complete_case_unadjusted","covariate_adjusted","overlap_and_adjusted"]
focus = {("cardiomyocyte","CYP2J2","PROGENy:Estrogen"),("fibroblast","CYP19A1","PROGENy:Estrogen")}
records, eligibility, focal = [],[],[]
info = []

def ranked(x):
    r = rankdata(x,axis=0,method="average")
    return r-r.mean(axis=0)

def permutation_p(x,y,rho):
    s = 0
    for _ in range(19999):
        s += abs(x @ rng.permutation(y)) >= abs(rho)-1e-12
    return (s+1)/20000

def permutation_screen(a,b,rho,asymptotic,q):
    """Conservative screen with Monte Carlo p values for selected small-stratum pairs."""
    n = a.shape[0]
    needs = asymptotic<=.05 if n<=20 else np.abs(rho)>.999999
    result = np.ones_like(asymptotic) if n<=20 else asymptotic.copy()
    if not needs.any():
        return result
    xi,yi = np.where(needs)
    observed = np.abs(rho[xi,yi])
    exceed = np.zeros(len(xi),dtype=np.int64)
    # Compute the same donor permutations for all candidate pairs in the stratum.
    for start in range(0,19999,256):
        block = min(256,19999-start)
        perm = np.vstack([screen_rng.permutation(n) for _ in range(block)])
        bp = b[perm,:]
        if q is not None:
            bp = bp - q @ (q.T @ bp)
            norm = np.linalg.norm(bp,axis=1,keepdims=True)
            bp = np.divide(bp,norm,out=np.zeros_like(bp),where=norm>1e-10)
        for offset in range(0,len(xi),128):
            w = slice(offset,offset+128)
            permuted = np.einsum("nc,bnc->bc",a[:,xi[w]],bp[:,:,yi[w]])
            exceed[w] += (np.abs(permuted)>=observed[w]-1e-12).sum(axis=0)
    result[xi,yi]=(exceed+1)/20000
    return result

for i in cohorts:
    p = source_root/f"Cohort_{i}"
    read = lambda name: pd.read_csv(p/name,sep="\t",dtype={"donor_id":str,"sample":str})
    m = read("Donor_Metadata.tsv").set_index("sample",verify_integrity=True)
    study, = m.study.unique()
    g = read("Donor_CYP_Expression.tsv")
    score = read("CYP_Free_Program_Scores.tsv")
    score["program"] = score.program.str.replace("GO:GO:","GO:",regex=False)
    base = score.pivot(index="sample",columns="program",values="score")
    cpm = g.pivot(index="sample",columns="gene",values="CPM")
    expr = g.pivot(index="sample",columns="gene",values="log_CPM")
    assert cpm.columns.equals(expr.columns)
    labels = review.loc[review.study.eq(study)].set_index("cluster",verify_integrity=True)
    job, = labels.cluster_job.unique()
    cols = ["cardiac_region","biosample_id","donor_id","disease","cluster"]
    nuclei = pd.read_csv(Path("qc")/f"FreshCluster_{study}_{job}"/"Nucleus_Clusters.tsv.gz",
                         sep="\t",dtype=str,usecols=cols,keep_default_na=False)
    nuclei["lineage"] = nuclei.cluster.map(labels.analysis_lineage)
    assert nuclei.lineage.notna().all()
    if study in {"Koenig","DoxCM"}:
        mapped = anatomy.loc[anatomy.study.eq(study)].set_index("biosample_id",verify_integrity=True)
        assert nuclei.biosample_id.map(mapped.condition).equals(nuclei.disease)
        nuclei["cardiac_region"] = nuclei.biosample_id.map(mapped.cardiac_region)
    keys = ["cardiac_region","lineage","donor_id","disease"]
    tally = nuclei.loc[nuclei.lineage.ne("unresolved")].groupby(keys+["cluster"]).size().unstack(fill_value=0)
    tally = tally.reindex(pd.MultiIndex.from_frame(m[keys]),fill_value=0)
    assert np.array_equal(tally.sum(axis=1).to_numpy(),m.nuclei.to_numpy())
    tally.index = m.index
    tally.columns = ["OCC:"+c for c in tally.columns]
    occupancy = tally.div(m.nuclei,axis=0)
    programs = base.join(occupancy,how="outer")
    m["phenotype"] = m.disease.map(condition_map)
    assert m.phenotype.notna().all(),m.disease.unique()
    fields = ["age","sex"] if study in {"Chaffin","Simonson_ICM"} else ["sex"] if study in {"Reichart","Simonson_ARVC"} else []
    for (region,lin,disease,pheno), group in m.groupby(["cardiac_region","lineage","disease","phenotype"]):
        assert group.donor_id.is_unique
        for mode in modes:
            adjusted = mode in ("covariate_adjusted","overlap_and_adjusted")
            complete = adjusted or mode=="complete_case_unadjusted"
            selected = group
            why = "tested"
            if mode in ("exclude_flagged_overlap","overlap_and_adjusted"):
                selected = selected.loc[selected.overlap_flag.eq("no")]
                if len(selected)==len(group):
                    why = "identical_to_all_donors" if mode=="exclude_flagged_overlap" else "identical_to_covariate_adjusted"
            if complete:
                if not fields: why = "source_covariates_unavailable"
                selected = selected.dropna(subset=fields)
            n = len(selected)
            covariates = "+".join(fields) if complete and fields else ""
            row = dict(study=study,region=region,lineage=lin,disease=disease,phenotype=pheno,mode=mode,
                       n_donors=n,covariates=covariates)
            if why=="tested" and n<8:
                why = "fewer_than_8_donors"
            if why!="tested":
                eligibility.append(dict(row,status=why,tests=0));continue
            if not selected.index.isin(expr.index).all():
                eligibility.append(dict(row,status="scores_unavailable",tests=0));continue
            nuisance = np.ones((n,1))
            if adjusted:
                if selected.sex.nunique()!=2:
                    eligibility.append(dict(row,status="sex_invariant_or_missing",tests=0));continue
                nuisance = np.column_stack([nuisance,selected.sex.eq("female").astype(float)])
                if "age" in fields:
                    nuisance = np.column_stack([nuisance,ranked(selected.age.to_numpy())])
            k = nuisance.shape[1]
            if np.linalg.matrix_rank(nuisance)!=k or n-k-1<4:
                eligibility.append(dict(row,status="rank_or_residual_df_failure",tests=0));continue
            x = expr.loc[selected.index]
            detect = cpm.loc[selected.index]
            x = x.loc[:,detect.ge(1).sum(axis=0).ge(3)&x.nunique().ge(3)]
            y = programs.loc[selected.index].dropna(axis=1)
            keep = y.nunique().ge(3)
            occ_cols = y.columns.str.startswith("OCC:")
            keep.loc[occ_cols] &= y.loc[:,occ_cols].gt(0).sum(axis=0).ge(3)
            y = y.loc[:,keep]
            if x.empty or y.empty:
                eligibility.append(dict(row,status="no_variable_expressed_gene_or_score",tests=0));continue
            a,b = ranked(x.to_numpy()),ranked(y.to_numpy())
            q = None
            if adjusted:
                q,_ = np.linalg.qr(nuisance)
                a,b = a-q@(q.T@a), b-q@(q.T@b)
            norm_a,norm_b = np.linalg.norm(a,axis=0),np.linalg.norm(b,axis=0)
            ax,bx = norm_a>1e-10,norm_b>1e-10
            a,b = a[:,ax]/norm_a[ax],b[:,bx]/norm_b[bx]
            xnames,ynames = x.columns[ax],y.columns[bx]
            rho = np.clip(a.T@b,-1,1)
            df = n-k-1
            pval = 2*tdist.sf(np.abs(rho)*np.sqrt(df/np.maximum(1-rho*rho,1e-16)),df)
            screened = permutation_screen(a,b,rho,pval,q)
            frame = pd.DataFrame(dict(row,gene=np.repeat(xnames,len(ynames)),program=np.tile(ynames,len(xnames)),
                                      rho=rho.ravel(),p_value=screened.ravel(),asymptotic_p=pval.ravel(),df=df))
            eligibility.append(dict(row,status="tested",tests=len(frame)))
            records.append(frame)
            if pheno=="DCM" and region=="LV" and study in {"Chaffin","Koenig","Reichart"}:
                for lineage,gene,program in sorted(focus):
                    if lin!=lineage or gene not in xnames or program not in ynames: continue
                    ai,bi=xnames.get_loc(gene),ynames.get_loc(program)
                    fp = pval[ai,bi] if adjusted else permutation_p(a[:,ai],b[:,bi],rho[ai,bi])
                    focal.append(dict(row,gene=gene,program=program,rho=rho[ai,bi],p_value=fp,
                                      p_method="partial_rank_t" if adjusted else "19999_permutations"))
    coverage = read("Program_Coverage.tsv")
    coverage["program"] = coverage.program.str.replace("GO:GO:","GO:",regex=False)
    info.append(coverage)

results = pd.concat(records,ignore_index=True)
results["family"] = results.program.str.split(":",n=1).str[0]
assert not results.duplicated(["study","region","lineage","disease","mode","gene","program"]).any()
for _, x in results.groupby(["mode","family"]):
    results.loc[x.index,"FDR_BH"] = multipletests(x.p_value,method="fdr_bh")[1]
results = results.sort_values(["mode","family","FDR_BH","p_value"])
results.to_csv(root/"CYP_Donor_Associations.tsv.gz",sep="\t",index=False,compression="gzip")
pd.DataFrame(eligibility).to_csv(root/"Donor_Association_Eligibility.tsv",sep="\t",index=False)
pd.concat(info,ignore_index=True).to_csv(root/"Program_Coverage.tsv",sep="\t",index=False)
focal = pd.DataFrame(focal)
for _, x in focal.groupby("mode"):
    focal.loc[x.index,"Focal_FDR_BH"] = multipletests(x.p_value,method="fdr_bh")[1]
focal.to_csv(root/"Prespecified_DCM_Associations.tsv",sep="\t",index=False)

dcm = {"Chaffin","Koenig","Reichart"}
arvc = {"Reichart","Simonson_ARVC"}
grouped = results.loc[(results.phenotype.eq("DCM")&results.region.eq("LV")&results.study.isin(dcm) |
                       results.phenotype.eq("ARVC")&results.region.eq("RV")&results.study.isin(arvc)) &
                       results.family.ne("OCC")]
concordance=[]
for key, x in grouped.groupby(["phenotype","region","lineage","mode","family","gene","program"]):
    assert x.study.is_unique
    if len(x)<2: continue
    # Directional partial-conjunction test for evidence in at least two cohorts.
    positive = np.where(x.rho>0,x.p_value/2,1-x.p_value/2)
    negative = np.where(x.rho<0,x.p_value/2,1-x.p_value/2)
    factor = 2*(len(x)-1)
    pp,pn = np.sort(positive)[1]*factor,np.sort(negative)[1]*factor
    concordance.append(dict(zip(["phenotype","region","lineage","mode","family","gene","program"],key),
        cohorts_tested=len(x),all_directions_agree=x.rho.gt(0).all() or x.rho.lt(0).all(),
        direction="positive" if pp<pn else "negative",p_at_least_two=min(1,pp,pn),
        estimates="; ".join(f"{r.study}: rho={r.rho:+.2f}, n={r.n_donors}, p={r.p_value:.3g}" for r in x.itertuples())))
c = pd.DataFrame(concordance)
for _, x in c.groupby(["mode","family"]):
    c.loc[x.index,"Replication_FDR_BH"] = multipletests(x.p_at_least_two,method="fdr_bh")[1]
c = c.sort_values(["mode","family","Replication_FDR_BH","p_at_least_two"])
c.to_csv(root/"CYP_Donor_Cross_Cohort_Evidence.tsv.gz",sep="\t",index=False,compression="gzip")
main = results.loc[results["mode"].eq("all_donors")]
main.loc[main.FDR_BH.lt(.05)].to_csv(root/"Main_Screen_FDR_below_005.tsv",sep="\t",index=False)
c.loc[c.Replication_FDR_BH.lt(.05)].to_csv(root/"Cross_Cohort_FDR_below_005.tsv",sep="\t",index=False)
print(f"All-donor screen: {len(main)} tests; {main.FDR_BH.lt(.05).sum()} at FDR < 0.05")
print(c.groupby("mode").Replication_FDR_BH.agg(tests="size",fdr_below_005=lambda q:q.lt(.05).sum()).to_string())
print(focal.to_string(index=False),flush=True)
(root/"Methods.md").write_text(
    "# Donor-level CYP associations\n\n"
    "Study, cardiac region, disease condition and reviewed broad cell lineage define every model. A single donor contributes "
    "one expression profile per stratum. Models require at least eight donors, at least three with CYP expression at 1 CPM "
    "or greater, and variable expression and scores. Normalization is TMM within study, region and lineage. "
    "The screen uses Spearman rank correlation. For strata with at most 20 donors, a preliminary rank t screen "
    "at p<=0.05 is followed by 19,999 permutations; other pairs are assigned p=1 conservatively. "
    "For larger strata the t approximation is used except for nearly perfect correlations, which are permuted. "
    "Adjusted small-stratum permutations reshuffle nuisance residuals and reproject them before correlation. "
    "Adjusted models residualize ranks "
    "against verified sex and, where available, ranked age. Complete-case unadjusted fits provide a matched sensitivity. "
    "Overlap-screen fits exclude flagged identifiers. Source age plus sex were available in Chaffin and Simonson ICM, "
    "and sex in Reichart and Simonson ARVC; no values were imputed. FDR is separate by score family and model. "
    "The two preselected estrogen-footprint tests use two-sided 19,999-permutation p values in unadjusted strata, "
    "and partial-rank t tests after covariate adjustment. BH covers tested focal comparisons separately by model.\n\n"
    "Pathway scores omit all canonical CYP genes and the TBXAS1 and PTGIS aliases. GO scores average standardized "
    "expression of at least ten targets, and PROGENy scores use signed weights from at least twenty expressed "
    "targets among the top one hundred. Cluster occupancy is a donor's fraction in an existing cluster within its "
    "broad lineage. These clusters are study-specific and are excluded from cross-cohort replication testing. "
    "Cross-cohort evidence is a directional partial-conjunction test for at least two cohorts, with Bonferroni "
    "control across directions and BH across pairs within family and model. Cohort effects are not pooled. "
    "CYP expression and inferred RNA programs can co-vary for several reasons, including within-lineage cell state and composition.\n"
)
names=["CYP_Donor_Associations.tsv.gz","Donor_Association_Eligibility.tsv","Program_Coverage.tsv",
       "Prespecified_DCM_Associations.tsv","CYP_Donor_Cross_Cohort_Evidence.tsv.gz","Main_Screen_FDR_below_005.tsv",
       "Cross_Cohort_FDR_below_005.tsv","Methods.md"]
archive=root/"HeartMap_Donor_CYP_Associations_2969410.zip"
with ZipFile(archive,"w",ZIP_DEFLATED) as z:
    for n in names: z.write(root/n,n)
    for n in ("HeartMap_Research_Question.md","donor_cyp_link_scores_v1.R","fit_donor_cyp_links_v2.py","run_donor_cyp_links_v1.sh","run_fit_donor_cyp_links_v2.sh"):
        z.write(n,n)
with ZipFile(archive,"a",ZIP_DEFLATED) as z:
    bulk=Path("external_validation/GSE116250")
    for n in ("Bulk_DCM_External_Screen.tsv","validate_bulk_cyp.py"):
        z.write(bulk/n,f"Published_Bulk_Check/{n}")
    z.writestr("Published_Bulk_Check/README.md",
        "GSE116250 published LV bulk RNA-seq RPKM. 37 DCM and 14 non-failing donor samples. "
        "Welch tests used log2(RPKM+0.1) for two prespecified CYP genes with BH correction. "
        "Bulk data cannot localize expression changes to a cell type. "
        "The external single-nucleus validation is separate. Source: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE116250\n")
with ZipFile(archive) as z: assert z.testzip() is None
(root/"COMPLETED.txt").write_text(f"{len(results)} association tests; {len(c)} cross-cohort tests.\n")
