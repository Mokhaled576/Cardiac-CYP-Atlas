"""Test donor-level ligand-receptor co-expression and within-lineage cluster occupancy."""
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t

source, out = map(Path,sys.argv[1:3])
out.mkdir(parents=True,exist_ok=False)
obs = pd.read_csv(source/"Pseudobulk_Donors.tsv",sep="\t",dtype={"donor_id":str})
features = pd.read_csv(source/"Pseudobulk_Features.tsv",sep="\t",dtype=str,keep_default_na=False)
contrasts = pd.read_csv(source/"Contrasts.tsv",sep="\t",dtype=str)
study, = contrasts.study.unique()
counts = pd.read_csv(source/"Pseudobulk_Counts.tsv.gz",sep="\t",dtype={"feature_id":str},index_col=0)
assert counts.index.tolist()==features.feature_id.tolist() and counts.columns.tolist()==obs["sample"].tolist()
assert not obs.duplicated(["cardiac_region","lineage","donor_id"]).any()
assert obs.nuclei.ge(20).all()
library = counts.sum(axis=0).to_numpy()
assert (library>0).all()
unique = features.symbol.ne("") & ~features.symbol.duplicated(keep=False)
counts = counts.loc[unique.to_numpy()]
counts.index = features.loc[unique,"symbol"]
resource = pd.read_csv("source_metadata/Regulatory_Resources/LIANA_omni_resource.csv",dtype=str,keep_default_na=False)
pairs = resource.loc[resource.resource.eq("consensus"),["source_genesymbol","target_genesymbol"]].drop_duplicates()
pairs.columns = ["ligand","receptor"]
complexes = sorted(set(pairs.ligand)|set(pairs.receptor))
matched = {x:x.split("_") for x in complexes if set(x.split("_")) <= set(counts.index)}
pairs = pairs.loc[pairs.ligand.isin(matched) & pairs.receptor.isin(matched)].reset_index(drop=True)
assert len(pairs)>1000
complex_names = sorted(set(pairs.ligand)|set(pairs.receptor))
complex_cpm = pd.DataFrame(np.vstack([
    counts.loc[matched[x]].min(axis=0).to_numpy()/library*1e6 for x in complex_names
]),index=complex_names,columns=counts.columns)
ligand_rows = complex_cpm.index.get_indexer(pairs.ligand)
receptor_rows = complex_cpm.index.get_indexer(pairs.receptor)
assert (ligand_rows>=0).all() and (receptor_rows>=0).all()

def effects(matrix, case):
    """Row-wise two-sided Welch test with confidence intervals."""
    n, m = int(case.sum()),int((~case).sum())
    a, b = matrix[:,case],matrix[:,~case]
    delta = a.mean(axis=1)-b.mean(axis=1)
    va, vb = a.var(axis=1,ddof=1)/n,b.var(axis=1,ddof=1)/m
    se2 = va+vb
    ok = np.isfinite(se2) & (se2>0)
    q = np.full(len(delta),np.nan)
    ci = np.full(len(delta),np.nan)
    df = se2[ok]**2/(va[ok]**2/(n-1)+vb[ok]**2/(m-1))
    q[ok] = 2*student_t.sf(np.abs(delta[ok])/np.sqrt(se2[ok]),df)
    ci[ok] = student_t.ppf(.975,df)*np.sqrt(se2[ok])
    return delta,q,delta-ci,delta+ci

lr_results,lr_audit = [],[]
for c in contrasts.itertuples(index=False):
    in_scope = obs.loc[obs.cardiac_region.eq(c.region) & obs.disease.isin([c.case,c.control])]
    groups = {k:t.set_index("donor_id",verify_integrity=True) for k,t in in_scope.groupby("lineage")}
    for sender,receiver in product(sorted(groups),repeat=2):
        s,r = groups[sender].align(groups[receiver],join="inner",axis=0)
        assert s.disease.equals(r.disease) and s.overlap_flag.equals(r.overlap_flag)
        is_case = s.disease.eq(c.case).to_numpy()
        for pop in ("all_donors","exclude_flagged_overlap"):
            select = np.ones(len(s),dtype=bool) if pop=="all_donors" else s.overlap_flag.eq("no").to_numpy()
            case = is_case[select]
            n,m = int(case.sum()),int((~case).sum())
            rec = dict(contrast=c.contrast,study=study,region=c.region,tier=c.tier,population=pop,
                       sender=sender,receiver=receiver,case_donors=n,control_donors=m)
            if pop!="all_donors" and select.all():
                lr_audit.append(dict(rec,status="identical_to_all_donors",testable_pairs=0));continue
            if min(n,m)<3:
                lr_audit.append(dict(rec,status="fewer_than_3_shared_donors_per_group",testable_pairs=0));continue
            l = complex_cpm.loc[:,s.loc[select,"sample"]].to_numpy()[ligand_rows]
            r_exp = complex_cpm.loc[:,r.loc[select,"sample"]].to_numpy()[receptor_rows]
            observed = (l>=1)&(r_exp>=1)
            keep = observed.sum(axis=1)>=max(3,int(np.ceil(.2*len(case))))
            if not keep.any():
                lr_audit.append(dict(rec,status="no_pairs_pass_expression_gate",testable_pairs=0));continue
            score = .5*(np.log2(1+l[keep])+np.log2(1+r_exp[keep]))
            delta,p,lo,hi = effects(score,case)
            result = pairs.loc[keep].copy()
            result = result.assign(**rec,score_delta=delta,p_value=p,ci_low=lo,ci_high=hi,
                case_pair_detection=observed[keep][:,case].mean(axis=1),
                control_pair_detection=observed[keep][:,~case].mean(axis=1))
            lr_results.append(result)
            lr_audit.append(dict(rec,status="modeled",testable_pairs=int(np.isfinite(p).sum())))

assert lr_results
pd.concat(lr_results,ignore_index=True).to_csv(out/"Ligand_Receptor_Effects.tsv.gz",sep="\t",index=False,compression="gzip",na_rep="NA")
pd.DataFrame(lr_audit).to_csv(out/"Ligand_Receptor_Eligibility.tsv",sep="\t",index=False)
pairs.to_csv(out/"Available_LIANA_Consensus_Pairs.tsv",sep="\t",index=False)

review = pd.read_csv("results/CYP_Landscape_2958286/Cluster_Landscape_Review.tsv",sep="\t",dtype=str)
review = review.loc[review.study.eq(study)].set_index("cluster",verify_integrity=True)
job, = review.cluster_job.unique()
cluster_path = Path("qc")/f"FreshCluster_{study}_{job}"/"Nucleus_Clusters.tsv.gz"
columns = ["donor_id","biosample_id","disease","cardiac_region","cluster"]
clusters = pd.read_csv(cluster_path,sep="\t",dtype=str,usecols=columns,keep_default_na=False)
clusters["lineage"] = clusters.cluster.map(review.analysis_lineage)
assert clusters.lineage.notna().all()
if study in {"Koenig","DoxCM"}:
    anatomy = pd.read_csv("results/Marker_Donor_Review_2957389/Published_Cohort_Anatomy.tsv",sep="\t",dtype=str)
    mapping = anatomy.loc[anatomy.study.eq(study)].set_index("biosample_id",verify_integrity=True)
    assert clusters.biosample_id.map(mapping.condition).equals(clusters.disease)
    clusters["cardiac_region"] = clusters.biosample_id.map(mapping.cardiac_region)
keys = ["cardiac_region","lineage","donor_id","disease"]
clusters = clusters.loc[clusters.lineage.ne("unresolved")]
tally = clusters.groupby(keys+["cluster"],sort=False).size().unstack(fill_value=0)
aligned = tally.reindex(pd.MultiIndex.from_frame(obs[keys]),fill_value=0)
assert np.array_equal(aligned.sum(axis=1).to_numpy(),obs.nuclei.to_numpy())
occ,occ_audit = [],[]
for c in contrasts.itertuples(index=False):
    for lineage in sorted(obs.loc[obs.cardiac_region.eq(c.region),"lineage"].unique()):
        within = obs.cardiac_region.eq(c.region)&obs.lineage.eq(lineage)&obs.disease.isin([c.case,c.control])
        for pop in ("all_donors","exclude_flagged_overlap"):
            select = within if pop=="all_donors" else within & obs.overlap_flag.eq("no")
            case = obs.loc[select,"disease"].eq(c.case).to_numpy()
            n,m = int(case.sum()),int((~case).sum())
            rec = dict(contrast=c.contrast,study=study,region=c.region,tier=c.tier,population=pop,
                lineage=lineage,case_donors=n,control_donors=m)
            if pop!="all_donors" and select.equals(within):
                occ_audit.append(dict(rec,status="identical_to_all_donors",testable_clusters=0));continue
            if min(n,m)<3:
                occ_audit.append(dict(rec,status="fewer_than_3_per_group",testable_clusters=0));continue
            x = aligned.iloc[select.to_numpy()].to_numpy().T
            y = obs.loc[select,"nuclei"].to_numpy()
            keep = ((x>0).sum(axis=1)>=3)&((x<y).sum(axis=1)>=3)
            if not keep.any():
                occ_audit.append(dict(rec,status="no_variable_supported_clusters",testable_clusters=0));continue
            delta,p,lo,hi = effects(np.log2((x[keep]+.5)/(y-x[keep]+.5)),case)
            occ.append(pd.DataFrame(dict(rec,cluster=aligned.columns[keep],log2_odds_delta=delta,p_value=p,
                ci_low=lo,ci_high=hi,median_case_fraction=np.median(x[keep][:,case]/y[case],axis=1),
                median_control_fraction=np.median(x[keep][:,~case]/y[~case],axis=1))))
            occ_audit.append(dict(rec,status="modeled",testable_clusters=int(np.isfinite(p).sum())))

assert occ
pd.concat(occ,ignore_index=True).to_csv(out/"Cluster_Occupancy_Effects.tsv",sep="\t",index=False,na_rep="NA")
pd.DataFrame(occ_audit).to_csv(out/"Cluster_Occupancy_Eligibility.tsv",sep="\t",index=False)
(out/"COMPLETED.txt").write_text(f"{study}: LIANA consensus and within-lineage occupancy complete.\n")
