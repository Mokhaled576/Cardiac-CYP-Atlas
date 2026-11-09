"""Attach only covariates explicitly recorded in the cohort's count metadata."""
import sys
from pathlib import Path

import anndata as ad
import h5py
import pandas as pd

src, out = map(Path,sys.argv[1:3])
out.mkdir(parents=True,exist_ok=False)
samples = pd.read_csv(src/"Pseudobulk_Donors.tsv",sep="\t",dtype={"donor_id":str})
study, = pd.read_csv(src/"Contrasts.tsv",sep="\t").study.unique()
h5, = Path("qc").glob(f"ClusteringInput_{study}_*/Source_Counts.h5ad")
with h5py.File(h5) as handle:
    obs = handle["obs"]
    fields = [field for field in ("donor_id","age","sex") if field in obs]
    data = pd.DataFrame({field:ad.io.read_elem(obs[field]) for field in fields}).drop_duplicates()
data["donor_id"] = data.donor_id.astype(str)
if "age" in data:
    data["age"] = pd.to_numeric(data.age,errors="coerce").where(lambda x:x.between(18,110))
else:
    data["age"] = pd.NA
if "sex" in data:
    data["sex"] = data.sex.astype("string").str.lower().where(lambda x:x.isin(["female","male"]))
else:
    data["sex"] = pd.NA
assert samples.donor_id.isin(data.donor_id).all()
records = []
for donor,g in data.groupby("donor_id"):
    record = {"donor_id":donor}
    for field in ("age","sex"):
        unique = g[field].dropna().unique()
        record[field] = unique[0] if len(unique)==1 else pd.NA
        record[field+"_conflict"] = len(unique)>1
    records.append(record)
samples.merge(pd.DataFrame(records),on="donor_id",validate="many_to_one").to_csv(
    out/"Source_Donor_Covariates.tsv",sep="\t",index=False,na_rep="NA")
