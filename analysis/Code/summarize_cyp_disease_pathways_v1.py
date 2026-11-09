"""Summarize GO pathway enrichment across the disease-contrast family."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

root = Path(sys.argv[1])
folders = [root / f"Cohort_{i}" for i in (0,1,2,4,7,8,9)]
names = ["CYP_Disease_Pathways.tsv","Pathway_Gene_Set_Provenance.tsv"]
for folder in folders:
    assert all((folder / n).stat().st_size > 0 for n in names)
pathways = pd.concat([pd.read_csv(p / names[0],sep="\t") for p in folders],ignore_index=True)
provenance = pd.concat([pd.read_csv(p / names[1],sep="\t") for p in folders],ignore_index=True).drop_duplicates()
assert len(provenance)==5 and set(pathways.GO_ID)==set(provenance.GO_ID)
assert not pathways.duplicated(["contrast","lineage","GO_ID"]).any()
assert pathways.study.nunique()==7
pathways["GO_BH_across_tier"] = np.nan
for _, t in pathways.loc[pathways.pval.notna()].groupby("tier"):
    pathways.loc[t.index,"GO_BH_across_tier"] = multipletests(t.pval,method="fdr_bh")[1]
pathways = pathways.sort_values(["tier","GO_BH_across_tier","pval","contrast"])
pathways.to_csv(root / names[0],sep="\t",index=False,na_rep="NA")
provenance.to_csv(root / names[1],sep="\t",index=False)
(root/"Pathway_Methods.txt").write_text(
    "fgseaMultilevel with signed square root of the genome-wide edgeR quasi-likelihood F statistic. "
    "Ambiguous duplicate symbols omitted. Five prespecified human GO sets used GOALL symbol mappings "
    "from the Fir Bioconductor annotation release. Filtered sets require 10 to 2500 ranked features. "
    "BH was applied separately across all tested pathway-lineage-contrast combinations per primary "
    "or secondary tier, omitting results with missing p-values. Enrichment uses all-donor fits. "
    "Positive NES favors higher-ranked genes in disease. Pathway tests are exploratory.\n")
print(f"{len(pathways)} GO enrichment rows, {pathways.contrast.nunique()} contrasts",flush=True)
print(pathways.loc[pathways.tier.eq("primary"),
   ["contrast","lineage","pathway","NES","pval","GO_BH_across_tier"]].head(12).to_string(index=False),flush=True)
