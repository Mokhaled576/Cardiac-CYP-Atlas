from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t,ttest_ind
from statsmodels.stats.multitest import multipletests
base=Path(__file__).parent
x=pd.read_csv(base/"GSE116250_rpkm.txt.gz",sep="\t")
case=[c for c in x.columns if c.startswith("DCM")]
control=[c for c in x.columns if c.startswith("NF")]
assert len(case)==37 and len(control)==14
effects=[]
for gene in ("CYP2J2","CYP19A1"):
    found=x.loc[x.Common_name.eq(gene)]
    assert len(found)==1
    a=found[case].to_numpy(dtype=float).ravel()
    b=found[control].to_numpy(dtype=float).ravel()
    assert np.isfinite(a).all() and np.isfinite(b).all() and min(a.min(),b.min())>=0
    al,bl=np.log2(a+0.1),np.log2(b+0.1)
    effect=al.mean()-bl.mean()
    va,vb=al.var(ddof=1)/len(a),bl.var(ddof=1)/len(b)
    df=(va+vb)**2/(va**2/(len(a)-1)+vb**2/(len(b)-1))
    half=t.ppf(.975,df)*np.sqrt(va+vb)
    effects.append(dict(gene=gene,case_n=len(a),control_n=len(b),case_median_RPKM=np.median(a),control_median_RPKM=np.median(b),
        log2_RPKM_plus_0_1_delta=effect,ci95_low=effect-half,ci95_high=effect+half,
        p_value=ttest_ind(al,bl,equal_var=False).pvalue,case_detected=(a>0).sum(),control_detected=(b>0).sum()))
out=pd.DataFrame(effects)
out["FDR_BH"]=multipletests(out.p_value,method="fdr_bh")[1]
out.to_csv(base/"Bulk_DCM_External_Screen.tsv",sep="\t",index=False)
print(out.to_string(index=False))
