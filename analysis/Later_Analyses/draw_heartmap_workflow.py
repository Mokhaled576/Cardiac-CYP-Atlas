from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

out = Path(__file__).parent
with open("external_validation/GSE298023/Validation_3016861/Leuven_DCM_CYP_Validation.tsv") as f:
    effects = {r["gene"]:r for r in csv.DictReader(f,delimiter="\t") if r["model"]=="batch_adjusted"}
with open("results/Donor_CYP_Links_2969410/Permutation_Validated/Prespecified_DCM_Associations.tsv") as f:
    associations = {r["study"]:r for r in csv.DictReader(f,delimiter="\t")
                    if r["gene"]=="CYP2J2" and r["mode"]=="all_donors"}
assert set(effects)=={"CYP2J2","CYP19A1"} and set(associations)=={"Chaffin","Koenig","Reichart"}

ink="#142538"
muted="#47596d"
border="#d6e1e7"
teal="#006d70"
teal_light="#eaf6f4"
blue="#285080"
blue_light="#edf3fa"
warm="#954a34"
warm_light="#faf2ed"
pale="#f3f6f8"
fig=plt.figure(figsize=(15,8.6),facecolor="white")
ax=fig.add_axes([0,0,1,1],xlim=(0,100),ylim=(0,100))
ax.axis("off")

def text(x,y,s,size=12,color=ink,weight="normal",ha="left",va="center",**kw):
    ax.text(x,y,s,fontsize=size,color=color,fontweight=weight,ha=ha,va=va,**kw)
def box(x,y,w,h,fill="white",edge=border,r=1.0):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle=f"round,pad=0.12,rounding_size={r}",
                                facecolor=fill,edgecolor=edge,linewidth=1,zorder=0))
def panel(letter,x,y,title):
    box(x,y-1.7,2.4,3.4,fill=ink,edge=ink,r=.42)
    text(x+1.2,y,letter,12,"white","bold",ha="center")
    text(x+3.15,y,title,15,ink,"bold")

text(4,95.8,"Cell-type-specific CYP changes in human dilated cardiomyopathy",25,ink,"bold")
text(4,91.2,"HeartMap workflow, independent single-nucleus validation and donor-level association testing",13,muted)

panel("A",4,84.8,"From the CYP atlas to external validation")
stages=[
    ("01","CYP landscape","11 cardiac cohorts\n57 CYP genes",blue,blue_light),
    ("02","Disease contrasts","13 region-matched\ncontrasts in 7 cohorts",blue,blue_light),
    ("03","Biological context","Enrichment, signaling,\nregulons, communication\nand cell states",blue,blue_light),
    ("04","Donor associations","CYP and program scores\nwithin cohort and cell type",warm,warm_light),
    ("05","External validation","Leuven single-nucleus RNA-seq\n6 DCM and 4 non-failing hearts",teal,teal_light),
]
for i,(num,title,body,color,fill) in enumerate(stages):
    x=4+i*18.85
    box(x,64.7,16.6,15.5,fill=fill,edge=border,r=.85)
    text(x+1.15,77.8,num,10.5,color,"bold")
    text(x+1.15,74.5,title,11.7,ink,"bold")
    text(x+1.15,69.5,body,9.4,muted,linespacing=1.43)
    if i<4:
        ax.add_patch(FancyArrowPatch((x+16.86,72.4),(x+18.55,72.4),
                     arrowstyle="-|>",mutation_scale=11,linewidth=1.2,color="#61758a"))

box(4,16,51.2,43.2,fill="white")
box(57.4,16,38.6,43.2,fill="white")
panel("B",5.5,55.7,"Both CYP increases replicated in Leuven")
text(6.1,51.8,"DCM vs non-failing LV  |  batch-adjusted donor pseudobulk",10.5,muted)

text(6.2,45.5,"Gene / cell type",9.4,muted,"bold")
text(49.5,45.5,"FDR",9.4,muted,"bold",ha="center")
for gene,y in [("CYP2J2",39.4),("CYP19A1",30.5)]:
    r=effects[gene]
    text(6.2,y+1.15,gene,13,teal,"bold",fontstyle="italic")
    text(6.2,y-1.6,r["lineage"].capitalize(),10.2,muted)
    text(49.5,y,format(float(r["FDR_BH"]),".3f"),11,teal,"bold",ha="center")
forest=fig.add_axes([.195,.256,.268,.186],facecolor="none")
forest.set_xlim(-.42,6)
forest.set_ylim(-.57,1.30)
forest.set_yticks([])
forest.set_xticks([0,1,2,3,4,5,6])
forest.tick_params(axis="x",labelsize=9,colors=muted,length=3)
for side in ("left","right","top"): forest.spines[side].set_visible(False)
forest.spines["bottom"].set_color(border)
forest.axvline(0,color="#758596",linestyle=(0,(3,3)),linewidth=1.0,zorder=0)
for gene,y in [("CYP2J2",.82),("CYP19A1",-.077)]:
    r=effects[gene]
    x,lo,hi=[float(r[c]) for c in ("log2FC","CI_low","CI_high")]
    forest.errorbar(x,y,xerr=[[x-lo],[hi-x]],fmt="o",color=teal,elinewidth=2.2,
                    markersize=7,capsize=4,zorder=2)
    forest.text(x,y-.26,f"+{x:.2f}",ha="center",va="center",fontsize=9,color=teal,fontweight="bold")
forest.set_xlabel("log$_2$ fold change (95% CI)",fontsize=10,color=muted,labelpad=5)
text(6.2,18.6,"Both findings held in the shared-batch analysis (3 DCM, 4 controls).",9.5,muted)

panel("C",59,55.7,"CYP2J2 and estrogen by cohort")
text(59.8,51.8,"CYP2J2 vs inferred estrogen score within LV DCM cardiomyocytes",10.1,muted)
for x,label in [(59.8,"Cohort"),(77.9,"Donors"),(85.0,"Spearman $\\rho$"),(92.8,"FDR")]:
    text(x,45.5,label,9.6,muted,"bold",ha="left" if x==59.8 else "center")
for y,k in [(41.3,"Chaffin"),(36.8,"Koenig"),(32.3,"Reichart")]:
    row=associations[k]
    ax.hlines(y-2.0,59.8,94.2,linewidth=.65,color=border)
    text(59.8,y,k,11.1,ink,"bold")
    text(77.9,y,row["n_donors"],11,muted,ha="center")
    text(85.0,y,f'{float(row["rho"]):+.2f}',11,ink,ha="center")
    text(92.8,y,f'{float(row["Focal_FDR_BH"]):.2f}',11,muted,ha="center")
box(59.8,18.6,34.2,8.4,fill=warm_light,edge=warm_light,r=.65)
text(61.0,24.2,"No association passed the focal FDR threshold.",10.1,warm,"bold")
text(61.0,21.4,"The RNA-based estrogen score is an inferred pathway readout.",9.3,muted)

box(4,5.4,92,6.8,fill=teal_light,edge=teal_light,r=.65)
text(5.4,8.9,"Main finding",11,teal,"bold")
text(16.9,8.9,"DCM shows higher CYP2J2 RNA in cardiomyocytes and higher CYP19A1 RNA in fibroblasts.",11.6,ink,"bold")
text(4,2.7,"LV, left ventricle. DCM, dilated cardiomyopathy. FDR, false discovery rate. External cohort, GSE298023.",8.7,muted)
fig.savefig(out/"HeartMap_Workflow_and_Main_Findings.png",dpi=300,facecolor="white")
