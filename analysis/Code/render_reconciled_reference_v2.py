from pathlib import Path
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

out=Path("results/Reconciled_Reference_2965412")
data=pd.read_csv(out/"Reference_Embedding.tsv.gz",sep="\t")
assert len(data)==33720 and data.study.nunique()==11
figure, axes=plt.subplots(1,2,figsize=(14,6.4),layout="constrained")
palette=plt.get_cmap("tab20").colors
for ax,key,title in zip(axes,("study","analysis_lineage"),("Source study","Marker-supported lineage")):
    groups=list(data.groupby(key,sort=True))
    assert len(groups)<=len(palette)
    colors=palette[::2]+palette[1::2]
    for (label,frame),color in zip(groups,colors):
        ax.scatter(frame.UMAP1,frame.UMAP2,s=1.35,color=[color],label=label.replace("_"," "),rasterized=True)
    ax.set_title(title)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.legend(fontsize=7,markerscale=3,loc="upper left",bbox_to_anchor=(1.01,1),frameon=False)
figure.suptitle("Cardiac reference after reconciliation of shared HCA specimens",fontsize=13,fontweight="bold")
temp=out/"Reconciled_Reference_UMAP.pdf.part"
figure.savefig(temp,format="pdf",dpi=180)
os.replace(temp,out/"Reconciled_Reference_UMAP.pdf")
