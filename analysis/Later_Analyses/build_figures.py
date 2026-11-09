from pathlib import Path
import shutil
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
P=Path("manuscript/TAC_style_v1/Figures")
R=Path("results")
L=R/"CYP_Landscape_2958286"; D=R/"CYP_Disease_2962292"; F=R/"Full_Analysis_2965831"; S=R/"Signaling_Activity_2967196"
A=R/"Donor_CYP_Links_2969410/Permutation_Validated"
V=Path("external_validation")
read=lambda p:pd.read_csv(p,sep='\t')
plt.rcParams.update({'font.size':10,'axes.titlesize':14,'axes.titleweight':'bold','axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
ink='#142538'; teal='#006d70'; blue='#285080'; warm='#a2442d'
def save(fig,n):
    for ext in ('png','pdf'): fig.savefig(P/f'Figure_{n:02}.{ext}',dpi=300,facecolor='white',bbox_inches='tight',pad_inches=.25)
    plt.close(fig)
def heat(ax,z,title,cmap='RdBu_r',lim=None,label='',dots=None):
    data=np.ma.masked_invalid(z.to_numpy(dtype=float))
    c=plt.get_cmap(cmap).copy();c.set_bad('#e8edf1')
    opts=dict(vmin=0,vmax=1) if lim is None else dict(norm=TwoSlopeNorm(vmin=-lim,vcenter=0,vmax=lim))
    im=ax.imshow(data,aspect='auto',cmap=c,**opts)
    ax.set_xticks(range(len(z.columns)),z.columns,rotation=42,ha='right')
    ax.set_yticks(range(len(z.index)),z.index)
    ax.set_title(title,loc='left',pad=12,color=ink)
    ax.spines[:].set_visible(False)
    if dots is not None:
        for i,j in zip(*np.where(dots.fillna(False).to_numpy(dtype=bool))):
            ax.plot(j,i,'o',ms=3.5,color='#17253b',mec='white',mew=.35)
    ax.figure.colorbar(im,ax=ax,fraction=.035,pad=.025,label=label)
def nice(s):
    return str(s).replace('_',' ').replace('Simonson ', 'Sim. ').replace('Reichart ','Reich. ')
def wide(d,index,col,val,order=None):
    t=d.pivot(index=index,columns=col,values=val)
    if order is not None: t=t.reindex(order)
    return t
shutil.copyfile("figures/Workflow_3016861/HeartMap_Workflow_and_Main_Findings.png",P/"Figure_01.png")
ref=R/"Reconciled_Reference_2965412"
x=read(ref/"Reference_Embedding.tsv.gz")
land=read(L/"CYP_Control_Regional_Recurrence.tsv")
fig,(a,b)=plt.subplots(1,2,figsize=(15,10),gridspec_kw={'width_ratios':[1,1.4]},layout='constrained')
colors=plt.get_cmap('tab20').colors
for j,(k,g) in enumerate(x.groupby('analysis_lineage')):
    a.scatter(g.UMAP1,g.UMAP2,s=1.1,color=colors[j],label=nice(k),rasterized=True)
a.set(title='A  Reconciled cardiac reference',xlabel='UMAP 1',ylabel='UMAP 2')
a.legend(frameon=False,fontsize=7,markerscale=4,ncol=2,loc='upper center',bbox_to_anchor=(.5,-.1))
z=land[land.region.eq('LV')].pivot(index='gene',columns='lineage',values='study_strata_fraction_detected_in_majority')
z=z.reindex(read(L/"CYP_Source_Coverage.tsv").gene.drop_duplicates().tolist());z.columns=[nice(k) for k in z.columns]
heat(b,z,'B  Control LV CYP detection','viridis',label='Fraction of evaluable study strata')
b.tick_params(axis='y',labelsize=6)
save(fig,2)

d=read(D/"CYP_Disease_Effects.tsv")
c=read(D/"CYP_Disease_Cross_Cohort_Concordance.tsv")
sel=c[(c.phenotype_region.eq("DCM_LV")) & c.direction_concordant & c.tier_fdr_below_005.eq(3)]
keys=list(dict.fromkeys(list(zip(sel.lineage,sel.cyp_gene))+[("fibroblast","CYP19A1")]))
order=[f'{g}  |  {nice(l)}' for l,g in keys]
q=d[(d.population.eq('all_donors')) & d.apply(lambda r:(r.lineage,r.cyp_gene) in keys,axis=1)].copy()
q['label']=[f'{g}  |  {nice(l)}' for l,g in zip(q.lineage,q.cyp_gene)]
contrasts=read(D/"CYP_Disease_Contrast_Summary.tsv").contrast.tolist()
z=wide(q,'label','contrast','logFC',order).reindex(columns=contrasts)
sig=wide(q.assign(sig=q.CYP_BH_across_tier.lt(.05)),'label','contrast','sig',order).reindex(columns=contrasts)
z.columns=[nice(k) for k in z.columns]
fig,ax=plt.subplots(figsize=(13,7),layout='constrained');heat(ax,z,'Recurrent LV DCM CYP signals across matched disease contrasts',lim=max(4,np.nanmax(abs(z.to_numpy()))),label='log2 fold change',dots=sig)
save(fig,3)

phen=['Chaffin_DCM_LV','Koenig_DCM_LV','Reichart_DCM_LV']
p=read(D/"CYP_Disease_Pathways.tsv")
s=read(S/"Signaling_Activity_Effects.tsv")
def panel_data(t,key,value,fdr,cols):
    t=t[t.contrast.isin(phen) & t.lineage.isin(cols)].copy()
    t['model']=t.contrast.str.split('_').str[0]+' | '+t.lineage
    return wide(t,key,'model',value),wide(t.assign(sig=t[fdr].lt(.05)),key,'model','sig')
fig,axes=plt.subplots(2,1,figsize=(12,10),layout='constrained',gridspec_kw={'height_ratios':[1,2]})
for a,t,k,v,qv,title,label in [(axes[0],p,'pathway','NES','GO_BH_across_tier','A  Prespecified GO enrichment','Normalized enrichment score'),
(axes[1],s[s.population.eq('all_donors')],'pathway','score_delta','FDR_BH','B  Inferred PROGENy signaling','Difference in score SD units')]:
    z,q=panel_data(t,k,v,qv,['cardiomyocyte','fibroblast'])
    heat(a,z,title,lim=max(1.5,np.nanmax(abs(z.to_numpy()))),label=label,dots=q)
save(fig,4)

t=read(F/"Regulatory_State_Effects.tsv.gz")
pick={'THRB':'THRB regulon','GO:0042775':'Mitochondrial ATP synthesis','GO:0030968':'ER unfolded protein response','GO:0001666':'Response to hypoxia'}
t=t[t.program.isin(pick) & t.population.eq('all_donors')].copy()
t['display']=t.program.map(pick)
z,q=panel_data(t,'display','score_delta','FDR_BH',['cardiomyocyte','fibroblast'])
fig,axes=plt.subplots(1,2,figsize=(14,5.5),layout='constrained',gridspec_kw={'width_ratios':[1.65,1]})
heat(axes[0],z,'A  Selected regulon and state-program scores',lim=max(1,np.nanmax(abs(z.to_numpy()))),label='Difference in mean score',dots=q)
x=read(F/"Analysis_Test_Counts.tsv");x=x[x.tier.eq('primary')&x.population.eq('all_donors')].copy()
x['display']=x.domain.fillna(x.analysis).replace({'state_program':'State programs','Ligand receptor':'Ligand-receptor','Cluster occupancy':'Cluster occupancy'})
pct=100*x.fdr_below_005/x.finite_tests
a=axes[1];pos=np.arange(len(x))
a.barh(pos,pct,color=[blue,blue,'#a96a49','#a96a49'])
a.set_yticks(pos,x.display);a.invert_yaxis();a.set_xlim(0,max(pct)*1.48)
for i,(u,n,m) in enumerate(zip(pct,x.fdr_below_005,x.finite_tests)): a.text(u+1,i,f'{n:,} / {m:,}',va='center',fontsize=9)
a.set(xlabel='Tests at FDR < 0.05 (%)',title='B  Exploratory test families')
save(fig,5)

t=read(A/"Prespecified_DCM_Associations.tsv")
t=t[t["mode"].eq('all_donors')].set_index(['gene','study'])
fig,ax=plt.subplots(figsize=(9,5.2),layout='constrained')
for i,g in enumerate(['CYP2J2','CYP19A1']):
    for j,st in enumerate(['Chaffin','Koenig','Reichart']):
        r=t.loc[(g,st)];y=5-i*3.5-j
        ax.scatter(r.rho,y,s=60,color=[teal,blue][i]);ax.text(1.08,y,f'n={r.n_donors}   q={r.Focal_FDR_BH:.2f}',va='center',fontsize=9)
ax.axvline(0,ls='--',color='#7b8995',lw=1)
ax.set_yticks([5,4,3,1.5,.5,-.5],[f'{g}  |  {st}' for g in ['CYP2J2','CYP19A1'] for st in ['Chaffin','Koenig','Reichart']])
ax.set(xlim=(-1,1.65),ylim=(-1.05,5.65),xticks=np.linspace(-1,1,5),xlabel='Spearman correlation with inferred estrogen score',
       title='Within-donor CYP and estrogen associations in LV DCM')
save(fig,6)

v=read(V/"GSE298023/Validation_3016861/Leuven_DCM_CYP_Validation.tsv")
w=read(V/"GSE116250/Bulk_DCM_External_Screen.tsv")
fig,axes=plt.subplots(1,2,figsize=(13,4.8),layout='constrained')
for a,rows,title,labels,key,lo,hi in [
 (axes[0],v[v.model.isin(['batch_adjusted','shared_batch_only'])],'A  Leuven single-nucleus validation',None,'log2FC','CI_low','CI_high'),
 (axes[1],w,'B  Independent bulk LV RNA-seq',w.gene.tolist(),'log2_RPKM_plus_0_1_delta','ci95_low','ci95_high')]:
    y=np.arange(len(rows))[::-1]
    if labels is None: labels=[f'{r.gene}  |  {nice(r.model)}' for r in rows.itertuples()]
    a.errorbar(rows[key],y,xerr=[rows[key]-rows[lo],rows[hi]-rows[key]],fmt='o',capsize=4,color=teal,lw=2)
    a.axvline(0,ls='--',color='#7b8995',lw=1)
    a.set_yticks(y,labels)
    a.set_ylim(-.55,len(rows)-.45);a.set_title(title,loc='left')
    a.set_xlabel('log2 fold change (95% CI)' if key=='log2FC' else 'Difference in log2(RPKM + 0.1) (95% CI)')
    for j,r in zip(y,rows.itertuples()): a.text(.97,j+.13,f'q={r.FDR_BH:.3g}',transform=a.get_yaxis_transform(),ha='right',fontsize=8,color=ink)
save(fig,7)
assert len(list(P.glob("Figure_*.png")))==7
