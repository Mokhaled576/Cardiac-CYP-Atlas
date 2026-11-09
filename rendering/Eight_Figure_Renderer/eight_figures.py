"""Generate the final cardiac CYP display items from the audited local exports."""
from pathlib import Path
import io
import re
import shutil
import subprocess
import zipfile

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from pypdf import PdfReader, PdfWriter
from reportlab.lib.colors import HexColor
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen.canvas import Canvas

import figures as a

ROOT=a.ROOT
OUT=ROOT/"eight_figure_revision"
MAIN=OUT/"Main_Figures"
EXT=OUT/"Extended_Data"
DATA=ROOT/"eight_source"
PREV=ROOT/"ten_figure_revision"
for path in (OUT,MAIN,EXT,DATA):
    path.mkdir(exist_ok=True)
a.FIG=MAIN
a.DATA=DATA
LIN=["cardiomyocyte","fibroblast","endothelial"]
CH=["LV","RV","LA","RA"]
SHORT=dict(zip(LIN,("CM","FB","EC")))
DCM=["Chaffin_DCM_LV","Koenig_DCM_LV","Reichart_DCM_LV"]
ARVC=["Reichart_ARVC_RV","Simonson_ARVC_RV"]
ALL=a.read(a.D/"CYP_Disease_Effects.tsv").query("population == 'all_donors'").copy()
REG=a.read(a.L/"CYP_Control_Regional_Recurrence.tsv")
REG=REG[REG.region.isin(CH)&REG.lineage.isin(LIN)].copy()
REG["chamber_lineage"]=REG.region+" "+REG.lineage
REG["log_cpm"]=np.log10(1+REG.median_of_study_median_cpm)
GENES=sorted(REG.gene.unique(),key=lambda x:[int(t) if t.isdigit() else t for t in re.split(r"(\d+)",x)])
KEY=["CYP2J2","CYP19A1","CYP1B1","CYP4B1","CYP4X1","CYP4Z1"]
COLS=[r+" "+l for r in CH for l in LIN]
SIG=a.read(a.S/"Signaling_Activity_Effects.tsv").query("population == 'all_donors'")
FOCAL=a.read(a.F/"Prespecified_DCM_Associations.tsv")
FOLDER=a.RESULT/"Full_Analysis_2965831"
assert len(GENES)==57 and len(ALL)==1234


def archive(name,df):
    a.source(name,df)


def axis(p,x,y,w,ticks,lo,hi,label):
    p.line(x,y,x+w,y)
    for t in ticks:
        u=x+(t-lo)/(hi-lo)*w
        p.line(u,y,u,y+3)
        p.text(u,y+12,f"{t:g}",6.3,align="c")
    p.text(x+w/2,y+23,label,7,align="c")


def matrix(p,df,genes,field,y,rh,legend,maximum):
    x,w=103,390
    for j,k in enumerate(COLS):
        region,lineage=k.split(" ",1)
        p.text(x+(j+.5)*w/12,y-7,SHORT[lineage],6.5,align="c")
        if j%3==1:
            p.text(x+(j+.5)*w/12,y-17,region,7.2,True,align="c")
    a.heat(p,df,[(g,) for g in genes],COLS,["gene"],"chamber_lineage",field,None,
           (x,y,w,len(genes)*rh),genes,maximum,legend,signed=False)


def umap(p,frame,x,y,size,group,xcol,ycol):
    names=sorted(frame[group].astype(str).unique())
    colors={k:a.COLORS[i%len(a.COLORS)] for i,k in enumerate(names)}
    px=int(np.ceil(size*650/72))
    im=Image.new("RGB",(px,px),"white")
    dr=ImageDraw.Draw(im)
    xv=frame[xcol].to_numpy()
    yv=frame[ycol].to_numpy()
    f=(px-18)/max(np.ptp(xv),np.ptp(yv))
    u=(xv-(xv.min()+xv.max())/2)*f+px/2
    v=((yv.max()+yv.min())/2-yv)*f+px/2
    for xx,yy,k in zip(u,v,frame[group].astype(str)):
        dr.ellipse((xx-1.7,yy-1.7,xx+1.7,yy+1.7),fill=colors[k])
    p.c.drawImage(ImageReader(im),x,p.h-y-size,size,size)
    for i,k in enumerate(names):
        yy=y+8+i*8.4
        p.dot(x+size+4,yy-2,1.6,colors[k])
        p.text(x+size+8,yy,k.replace("_"," "),5.8)
    p.text(x+size/2,y+size+10,"UMAP 1 · UMAP 2",6.3,align="c")


def fig1():
    qc=a.read(a.L/"CYP_Landscape_QC.tsv")
    emb=a.read(a.R/"Reference_Embedding.tsv.gz")
    assert len(emb)==33720
    p=a.Page(1,424)
    p.panel(8,14,"a","Cohort assembly and lineage review")
    for j,(top,bottom) in enumerate([
        ("11 source cohorts","3,200,270 clustered nuclei"),
        ("239 candidate clusters","172 passed lineage review"),
        ("2,778,472 retained nuclei","67 unresolved clusters held out"),
    ]):
        x=12+168*j
        p.rect(x,26,157,38,"#e8f2f4")
        p.text(x+7,42,top,8,True)
        p.text(x+7,55,bottom,7)
        if j<2:
            p.text(x+160,47,"→",8)
    p.panel(8,86,"b","Retained nuclei and cluster support by cohort")
    for i,r in enumerate(qc.itertuples()):
        y=101+10*i
        p.text(90,y+6,r.study.replace("_"," "),6.7,align="r")
        p.rect(97,y,r.clustered_nuclei/2800,6.7,"#dfe4ea")
        p.rect(97,y,r.marker_supported_nuclei/2800,6.7,a.TEAL)
        p.text(434,y+6,f"{r.supported_clusters}/{r.clusters} clusters",6.4)
    p.text(98,222,"Teal: retained    Gray: unresolved for lineage-level inference",6.6)
    p.panel(8,246,"c","Harmony reference by lineage")
    p.panel(268,246,"d","Same reference by cohort")
    umap(p,emb,12,257,143,"analysis_lineage","UMAP1","UMAP2")
    umap(p,emb,272,257,143,"study","UMAP1","UMAP2")
    p.finish()
    archive("Fig1a_b_cohort_retention",qc)
    archive("Fig1c_d_reference_coordinates",emb[["study","analysis_lineage","UMAP1","UMAP2"]])


def fig2():
    donor=a.read(a.L/"CYP_Donor_Landscape.tsv.gz")
    controls={"NF","Control","Donor","normal","PATO:0000461"}
    sub=donor[donor.cardiac_region.eq("LV")&donor.disease.isin(controls)&donor.donor_eligible&
              donor.lineage.isin(LIN)&donor.gene.isin(["CYP2J2","CYP4X1","CYP1B1"])&donor.study.ne("Kuppe")].copy()
    assert sub.study.nunique()==7
    sub["log_cpm"]=np.log10(1+sub.cpm)
    archive("Fig2a_b_regional_baseline",REG)
    archive("Fig2c_d_LV_donor_profiles",sub.drop(columns="donor_id"))
    p=a.Page(2,681)
    p.panel(8,13,"a","Abundance across the 57-gene CYP roster")
    matrix(p,REG,GENES,"log_cpm",38,6.0,"log₁₀(1 + median study-median CPM)",3)
    p.panel(8,417,"b","Non-failing LV recurrence across study strata")
    for j,l in enumerate(LIN):
        x=17+j*167
        p.text(x,431,l,7.1,True)
        for i,g in enumerate(KEY):
            row=REG[REG.region.eq("LV")&REG.lineage.eq(l)&REG.gene.eq(g)].iloc[0]
            y=438+i*11
            p.text(x,y+6,g,6)
            p.rect(x+43,y,61,7,"#e5eaee")
            p.rect(x+43,y,61*row.study_strata_fraction_detected_in_majority,7,a.TEAL)
            p.text(x+108,y+6,f"{row.study_strata_detected_in_majority}/{row.study_strata}",6)
    p.panel(8,520,"c","LV donor-level abundance distributions")
    p.panel(274,520,"d","Individual LV donor profiles")
    for j,g in enumerate(["CYP2J2","CYP4X1","CYP1B1"]):
        y=551+j*35
        p.text(8,y+1,g,6.5)
        p.text(274,y+1,g,6.5)
        for k,l in enumerate(LIN):
            s=sub[sub.gene.eq(g)&sub.lineage.eq(l)]
            vals=s.log_cpm.to_numpy()
            if len(vals)==0:
                continue
            yy=y-7+k*7
            grid=np.linspace(0,3,64)
            den=np.exp(-.5*((grid[:,None]-vals[None,:])/.16)**2).sum(axis=1)
            amp=2.7*den/max(den)
            path=p.c.beginPath()
            p.c.setFillColor(HexColor(a.COLORS[k]))
            path.moveTo(66,p.h-yy)
            for u,z in zip(grid,amp):
                path.lineTo(66+55*u,p.h-yy+z)
            for u,z in zip(grid[::-1],amp[::-1]):
                path.lineTo(66+55*u,p.h-yy-z)
            path.close()
            p.c.drawPath(path,stroke=0,fill=1)
            for i,v in enumerate(vals):
                p.dot(334+52*v,yy+((i*13)%9-4)*.30,1,a.COLORS[k])
    axis(p,66,644,165,[0,1,2,3],0,3,"log₁₀(1 + CPM)")
    axis(p,334,644,156,[0,1,2,3],0,3,"log₁₀(1 + CPM)")
    for j,l in enumerate(LIN):
        p.dot(125+110*j,674,1.8,a.COLORS[j])
        p.text(131+110*j,676,l,6.4)
    p.finish()


def fig3():
    q="CYP_BH_across_tier"
    c=a.read(a.D/"CYP_Disease_Contrast_Summary.tsv")
    hits=ALL[ALL.lineage.eq("fibroblast")&ALL.cyp_gene.eq("CYP4B1")]
    hits=hits.set_index("contrast").loc[[k for k in a.CONTRASTS if k in set(hits.contrast)]].reset_index()
    highlight=ALL[(ALL.contrast.eq("AS_IVS")&ALL.lineage.isin(LIN)&ALL.cyp_gene.isin(["CYP1B1","CYP46A1"]))|
                   (ALL.contrast.eq("DoxCM_vs_NICM_LV")&ALL.lineage.eq("cardiomyocyte")&ALL.cyp_gene.eq("CYP19A1"))]
    archive("Fig3a_all_pseudobulk_tests",ALL)
    archive("Fig3b_fibroblast_CYP4B1",hits)
    archive("Fig3c_etiology_highlights",highlight)
    p=a.Page(3,470)
    p.panel(8,13,"a","Matched pseudobulk tests across 13 contrasts")
    for i,r in enumerate(ALL.itertuples()):
        u=44+(min(max(r.logFC,-8),8)+8)*14.7
        yy=111-min(-np.log10(max(getattr(r,q),1e-13)),13)*5.9
        p.dot(u,yy,1.15,(a.TEAL if getattr(r,q)<.05 else "#a4aeb9") if r.tier=="primary" else "#a05d32")
    p.line(44,34,44,114)
    for tick in (0,4,8,12):
        y=111-tick*5.9
        p.line(41,y,44,y)
        p.text(39,y+2,tick,5.6,align="r")
    axis(p,44,114,235,[-8,-4,0,4,8],-8,8,"log₂ fold change")
    p.text(20,31,"−log₁₀(q)",6)
    p.text(48,24,"Teal: primary q < 0.05   Gray: primary q ≥ 0.05   Brown: secondary",5.8)
    for i,r in enumerate(c.itertuples()):
        y=25+8*i
        p.text(394,y+4,r.contrast.replace("_"," "),5.5,align="r")
        p.rect(400,y,76*r.CYP_tests/c.CYP_tests.max(),5,"#e2e6ec")
        p.rect(400,y,76*r.tier_FDR_below_005/c.CYP_tests.max(),5,a.TEAL)
    p.panel(8,153,"b","Fibroblast CYP4B1 effects by cohort and etiology")
    a.dotplot(p,154,170,230,[f"{r.study} · {a.CONTRASTS[r.contrast][0]} {r.region}" for r in hits.itertuples()],
              hits.logFC,(-6,3),[-6,-4,-2,0,2],"Case vs comparator log₂ fold change",
              [f"q = {v:.2g}" for v in hits[q]])
    p.panel(8,395,"c","Selected etiology-specific contrasts")
    rows=highlight[(highlight.contrast.ne("AS_IVS")) |
                   (highlight.cyp_gene.eq("CYP1B1")&highlight.lineage.eq("cardiomyocyte"))|
                   (highlight.cyp_gene.eq("CYP46A1")&highlight.lineage.eq("fibroblast"))]
    for i,r in enumerate(rows.itertuples()):
        y=418+15*i
        p.text(14,y,f"{a.CONTRASTS[r.contrast][0]} · {r.lineage} · {r.cyp_gene}",7.3)
        p.text(311,y,f"log₂ FC {r.logFC:+.2f}",7)
        p.text(412,y,f"q {getattr(r,q):.2g}",7)
    p.finish()


def fig4():
    p=a.Page(4,380)
    con=a.read(a.D/"CYP_Disease_Cross_Cohort_Concordance.tsv")
    for x,letter,title,cohorts,condition in [
        (8,"a","LV DCM cohort effect concordance",DCM,"DCM_LV"),
        (267,"b","RV ARVC cohort effect concordance",ARVC,"ARVC_RV")]:
        p.panel(x,13,letter,title)
        q=ALL[ALL.contrast.isin(cohorts)].pivot(index=["lineage","cyp_gene"],columns="contrast",values="logFC").dropna()
        g=con[con.phenotype_region.eq(condition)&con.tested_studies.eq(len(cohorts))]
        supported=g[g.direction_concordant&g.tier_fdr_below_005.eq(len(cohorts))]
        supported=set(zip(supported.lineage,supported.cyp_gene))
        assert len(q)==(70 if letter=="a" else 55)
        lim=max(5,int(np.ceil(np.abs(q.to_numpy()).max())))
        left=x+25
        scale=lambda z:left+(z+lim)/(2*lim)*185
        ymap=lambda z:150-(z+lim)/(2*lim)*121
        p.line(scale(-lim),ymap(-lim),scale(lim),ymap(lim),"#9aa5b2",.5,True)
        p.line(scale(0),29,scale(0),150,"#cad1d8")
        p.line(left,ymap(0),left+185,ymap(0),"#cad1d8")
        for key,r in q.iterrows():
            for j,other in enumerate(cohorts[1:]):
                p.dot(scale(r[cohorts[0]]),ymap(r[other]),2 if key in supported else 1.25,
                      [a.TEAL,"#b25236"][j] if key in supported else ["#a4bcc6","#d7aca1"][j])
        p.line(left,29,left,150)
        for t in (-lim,0,lim):
            p.line(left-3,ymap(t),left,ymap(t))
            p.text(left-5,ymap(t)+2,t,5.6,align="r")
        axis(p,left,150,185,[-lim,0,lim],-lim,lim,f"{a.CONTRASTS[cohorts[0]][2]} log₂ FC")
        p.text(left,185,"Y: Koenig (teal) · Reichart (brown)" if letter=="a" else "Y: Simonson ARVC",6.6)
        p.text(left,197,f"{len(supported)} of {len(q)} pass direction and FDR criteria",6.4)
        archive(f"Fig4{letter}_study_effect_concordance",q.reset_index())
    p.panel(8,221,"c","Twelve recurrent LV DCM gene–lineage effects")
    selected=con[con.phenotype_region.eq("DCM_LV")&con.direction_concordant&
                 con.tier_fdr_below_005.eq(3)]
    pairs=list(zip(selected.lineage,selected.cyp_gene))
    z=ALL[ALL.contrast.isin(DCM)&ALL.set_index(["lineage","cyp_gene"]).index.isin(pairs)]
    for j,k in enumerate(DCM):
        p.text(175+(j+.5)*310/3,237,a.CONTRASTS[k][2],6.7,align="c")
    a.heat(p,z,pairs,DCM,["lineage","cyp_gene"],"contrast","logFC",None,(175,244,310,105),
           [f"{l.replace('_',' ')} · {g}" for l,g in pairs],5,"Study-specific log₂ fold change")
    p.finish()
    archive("Fig4c_12_recurrent_DCM_pairs",z)


def fig5():
    go=a.read(a.D/"CYP_Disease_Pathways.tsv")
    go_names=sorted(go.pathway.unique())
    assert len(go)==533 and go.GO_BH_across_tier.lt(.05).sum()==203
    p=a.Page(5,445)
    p.panel(8,13,"a","Ranked GO enrichment across selected lineages")
    x,w=177,323
    a.contrast_header(p,x,25,w)
    for i,l in enumerate(LIN):
        for j,g in enumerate(go_names):
            y=52+(i*len(go_names)+j)*8
            p.text(x-5,y+2,f"{SHORT[l]} · {g}",5.8,align="r")
            for k,c in enumerate(a.CONTRASTS):
                row=go[go.lineage.eq(l)&go.pathway.eq(g)&go.contrast.eq(c)]
                if len(row):
                    r=row.iloc[0]
                    p.dot(x+(k+.5)*w/13,y,min(3.6,1.15+abs(r.NES)),a.RED if r.NES>0 else a.BLUE,
                          a.INK if r.GO_BH_across_tier<.05 else None)
    p.text(178,175,"Radius: |NES|    Blue: negative    Red: positive    Outline: q < 0.05",6.3)
    p.panel(8,195,"b","PROGENy footprint shifts")
    pathways=["Estrogen","WNT","TNFa","Hypoxia","TGFb"]
    pairs=[(l,k) for l in LIN for k in pathways]
    a.contrast_header(p,x,206,w)
    s=SIG[SIG.lineage.isin(LIN)&SIG.pathway.isin(pathways)]
    a.heat(p,s,pairs,list(a.CONTRASTS),["lineage","pathway"],"contrast","score_delta","FDR_BH",
           (x,230,w,91),[f"{SHORT[l]} · {k}" for l,k in pairs],2.5,"Within-model SD")
    p.panel(8,355,"c","Inferred cardiomyocyte estrogen response in LV DCM")
    z=SIG[SIG.contrast.isin(DCM)&SIG.lineage.eq("cardiomyocyte")&SIG.pathway.eq("Estrogen")]
    z=z.set_index("contrast").loc[DCM].reset_index()
    a.dotplot(p,124,371,248,z.study,z.score_delta,(-2.5,.25),[-2.5,-2,-1,0],
              "Case vs comparator (within-model SD)",[f"q = {v:.2g}" for v in z.FDR_BH],
              list(zip(z["CI.L"],z["CI.R"])))
    p.finish()
    archive("Fig5a_GO_enrichment",go.drop(columns="leadingEdge"))
    archive("Fig5b_PROGENy_footprints",s)
    archive("Fig5c_estrogen_focal",z)


def fig7():
    counts=a.read(PREV/"new_source_data/Fig7c_cross_cohort_test_summary.tsv")
    val=a.read(ROOT/"source_data/Fig6b_Leuven.tsv")
    bulk=a.read(ROOT/"source_data/Fig6c_bulk.tsv")
    p=a.Page(7,460)
    p.panel(8,13,"a","Donor association and recurrence screens")
    for y,head,number in [(32,"All-donor CYP–program and occupancy tests","109,637"),
                          (47,"Cross-cohort tests across eligible analysis modes",f"{counts.tested.sum():,}"),
                          (62,"Cross-cohort tests with FDR < 0.05","0")]:
        p.text(16,y,head,7.3)
        p.text(476,y,number,7.5,True,align="r")
    p.panel(8,84,"b","CYP2J2 and CYP19A1 vs CYP-excluded estrogen scores")
    f=FOCAL[FOCAL["mode"].eq("all_donors")].sort_values(["gene","study"],ascending=[False,True])
    a.dotplot(p,141,101,239,[f"{r.gene} · {r.study}" for r in f.itertuples()],f.rho,(-1,1),
              [-1,-.5,0,.5,1],"Cohort Spearman correlation",
              [f"n = {r.n_donors}; q = {r.Focal_FDR_BH:.2g}" for r in f.itertuples()])
    p.panel(8,224,"c","Leuven single-nucleus LV validation")
    a.dotplot(p,141,241,239,[f"{r.gene} · {'shared batch' if r.model=='shared_batch_only' else 'batch adjusted'}" for r in val.itertuples()],
              val.log2FC,(-1,6.5),[-1,0,2,4,6],"Donor pseudobulk log₂ fold change",
              [f"n = {r.case_n}/{r.control_n}; q = {r.FDR_BH:.2g}" for r in val.itertuples()],
              list(zip(val.CI_low,val.CI_high)))
    p.panel(8,339,"d","Independent bulk LV RNA validation")
    a.dotplot(p,141,356,239,bulk.gene,bulk.log2_RPKM_plus_0_1_delta,(-1.3,1.3),[-1,-.5,0,.5,1],
              "DCM minus non-failing mean log₂(RPKM + 0.1)",
              [f"n = {r.case_n}/{r.control_n}; q = {r.FDR_BH:.2g}" for r in bulk.itertuples()],
              list(zip(bulk.ci95_low,bulk.ci95_high)))
    p.finish()
    for key,t in (("a",counts),("b",FOCAL),("c",val),("d",bulk)):
        archive(f"Fig7{key}_association_or_validation",t)


def fig8():
    p=a.Page(8,445)
    p.panel(8,13,"a","Baseline expression in four cardiac chambers")
    matrix(p,REG,KEY,"log_cpm",42,9,"log₁₀(1 + median study-median CPM)",3)
    p.panel(8,132,"b","Separate Reichart LV and RV effects")
    v=["Reichart_DCM_LV","Reichart_DCM_RV","Reichart_ARVC_LV","Reichart_ARVC_RV"]
    pairs=[(l,g) for l,g in zip([LIN[0],LIN[1],LIN[2],LIN[1],LIN[1],LIN[2]],KEY)]
    q=ALL[ALL.contrast.isin(v)&ALL.set_index(["lineage","cyp_gene"]).index.isin(pairs)]
    for j,n in enumerate(v):
        p.text(155+(j+.5)*329/4,151,n.removeprefix("Reichart_").replace("_"," "),6.6,align="c")
    a.heat(p,q,pairs,v,["lineage","cyp_gene"],"contrast","logFC","CYP_BH_across_tier",
           (155,157,329,69),[f"{l} · {g}" for l,g in pairs],5,"Case vs comparator log₂ fold change")
    p.panel(8,267,"c","Evidence map and functional hypotheses")
    for i,(name,l1,l2) in enumerate([
        ("Cardiomyocytes","LV DCM: higher CYP2J2 and CYP11A1","RV ARVC: lower CYP2E1 and CYP4F12"),
        ("Fibroblasts","LV DCM: lower CYP4X1 and CYP4Z1","Higher CYP19A1 in evaluable DCM cohorts"),
        ("Endothelial cells","LV DCM: lower CYP1B1, CYP4X1, CYP4Z1","RV ARVC: lower CYP1B1 and CYP4Z1"),
    ]):
        x=12+i*168
        p.rect(x,281,161,48,["#e6f0f7","#e6f3f0","#f2edf4"][i])
        p.text(x+6,294,name,7.6,True)
        p.text(x+6,308,l1,6.2)
        p.text(x+6,321,l2,6.2)
    for y,head,l1,l2,color in [
        (337,"Inferred RNA programs","THRB and mitochondrial ATP synthesis scores were lower in LV DCM.",
         "CYP–estrogen correlations did not survive cross-cohort FDR correction.","#eef0f3"),
        (387,"Functional questions","Does cardiomyocyte CYP2J2 track arachidonate epoxide production?",
         "Does fibroblast CYP19A1 track local aromatase conversion?","#fcf2e3"),
    ]:
        p.rect(12,y,497,43,color)
        p.text(19,y+12,head,7.7,True)
        p.text(19,y+25,l1,6.9)
        p.text(19,y+37,l2,6.9)
    p.finish()
    archive("Fig8a_chamber_baseline",REG[REG.gene.isin(KEY)])
    archive("Fig8b_Reichart_chamber_disease_effects",q)


def extended():
    a.FIG=EXT
    source=PREV/"Cardiac_CYP_Updated_Source_Data.zip"
    with zipfile.ZipFile(source) as z:
        names=[n for n in z.namelist() if n.endswith("_ab_coordinates.tsv.gz")]
        p=a.Page("ED1",400)
        for i,key in enumerate(["Chaffin","Kuppe","AS","Reichart"]):
            name=next(n for n in names if re.search(rf"_{key}_ab_coordinates",n))
            df=pd.read_csv(io.BytesIO(z.read(name)),sep="\t",compression="gzip")
            col="marker_candidate"
            if col not in df:
                col=next(c for c in df if "candidate" in c)
            x=8+260*(i%2); y=13+195*(i//2)
            p.panel(x,y,"abcd"[i],f"{key} candidate lineage assignment")
            umap(p,df,x+4,y+14,144,col,"umap_1","umap_2")
        p.finish()

    p=a.Page("ED2",445)
    p.panel(8,13,"a","Complete four-chamber CYP detection recurrence")
    matrix(p,REG,GENES,"study_strata_fraction_detected_in_majority",40,6.45,
           "Fraction of strata with majority-donor detection",1)
    p.finish()
    archive("ED2_four_chamber_CYP_recurrence",REG)

    p=a.Page("ED3",315)
    c=a.read(FOLDER/"Covariate_CYP_Effects.tsv")
    rows=c[c.contrast.isin(DCM)&((c.lineage.eq("cardiomyocyte")&c.cyp_gene.eq("CYP2J2"))|
         (c.lineage.eq("fibroblast")&c.cyp_gene.eq("CYP19A1")))].copy()
    p.panel(8,13,"a","Covariate sensitivity in cohorts with available complete cases")
    labels=[f"{r.study} · {r.cyp_gene} · {r.mode.replace('_',' ')}" for r in rows.itertuples()]
    a.dotplot(p,218,35,163,labels,rows.logFC,(-1,3.6),[-1,0,1,2,3],
              "Case vs comparator log₂ fold change",[f"q = {r.FDR_BH:.2g}" for r in rows.itertuples()])
    p.finish()
    archive("ED3_covariate_sensitivity",rows)

    p=a.Page("ED4",285)
    overlap=a.read(a.D/"CYP_Disease_Overlap_Sensitivity.tsv")
    rows=overlap[overlap.contrast.isin(DCM)&((overlap.lineage.eq("cardiomyocyte")&overlap.cyp_gene.eq("CYP2J2"))|
         (overlap.lineage.eq("fibroblast")&overlap.cyp_gene.eq("CYP19A1")))].copy()
    p.panel(8,13,"a","Flagged-donor exclusion and the two focal CYP effects")
    p.text(170,35,"All donors",7.3,True)
    p.text(339,35,"Overlap-exclusion model",7.3,True)
    for i,r in enumerate(rows.itertuples()):
        y=57+32*i
        p.text(12,y,f"{r.study} · {r.cyp_gene}",7.3)
        p.text(170,y,f"FC {r.logFC:+.2f}; q {r.CYP_BH_across_tier:.2g}",7.2)
        p.text(170,y+11,f"n = {r.case_donors}/{r.control_donors}",6.4)
        if pd.notna(r.logFC_overlap_screen):
            p.text(339,y,f"FC {r.logFC_overlap_screen:+.2f}; q {r.CYP_BH_across_tier_overlap_screen:.2g}",7.2)
            p.text(339,y+11,f"n = {r.case_donors_overlap_screen:g}/{r.control_donors_overlap_screen:g}",6.4)
        else:
            p.text(339,y,r.status_overlap_screen.replace("_"," "),7)
    p.finish()
    archive("ED4_overlap_sensitivity",rows)

    p=a.Page("ED5",412)
    p.panel(8,13,"a","Expanded PROGENy footprint contrasts across 14 pathways")
    s=SIG[SIG.lineage.isin(LIN)]
    pairs=[(l,k) for l in LIN for k in sorted(s.pathway.unique())]
    a.contrast_header(p,143,25,363)
    a.heat(p,s,pairs,list(a.CONTRASTS),["lineage","pathway"],"contrast","score_delta","FDR_BH",
           (143,49,363,328),[f"{SHORT[l]} · {k}" for l,k in pairs],2.5,"Within-model SD")
    p.finish()
    archive("ED5_all_PROGENy_selected_lineages",s)

    p=a.Page("ED6",286)
    lr=a.read(FOLDER/"Top_Ligand_Receptor.tsv")
    v=lr[lr.contrast.isin(DCM)&lr.population.eq("all_donors")].nsmallest(12,"FDR_BH")
    p.panel(8,13,"a","Selected LIANA-supported ligand–receptor RNA scores in LV DCM")
    labels=[f"{r.study} · {r.ligand}/{r.receptor} · {r.sender.replace('_',' ')} → {r.receiver.replace('_',' ')}"
            for r in v.itertuples()]
    lim=(min(-1,int(np.floor(v.ci_low.min()))),max(1,int(np.ceil(v.ci_high.max()))))
    a.dotplot(p,246,34,133,labels,v.score_delta,lim,list(range(lim[0],lim[1]+1)),
              "Paired RNA score difference",[f"q = {r.FDR_BH:.2g}" for r in v.itertuples()],
              list(zip(v.ci_low,v.ci_high)))
    p.text(10,258,"Scores summarize co-expression. Direction labels indicate sender and receiver lineages.",6.6)
    p.finish()
    archive("ED6_selected_LIANA_scores",v)
    for p in EXT.glob("Figure_ED*.pdf"):
        p.replace(p.with_name(p.name.replace("Figure_ED","Extended_Data_Fig_")))
    a.FIG=MAIN


def main():
    for f in (fig1,fig2,fig3,fig4,fig5):
        f()
    fig6()
    fig7()
    fig8()
    extended()
    for f in list(MAIN.glob("*.pdf"))+list(EXT.glob("*.pdf")):
        subprocess.run(["pdftoppm","-r","600","-png","-singlefile",str(f),str(f.with_suffix(""))],check=True)


def fig6():
    focal=a.read(ROOT/"atlas_manuscript/New_Main_Figures/Figure_8_Source_Data.tsv")
    state=a.read(PREV/"new_source_data/Fig6c_other_cell_state_programs.tsv")
    p=a.Page(6,362)
    studies=["Chaffin","Koenig","Reichart"]
    lineages=["cardiomyocyte","endothelial","fibroblast"]
    for i,(term,title) in enumerate([
        ("THRB","THRB regulon"),
        ("mitochondrial ATP synthesis coupled electron transport","Mitochondrial ATP synthesis"),
    ]):
        x=12+253*i
        p.panel(x,15,"ab"[i],title)
        rows=focal[focal.term.eq(term)].set_index(["lineage","study"]).loc[
            pd.MultiIndex.from_product([lineages,studies])].reset_index()
        labels=[f"{SHORT[l]} / {s}" for l,s in zip(rows.lineage,rows.study)]
        a.dotplot(p,x+77,36,122,labels,rows.score_delta,(-1.4,.15),[-1.4,-1,-.5,0],
                  "Case vs comparator score difference",
                  ["*" if q<.05 else "" for q in rows.FDR_BH],list(zip(rows["CI.L"],rows["CI.R"])))
    p.text(12,219,"* Tier-wide FDR < 0.05. CM, cardiomyocyte. EC, endothelial. FB, fibroblast.",6.4)
    p.panel(8,241,"c","Additional RNA-derived state programs in LV DCM")
    labels={"angiogenesis":"Angiogenesis","response to hypoxia":"Hypoxia",
            "response to oxidative stress":"Oxidative stress"}
    terms=list(state.term.drop_duplicates())
    pairs=[(l,t) for l in lineages for t in terms]
    for j,s in enumerate(studies):
        p.text(172+(j+.5)*326/3,261,s,6.5,align="c")
    a.heat(p,state,pairs,studies,["lineage","term"],"study","score_delta","FDR_BH",
           (172,267,326,64),[f"{SHORT[l]} / {labels.get(t,t)}" for l,t in pairs],
           .5,"Case vs comparator score difference")
    p.finish()


if __name__=="__main__":
    main()
