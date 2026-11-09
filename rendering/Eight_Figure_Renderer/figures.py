from pathlib import Path
import os
import subprocess

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from reportlab.lib.colors import HexColor, Color
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.pdfgen.canvas import Canvas

ROOT = Path(__file__).resolve().parent
PACKAGE = ROOT/"package"
FIG = PACKAGE/"Main_Figures"
DATA = ROOT/"source_data"
LOCAL = Path(os.environ.get("HEARTMAP_FIR_RESULTS","D:/HeartMap_CYP_Fresh/Fir_Results"))
RESULT = LOCAL/"results"
L = RESULT/"CYP_Landscape_2958286"
D = RESULT/"CYP_Disease_2962292"
S = RESULT/"Signaling_Activity_2967196"
R = RESULT/"Reconciled_Reference_2965412"
F = RESULT/"Donor_CYP_Links_2969410/Permutation_Validated"
OLD = ROOT.parent/"heartmap-tac"
CONTRASTS = {
    "AS_IVS": ("AS", "IVS", "AS"),
    "Chaffin_DCM_LV": ("DCM", "LV", "Chaffin"),
    "Koenig_DCM_LV": ("DCM", "LV apex", "Koenig"),
    "Reichart_DCM_LV": ("DCM", "LV", "Reichart"),
    "Chaffin_HCM_LV": ("HCM", "LV", "Chaffin"),
    "Simonson_ICM_LV": ("ICM", "LV", "Simonson"),
    "DoxCM_LV": ("DoxCM", "LV", "Guo"),
    "NICM_LV": ("NICM", "LV", "Guo"),
    "DoxCM_vs_NICM_LV": ("Dox/NICM", "LV", "Guo"),
    "Reichart_DCM_RV": ("DCM", "RV", "Reichart"),
    "Reichart_ARVC_RV": ("ARVC", "RV", "Reichart"),
    "Simonson_ARVC_RV": ("ARVC", "RV", "Simonson"),
    "Reichart_ARVC_LV": ("ARVC", "LV", "Reichart"),
}
CM, FB, EC = "cardiomyocyte", "fibroblast", "endothelial"
BLUE, RED, TEAL, INK, GRAY = "#205fa1", "#b24131", "#117477", "#24282e", "#e0e4e8"
MID = "#fbf7f0"
COLORS = ["#3c78a4","#d17845","#44926f","#91629b","#b49736","#c26482",
          "#4c999b","#7f6b46","#a44e4a","#777e8e","#5992bc","#7a8c46","#815d88","#ad7d5c"]
for font, file in [("Arial","arial.ttf"),("ArialBold","arialbd.ttf")]:
    pdfmetrics.registerFont(TTFont(font,str(Path("C:/Windows/Fonts")/file)))


def read(p):
    return pd.read_csv(p,sep="\t")


def source(name, data):
    data.to_csv(DATA/f"{name}.tsv",sep="\t",index=False,na_rep="NA")


def interpolate(a, b, fraction):
    a, b = HexColor(a), HexColor(b)
    return Color(*[x+(y-x)*fraction for x,y in zip(a.rgb(),b.rgb())])


class Page:
    def __init__(self, num, height):
        self.h = height
        self.pdf = FIG/f"Figure_{num}.pdf"
        self.c = Canvas(str(self.pdf),pagesize=(183*mm,height),pageCompression=1)
        self.c.setTitle(f"Figure {num} | Human cardiac cytochrome P450 atlas")

    def text(self, x, y, text, size=6.3, bold=False, align="l", color=INK):
        self.c.setFillColor(HexColor(color))
        self.c.setFont("ArialBold" if bold else "Arial",max(6,size))
        {"l":self.c.drawString,"r":self.c.drawRightString,"c":self.c.drawCentredString}[align](x,self.h-y,str(text))

    def panel(self, x, y, letter, title):
        self.text(x,y,letter,8.5,True)
        self.text(x+12,y,title,7.4,True)

    def rect(self,x,y,w,h,color):
        self.c.setFillColor(HexColor(color) if isinstance(color,str) else color)
        self.c.rect(x,self.h-y-h,w,h,stroke=0,fill=1)

    def line(self,x,y,xx,yy,color=INK,width=.45,dash=False):
        self.c.setStrokeColor(HexColor(color))
        self.c.setLineWidth(width)
        self.c.setDash(2,2) if dash else self.c.setDash()
        self.c.line(x,self.h-y,xx,self.h-yy)
        self.c.setDash()

    def dot(self,x,y,r,color=INK,stroke=None):
        self.c.setFillColor(HexColor(color))
        if stroke:
            self.c.setStrokeColor(HexColor(stroke))
            self.c.setLineWidth(.45)
        self.c.circle(x,self.h-y,r,fill=1,stroke=bool(stroke))

    def finish(self):
        self.c.showPage()
        self.c.save()


def heat(page, frame, rows, cols, row_keys, col_key, effect, q, box, labels, limit, scale, signed=True):
    x,y,w,h = box
    cw,rh = w/len(cols),h/len(rows)
    assert not frame.duplicated([*row_keys,col_key]).any()
    values = frame.set_index([*row_keys,col_key]).to_dict("index")
    for i,(row,label) in enumerate(zip(rows,labels)):
        page.text(x-4,y+(i+.69)*rh,label,6,align="r")
        for j,col in enumerate(cols):
            cell = values.get((*row,col))
            val = cell[effect] if cell else np.nan
            color = GRAY if pd.isna(val) else interpolate(MID,(BLUE if val<0 else RED) if signed else TEAL,min(1,abs(val)/limit))
            page.rect(x+j*cw,y+i*rh,cw-.36,rh-.36,color)
            if cell and q and pd.notna(cell[q]) and cell[q]<.05:
                page.dot(x+(j+.5)*cw,y+(i+.5)*rh,1.35,INK,"#ffffff")
    sy = y+h+8
    low = -limit if signed else 0
    for j,v in enumerate(np.linspace(low,limit,90)):
        page.rect(x+j*.7,sy,.8,4,interpolate(MID,(BLUE if v<0 else RED) if signed else TEAL,abs(v)/limit))
    page.text(x,sy+12,f"{low:g}",5.5)
    page.text(x+63,sy+12,f"{limit:g}",5.5,align="r")
    page.text(x+71,sy+4,scale,6)
    page.rect(x+w-78,sy,6,5,GRAY)
    page.text(x+w-68,sy+4,"Unavailable",6)


def contrast_header(page,x,y,w):
    short_study={"Chaffin":"Chaf.","Koenig":"Koen.","Reichart":"Reich.","Simonson":"Sim."}
    for i,(case,region,study) in enumerate(CONTRASTS.values()):
        xx=x+(i+.5)*w/len(CONTRASTS)
        for j,s in enumerate(("Dox/N" if case=="Dox/NICM" else case,region,short_study.get(study,study))):
            page.text(xx,y+j*7,s,6,j==0,align="c")


def umap(page,frame,column,x,y):
    labels=sorted(frame[column].unique())
    rgb={k:tuple(bytes.fromhex(COLORS[i%len(COLORS)][1:])) for i,k in enumerate(labels)}
    size=760
    im=Image.new("RGB",(size,size),"white")
    d=ImageDraw.Draw(im)
    a=frame.UMAP1.to_numpy()
    b=frame.UMAP2.to_numpy()
    scale=(size-12)/max(np.ptp(a),np.ptp(b))
    cx=(a.min()+a.max())/2
    cy=(b.min()+b.max())/2
    for xx,yy,k in zip((a-cx)*scale+size/2,(cy-b)*scale+size/2,frame[column]):
        d.ellipse((xx-1,yy-1,xx+1,yy+1),fill=rgb[k])
    page.c.drawImage(ImageReader(im),x,page.h-y-118,118,118)
    for i,k in enumerate(labels):
        yy=y+8+i*8
        page.dot(x+123,yy-2,2,COLORS[i%len(COLORS)])
        page.text(x+129,yy,str(k).replace("_"," "),5.5)
    page.text(x+58,y+126,"UMAP 1",5.8,align="c")
    page.c.saveState()
    page.c.translate(x+3,page.h-y-61)
    page.c.rotate(90)
    page.c.setFont("Arial",5.8)
    page.c.setFillColor(HexColor(INK))
    page.c.drawCentredString(0,0,"UMAP 2")
    page.c.restoreState()


def overview():
    page=Page(1,304)
    qc=read(L/"CYP_Landscape_QC.tsv")
    coverage=read(R/"Cohort_Coverage.tsv")
    emb=read(R/"Reference_Embedding.tsv.gz")
    assert len(emb)==33720
    emb=emb[["study","analysis_lineage","UMAP1","UMAP2"]]
    page.panel(8,12,"a","Nuclei and lineage support across source cohorts")
    page.panel(263,12,"b","Reference sampling")
    for i,r in enumerate(qc.itertuples()):
        y=30+i*10
        page.text(10,y+4,r.study.replace("_"," "),6.1)
        page.rect(80,y-1,r.clustered_nuclei/6500,6,"#d9dfe6")
        page.rect(80,y-1,r.marker_supported_nuclei/6500,6,TEAL)
        page.text(246,y+4,f"{r.supported_clusters}/{r.clusters}",5.8,align="r")
        ss=coverage.loc[coverage.study.eq(r.study),"sampled_nuclei"].iat[0]
        page.rect(267,y-1,ss/48,6,BLUE)
        page.text(467,y+4,f"{ss:,}",5.8,align="r")
    page.rect(80,143,7,5,TEAL); page.text(91,148,"Lineage supported",5.7)
    page.rect(153,143,7,5,"#d9dfe6"); page.text(164,148,"Held",5.7)
    page.text(246,148,"Clusters passing",5.5,align="r")
    page.panel(8,167,"c","Harmony reference by broad lineage")
    page.panel(263,167,"d","Same embedding by source cohort")
    umap(page,emb,"analysis_lineage",14,172)
    umap(page,emb,"study",269,172)
    page.finish()
    source("Fig1a_lineage_coverage",qc)
    source("Fig1b_reference_coverage",coverage)
    source("Fig1cd_reference_coordinates",emb)


def baseline():
    page=Page(2,358)
    df=read(L/"CYP_Control_Regional_Recurrence.tsv")
    genes=["CYP1B1","CYP2J2","CYP2E1","CYP4B1","CYP4F2","CYP4X1","CYP4Z1","CYP5A1","CYP7B1","CYP8A1","CYP11A1","CYP19A1","CYP46A1"]
    df=df[df.gene.isin(genes) & df.region.isin(["LV","RV","LA","RA"]) & df.lineage.isin([CM,FB,EC])].copy()
    df["column"]=df.region+"/"+df.lineage
    df["log_cpm"]=np.log10(1+df.median_of_study_median_cpm)
    cols=[f"{r}/{l}" for r in ["LV","RV","LA","RA"] for l in [CM,FB,EC]]
    names=[{"CYP5A1":"CYP5A1 / TBXAS1","CYP8A1":"CYP8A1 / PTGIS"}.get(g,g) for g in genes]
    for letter,title,y,field,limit,scale in [
        ("a","Transcript abundance",41,"log_cpm",3,"log₁₀(1 + median of study median CPM)"),
        ("b","Recurrence across study strata",217,"study_strata_fraction_detected_in_majority",1,"Fraction of strata with majority-donor detection")
    ]:
        page.panel(8,y-30,letter,title)
        for i,c in enumerate(cols):
            r,l=c.split("/")
            xx=84+(i+.5)*420/12
            page.text(xx,y-16,r,6.2,True,align="c")
            page.text(xx,y-7,{CM:"CM",FB:"FB",EC:"EC"}[l],6,align="c")
        heat(page,df,[(g,) for g in genes],cols,["gene"],"column",field,None,(84,y,420,110),
             names,limit,scale,signed=False)
    page.finish()
    source("Fig2ab_baseline",df)


def contrasts(effect):
    page=Page(3,307)
    coverage=read(D/"CYP_Disease_Contrast_Summary.tsv").set_index("contrast").loc[list(CONTRASTS)].reset_index()
    page.panel(8,12,"a","CYP tests across region-matched comparisons")
    x,w=101,405
    contrast_header(page,x,27,w)
    for i,r in enumerate(coverage.itertuples()):
        xx=x+i*w/len(CONTRASTS)+8
        page.rect(xx,83-r.CYP_tests*.25,15,r.CYP_tests*.25,"#d9dfe6")
        page.rect(xx,83-r.tier_FDR_below_005*.25,15,r.tier_FDR_below_005*.25,TEAL)
        page.text(xx+7.5,92,str(r.CYP_tests),5.5,align="c")
    page.text(9,63,"Tests (grey)",6)
    page.text(9,73,"FDR < 0.05 (teal)",6)
    page.panel(8,110,"b","Direction and magnitude of selected CYP changes")
    pairs=[(FB,g) for g in ["CYP4B1","CYP4X1","CYP4Z1","CYP19A1","CYP46A1","CYP8A1"]]
    pairs += [(CM,g) for g in ["CYP1B1","CYP2J2","CYP4B1","CYP11A1","CYP19A1","CYP2E1"]]
    pairs += [(EC,g) for g in ["CYP4B1","CYP4X1","CYP4Z1","CYP8A1"]]+[("myeloid","CYP5A1")]
    contrast_header(page,x,124,w)
    subset=effect[effect.set_index(["lineage","cyp_gene"]).index.isin(pairs)]
    heat(page,subset,pairs,list(CONTRASTS),["lineage","cyp_gene"],"contrast","logFC","CYP_BH_across_tier",
         (x,147,w,130),[f"{ {CM:'CM',FB:'FB',EC:'EC','myeloid':'Myeloid'}[l]} · {g}" for l,g in pairs],5,
         "Case versus comparator log₂ fold change")
    page.finish()
    source("Fig3a_test_coverage",coverage)
    source("Fig3b_CYP_effects",subset)


def pathways():
    page=Page(4,366)
    go=read(D/"CYP_Disease_Pathways.tsv")
    sig=read(S/"Signaling_Activity_Effects.tsv").query("population == 'all_donors'")
    names={"arachidonate metabolic process":"Arachidonate metabolism","icosanoid metabolic process":"Eicosanoid metabolism",
           "inflammatory response":"Inflammation","response to oxidative stress":"Oxidative stress","apoptotic process":"Apoptosis"}
    for letter,title,y,frame,pairs,field,q,limit,scale in [
        ("a","Ranked gene set enrichment",46,go,[(l,p) for l in (CM,FB,EC) for p in names],"NES","GO_BH_across_tier",2.5,"Normalized enrichment score"),
        ("b","Inferred PROGENy signaling",226,sig,[(l,p) for l in (CM,FB) for p in ("Estrogen","WNT","Androgen","TNFa","TGFb","NFkB")],"score_delta","FDR_BH",2,"Case versus comparator score change (SD)")
    ]:
        page.panel(8,y-34,letter,title)
        contrast_header(page,113,y-22,393)
        subset=frame[frame.set_index(["lineage","pathway"]).index.isin(pairs)]
        heat(page,subset,pairs,list(CONTRASTS),["lineage","pathway"],"contrast",field,q,(113,y,393,112),
             [f"{ {CM:'CM',FB:'FB',EC:'EC'}[l]} · {names.get(p,p)}" for l,p in pairs],limit,scale)
        source(f"Fig4{letter}_pathways",subset)
    page.finish()


def dotplot(page,x,y,w,labels,values,limits,ticks,axis,meta,intervals=None):
    scale=lambda v:x+w*(v-limits[0])/(limits[1]-limits[0])
    dy=14
    height=dy*len(labels)
    page.line(scale(0),y-6,scale(0),y+height-4,"#818893",.5,True)
    for i,(name,val,tag) in enumerate(zip(labels,values,meta)):
        yy=y+i*dy
        page.text(x-6,yy+2,name,6.1,align="r")
        page.text(x+w+10,yy+2,tag,6)
        if pd.notna(val):
            if intervals is not None:
                lo,hi=intervals[i]
                page.line(scale(lo),yy,scale(hi),yy,TEAL,.75)
                page.line(scale(lo),yy-2,scale(lo),yy+2,TEAL,.65)
                page.line(scale(hi),yy-2,scale(hi),yy+2,TEAL,.65)
            page.dot(scale(val),yy,2.4,TEAL)
    end=y+height-3
    page.line(x,end,x+w,end)
    for t in ticks:
        page.line(scale(t),end,scale(t),end+3)
        page.text(scale(t),end+10,f"{t:g}",5.8,align="c")
    page.text(x+w/2,end+20,axis,6.1,align="c")


def donor():
    page=Page(5,189)
    f=read(F/"Prespecified_DCM_Associations.tsv")
    for i,(gene,lineage) in enumerate([("CYP2J2",CM),("CYP19A1",FB)]):
        y=12+i*87
        page.panel(8,y,"ab"[i],f"{gene} in {lineage}s")
        df=f[f.gene.eq(gene)&f["mode"].eq("all_donors")].sort_values("study")
        dotplot(page,90,y+16,300,df.study,df.rho,(-1,1),[-1,-.5,0,.5,1],"Spearman correlation",
                [f"n = {r.n_donors}; q = {r.Focal_FDR_BH:.2f}" for r in df.itertuples()])
    page.finish()
    source("Fig5ab_focal_associations",f)


def validation(effect):
    page=Page(6,330)
    subset=effect[effect.contrast.isin(["Chaffin_DCM_LV","Koenig_DCM_LV","Reichart_DCM_LV"]) &
                  effect.set_index(["lineage","cyp_gene"]).index.isin([(CM,"CYP2J2"),(FB,"CYP19A1")])]
    idx=pd.MultiIndex.from_product([["CYP2J2","CYP19A1"],["Chaffin","Koenig","Reichart"]],names=["cyp_gene","study"])
    a=subset.set_index(["cyp_gene","study"]).reindex(idx).reset_index()
    page.panel(8,12,"a","Discovery estimates in LV DCM")
    dotplot(page,122,31,248,[f"{r.cyp_gene} · {r.study}" for r in a.itertuples()],a.logFC,(-1,3),
            [-1,0,1,2,3],"Donor pseudobulk log₂ fold change",
            [f"n = {int(r.case_donors)}/{int(r.control_donors)}; q = {r.CYP_BH_across_tier:.2g}" if pd.notna(r.logFC)
             else "No eligible estimate" for r in a.itertuples()])
    v=read(LOCAL/"external_validation/GSE298023/Validation_3016861/Leuven_DCM_CYP_Validation.tsv")
    v=v[v.model.isin(["batch_adjusted","shared_batch_only"])].sort_values(["gene","model"],ascending=[False,True])
    assert len(v)==4
    page.panel(8,145,"b","Leuven single-nucleus validation")
    dotplot(page,122,163,248,[f"{r.gene} · {'adjusted' if r.model=='batch_adjusted' else 'shared batch'}"
                              for r in v.itertuples()],v.log2FC,(-1,6),[-1,0,1,2,3,4,5,6],
            "Log₂ fold change (95% CI)",[f"n = {r.case_n}/{r.control_n}; q = {r.FDR_BH:.2g}" for r in v.itertuples()],
            list(zip(v.CI_low,v.CI_high)))
    bulk_path=ROOT/"Bulk_DCM_External_Screen.tsv"
    b=read(bulk_path if bulk_path.exists() else OLD/"External/Bulk_DCM_External_Screen.tsv")
    page.panel(8,248,"c","Independent bulk LV RNA sequencing")
    dotplot(page,122,265,248,b.gene,b.log2_RPKM_plus_0_1_delta,(-1.25,1.25),[-1,-.5,0,.5,1],
            "Difference in mean log₂(RPKM + 0.1) (95% CI)",
            [f"n = {r.case_n}/{r.control_n}; q = {r.FDR_BH:.2g}" for r in b.itertuples()],
            list(zip(b.ci95_low,b.ci95_high)))
    page.finish()
    source("Fig6a_discovery_DCM",a)
    source("Fig6b_Leuven",v)
    source("Fig6c_bulk",b)


def build():
    FIG.mkdir(parents=True,exist_ok=True)
    DATA.mkdir(exist_ok=True)
    effect=read(D/"CYP_Disease_Effects.tsv").query("population == 'all_donors'")
    assert len(effect)==1234
    overview(); baseline(); contrasts(effect); pathways(); donor(); validation(effect)
    for pdf in sorted(FIG.glob("Figure_*.pdf")):
        subprocess.run(["pdftoppm","-r","300","-png","-singlefile",str(pdf),str(pdf.with_suffix(""))],check=True)


if __name__=="__main__":
    build()
