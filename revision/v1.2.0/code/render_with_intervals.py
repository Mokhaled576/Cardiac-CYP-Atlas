from pathlib import Path
import argparse, io, math, re, zipfile
import numpy as np
import pandas as pd
from reportlab.graphics.shapes import Drawing, Rect, Line, Circle, String, Polygon, Group, Path as ShapePath
from reportlab.graphics import renderPDF, renderSVG
from reportlab.lib.colors import HexColor
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

parser=argparse.ArgumentParser(description="Render the 8 main figures and 6 Extended Data figures from archived Fresh CYP results.")
parser.add_argument("output", type=Path)
parser.add_argument("--source-data", type=Path, required=True, help="Path to Figure_Source_Data.zip from the submission package")
parser.add_argument("--arial-dir", type=Path, default=Path("C:/Windows/Fonts"), help="Directory containing arial.ttf and arialbd.ttf")
names=[f"fig{i}" for i in range(1,9)]+[f"ed{i}" for i in range(1,7)]
parser.add_argument("--figures", nargs="+", choices=names, help="Optional subset of figures")
parser.add_argument("--ci",type=Path,default=Path(__file__).resolve().parents[1]/"refit"/"CYP_Discovery_With_CI.tsv")
args=parser.parse_args()
package=args.output
z=zipfile.ZipFile(args.source_data)
package.mkdir(parents=True,exist_ok=True)
def get(name,prior=False):
    prefix='Previous_Display_and_Atlas/' if prior else 'Eight_Figure_Displays/'
    return pd.read_csv(io.BytesIO(z.read(prefix+name+'.tsv')),sep='\t',low_memory=False)
pdfmetrics.registerFont(TTFont("Arial",str(args.arial_dir/"arial.ttf")))
pdfmetrics.registerFont(TTFont("Arial-Bold",str(args.arial_dir/"arialbd.ttf")))

ink='#172532'; muted='#435668'; grid='#DCE5EB'; bg='#F2F7FA'; grey='#78848E'
cm='cardiomyocyte'; fb='fibroblast'; ec='endothelial'
lin={cm:'#0072B2',fb:'#D55E00',ec:'#008467','pericyte':'#7556A4','smooth_muscle':'#CC79A7',
     'myeloid':'#00838B','T_NK':'#253B75','B_cell':'#253B75','adipocyte':'#80632C',
     'mast_cell':'#9A4925','neural':'#586D60','mesothelial':'#A58042','lymphatic_endothelial':'#339987','unresolved':'#8B949C'}
ls={cm:'CM',fb:'FB',ec:'EC','pericyte':'Pericyte','smooth_muscle':'vSMC','T_NK':'T / NK',
    'B_cell':'B cell','lymphatic_endothelial':'Lymph. EC','mast_cell':'Mast','neural':'Neural',
    'adipocyte':'Adipocyte','mesothelial':'Mesothelial','myeloid':'Myeloid'}
disease={'DCM':'#0072B2','ARVC':'#D55E00','ICM':'#A6538C','HCM':'#007D5B','DoxCM':'#8A6B00',
         'NICM':'#62499A','AS':'#3E86A8','Control':'#666666'}
studies={'Chaffin':'Chaffin','Koenig':'Koenig','Reichart':'Reichart','AS':'Lindman',
         'DoxCM':'Guo','Litvinukova':'Litvinukova','Simonson_ICM':'Simonson ICM','Simonson_ARVC':'Simonson ARVC'}
dc=['Chaffin_DCM_LV','Koenig_DCM_LV','Reichart_DCM_LV']
ar=['Reichart_ARVC_RV','Simonson_ARVC_RV']
dn=['Chaffin','Koenig','Reichart']
regs=['LV','RV','LA','RA']
q='CYP_BH_across_tier'
ch={
'AS_IVS':'Lindman AS', 'Chaffin_DCM_LV':'Chaffin DCM','Chaffin_HCM_LV':'Chaffin HCM',
'Koenig_DCM_LV':'Koenig DCM','Reichart_DCM_LV':'Reichart DCM LV','Reichart_DCM_RV':'Reichart DCM RV',
'Reichart_ARVC_RV':'Reichart ARVC RV','Reichart_ARVC_LV':'Reichart ARVC LV','Simonson_ARVC_RV':'Simonson ARVC',
'Simonson_ICM_LV':'Simonson ICM','DoxCM_LV':'Guo DoxCM','NICM_LV':'Guo NICM','DoxCM_vs_NICM_LV':'Guo Dox / NICM'}
kind={k:next((d for d in ('ARVC','DoxCM','NICM','DCM','ICM','HCM','AS') if d in k),'Control') for k in ch}
short={k:v.replace('Reichart','Reich.').replace('Simonson','Sim.').replace('Chaffin','Chaf.') for k,v in ch.items()}
blue='#2865A6'; red='#C24B2A'

def direction(value, supported):
    return (blue if value<0 else red) if supported else '#89959F'
def blend(a,b,t):
    a,b=a.lstrip('#'),b.lstrip('#')
    t=max(0,min(1,float(t)))
    return '#'+''.join(f'{round((1-t)*int(a[i:i+2],16)+t*int(b[i:i+2],16)):02x}' for i in (0,2,4))
def signed(x,v):
    return blend('#F7FAFC',blue if x<0 else red,abs(x)/v) if pd.notna(x) else '#E2E7EC'
def logfmt(x): return f'{x:g}'
def noteq(x): return f'{x:.2g}'

class Page:
    W=510.236; H=415
    def __init__(self,n,title,rows,heights,notes):
        self.n=n
        self.H=12+sum(heights)+7*(len(heights)-1)
        self.D=Drawing(self.W,self.H)
        self.box={}
        self.boundaries=[]
        self.rect(0,0,self.W,self.H,'#FFFFFF')
        y=6
        for (letters,ratios),h in zip(rows,heights,strict=True):
            x=11; available=self.W-22-10*(len(letters)-1)
            for letter,ratio in zip(letters,ratios,strict=True):
                w=available*ratio/sum(ratios);self.box[letter]=(x,y,w,h);x+=w+10
            y+=h+7
    def txt(self,x,y,s,size=6.6,bold=False,anchor='start',fill=ink):
        font='Arial-Bold' if bold else 'Arial'
        self.D.add(String(x,self.H-y,str(s),fontName=font,fontSize=size,textAnchor=anchor,fillColor=HexColor(fill)))
    def rect(self,x,y,w,h,c,stroke=None,sw=.35):
        self.D.add(Rect(x,self.H-y-h,w,h,fillColor=HexColor(c) if c else None,
                        strokeColor=HexColor(stroke) if stroke else None,strokeWidth=sw))
    def line(self,x,y,u,v,c=ink,sw=.5):
        if c==grid:sw=.25
        self.D.add(Line(x,self.H-y,u,self.H-v,strokeColor=HexColor(c),strokeWidth=sw))
    def dot(self,x,y,r,c,filled=True):
        self.D.add(Circle(x,self.H-y,r,fillColor=HexColor(c if filled else '#FFFFFF'),strokeColor=HexColor(c),strokeWidth=.45))
    def cloud(self,x,y,rad,c):
        self.D.add(Circle(x,self.H-y,rad,fillColor=HexColor(c),strokeColor=None))
    def clouds(self,xy,r,c):
        path=ShapePath(fillColor=HexColor(c),strokeColor=None)
        for x,y in xy:
            y=self.H-y
            path.moveTo(x-r,y-r);path.lineTo(x+r,y-r)
            path.lineTo(x+r,y+r);path.lineTo(x-r,y+r);path.closePath()
        self.D.add(path)
    def arrow(self,x,y,u,c=muted):
        self.line(x,y,u-4,y,c,.7)
        self.D.add(Polygon([u,self.H-y,u-5,self.H-y+2.5,u-5,self.H-y-2.5],
                           fillColor=HexColor(c),strokeColor=None))
    def panel(self,k,title):
        x,y,w,h=self.box[k]
        self.boundaries.append((k,len(self.D.contents)))
        self.txt(x,y+9,k,9.5,True)
        self.txt(x+13,y+9,title,7,True)
        return x,y+14,w,h-14
    def flow(self,title,cards):
        x,y,w,h=self.panel('A',title)
        g=14;cw=(w-g*(len(cards)-1))/len(cards)
        for i,(label,detail,color) in enumerate(cards):
            xx=x+i*(cw+g)
            self.rect(xx,y+2,cw,h-6,bg)
            self.rect(xx,y+2,2.5,h-6,color)
            self.txt(xx+8,y+13,label,7,True)
            self.txt(xx+8,y+23,detail,6.3,fill=muted)
            if i<len(cards)-1:self.arrow(xx+cw+1,y+13,xx+cw+g-1)
    def save(self):
        assert list(self.box)==list('ABCDEFGHI'[:len(self.box)])
        assert 8<=len(self.box)<=12
        name=dest/f'Figure_{self.n:02d}'
        renderPDF.drawToFile(self.D,str(name.with_suffix('.pdf')))
        renderSVG.drawToFile(self.D,str(name.with_suffix('.svg')))

def layout(n,title,notes,rows=None,heights=None):
    rows=rows or [('A',[1]),('BCD',[1,1,1]),('EF',[1,1]),('GH',[1,1])]
    return Page(n,title,rows,heights or [48,105,115,112],notes)
def axes_x(p,x,y,w,lim,ticks,label=None):
    lo,hi=lim
    p.line(x,y,x+w,y,muted)
    for t in ticks:
        xx=x+(t-lo)/(hi-lo)*w
        p.line(xx,y,xx,y+2.5,muted)
        p.txt(xx,y+10,logfmt(t),6.1,anchor='middle',fill=muted)
    if label:p.txt(x+w/2,y+19,label,6.5,anchor='middle',fill=muted)
def tickset(lo,hi):
    return list(np.linspace(lo,hi,5).round(4))
def estimates(p,k,title,records,lim,ticks,label,lwidth=.42):
    # (row label, estimate, lower, upper, color, supported)
    x,y,w,h=p.panel(k,title)
    left=x+w*lwidth; width=w*(1-lwidth)-7; top=y+2; bottom=y+h-19
    def X(t):return left+(float(t)-lim[0])/(lim[1]-lim[0])*width
    if lim[0]<0<lim[1]:p.line(X(0),top,X(0),bottom,grey,.58)
    dy=(bottom-top)/max(1,len(records))
    for i,(name,value,lo,hi,color,hit) in enumerate(records):
        yy=top+dy*(i+.5)
        if not hit:color='#89959F'
        p.txt(x+1,yy+2,name,6.25)
        p.line(left,yy,left+width,yy,grid,.34)
        if pd.isna(value):
            p.txt(left+width/2,yy+2,'Unavailable',6.1,anchor='middle',fill=muted)
            continue
        assert lim[0]<=value<=lim[1],(title,value,lim)
        if lo is not None and pd.notna(lo):
            assert lim[0]<=lo<=hi<=lim[1],(title,lo,hi,lim)
            p.line(X(lo),yy,X(hi),yy,color,1)
            for t in (lo,hi):p.line(X(t),yy-2,X(t),yy+2,color,.7)
        p.dot(X(value),yy,2.3,color,bool(hit))
    axes_x(p,left,bottom,width,lim,ticks,label)
def bars(p,k,title,records,lim,ticks,label,lwidth=.4):
    x,y,w,h=p.panel(k,title)
    left=x+w*lwidth; width=w*(1-lwidth)-7; top=y+1;bottom=y+h-19
    dy=(bottom-top)/len(records)
    for i,(name,val,col) in enumerate(records):
        yy=top+dy*(i+.5)
        p.txt(x,yy+2,name,6.2)
        p.rect(left,yy-dy*.24,val/lim*width,dy*.48,col)
    axes_x(p,left,bottom,width,(0,lim),ticks,label)
def tiles(p,k,title,rows,cols,lookup,limit,*,focal=False,lwidth=.43,labels=None,mask_nonsignificant=True):
    x,y,w,h=p.panel(k,title)
    xx=x+w*lwidth;cw=(w*(1-lwidth)-4)/len(cols)
    top=y+10;dy=(h-12)/len(rows)
    for j,col in enumerate(cols):p.txt(xx+(j+.5)*cw,y+7,(labels or {}).get(col,col),6.15,anchor='middle',fill=muted)
    for i,(key,name) in enumerate(rows):
        yy=top+i*dy;p.txt(x+1,yy+dy*.5+2,name,min(6.35,dy*.76))
        for j,col in enumerate(cols):
            value,hit=lookup.get((key,col),(np.nan,False)); cc=signed(value,limit) if not focal else blend('#F3F8FB',blue,value/limit) if pd.notna(value) else '#E2E7EC'
            if mask_nonsignificant and pd.notna(value) and not hit:cc='#E7EBEE'
            left=xx+j*cw;p.rect(left,yy,cw-.65,dy-.7,cc)
            if pd.notna(value):
                s=(f'{value:+.2f}' if abs(value)<0.10 else f'{value:+.1f}') if not focal else f'{value:.1f}'
                if cw>=15:p.txt(left+(cw-.6)/2,yy+dy*.52+2,s+('*' if hit else ''),min(6.0,dy*.70),anchor='middle',fill='#FFFFFF' if abs(value)/limit>.72 and (hit or not mask_nonsignificant) else ink)
            else:p.txt(left+cw/2,yy+dy*.52+2,'NA',5.5,anchor='middle',fill=muted)
def marker(p,x,y,r,c,shape):
    if shape=='diamond':
        p.D.add(Polygon([x,p.H-y+r,x+r,p.H-y,x,p.H-y-r,x-r,p.H-y],
                        fillColor=HexColor(c),strokeColor=None))
    elif shape=='square':p.rect(x-r,y-r,2*r,2*r,c)
    else:p.cloud(x,y,r,c)

def dots(p,k,title,rs,limx,limy,xticks,yticks,xlab,ylab,diagonal=False):
    x,y,w,h=p.panel(k,title)
    xx=x+26;yy=y+3;ww=w-34;hh=h-25
    def X(t):return xx+(float(t)-limx[0])/(limx[1]-limx[0])*ww
    def Y(t):return yy+hh-(float(t)-limy[0])/(limy[1]-limy[0])*hh
    for t in yticks:
        p.line(xx,Y(t),xx+ww,Y(t),grid,.4);p.txt(xx-3,Y(t)+2,logfmt(t),6.0,anchor='end',fill=muted)
    p.line(xx,yy,xx,yy+hh,muted)
    if limx[0]<0<limx[1]:p.line(X(0),yy,X(0),yy+hh,grey,.5)
    if limy[0]<0<limy[1]:p.line(xx,Y(0),xx+ww,Y(0),grey,.5)
    if ylab=='-log10(q)':
        threshold=Y(-math.log10(.05))
        for offset in range(0,int(ww)-2,6):p.line(xx+offset,threshold,xx+offset+3,threshold,grey,.55)
    if diagonal:p.line(X(limx[0]),Y(limy[0]),X(limx[1]),Y(limy[1]),grey,.5)
    for a,b,c,hit,*shape in rs:
        if not(pd.notna(a) and pd.notna(b)):continue
        a=float(np.clip(a,*limx));b=float(np.clip(b,*limy))
        marker(p,X(a),Y(b),1.8 if hit else 1.15,c,shape[0] if shape else 'circle')
    axes_x(p,xx,yy+hh,ww,limx,xticks,xlab)
    if ylab:p.txt(xx+2,yy+7,ylab,6,fill=muted)
def forestrow(name,r,col,field='score_delta',lo='CI.L',hi='CI.R',qcol='FDR_BH'):
    return (name,float(r[field]),float(r[lo]),float(r[hi]),col,float(r[qcol])<.05)
def pointrow(name,r,field='logFC',qcol=q,col=None):
    return (name,float(r[field]),None,None,col or lin.get(r['lineage'],grey),float(r[qcol])<.05)
def lookup(tab,rowfields,column,value,fdr=None):
    return {(tuple(r[f] for f in rowfields),r[column]):(float(r[value]), bool(r[fdr]<.05) if fdr else False) for _,r in tab.iterrows()}
def keyrow(tab,keys):
    mask=pd.Series(True,index=tab.index)
    for k,v in keys.items():mask &= tab[k].eq(v)
    part=tab.loc[mask]
    if part.empty:return None
    assert len(part)==1,(keys,len(part))
    return part.iloc[0]
allfx=get('Fig3a_all_pseudobulk_tests')
ci=pd.read_csv(args.ci,sep='\t')
keys=['contrast','lineage','cyp_gene']
intervals=ci[keys+['logFC_unshrunk','CI_95_low','CI_95_high']]
assert len(intervals)==1234 and not intervals.duplicated(keys).any()
def with_ci(d):
    result=d.merge(intervals,on=keys,how='left',validate='one_to_one')
    assert result[['logFC_unshrunk','CI_95_low','CI_95_high']].notna().all().all()
    return result
allfx=with_ci(allfx)
def ci_row(label,r):
    return forestrow(label,r,direction(r.logFC,r[q]<.05),'logFC_unshrunk','CI_95_low','CI_95_high',q)
def ci_axis(records):
    low=min([0]+[r[2] for r in records if r[2] is not None])
    high=max([0]+[r[3] for r in records if r[3] is not None])
    step=max(1,math.ceil((high-low)/5))
    ends=(step*(math.floor(low/step)-1),step*(math.ceil(high/step)+1))
    return ends,list(range(ends[0],ends[1]+1,step))

base=get('Fig2a_b_regional_baseline')
assert len(allfx)==1234 and allfx.contrast.nunique()==13 and base.gene.nunique()==57

def fig1():
    p=layout(1,'Cohort assembly and lineage review',[
        'B–C: colored bars show retained counts and light bars show total counts. H–I: sampled Harmony reference.',
        'The 3,200,270 denominator counts clustered, quality-screened nuclei. Donors can recur across source studies.'],
        [('A',[1]),('BCD',[1.15,1.05,.8]),('EFG',[1,1,1]),('HI',[1,1])],[48,111,91,130])
    p.flow('From source cohorts to the cardiac atlas',[
        ('11 source cohorts','3,200,270 clustered nuclei',blue),('Lineage review','172 of 239 clusters supported',lin[fb]),
        ('2,778,472 retained nuclei','7 cohorts supply 13 matched contrasts',lin[ec])])
    qc=get('Fig1a_b_cohort_retention').sort_values('clustered_nuclei',ascending=False)
    rc=get('Fig1b_reference_coverage',True).set_index('study')
    names=lambda s:studies.get(s,s)
    for k,title,field,lim,ticks,axis,c in [
        ('B','Clustered and retained nuclei','clustered_nuclei',900,[0,300,600,900],'Nuclei (thousands)',blue),
        ('C','Cluster review by cohort','clusters',30,[0,10,20,30],'Clusters',lin[ec]),
        ('D','Lineage retention','retained',100,[0,50,100],'Retained nuclei (%)',lin[fb])]:
        x,y,w,h=p.panel(k,title)
        left=x+w*.44;ww=w*.56-5;bottom=y+h-19;dy=(bottom-y)/len(qc)
        for i,r in enumerate(qc.itertuples(index=False)):
            yy=y+i*dy+dy*.5
            p.txt(x,yy+1.9,names(r.study),6.05)
            v={'B':r.clustered_nuclei/1000,'C':r.clusters,'D':100*r.marker_supported_nuclei/r.clustered_nuclei}[k]
            p.rect(left,yy-2.3,ww*v/lim,4.6,c if k=='D' else '#E1EAF0')
            if k in 'BC':
                hit=r.marker_supported_nuclei/1000 if k=='B' else r.supported_clusters
                p.rect(left,yy-2.3,ww*hit/lim,4.6,c)
        axes_x(p,left,bottom,ww,(0,lim),ticks,axis)
    cov=[(names(s),int(rc.loc[s,'sampled_nuclei']),blue) for s in qc.study]
    bars(p,'E','Sampled reference by cohort',cov,9000,[0,3000,6000,9000],'Nuclei',.42)
    x,y,w,h=p.panel('F','Duplicate-library reconciliation')
    for i,(label,detail) in enumerate([('23 biospecimen matches','Litvinukova and Reichart'),('81,753 nuclei removed','From the Reichart reference input'),('33,720 sampled nuclei','After reconciliation and subsampling')]):
        yy=y+i*22
        p.txt(x,yy+8,label,8.2,True,fill=blue);p.txt(x,yy+17,detail,6.5,fill=muted)
    x,y,w,h=p.panel('G','Matched model coverage')
    for i,(label,detail) in enumerate([('7 modeling cohorts','4 cohorts contribute atlas data only'),('13 contrasts','10 primary and 3 secondary'),('Donor-level inference','Eligibility varies by lineage model')]):
        yy=y+i*22
        p.txt(x,yy+8,label,8.2,True,fill=lin[ec]);p.txt(x,yy+17,detail,6.5,fill=muted)
    ref=get('Fig1c_d_reference_coordinates')
    assert len(ref)==33720
    colors=['#0072B2','#D55E00','#008467','#7256A2','#A56C00','#4975A0','#AF4F77','#626E77','#008C9C','#896324','#778139']
    assert len(ref.study.unique())==len(colors)
    for k,title,col,palette in [('H','Harmony reference by lineage','analysis_lineage',lin),
                                ('I','Harmony reference by study','study',dict(zip(sorted(ref.study.unique()),colors,strict=True)))]:
        x,y,w,h=p.panel(k,title)
        minx,maxx=ref.UMAP1.min(),ref.UMAP1.max(); miny,maxy=ref.UMAP2.min(),ref.UMAP2.max()
        scale=min((w*.64-8)/(maxx-minx),(h-11)/(maxy-miny))
        xx=x+5+(w*.62-(maxx-minx)*scale)/2; yy=y+1
        for group,v in ref.groupby(col,sort=True):
            color=palette[group]
            p.clouds(zip(xx+(v.UMAP1-minx)*scale,yy+(maxy-v.UMAP2)*scale,strict=True),.32,color)
        for i,group in enumerate(sorted(ref[col].unique())):
            px=x+w*.67;py=y+5+i*7.2
            p.dot(px,py-1.4,1.6,palette[group])
            p.txt(px+4.5,py,(ls if col=='analysis_lineage' else studies).get(group,group),6.0)
        p.txt(x+w*.33,y+h-1,'UMAP 1 / UMAP 2',6.4,anchor='middle',fill=muted)
    p.save()

def fig2():
    p=layout(2,'Baseline CYP abundance across cardiac lineages',[
        'Abundance is the median of study median CPM. Detection requires RNA in at least half of eligible donors in a stratum.',
        'Panel B color scale is log10(1 + CPM), white at 0 and blue at 2.6. Gray denotes unavailable data.'],
        [('A',[1]),('B',[1]),('CDE',[1,1,1]),('FGH',[1,1,1])],[48,133,101,98])
    p.flow('Quantifying baseline expression',[
        ('Non-failing donor profiles','Original counts in region and lineage',blue),('Donor-level CPM','One estimate per donor profile',lin[fb]),
        ('Study summaries','Abundance and detection recurrence',lin[ec])])
    x,y,w,h=p.panel('B','Four-chamber baseline for all 57 CYP genes')
    genes=sorted(base.gene.unique(),key=lambda s:[int(t) if t.isdigit() else t for t in re.split(r'(\d+)',s)])
    lut={(r.gene,r.region,r.lineage):r.log_cpm for r in base.itertuples(index=False)}
    for i in range(3):
        xx=x+i*(w/3); cw=(w/3-38)/12
        for j,region in enumerate(regs):
            p.txt(xx+34+(j*3+1.5)*cw,y+7,region,6.1,True,anchor='middle')
            for l,lineage in enumerate((cm,fb,ec)):
                p.txt(xx+34+(j*3+l+.5)*cw,y+13,ls[lineage],5.35,anchor='middle',fill=lin[lineage])
        dy=(h-16)/19
        for a,g in enumerate(genes[i*19:(i+1)*19]):
            yy=y+16+a*dy;p.txt(xx,yy+dy*.83,g,5.5)
            for j,region in enumerate(regs):
                for l,lineage in enumerate((cm,fb,ec)):
                    val=lut.get((g,region,lineage),np.nan)
                    p.rect(xx+34+(j*3+l)*cw,yy,cw-.2,dy-.28,blend('#F8FBFD',blue,val/2.6) if pd.notna(val) else '#DDE3E8')
    prof=get('Fig2c_d_LV_donor_profiles')
    for k,g in zip('CDE',('CYP2J2','CYP4X1','CYP1B1'),strict=True):
        x,y,w,h=p.panel(k,g+' in non-failing LV')
        px=x+25; ww=w-33; top=y+4; bottom=y+h-19
        lim=(0,3.5)
        for i,lineage in enumerate((cm,fb,ec)):
            vv=prof.loc[prof.gene.eq(g)&prof.lineage.eq(lineage),'log_cpm'].dropna().tolist()
            yy=top+(bottom-top)*(i+.5)/3
            p.txt(x,yy+2,ls[lineage],6.5,fill=lin[lineage])
            p.line(px,yy,px+ww,yy,grid)
            for j,v in enumerate(vv):
                p.cloud(px+ww*min(v,lim[1])/lim[1],yy+((j*37%23)-11)*.34,.89,lin[lineage])
            if vv:
                mx=px+ww*np.median(vv)/lim[1]
                p.line(mx,yy-5,mx,yy+5,ink,1.1)
        axes_x(p,px,bottom,ww,lim,[0,1,2,3],'log10(1 + donor CPM)')
    lv=base.query("region=='LV'")
    chosen=[('CYP2J2',cm),('CYP4X1',fb),('CYP4X1',ec),('CYP1B1',fb),('CYP1B1',ec),('CYP19A1',fb)]
    x,y,w,h=p.panel('F','LV detection across study strata')
    for i,(g,l) in enumerate(chosen):
        r=keyrow(lv,{'gene':g,'lineage':l}); yy=y+6+i*11
        p.txt(x,yy,g+'  '+ls[l],6.45)
        n=int(r.study_strata); nn=int(r.study_strata_detected_in_majority)
        for j in range(n):p.dot(x+53+j*7.5,yy-1.6,2.25,lin[l],j<nn)
        p.txt(x+w-2,yy,f'{nn}/{n}',6.4,anchor='end')
    for k,g in [('G','CYP2J2'),('H','CYP4X1')]:
        vals=[(ls[l],float(keyrow(lv,{'gene':g,'lineage':l}).median_of_study_median_cpm),lin[l]) for l in (cm,fb,ec)]
        bars(p,k,g+' LV abundance',vals,100,[0,25,50,75,100],'Median study median CPM',.20)
    p.save()

def fig3():
    p=layout(3,'Matched disease effects across etiologies',[
        'Discovery effects are edgeR log2 fold changes. Color identifies the case diagnosis. Gray dots have tier-wide q >= 0.05.',
        'Dashed line: q = 0.05. Forest points and approximate pointwise 95% intervals use unshrunk log2 fold changes.'],
        [('A',[1]),('BCD',[1,1,1]),('EFG',[1,1,1.22]),('HI',[1,1])],[48,105,117,110])
    p.flow('Donor-level differential expression',[
        ('Original RNA counts','Sum by donor, region and lineage',blue),('Matched edgeR models','Fit within study and sampled region',lin[fb]),
        ('Tier-wide FDR','Primary and secondary tests corrected separately',lin[ec])])
    color_legend={'B':[('Higher',red),('Lower',blue),('q ≥ 0.05','#89959F')]}
    shapes={'AS_IVS':'diamond','NICM_LV':'square','DoxCM_vs_NICM_LV':'diamond'}
    shape_keys={'E':[('HCM','circle'),('AS','diamond')],'F':[('DoxCM','circle'),('NICM','square'),('Dox / NICM','diamond')]}
    for k,title,keys in [('B','Dilated cardiomyopathy',dc+['Reichart_DCM_RV']),('C','Arrhythmogenic cardiomyopathy',ar+['Reichart_ARVC_LV']),
        ('D','Ischemic cardiomyopathy',['Simonson_ICM_LV']),('E','HCM and aortic stenosis',['Chaffin_HCM_LV','AS_IVS']),
        ('F','DoxCM and NICM',['DoxCM_LV','NICM_LV','DoxCM_vs_NICM_LV'])]:
        d=allfx[allfx.contrast.isin(keys)]
        rs=[(r.logFC,-np.log10(max(r[q],1e-300)),direction(r.logFC,r[q]<.05),r[q]<.05,shapes.get(r.contrast,'circle')) for _,r in d.iterrows()]
        xstep=max(1,math.ceil(max(abs(r[0]) for r in rs)/2))
        ystep=max(1,math.ceil(max(r[1] for r in rs)/4))
        dots(p,k,title,rs,(-2.15*xstep,2.15*xstep),(0,4.25*ystep),[-2*xstep,-xstep,0,xstep,2*xstep],[i*ystep for i in range(5)],'log2 fold change','-log10(q)')
        if k in color_legend:
            x,y,w,h=p.box[k]
            for i,(label,d) in enumerate(color_legend[k]):
                p.dot(x+w-42,y+21+i*8,1.7,d);p.txt(x+w-37,y+23+i*8,label,6.0)
        if k in shape_keys:
            x,y,w,h=p.box[k]
            for i,(label,shape) in enumerate(shape_keys[k]):
                marker(p,x+w-48,y+21+i*8,1.6,muted,shape);p.txt(x+w-43,y+23+i*8,label,6.0)
    d=with_ci(get('Fig3b_fibroblast_CYP4B1')).sort_values('logFC')
    rec=[ci_row(short[r.contrast],r) for _,r in d.iterrows()]
    limits,ticks=ci_axis(rec)
    estimates(p,'G','Fibroblast CYP4B1',rec,limits,ticks,'log2 fold change',.55)
    d=with_ci(get('Fig3c_etiology_highlights'))
    rec=[ci_row(('AS ' if r.contrast=='AS_IVS' else 'Dox / NICM ')+r.cyp_gene+' '+ls[r.lineage],r) for _,r in d.iterrows()]
    limits,ticks=ci_axis(rec)
    estimates(p,'H','Selected septal and DoxCM effects',rec,limits,ticks,'log2 fold change',.48)
    x,y,w,h=p.panel('I','Eligible tests by contrast')
    left=x+76; ww=w-82;bottom=y+h-19;dy=(bottom-y)/13
    for i,(contrast,d) in enumerate(allfx.groupby('contrast',sort=False)):
        yy=y+(i+.5)*dy
        p.txt(x,yy+1.65,short[contrast],5.65)
        p.rect(left,yy-1.5,len(d)/140*ww,3,'#DFE8EE')
        p.rect(left,yy-1.5,(d[q]<.05).sum()/140*ww,3,disease[kind[contrast]])
    axes_x(p,left,bottom,ww,(0,140),[0,35,70,105,140],'Eligible (light) and q < 0.05 (color)')
    p.save()

def fig4():
    p=layout(4,'Recurrent CYP changes within diagnoses',[
        'Recurrence requires an eligible result, concordant direction and tier-wide q < 0.05 in each study.',
        'Heatmap cells show log2 fold change. * q < 0.05. Panel H shows unshrunk effects with approximate pointwise 95% intervals.'])
    p.flow('Comparing effects across cohorts',[
        ('Complete case set','Retain pairs tested in every cohort',blue),('Directional agreement','Compare study-specific effect signs',lin[fb]),
        ('FDR support','Require q < 0.05 in every cohort',lin[ec])])
    dcm=get('Fig4a_study_effect_concordance'); arvc=get('Fig4b_study_effect_concordance'); rr=get('Fig4c_12_recurrent_DCM_pairs')
    pairs_dc=set(zip(rr.lineage,rr.cyp_gene))
    av=allfx[allfx.contrast.isin(ar)]
    pairs_ar=set()
    for key,g in av.groupby(['lineage','cyp_gene']):
        if len(g)==2 and len(set(np.sign(g.logFC)))==1 and (g[q]<.05).all():pairs_ar.add(key)
    assert len(pairs_dc)==12 and len(pairs_ar)==4
    for k,title,d,s,t,hits in [('B','Chaffin vs Koenig',dcm,dc[0],dc[1],pairs_dc),('C','Chaffin vs Reichart',dcm,dc[0],dc[2],pairs_dc),
                                ('D','Reichart vs Simonson ARVC',arvc,ar[0],ar[1],pairs_ar)]:
        bound=6 if max(d[[s,t]].abs().max())<=6 else 8
        pts=[(r[s],r[t],direction(r[s],(r.lineage,r.cyp_gene) in hits),(r.lineage,r.cyp_gene) in hits) for _,r in d.iterrows()]
        dots(p,k,title,pts,(-bound,bound),(-bound,bound),[-bound,0,bound],[-bound,0,bound],s.split('_')[0]+' log2 FC',t.split('_')[0],True)
    tiles(p,'E','Twelve recurrent LV DCM pairs',[(a,ls.get(a[0],a[0])+' '+a[1]) for a in sorted(pairs_dc)],dc,lookup(rr,['lineage','cyp_gene'],'contrast','logFC',q),
          3.6,lwidth=.36,labels=dict(zip(dc,dn)))
    x,y,w,h=p.panel('F','Concordance and support')
    for i,(name,nums,color) in enumerate([('LV DCM',[70,36,12],disease['DCM']),('RV ARVC',[55,30,4],disease['ARVC'])]):
        yy=y+i*43
        p.txt(x,yy+10,name,7.3,True,fill=color)
        for j,(v,label) in enumerate(zip(nums,('Eligible','Same sign','q-supported'),strict=True)):
            xx=x+2+j*w/3;p.txt(xx,yy+24,v,10.7,True,fill=color);p.txt(xx,yy+34,label,6.4)
            if j<2:p.arrow(xx+24,yy+21,xx+w/3-5,color)
    tiles(p,'G','Four recurrent RV ARVC pairs',[(a,ls.get(a[0],a[0])+' '+a[1]) for a in sorted(pairs_ar)],ar,lookup(av,['lineage','cyp_gene'],'contrast','logFC',q),
          4.7,lwidth=.37,labels=dict(zip(ar,['Reichart','Simonson'])))
    chosen=[(cm,'CYP2J2'),(cm,'CYP11A1'),(fb,'CYP19A1')]
    rec=[]
    for l,gene in chosen:
        for s,ss in zip(dc,('Chaf.','Koen.','Reich.'),strict=True):
            r=keyrow(allfx,{'contrast':s,'lineage':l,'cyp_gene':gene})
            label=f'{gene}  {ls[l]}  {ss}'
            rec.append(ci_row(label,r) if r is not None else (label,np.nan,None,None,lin[l],False))
    limits,ticks=ci_axis(rec)
    estimates(p,'H','Focal cardiomyocyte and fibroblast effects',rec,limits,ticks,'log2 fold change',.40)
    p.save()

def fig5():
    p=layout(5,'Pathway enrichment and inferred signaling',[
        'GO tiles show normalized enrichment scores (NES). PROGENy effects are in within-model SD. Bars are 95% intervals.',
        '* or filled points indicate tier-wide q < 0.05. RNA footprints do not measure estrogen levels or pathway flux.'])
    p.flow('Two complementary RNA readouts',[
        ('Matched RNA statistics','Rank the signed edgeR statistics',blue),('GO enrichment','Test five prespecified pathways',lin[fb]),
        ('PROGENy footprints','Fit donor-level scores within contrast',lin[ec])])
    go=get('Fig5a_GO_enrichment'); pr=get('Fig5b_PROGENy_footprints'); focal=get('Fig5c_estrogen_focal')
    terms=[('icosanoid metabolic process','Eicosanoid'),('arachidonate metabolic process','Arachidonate'),
        ('inflammatory response','Inflammation'),('response to oxidative stress','Oxidative stress'),('apoptotic process','Apoptosis')]
    for k,l in zip('BCD',(cm,fb,ec),strict=True):
        d=go[go.lineage.eq(l)&go.contrast.isin(dc)]
        tiles(p,k,ls[l]+' GO enrichment in LV DCM',[((a,),b) for a,b in terms],dc,lookup(d,['pathway'],'contrast','NES','GO_BH_across_tier'),
              2.7,lwidth=.40,labels=dict(zip(dc,['Chaf.','Koen.','Reich.'])))
    d=go[go.pathway.eq('icosanoid metabolic process')&go.lineage.isin((cm,fb,ec))]
    keys=['Chaffin_DCM_LV','Koenig_DCM_LV','Reichart_DCM_LV','Simonson_ICM_LV','Simonson_ARVC_RV','AS_IVS','DoxCM_vs_NICM_LV']
    tiles(p,'E','Eicosanoid enrichment across etiologies',[((k,),short[k]) for k in keys],(cm,fb,ec),
          lookup(d,['contrast'],'lineage','NES','GO_BH_across_tier'),2.7,lwidth=.45,labels=ls)
    estimates(p,'F','Cardiomyocyte estrogen footprint',[forestrow(s,keyrow(focal,{'contrast':k}),lin[cm]) for k,s in zip(dc,dn,strict=True)],
              (-2.5,.5),[-2.5,-2,-1.5,-1,-.5,0,.5],'Score difference (SD)',.35)
    for k,l,path,title in [('G',fb,'WNT','Fibroblast WNT footprint'),('H',cm,'Estrogen','Estrogen footprint across etiologies')]:
        tab=pr[pr.pathway.eq(path)&pr.lineage.eq(l)].drop_duplicates('contrast')
        chosen=[keyrow(tab,{'contrast':c}) for c in ['Chaffin_DCM_LV','Chaffin_HCM_LV','Simonson_ICM_LV','AS_IVS','DoxCM_LV','Simonson_ARVC_RV']]
        rec=[forestrow(short[r.contrast],r,lin[l]) for r in chosen if r is not None]
        assert rec
        lo=min(-3,math.floor(min(r[2] for r in rec))); hi=max(3,math.ceil(max(r[3] for r in rec)))
        estimates(p,k,title,rec,(lo,hi),tickset(lo,hi),'Score difference (SD)',.39)
    p.save()

def fig6():
    p=layout(6,'Regulons and cell-state programs in LV DCM',[
        'Intervals are 95% donor-model CIs. Filled points and * indicate tier-wide q < 0.05 in the relevant test family.',
        'LIANA scores compare candidate sender-to-receiver RNA co-expression. Peri, pericyte. Myel, myeloid.'])
    p.flow('From donor RNA to inferred transcriptional scores',[
        ('Donor pseudobulk','TMM-normalized RNA profiles',blue),('Signed regulon targets','Infer THRB footprint from target RNA',lin[fb]),
        ('GO state programs','Score ATP, stress and vascular programs',lin[ec])])
    tab=get('Figure_6_Source_Data',True)
    for k,l in zip('BCD',(cm,fb,ec),strict=True):
        rec=[forestrow(s,keyrow(tab,{'lineage':l,'contrast':c,'program':'THRB'}),lin[l]) for c,s in zip(dc,dn,strict=True)]
        estimates(p,k,ls[l]+' THRB regulon',rec,(-.75,.25),[-.75,-.5,-.25,0,.25],'Score difference',.33)
    atp=tab[tab.program.eq('GO:0042775')&tab.contrast.isin(dc)]
    rec=[forestrow(ls[l]+' '+s,keyrow(atp,{'lineage':l,'contrast':c}),lin[l]) for l in (cm,fb,ec) for c,s in zip(dc,['Chaf.','Koen.','Reich.'],strict=True)]
    lo=min(-1.5,math.floor(min(r[2] for r in rec)*2)/2)
    estimates(p,'E','Mitochondrial ATP synthesis program',rec,(lo,.25),[x/2 for x in range(math.ceil(lo*2),1)],'Score difference',.35)
    other=get('Fig6c_other_cell_state_programs',True)
    terms=[('response to oxidative stress','Oxidative stress'),('response to hypoxia','Hypoxia'),('angiogenesis','Angiogenesis')]
    tiles(p,'F','Additional RNA state programs',[((l,t),ls[l]+' '+label) for l in (cm,fb,ec) for t,label in terms],dc,
          lookup(other,['lineage','term'],'contrast','score_delta','FDR_BH'),.55,lwidth=.44,labels=dict(zip(dc,dn)))
    pair=get('ED6_selected_LIANA_scores')
    lr={**ls,'pericyte':'Peri','myeloid':'Myel'}
    for k,c,s in [('G','Koenig_DCM_LV','Koenig'),('H','Reichart_DCM_LV','Reichart')]:
        table=pair[pair.contrast.eq(c)]
        rec=[forestrow(f'{r.ligand}/{r.receptor}  {lr.get(r.sender,r.sender)} → {lr.get(r.receiver,r.receiver)}',r,lin.get(r.sender,grey),lo='ci_low',hi='ci_high') for _,r in table.iterrows()]
        assert rec
        low=min(-4,math.floor(min(r[2] for r in rec)));high=max(4,math.ceil(max(r[3] for r in rec)))
        estimates(p,k,s+' ligand–receptor RNA scores',rec,(low,high),tickset(low,high),'Paired RNA score difference',.53)
    p.save()

def fig7():
    p=layout(7,'External RNA validation and donor associations',[
        'Focal correlations are within LV DCM. Leuven intervals use the donor models, and bulk intervals use Welch tests.',
        'Filled points indicate q < 0.05. The Leuven shared-batch subset reuses all four non-failing hearts.'],
        [('A',[1]),('BCD',[1.18,1,1]),('EF',[1,1]),('GH',[1,1])],[48,105,115,112])
    p.flow('External checks at two measurement scales',[
        ('Discovery in LV DCM','Cardiomyocyte CYP2J2 and fibroblast CYP19A1',blue),
        ('Leuven single-nucleus RNA','Lineage-specific donor pseudobulk',lin[fb]),
        ('GSE116250 bulk LV RNA','Whole-tissue expression comparison',lin[ec])])
    x,y,w,h=p.panel('B','Donor association screens')
    screen=get('Fig7a_association_or_validation')
    assert int(screen.below_005.sum())==0
    for i,(big,small) in enumerate([('109,637','Within-stratum all-donor tests'),
                (f'{int(screen.tested.sum()):,}','Cross-cohort entries across modes'),('0','Cross-cohort hits at q < 0.05')]):
        yy=y+8+i*22;p.txt(x,yy,big,9,True,fill=blue);p.txt(x,yy+9,small,6.3,fill=muted)
    assoc=get('Fig7b_association_or_validation')
    for k,l,g in [('C',cm,'CYP2J2'),('D',fb,'CYP19A1')]:
        a=assoc[assoc['mode'].eq('all_donors')&assoc.lineage.eq(l)]
        rows=[pointrow(r.study+' n='+str(r.n_donors),r,field='rho',qcol='Focal_FDR_BH') for _,r in a.iterrows()]
        assert len(rows)==3
        estimates(p,k,g+' and estrogen footprint',rows,(-1,1),[-1,-.5,0,.5,1],'Spearman rho',.43)
    val=get('Fig7c_association_or_validation');bulk=get('Fig7d_association_or_validation')
    for k,l,g in [('E',cm,'CYP2J2'),('F',fb,'CYP19A1')]:
        d=val[val.gene.eq(g)]
        rows=[forestrow(('Batch adjusted' if r.model=='batch_adjusted' else 'Shared batch')+f'  {r.case_n}/{r.control_n}',r,
             lin[l],field='log2FC',lo='CI_low',hi='CI_high') for _,r in d.iterrows()]
        estimates(p,k,'Leuven '+ls[l]+' '+g,rows,(-1,6),[-1,0,1,2,3,4,5,6],'Pseudobulk log2 fold change',.49)
    for k,g in [('G','CYP2J2'),('H','CYP19A1')]:
        r=keyrow(bulk,{'gene':g})
        rows=[forestrow('DCM / non-failing  37/14',r,'#465F74',field='log2_RPKM_plus_0_1_delta',lo='ci95_low',hi='ci95_high')]
        estimates(p,k,'GSE116250 bulk '+g,rows,(-1.5,1.5),[-1.5,-1,-.5,0,.5,1,1.5],'Mean difference in log2(RPKM + 0.1)',.43)
    p.save()

def fig8():
    p=layout(8,'Chamber context and evidence synthesis',[
        'Baseline n counts eligible study strata. Separate LV and RV fits do not estimate a disease-by-chamber interaction.',
        'Heatmap tiles show log2 fold change. * q < 0.05. Functional readouts in H are proposed measurements.'])
    p.flow('From chamber-resolved RNA to testable hypotheses',[
        ('Regional baseline','Compare eligible non-failing strata',blue),('Matched chamber effects','Separate LV and RV case–control models',lin[fb]),
        ('Functional hypotheses','Measure products in the implicated lineages',lin[ec])])
    for k,l,g in [('B',cm,'CYP2J2'),('C',fb,'CYP19A1'),('D',ec,'CYP1B1')]:
        r=[keyrow(base,{'region':c,'lineage':l,'gene':g}) for c in regs]
        v=[(c+'  n='+str(int(row.study_strata)),float(row.median_of_study_median_cpm),lin[l]) for c,row in zip(regs,r,strict=True)]
        maximum=max(a[1] for a in v)
        unit=10**math.floor(math.log10(maximum))
        top=math.ceil(maximum/unit)*unit
        bars(p,k,ls[l]+' '+g+' baseline',v,top,tickset(0,top),'Median study median CPM',.30)
    effects=allfx[allfx.contrast.isin(['Reichart_DCM_LV','Reichart_DCM_RV','Reichart_ARVC_LV','Reichart_ARVC_RV'])]
    rows=[(a,ls[a[0]]+' '+a[1]) for a in [(cm,'CYP2J2'),(fb,'CYP19A1'),(ec,'CYP1B1'),(fb,'CYP4X1'),(ec,'CYP4X1'),(fb,'CYP4B1')]]
    for k,d in [('E','DCM'),('F','ARVC')]:
        cs=[f'Reichart_{d}_LV',f'Reichart_{d}_RV']
        tiles(p,k,'Reichart '+d+' chamber estimates',rows,cs,lookup(effects,['lineage','cyp_gene'],'contrast','logFC',q),
              3.5,lwidth=.42,labels=dict(zip(cs,['LV','RV'])),mask_nonsignificant=True)
    x,y,w,h=p.panel('G','What the RNA evidence supports')
    for i,(big,small,color) in enumerate([
        ('LV DCM recurrence','12 CYP–lineage pairs across three cohorts',disease['DCM']),
        ('RV ARVC recurrence','4 CYP–lineage pairs across two cohorts',disease['ARVC']),
        ('External RNA','Leuven supports both focal lineage increases',lin[ec])]):
        yy=y+6+i*24;p.txt(x,yy,big,8,True,fill=color);p.txt(x,yy+10,small,6.8)
    x,y,w,h=p.panel('H','Direct functional tests')
    for i,(name,detail,result,color) in enumerate([
        ('Cardiomyocyte CYP2J2','Measure CYP2J2 protein and','arachidonate-derived epoxides and diols.',lin[cm]),
        ('Fibroblast CYP19A1','Measure aromatase protein and','labeled androgen-to-estrogen conversion.',lin[fb])]):
        yy=y+2+i*44;p.rect(x,yy,w,36,bg);p.rect(x,yy,2.5,36,color)
        p.txt(x+8,yy+10,name,7.7,True,fill=color);p.txt(x+8,yy+20,detail,6.6);p.txt(x+8,yy+29,result,6.6)
    p.save()

def extpage(n,rows,heights):
    p=Page(n,'',rows,heights,[])
    p.extended=True
    return p

def ed1():
    p=extpage(1,[('AB',[1,1]),('CD',[1,1])],[143,143])
    for k,n,title in [
       ('A','02_Chaffin','Chaffin'),('B','06_Kuppe','Kuppe'),
       ('C','01_AS','Lindman aortic stenosis'),('D','08_Reichart','Reichart')]:
        path='Previous_Display_and_Atlas/SupplementaryFig'+n+'_ab_coordinates.tsv.gz'
        ref=pd.read_csv(io.BytesIO(z.read(path)),sep='\t',compression='gzip')
        x,y,w,h=p.panel(k,title)
        minx,maxx=ref.umap_1.min(),ref.umap_1.max()
        miny,maxy=ref.umap_2.min(),ref.umap_2.max()
        scale=min((w*.66-8)/(maxx-minx),(h-14)/(maxy-miny))
        xx=x+3+(w*.64-(maxx-minx)*scale)/2
        yy=y+1+(h-12-(maxy-miny)*scale)/2
        for i,(l,d) in enumerate(ref.groupby('marker_candidate',sort=True)):
            color=lin.get(l,grey)
            p.clouds(zip(xx+(d.umap_1-minx)*scale,yy+(maxy-d.umap_2)*scale,strict=True),.28,color)
            px=x+w*.68;py=y+8+i*8.1
            p.dot(px,py-1.7,1.8,color)
            p.txt(px+4.8,py,ls.get(l,l.replace('_',' ')),6.3)
        p.txt(x+w*.34,y+h-2,'UMAP 1 / UMAP 2',6.5,anchor='middle',fill=muted)
    p.save()

def color_key(p,x,y,lo,hi,name,fn,width=72):
    for j in range(24):
        v=lo+(hi-lo)*j/23
        p.rect(x+j*width/24,y,width/24+.1,5,fn(v))
    p.txt(x,y+14,f'{lo:g}',6.2)
    p.txt(x+width,y+14,f'{hi:g}',6.2,anchor='end')
    p.txt(x+width+12,y+5,name,6.3)

def ed2():
    p=extpage(2,[('A',[1])],[278])
    x,y,w,h=p.panel('A','Detection recurrence')
    df=get('ED2_four_chamber_CYP_recurrence')
    genes=sorted(df.gene.unique(),key=lambda s:[int(t) if t.isdigit() else t for t in re.split(r'(\d+)',s)])
    lut={(r.gene,r.region,r.lineage):r.study_strata_fraction_detected_in_majority for r in df.itertuples(index=False)}
    heat=lambda v:blend('#F8FBFD',blue,v)
    assert len(genes)==57
    for side,gg in enumerate((genes[:29],genes[29:])):
        xx=x+side*(w/2)
        cw=(w/2-46)/12
        for j,region in enumerate(regs):
            p.txt(xx+37+(j*3+1.5)*cw,y+8,region,6.5,True,anchor='middle')
            for a,l in enumerate((cm,fb,ec)):
                p.txt(xx+37+(j*3+a+.5)*cw,y+16,ls[l],6.1,anchor='middle',fill=lin[l])
        for i,g in enumerate(gg):
            yy=y+20+i*7.6
            p.txt(xx,yy+6,g,6.15)
            for j,r in enumerate(regs):
                for a,l in enumerate((cm,fb,ec)):
                    v=lut.get((g,r,l),np.nan)
                    p.rect(xx+37+(j*3+a)*cw,yy,cw-.4,7.2,heat(v) if pd.notna(v) else '#DCE3E9')
    color_key(p,x,y+245,0,1,'Fraction of eligible study strata',heat)
    p.rect(x+268,y+245,9,5,'#DCE3E9')
    p.txt(x+282,y+250,'Unavailable',6.3)
    p.save()

def ed3():
    p=extpage(3,[('A',[1])],[156])
    d=get('ED3_covariate_sensitivity')
    rec=[]
    for _,r in d.iterrows():
        adjusted=r['mode']=='complete_case_adjusted'
        label=f"{r.study}  {r.cyp_gene}  {'Adjusted' if adjusted else 'Unadjusted'}  n={r.case_donors}/{r.control_donors}  q={noteq(r.FDR_BH)}"
        rec.append(pointrow(label,r,qcol='FDR_BH',col=direction(r.logFC,r.FDR_BH<.05)))
    estimates(p,'A','Complete-case models',rec,(-.5,2.5),[-.5,0,.5,1,1.5,2,2.5],'log2 fold change',.49)
    p.save()

def ed4():
    p=extpage(4,[('A',[1])],[170])
    x,y,w,h=p.panel('A','Flagged-donor exclusion')
    d=get('ED4_overlap_sensitivity')
    p.txt(x+130,y+8,'All donors',7,True)
    p.txt(x+305,y+8,'Overlap-exclusion model',7,True)
    for i,r in enumerate(d.itertuples(index=False)):
        yy=y+15+i*26
        p.line(x,yy,x+w-4,yy,grid)
        p.txt(x,yy+11,r.study+'  '+r.cyp_gene,7)
        p.txt(x+130,yy+11,f'{r.logFC:+.2f}  |  q={noteq(r.CYP_BH_across_tier)}',7,fill=direction(r.logFC,r.CYP_BH_across_tier<.05))
        p.txt(x+130,yy+20,f'n={r.case_donors}/{r.control_donors}',6.4)
        if pd.notna(r.logFC_overlap_screen):
            p.txt(x+305,yy+11,f'{r.logFC_overlap_screen:+.2f}  |  q={noteq(r.CYP_BH_across_tier_overlap_screen)}',
                  7,fill=direction(r.logFC_overlap_screen,r.CYP_BH_across_tier_overlap_screen<.05))
            p.txt(x+305,yy+20,f'n={r.case_donors_overlap_screen:g}/{r.control_donors_overlap_screen:g}',6.4)
        else:p.txt(x+305,yy+11,'Identical to all donors',6.8,fill=muted)
    p.save()

def ed5():
    p=extpage(5,[('A',[1])],[386])
    x,y,w,h=p.panel('A','PROGENy response footprints')
    d=get('ED5_all_PROGENy_selected_lineages').query("population=='all_donors'")
    cols=['AS_IVS','Chaffin_DCM_LV','Koenig_DCM_LV','Reichart_DCM_LV','Chaffin_HCM_LV',
          'Simonson_ICM_LV','DoxCM_LV','NICM_LV','DoxCM_vs_NICM_LV','Reichart_DCM_RV',
          'Reichart_ARVC_RV','Simonson_ARVC_RV','Reichart_ARVC_LV']
    names=[('AS','IVS','Lind.'),('DCM','LV','Chaf.'),('DCM','LV apex','Koen.'),('DCM','LV','Reich.'),
           ('HCM','LV','Chaf.'),('ICM','LV','Sim.'),('DoxCM','LV','Guo'),('NICM','LV','Guo'),
           ('Dox/N','LV','Guo'),('DCM','RV','Reich.'),('ARVC','RV','Reich.'),('ARVC','RV','Sim.'),('ARVC','LV','Reich.')]
    lt=lookup(d,['lineage','pathway'],'contrast','score_delta','FDR_BH')
    pw=sorted(d.pathway.unique())
    assert len(pw)==14
    cw=(w-97)/13
    for j,labels in enumerate(names):
        for i,label in enumerate(labels):p.txt(x+93+(j+.5)*cw,y+6+i*7,label,6.1,anchor='middle')
    for i,l in enumerate((cm,fb,ec)):
        yy=y+26+i*110
        for k,t in enumerate(pw):
            a=yy+k*7.55
            p.txt(x,a+6,ls[l]+'  '+t,6.2,fill=lin[l])
            for j,col in enumerate(cols):
                v,hit=lt.get(((l,t),col),(np.nan,False))
                xx=x+93+j*cw
                p.rect(xx,a,cw-.45,7.1,'#E7EBEE' if pd.notna(v) and not hit else signed(v,2.5))
                if hit:p.dot(xx+cw/2,a+3.5,.85,ink)
    color_key(p,x,y+358,-2.5,2.5,'Within-model SD',lambda v:signed(v,2.5))
    p.rect(x+218,y+358,8,5,'#E2E7EC');p.txt(x+230,y+363,'Unavailable',6.3)
    p.rect(x+285,y+358,8,5,'#E7EBEE');p.txt(x+297,y+363,'q ≥ 0.05',6.3)
    p.dot(x+352,y+361,1.3,ink);p.txt(x+358,y+363,'q < 0.05',6.3)
    p.save()

def ed6():
    p=extpage(6,[('A',[1])],[207])
    d=get('ED6_selected_LIANA_scores').sort_values('FDR_BH',kind='stable')
    lr={**ls,'pericyte':'Peri','myeloid':'Myel'}
    rec=[]
    for _,r in d.iterrows():
        label=f'{r.study}  {r.ligand}/{r.receptor}  {lr.get(r.sender,r.sender)} → {lr.get(r.receiver,r.receiver)}  q={noteq(r.FDR_BH)}'
        rec.append(forestrow(label,r,lin.get(r.sender,grey),lo='ci_low',hi='ci_high'))
    low=min(-4,math.floor(min(r[2] for r in rec)))
    high=max(3,math.ceil(max(r[3] for r in rec)))
    estimates(p,'A','Selected ligand–receptor RNA scores',rec,(low,high),list(range(low,high+1)),'Paired RNA score difference',.51)
    p.save()

def svg_layers(file,panel_letters):
    from xml.etree import ElementTree as ET
    url='http://www.w3.org/2000/svg'
    ET.register_namespace('',url)
    tree=ET.parse(file)
    root=tree.getroot()
    labels=set(panel_letters)
    for group in root.iter(f'{{{url}}}g'):
        for text in group.findall(f'{{{url}}}text'):
            if text.text in labels:
                group.set('id','panel-'+text.text)
                labels.remove(text.text)
                break
    assert not labels,labels
    tree.write(file,encoding='utf-8',xml_declaration=True)

def render(p):
    is_ed=getattr(p,'extended',False)
    title=('Extended_Data_Fig_' if is_ed else 'Figure_')+f'{p.n:02d}'
    parent=package/('Extended_Data/Figures' if is_ed else 'Main_Figures')
    parent.mkdir(parents=True,exist_ok=True)
    panel_root=package/'Individual_Panels'/('Extended_Data' if is_ed else 'Main_Figures')/title
    panel_root.mkdir(parents=True,exist_ok=True)
    order=p.boundaries+[(None,len(p.D.contents))]
    assert [k for k,_ in order[:-1]]==list(p.box)
    layers=Drawing(p.W,p.H)
    layers.add(p.D.contents[0])
    for (letter,i),(_,j) in zip(order,order[1:]):
        elements=p.D.contents[i:j]
        layers.add(Group(*elements))
        x,y,w,h=p.box[letter]
        solo=Drawing(w+8,h+8)
        solo.add(Rect(0,0,w+8,h+8,fillColor=HexColor('#FFFFFF'),strokeColor=None))
        group=Group(*elements)
        group.translate(4-x,4-(p.H-y-h))
        solo.add(group)
        save_pair(solo,panel_root/(title+'_'+letter),[letter])
    save_pair(layers,parent/title,list(p.box))

def save_pair(drawing,path,letters):
    pdf=path.with_suffix('.pdf')
    svg=path.with_suffix('.svg')
    renderPDF.drawToFile(drawing,str(pdf))
    renderSVG.drawToFile(drawing,str(svg))
    svg_layers(svg,letters)

Page.save=render
for name in (args.figures or names):
    print('Rendering',name,flush=True)
    globals()[name]()
