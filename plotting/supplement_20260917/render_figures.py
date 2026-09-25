import os
"""Draw S5/S7 as native PDF vectors from current CSV-derived records only.
No prior PDF, crop coordinates, raster panels, or historical panel mapping is read.
"""
from common import *
import math
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4,landscape
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import Paragraph
from reportlab.lib.styles import ParagraphStyle
from xml.sax.saxutils import escape

pdfmetrics.registerFont(TTFont('Arial',os.environ['ICI_FONT_REGULAR']))
pdfmetrics.registerFont(TTFont('Arial-Bold',os.environ['ICI_FONT_BOLD']))
W,H=landscape(A4)
S=jr(I/'plot_data/plot_spec.json')
T=dict(font='Arial',heading_pt=12,panel_letter_pt=14,panel_title_pt=9,axis_pt=8,tick_pt=8,legend_pt=8,caption_pt=9.5,caption_leading_pt=11.8,axis_width_pt=0.65,median_width_pt=1.6,zero_width_pt=0.5,point_radius_pt=2.4,page_margin_pt=32,
 condition_colors={'Healthy':'#9AA4AD','ICI colitis':'#B84B46','UC':'#2E8079'},family_colors={'Stem':'#2C79A6','Transit_Amplifying':'#C14A43','Colonocyte':'#795AA6','Goblet':'#C59633'},donor_shapes={'DNW14':'circle','DNW15':'triangle','DNW9':'square'},NE_color='#D4D9DE',heatmap_low='#356D95',heatmap_mid='#FFFFFF',heatmap_high='#C14A43',NES_limit=S['NES_limit'])
jw(I/'visual_specification.json',T)
ledger=[];geometries=[];current_figure='';current_page=0
def rgb(s):return tuple(int(s[i:i+2],16)/255 for i in (1,3,5))
def text(c,x,y,s,size=None,bold=False,align='left',color='#25333C'):
 c.setFillColorRGB(*rgb(color));c.setFont('Arial-Bold' if bold else 'Arial',size or T['axis_pt'])
 getattr(c,{'left':'drawString','center':'drawCentredString','right':'drawRightString'}[align])(x,y,str(s))
def line(c,x1,y1,x2,y2,width=0.65,color='#252A2E',dash=None):
 c.setStrokeColorRGB(*rgb(color));c.setLineWidth(width);c.setDash(dash or []);c.line(x1,y1,x2,y2);c.setDash([])
def page(c,fig,pagenum,first=False):
 global current_figure,current_page
 current_figure=fig;current_page=pagenum
 title={'S5':'Targeted donor-adjusted organoid responses','S7':'Oxford and Xenium spatial measurements'}[fig]
 text(c,32,H-30,f'Supplementary Figure {fig}. {title}' if first else f'Supplementary Figure {fig} (continued)',12,True,color='#000000')
 text(c,W/2,12,str(pagenum),8,align='center',color='#000000')
def panel_title(c,panel,title,x,y):
 text(c,x-45,y,panel,14,True);text(c,x,y,title,9,True)
def record(r,kind,**more):ledger.append(dict(figure=current_figure,page=current_page,kind=kind,record=r,**more))
def point(c,x,y,color,shape='circle',r=None):
 r=r or T['point_radius_pt'];c.setFillColorRGB(*rgb(color));c.setStrokeColorRGB(*rgb(color));c.setLineWidth(0)
 if shape=='circle':c.circle(x,y,r,stroke=0,fill=1)
 elif shape=='square':c.rect(x-r,y-r,2*r,2*r,stroke=0,fill=1)
 elif shape=='triangle':
  p=c.beginPath();p.moveTo(x,y+r);p.lineTo(x-r,y-r);p.lineTo(x+r,y-r);p.close();c.drawPath(p,stroke=0,fill=1)
 else:raise ValueError(shape)
def heatcolor(v):
 limit=float(S['NES_limit']);f=abs(v)/limit;assert f<=1
 end=rgb(T['heatmap_high'] if v>=0 else T['heatmap_low']);mid=rgb(T['heatmap_mid'])
 return tuple(a+(b-a)*f for a,b in zip(mid,end))
def heatmap(c,rows,x,top,cw_=62,ch=20,reference=False):
 for r in rows:
  left=x+r['display_col']*cw_;bottom=top-(r['display_row']+1)*ch
  color=rgb(T['NE_color']) if r['status']!='ESTIMATED' else heatcolor(float(r['NES']))
  c.setFillColorRGB(*color);c.setStrokeColorRGB(1,1,1);c.setLineWidth(.5);c.rect(left,bottom,cw_,ch,stroke=1,fill=1)
  txt=r['display_text'];tc='#FFFFFF' if min(color)<.3 and sum(color)/3<.48 else '#17242D'
  text(c,left+cw_/2,bottom+ch/2-2.7,txt,8,align='center',color=tc)
  record(r,'heatmap_cell',box=[left,bottom,cw_,ch],display=txt)
 nrows=max(r['display_row'] for r in rows)+1
 for i in range(nrows):
  label=S['modules'][i][1] if reference else S['ligands'][i].replace('_',' ').replace('No Cytomix','No cytomix')
  yy=top-(i+.5)*ch-2.7
  if reference and label=='IFN-gene-excluded residual':
   text(c,x-8,yy+5,'IFN-gene-excluded',8,align='right');text(c,x-8,yy-5,'residual',8,align='right')
  else:text(c,x-8,yy,label,8,align='right')
 for i,(_,label) in enumerate(S['families']):text(c,x+(i+.5)*cw_,top-nrows*ch-14,label,8,align='center')
 if not reference:
  for boundary in [1,nrows-1]:line(c,x,top-boundary*ch,x+4*cw_,top-boundary*ch,.65,'#66737E')
def colorbar(c,x,y,w=230):
 limit=float(S['NES_limit'])
 for i in range(180):
  c.setFillColorRGB(*heatcolor(-limit+2*limit*i/179));c.rect(x+w*i/180,y,w/180+.05,8,stroke=0,fill=1)
 for v in [-limit,-limit/2,0,limit/2,limit]:
  xx=x+(v+limit)/(2*limit)*w;text(c,xx,y-12,f'{v:g}',8,align='center')
 text(c,x-10,y,'NES',8,align='right')
def nice_ticks(vals,include_zero=False):
 lo=min(vals);hi=max(vals)
 if include_zero:lo=min(0,lo);hi=max(0,hi)
 span=hi-lo
 if span==0:span=max(abs(hi)*.1,.01)
 lo-=span*.09;hi+=span*.09
 rough=(hi-lo)/4;power=10**math.floor(math.log10(rough));f=rough/power
 step=min([1,2,2.5,5,10],key=lambda a:abs(a-f))*power
 ticks=[];v=math.ceil((lo-1e-12)/step)*step
 while v<=hi+step*1e-8:ticks.append(0.0 if abs(v)<step*1e-8 else v);v+=step
 decimals=max(0,-math.floor(math.log10(step)))+(1 if abs(step/power-2.5)<1e-8 else 0)
 return lo,hi,ticks,decimals
def axis(c,box,limits,direction,category_labels=None,xlabel=None,ylabel=None):
 x,y,w,h=box;lo,hi,ticks,decimals=limits
 line(c,x,y,x+w,y);line(c,x,y,x,y+h)
 for val in ticks:
  if direction=='x':
   xx=x+(val-lo)/(hi-lo)*w;line(c,xx,y,xx,y-3,.5);text(c,xx,y-14,f'{val:.{decimals}f}',8,align='center')
  else:
   yy=y+(val-lo)/(hi-lo)*h;line(c,x-3,yy,x,yy,.5);text(c,x-7,yy-2.7,f'{val:.{decimals}f}',8,align='right')
 if xlabel:text(c,x+w/2,y-29,xlabel,8,align='center')
 if ylabel:
  c.saveState();c.translate(x-47,y+h/2);c.rotate(90);text(c,0,0,ylabel,8,align='center');c.restoreState()
def paired(c,panel,title,x,y,w=280,h=172):
 rows=[r for r in S['S5_donors'] if r['panel']==panel];lim=nice_ticks([float(r['change']) for r in rows],True)
 panel_title(c,panel,title,x,y+h+15);axis(c,[x,y,w,h],lim,'x',xlabel='Treatment minus control raw rank')
 ymin=.5;ymax=len(S['ligands'])+.5
 fx=lambda v:x+(float(v)-lim[0])/(lim[1]-lim[0])*w
 fy=lambda v:y+(float(v)-ymin)/(ymax-ymin)*h
 line(c,fx(0),y,fx(0),y+h,.5,'#A3ADB7',[3,3])
 for i,ligand in enumerate(S['ligands']):
  yy=len(S['ligands'])-i;text(c,x-6,fy(yy)-2.7,ligand.replace('No_Cytomix','No cytomix'),8,align='right')
 for v in [len(S['ligands'])-.5,1.5]:line(c,x,fy(v),x+w,fy(v),.5,'#CAD1D8')
 for r in rows:
  px,py=fx(r['change']),fy(r['display_y']);co=T['family_colors'][r['celltype']];sh=T['donor_shapes'][r['donor']]
  point(c,px,py,co,sh);record(r,'donor_point',x=px,y=py,radius=T['point_radius_pt'],color=co,shape=sh,data_x=r['change'],data_y=r['display_y'])
 geometries.append(dict(figure='S5',panel=panel,page=current_page,box=[x,y,w,h],x_range=lim[:2],y_range=[ymin,ymax]))
def caption(c,fig,top):
 txt=(I/'plot_data'/f'{fig}_caption.txt').read_text(encoding='utf8')
 style=ParagraphStyle('caption',fontName='Arial',fontSize=9.5,leading=11.8,textColor='#000000')
 p=Paragraph(escape(txt),style);ww,hh=p.wrap(W-64,top-30)
 assert top-hh>=30,(fig,top,hh)
 p.drawOn(c,32,top-hh);record({'figure':fig},'caption',box=[32,top-hh,ww,hh],text=txt)
def condition_panel(c,panel,title,x,y,w=300,h=120,oxford=True):
 rows=[r for r in S['S7_Oxford' if oxford else 'S7_epithelium'] if r['panel']==panel]
 lim=nice_ticks([float(r['value']) for r in rows]);panel_title(c,panel,title,x,y+h+29)
 ann=next(r for r in S['S7_annotations'] if r['panel']==panel);text(c,x,y+h+14,ann['annotation'],8)
 record(ann,'annotation',text=ann['annotation'])
 axis(c,[x,y,w,h],lim,'y',ylabel='Disease-oriented mean gene rank' if oxford else 'Disease-oriented panel rank')
 fy=lambda v:y+(float(v)-lim[0])/(lim[1]-lim[0])*h
 for k,cond in enumerate(S['conditions']):
  group=sorted([r for r in rows if r['display_condition']==cond],key=lambda r:r['patient'])
  center=x+(k+.5)*w/3
  for j,r in enumerate(group):
   dx=(j-(len(group)-1)/2)*3.2;px=center+dx;py=fy(r['value']);co=T['condition_colors'][cond]
   point(c,px,py,co);record(r,'donor_point',x=px,y=py,radius=T['point_radius_pt'],color=co,shape='circle',data_y=r['value'],categorical_offset_pt=dx)
  med=next(r for r in S['S7_medians'] if r['panel']==panel and r['condition']==cond)
  yy=fy(med['value']);line(c,center-15,yy,center+15,yy,1.6,T['condition_colors'][cond])
  record(med,'median',x0=center-15,x1=center+15,y=yy,data_y=med['value'],color=T['condition_colors'][cond])
  text(c,center,y-14,cond,8,align='center');text(c,center,y-25,f"n={med['n']}",8,align='center')
 geometries.append(dict(figure='S7',panel=panel,page=current_page,box=[x,y,w,h],y_range=lim[:2]))
def gene_panel(c,panel,gene,x,y,w=660,h=102):
 rows=[r for r in S['S7_genes'] if r['panel']==panel];lim=nice_ticks([float(r['value']) for r in rows],True)
 panel_title(c,panel,gene,x,y+h+13);axis(c,[x,y,w,h],lim,'y',ylabel='Mean molecules / segmented cell')
 fy=lambda v:y+(float(v)-lim[0])/(lim[1]-lim[0])*h
 for k,(comp,label) in enumerate(S['compartments']):
  cx=x+(k+.5)*w/len(S['compartments']);text(c,cx,y-14,label,8,align='center')
  for l,cond in enumerate(S['conditions']):
   group=sorted([r for r in rows if r['compartment']==comp and r['display_condition']==cond],key=lambda r:r['patient'])
   for j,r in enumerate(group):
    dx=(l-1)*12+(j-(len(group)-1)/2)*3.2;px=cx+dx;py=fy(r['value']);co=T['condition_colors'][cond]
    point(c,px,py,co);record(r,'donor_point',x=px,y=py,radius=T['point_radius_pt'],shape='circle',color=co,data_y=r['value'],categorical_offset_pt=dx)
 geometries.append(dict(figure='S7',panel=panel,page=current_page,box=[x,y,w,h],y_range=lim[:2]))

out=I/'new_figures';out.mkdir(parents=True,exist_ok=True)
c=canvas.Canvas(str(out/'S5_source_redraw.pdf'),pagesize=(W,H),pageCompression=1,invariant=1)
page(c,'S5',9,True)
for i,(mod,title) in enumerate(S['modules']):
 panel=chr(ord('a')+i);x=[132,526][i%2];top=[514,278][i//2]
 panel_title(c,panel,title,x,top+16);heatmap(c,[r for r in S['S5_heatmaps'] if r['panel']==panel],x,top)
colorbar(c,307,61);c.showPage()
page(c,'S5',10)
for i,(mod,title) in enumerate(S['modules']):paired(c,chr(ord('e')+i),title,[116,520][i%2],[348,112][i//2])
x=151;y=39;text(c,x-13,y-2.7,'Family',8,align='right')
for cell,label in S['families']:
 point(c,x,y,T['family_colors'][cell],r=2.1);text(c,x+7,y-2.7,label,8);x+={'Stem':61,'TA':44,'Colonocyte':85,'Goblet':69}[label]
text(c,x+2,y-2.7,'Donor',8);x+=38
for donor in S['donors']:
 point(c,x,y,'#000000',T['donor_shapes'][donor],r=2.1);text(c,x+7,y-2.7,donor,8);x+=63
c.showPage();page(c,'S5',11);panel_title(c,'i','No-cytomix assay reference',278,513)
heatmap(c,S['S5_reference'],278,482,cw_=90,ch=32,reference=True);colorbar(c,320,305);caption(c,'S5',247);c.showPage();c.save()

c=canvas.Canvas(str(out/'S7_source_redraw.pdf'),pagesize=(W,H),pageCompression=1,invariant=1)
page(c,'S7',14,True)
for i,(mod,title) in enumerate(S['modules']):condition_panel(c,chr(ord('a')+i),title,[110,506][i%2],[390,195][i//2])
c.showPage();page(c,'S7',15)
for i,gene in enumerate(S['genes']):gene_panel(c,chr(ord('e')+i),gene,120,[406,249,92][i])
x=300
for cond in S['conditions']:
 point(c,x,41,T['condition_colors'][cond]);text(c,x+8,38.3,cond,8);x+=100
c.showPage();page(c,'S7',16)
for i,(mod,title) in enumerate(S['modules']):condition_panel(c,chr(ord('h')+i),title,[110,506][i%2],[390,195][i//2],oxford=False)
caption(c,'S7',157);c.showPage();c.save()
jw(I/'plot_data/render_ledger.json',ledger);jw(I/'plot_data/panel_geometry.json',geometries)
print('Native-vector redraw complete:',len(ledger),'logged drawing records')
