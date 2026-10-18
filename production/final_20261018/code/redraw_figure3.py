import os
"""Four-panel display; all saved comparisons and extrema are retained."""
from plot_main import *
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Rectangle

f=fignew(170,42,'a','Coverage across the planned ligand-family atlas')
ax=axis(f,[.17,.32,.57,.48]);counts=[]
for i,fa in enumerate(FAM):
 z=C[C.family==fa];n=int(z.eligible.sum());ne=len(z)-n
 ax.barh(i,n,height=.5,color=FC[fa]);ax.barh(i,ne,left=n,height=.5,color=P['ne'])
 ax.text(len(z)+2,i,f'{n} eligible  /  {ne} NE',va='center',fontsize=7.5);counts.append((fa,n,ne))
ax.set_yticks(range(4),[FN[x] for x in FAM]);ax.invert_yaxis();ax.set_xlim(0,106);ax.set_xlabel('Ligand-family combinations (n)');ax.tick_params(axis='y',length=0)
f.text(.1,.035,f'{len(C)} planned; {int(C.eligible.sum())} eligible; {int((~C.eligible).sum())} NE. No-cytomix is separate.',fontsize=7)
save(f,'F3a',C,'Complete planned atlas counts; family colours and NE retained.')

cc=read('classification_counts_by_method.tsv');z=cc[cc.family=='ALL']
f=fignew(170,51,'b','Joint direction requires agreement in all three donors')
ax=axis(f,[.15,.32,.37,.45]);ax2=axis(f,[.57,.23,.40,.59]);ax2.axis('off')
for i,me in enumerate(METHODS):
 left=0
 for cl in CLASSES:
  n=int(z[(z.method==me)&(z.classification==cl)].n_combinations.iloc[0]);ax.barh(i,n,left=left,color=STYLE['class_colors'][cl],height=.55);left+=n
ax.set_yticks(range(3),['Raw rank','Relative rank','Non-rank']);ax.invert_yaxis();ax.set_xlim(0,len(C));ax.set_xlabel('Combinations (n)');ax.tick_params(axis='y',length=0)
for j,cl in enumerate(CLASSES):
 label={'BOTH':'Both','UP_ONLY':'UP\nonly','DEPLETED_ONLY':'Depleted\nonly','OTHER_OBSERVED':'Other','NE':'NE'}[cl]
 ax2.text((j+.5)/5,.97,label,ha='center',va='top',fontsize=7,color=P['ink'])
 ax2.plot([j/5+.03,(j+1)/5-.03],[.68,.68],lw=3,color=STYLE['class_colors'][cl])
 for i,me in enumerate(METHODS):
  n=int(z[(z.method==me)&(z.classification==cl)].n_combinations.iloc[0]);ax2.text((j+.5)/5,.53-i*.22,str(n),ha='center',va='center',fontsize=7.5)
ax2.set_xlim(0,1);ax2.set_ylim(0,1)
f.text(.1,.06,'Both: disease-up decreases AND depleted increases in each donor. Counts are descriptive.',fontsize=7)
save(f,'F3b',z,'All saved classifications, with small categories in a separate aligned count table.')

for me,pid,letter,oldpid in [('raw','F3c','c','F3b'),('nonrank','F3d','d','F3c')]:
 inp=Path(os.environ['ICI_DISPLAY_INPUT_ROOT'])/'archive/v1_panel_data'/f'{oldpid}.tsv'
 d=pd.read_csv(inp,sep='\t',float_precision='round_trip');q=d[['ligand','celltype','display_x','display_y']].drop_duplicates();assert len(q)==C.eligible.sum()
 f=fignew(85,103,letter,('Raw-rank' if me=='raw' else 'Non-rank')+' mean responses')
 ax=axis(f,[.20,.31,.73,.40]);full=axis(f,[.61,.77,.32,.14])
 # Display window only; full unmodified range remains in the separate overview.
 limits=(-.012,.018,-.022,.016) if me=='raw' else (-.20,.30,-.30,.25)
 for fa in FAM:
  s=q[q.celltype==fa]
  for a,size in [(ax,11),(full,3)]:a.scatter(s.display_x,s.display_y,s=size,color=FC[fa],alpha=.65,edgecolors='none')
 for a in [ax,full]:
  a.axhline(0,c=P['muted'],lw=.5);zero(a);a.xaxis.set_major_locator(MaxNLocator(3));a.yaxis.set_major_locator(MaxNLocator(3));a.tick_params(length=2,labelsize=7)
 ax.set_xlim(limits[0:2]);ax.set_ylim(limits[2:4]);ax.set_xlabel('- Disease-up change\n'+('Rank units' if me=='raw' else 'Mean member log2 units'));ax.set_ylabel('Depleted change')
 full.add_patch(Rectangle((limits[0],limits[2]),limits[1]-limits[0],limits[3]-limits[2],fill=False,lw=.7,edgecolor=P['ink']))
 full.set_title('Full range',fontsize=7,pad=3)
 bmp=q[(q.ligand=='BMP4')&(q.celltype=='Colonocyte')].iloc[0]
 full.annotate('BMP4',(bmp.display_x,bmp.display_y),xytext=(.02,.75),textcoords='axes fraction',fontsize=7,ha='left',arrowprops={'arrowstyle':'-','lw':.5,'color':P['ink']})
 for lig,fa,ofs in [('TNFSF12','Transit_Amplifying',(10,13))]:
  r=q[(q.ligand==lig)&(q.celltype==fa)].iloc[0];ax.scatter(r.display_x,r.display_y,s=29,facecolors='none',edgecolors=FC[fa],lw=.8);ax.annotate('TNFSF12 / TA',(r.display_x,r.display_y),xytext=ofs,textcoords='offset points',fontsize=7,arrowprops={'arrowstyle':'-','lw':.4})
 n=int(((q.display_x>0)&(q.display_y>0)).sum());total=len(q)
 f.text(.07,.84,f'All {total} means\nin overview.\nBox = enlarged range.',fontsize=7,linespacing=1.3)
 f.text(.08,.12,f'Upper-right: {n}/{total} means.\nMean direction differs from three-donor concordance.',fontsize=7,linespacing=1.3)
 hs=[Line2D([],[],ls='',marker='o',markersize=3,color=FC[fa],label=FN[fa]) for fa in FAM]
 f.legend(handles=hs,loc='lower center',bbox_to_anchor=(.52,.005),ncol=2,frameon=False,fontsize=7,handletextpad=.3,columnspacing=.6)
 save(f,pid,d,'Unmodified saved means; enlarged window and full-range overview preserve all extrema, including BMP4.')
layout=json.loads((R/'config/layout.json').read_text());layout['3']={'size':[170,196],'panels':[['F3a',0,0,170,42],['F3b',0,42,170,51],['F3c',0,93,85,103],['F3d',85,93,85,103]]};(R/'config/layout.json').write_text(json.dumps(layout,indent=2))
pd.DataFrame(PANELS).to_csv(R/'qa/figure3_source_map.tsv',sep='\t',index=False)
print('Figure 3: four panels; full saved range retained.')
