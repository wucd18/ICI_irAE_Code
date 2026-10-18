"""Bounded, source-supported narrative edits to the exact reviewed ZIP baseline."""
from pathlib import Path
from lxml import etree as E
import zipfile, json, copy, re, hashlib, os
F=Path(os.environ.get('ICI_PRODUCTION_ROOT',Path(__file__).resolve().parents[1]))
O=Path(os.environ.get('ICI_SUBMISSION_ROOT',F/'submission_final'))
I=Path(os.environ.get('ICI_SUBMISSION_BASELINE',F/'inputs/submission_final'))
NS={'w':'http://schemas.openxmlformats.org/wordprocessingml/2006/main','r':'http://schemas.openxmlformats.org/officeDocument/2006/relationships','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
W='{'+NS['w']+'}'
def load(p):
 with zipfile.ZipFile(p) as z:d={n:z.read(n) for n in z.namelist()}
 r=E.fromstring(d['word/document.xml']);return d,r,r.findall('.//w:p',NS)
def text(p):return ''.join(p.xpath('.//w:t/text()',namespaces=NS))
def settext(p,t):
 prop=p.find('w:pPr',NS);rr=p.find('.//w:rPr',NS);rr=copy.deepcopy(rr) if rr is not None else None
 for c in list(p):
  if c is not prop:p.remove(c)
 r=E.SubElement(p,W+'r')
 if rr is not None:r.append(rr)
 q=E.SubElement(r,W+'t');q.text=t;q.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
def save(d,r,p):
 d['word/document.xml']=E.tostring(r,xml_declaration=True,encoding='UTF-8',standalone=True)
 with zipfile.ZipFile(p,'w',zipfile.ZIP_DEFLATED) as z:
  for n,b in d.items():z.writestr(n,b)
d,r,ps=load(I/'Manuscript.docx');original=[text(p) for p in ps]
edits={
13:"Recurrent epithelial transcriptional programmes provide a common readout for comparing patient tissue with experimental perturbations in immune checkpoint inhibitor (ICI) colitis. In this secondary analysis of published data, we defined disease-up and disease-depleted programmes by patient-level differential expression in GSE253720 and assessed the same memberships in external cohorts. We reanalysed organoid counts using donor-adjusted models and evaluated fixed-programme responses across 79 ligands with matched cytomix controls. The 84-gene disease-up and 275-gene disease-depleted programmes retained their expected enrichment directions in four GSE206300 epithelial families and GSE210037 whole-tissue sample pseudobulks. Oxford donor-level comparisons did not pass the four-programme false discovery rate threshold. Of 316 planned ligand-family combinations, 278 were estimable and 38 were non-estimable. Concurrent disease-up decreases and depleted increases in all three donors occurred in 8, 12 and 7 combinations under raw-rank, matched-background-relative and non-rank representations, respectively. IFN-γ in stem cells showed a disease-direction response, whereas TNFSF12 in transit-amplifying cells retained depleted increases but had representation-dependent disease-up changes. BMP4-associated colonocyte responses accompanied increased recovered mature-state fractions; fixed state weights attenuated the depleted change without reversing its direction. Depleted responses correlated with IBD2 membership readouts after shared genes were removed. Joint evaluation of the two disease readouts distinguishes depleted-programme increases from concurrent reversal of both programme scores and relates these responses to inflammatory and maturation context.",
18:"A published human intestinal-organoid atlas linked inflammatory bowel disease (IBD) programmes to responses to 79 secreted ligands [14]. Fixed ICI-colitis programmes offer a way to compare that perturbation resource with patient disease associations. Joint evaluation retains both disease-up and disease-depleted readouts and distinguishes an isolated change in one programme from opposing changes in both.",
19:"We followed fixed epithelial memberships from a discovery cohort through external patient assessment and all eligible ligand-family combinations in the organoid atlas. We then examined inflammatory references, recovered maturation states and overlap with published IBD memberships to interpret the responses.",
22:"This secondary analysis of published human data links patient-level discovery, external assessment and paired-donor organoid perturbations (Fig. 1; Table 1). Patients, released tissue samples and organoid donors are the biological replicates; cells, spots and repeated sections are nested measurements. Available public cohorts determined sample size. No new participants were recruited. Programmes defined in the discovery cohort were held fixed for external assessment and organoid analysis.",
23:"GSE253720 [11] supplies discovery epithelium; GSE206300 [12] supplies external cell-family assessment; GSE210037 [10] supplies whole-tissue sample-pseudobulk assessment; and GSE313368 [14] supplies matched organoid perturbations. Oxford and Xenium data provide additional tissue assessments described below.",
44:"All-ligand perturbations; 3 paired donors; 278/316 eligible combinations; no-cytomix separate",
53:"Cohort-wide coverage and analysis-specific coverage are distinct. External epithelial families can share patients and must not be added as independent participants. The organoid no-cytomix reference is separate from the ligand-family denominator.",
65:"2.4. All-ligand organoid models and enrichment",
66:"The published screen used differentiated colon organoids from three donors [14]. On day 5 after seeding, ligands were added together with low-dose cytomix (TNF 10 ng/mL, IFN-γ 10 ng/mL and IL-1β 1 ng/mL) in medium containing 0.5% Matrigel. After 8 h, cytomix was removed by two Advanced DMEM/F12 washes; ligand-containing medium was restored and refreshed every 1–2 days until scRNA-seq collection on day 8. Doses were selected relative to vendor EC50 information and preliminary responses: most ligands used 100 ng/mL, IFN-γ and TNF used 50 ng/mL, and selected factors, including TNFSF12, used 250 ng/mL. Each sequencing channel included no-cytomix and cytomix-only controls. Table S4 records mapped doses; unresolved mappings remain NA. IFNG and TNFA in source-labelled figures and tables denote IFN-γ and TNF, respectively. We aggregated released RNA counts by donor, family and condition and matched treatments to cytomix-only Sato controls within experimental pools. Pools required ≥20 cells per arm and were combined within donor. Released family labels and their documented consolidation were retained. Eligible contrasts required three paired donors. DESeq2 used the design ~ donor + treatment [32], retaining genes with total count ≥10 across six pseudobulks. Independent filtering and automatic outlier replacement were disabled. Finite, unshrunk treatment Wald statistics supplied the ranking. No-cytomix cultures were compared with matched cytomix controls using the same rules.",
67:"Preranked enrichment used fgseaSimple with 10,000 gene-set permutations, seed 20260826 and set sizes of 10–500. Atlas-wide FDR adjustment used separate planned test families: 948 ligand-by-family-by-programme tests for disease-up, depleted and IFN; 316 residual-programme tests; and 16 no-cytomix reference tests. Unavailable estimates remained in the planned denominators. GSEA assesses enrichment of gene-level treatment statistics; its q values apply to enrichment tests, separately from donor-score changes, composition changes and directional classifications.",
68:"2.5. Programme-score representations and directional classes",
69:"Raw-rank changes are treatment-minus-matched-control differences in the mean within-sample expression rank of available fixed members. Depleted scores retain their raw direction in perturbation comparisons. Background matching used unique cytomix controls from eligible comparisons, deduplicated by donor, family and pool. Within each family, control counts and detection counts were pooled across these pools for each donor; log2(CPM + 1) expression and cell-detection fractions were then averaged equally across donors. Genes were assigned jointly to expression and detection deciles, with tied quantile boundaries collapsed. A stratum with fewer nonmember candidates than required target members was merged with the nearest stratum in two-dimensional bin coordinates, using a deterministic tie-break. Each of 1,000 background sets matched target size and stratum counts, sampled without replacement within a set and excluded target members. Family-and-membership-specific seeds derived from SHA-256 and seed 20260826 made the libraries reusable across conditions. Background-relative changes subtract the median background-set change from the raw programme change. Continuous matching imbalance is reported in Table S6.",
70:"Non-rank changes are the mean, across available fixed members, of treatment-minus-control log2(size-factor-normalized count + 1). Size factors were calculated separately for each six-pseudobulk comparison using the DESeq2 ratio formula: genes with total count ≥10 and positive counts in all six samples supplied geometric-mean references, and each sample factor was the median count-to-reference ratio. The factors were applied to the complete measured gene universe for programme scoring. This calculation matched the corresponding model size-factor extrema. Raw ranks, background-relative ranks and non-rank expression changes retain their distinct units; the non-rank main-figure display is a presentation choice within this retrospective analysis.",
71:"A descriptive Both class requires a negative disease-up change and a positive depleted change in each of the three matched donors under the specified representation. UP only requires the disease-up decrease alone, depleted only requires the depleted increase alone, and other observed denotes neither all-donor condition. Insufficient three-donor coverage is NE. The three representations use the same experiments and donors. TNFSF12 was an earlier candidate of interest; its assessment was extended to the full ligand range. IFN-γ, TNF, TNFSF12 and BMP4 were selected after inspection as explanatory examples. Original and additional Both combinations are reported in Table S5 and Supplementary Fig. S6.",
73:"Pool-specific analyses identified every donor-family-pool control and its reuse. Each ligand was compared with its matched cytomix-only control. Independent-pool and leave-one-pool-out analyses applied the same minimum-donor requirements and recorded NE reasons. Pool-wide ligand medians provided an additional reference alongside the matched control. Structural TGF-beta-superfamily summaries were descriptive because donors, controls and pool assignments were shared and the experiment was not fully crossed.",
75:"Composition analyses used released within-family substates and mean-cell ln(1 + 10,000 counts/library) expression, a scale distinct from the pseudobulk non-rank log2 change. Observed mixture weights and fixed control-state weights were compared on the same common-support population. Support required at least 20 cells per arm/substate and at least 80% of cells in each arm. Unsupported comparisons were NE. Stem and TA had no finer released substate labels.",
76:"The member-level BMP4 display includes every available fixed member in each donor. IBD1 and IBD2 scores were calculated from published gene memberships [14]. Within-family Spearman correlations compared common ligands using either three-donor means or the explicitly identified donor level. Symmetric overlap removal excluded shared genes from both sets and changed their sizes. Correlations were descriptive, without additional P values.",
78:"Analyses used R 4.4.1 with edgeR 4.2.2, DESeq2 1.44.0 and fgsea 1.30.0. Python 3.12.0 with NumPy 2.0.2 and pandas 2.2.3 was used for donor-score, matched-background, pool and composition summaries; Matplotlib was used for figures. FDR adjustment retained the specified analysis families. Negative, discordant, borderline and non-estimable results are reported throughout; NE denotes an unavailable estimate.",
86:"Oxford donor-level estimates were directionally compatible with increased disease-oriented scores, but none of the four programme comparisons passed the four-programme FDR threshold (q range 0.067–0.117; Fig. 2e; Supplementary Fig. S3). UC observations were descriptive. Xenium measured part of each programme in donors overlapping the discovery cohort. Observed IFN-score ranges had limited case–control overlap in the conditional analyses (Supplementary Fig. S2; Table S9).",
90:"Raw changes commonly moved both disease programmes downward in the same donor. Matched random-gene backgrounds also tended to move downward, changing the background-relative directions and classifications. Residual expression/detection imbalance remained visible. Pool deletion produced classification changes and loss of estimability, including BMP6–Goblet when only two donors remained (Supplementary Fig. S5; Table S6).",
95:"BMP4–Colonocyte, an original dual-direction combination with maturation-substate information, showed UP decreases and depleted increases in all three donors across pooled-score representations (Fig. 5a). The member-level display includes all 84 available UP members and 260 of 275 fixed depleted members per donor (Fig. 5b; Table S7). Member medians and direction proportions summarize these log-normalized gene changes.",
96:"The recovered mature-colonocyte fraction increased by 54.9 to 65.5 percentage points across donors (Fig. 5c). On the same common-support cells, the depleted cell ln-normalized change ranged from 0.188 to 0.213 before standardization and 0.0518 to 0.065 after fixing control-state weights. Its direction remained positive, and UP remained negative (Fig. 5d,e). Other evaluated combinations and unsupported states are reported in Supplementary Fig. S7 and Table S7; Stem/TA lacked finer substate annotation, and insufficient common support remained NE.",
99:"Disease-up relationships to IFN and IBD readouts varied with family and representation. Background-relative and individual-donor correlations, common-ligand counts, remaining set sizes and all overlap records are reported in Table S8. The other three family scatterplots appear in Supplementary Fig. S7g–l, complementing the TA displays in Fig. 6c,d.",
101:"Fixed patient gene sets connect recurrent disease associations with responses to experimental ligands. Transfer across patient cohorts supports recurrence of the transcriptional patterns, while heterogeneous organoid responses show the value of evaluating disease-up and depleted readouts jointly and alongside inflammatory and maturation context.",
102:"The three scoring representations preserve different quantities. Raw ranks describe relative position within an expression profile; background-relative ranks describe movement relative to matched gene backgrounds; non-rank changes retain normalized expression magnitude. Their comparison identifies shared directional patterns and disagreements, including changes close to zero, in the same three donors.",
103:"Inflammatory reference perturbations help interpret programme direction. IFN-γ–Stem showed a donor-concordant disease-direction pattern across scoring representations, whereas TNF responses depended on method and family. TNFSF12–TA retained depleted increases with unstable concurrent UP reduction. Prior TWEAK/Fn14 studies also report context-dependent responses in acute injury, chronic colitis, fibroblast signalling and tumour settings [33, 34, 35, 36].",
104:"BMP4-associated colonocyte changes coincided with a larger recovered mature-state fraction. Fixing control-state weights attenuated the depleted change on common support while preserving its sign, linking the programme response to maturation-associated composition. Depleted–IBD2 response correlations persisted after shared genes were removed and across scoring representations.",
105:"Several limitations shape interpretation. Discovery scores are in-sample descriptions; external cohorts are small; IFN-adjusted models have limited observed support; Oxford comparisons did not pass FDR; Xenium has partial coverage and overlapping discovery donors; and GSE210037 lacks coordinates for localisation. The organoid experiment used three healthy donors in an IBD-like inflammatory background, with shared controls and incomplete crossing of ligands and pools. Background matching retained imbalance, and pool deletion sometimes yielded NE. Finer Stem/TA substates and sufficient composition support were unavailable in some comparisons. Scoring representations and post hoc examples are not independent validations. Recovered proportions and fixed-weight comparisons cannot distinguish cell-state transition, selective recovery and cell-intrinsic regulation or establish a causal mediated fraction. IBD membership scores differ from the original atlas activity scores, and low correlations alone do not establish ICI specificity. These transcriptional associations require functional testing before claims about epithelial protection or therapeutic safety.",
106:"Further perturbation experiments that pair fixed transcriptional readouts with barrier, viability and antitumour measurements would clarify their functional relevance.",
108:"Fixed ICI-colitis epithelial programmes connect recurrent patient associations to distinct inflammatory and maturation-associated responses in human organoids. IFN-γ induced disease-concordant changes in stem cells, whereas BMP4 produced opposing programme changes in colonocytes accompanied by a greater recovered mature-state fraction. TNFSF12 increased the depleted programme in TA cells without a consistent accompanying disease-up decrease across scoring representations. Joint evaluation of the two readouts therefore distinguishes depleted-programme increases from concurrent reversal of both programme scores.",
127:"Supplementary_Figures.pdf contains Figures S1–S7 with panel legends. Supplementary_Tables.xlsx contains nine editable result worksheets, Tables S1–S9, with comparison definitions, units, sample counts, statistical scope and table notes."
}
replacements={
56:[('Saved edgeR','edgeR')],
57:[(' Historical BALF comparisons retain their original scope in the project archive.','')],
58:[('descriptive case–control separation and repeatability, not diagnostic performance','descriptive case–control separation')],
59:[('scores reconstructed from saved pseudobulks','pseudobulk-derived scores'),('saved omit-one','leave-one-patient-out'),(' Individual iterations and complete model records remain in the project archive.','')],
64:[(' The source patient mappings and memberships were not changed during integration.','')],
81:[(' These discovery displays describe the data used to define the programmes and are not independent validation.','')],
85:[('All saved percentile intervals','All percentile intervals'),(' These estimates quantify separation in available cohorts rather than clinical diagnostic performance.','')],
88:[('saved GSEA','GSEA')],
89:[('These are method-specific descriptive direction counts, not protective hit rates.','These are method-specific descriptive direction counts.'),(' Pooled-score agreement does not imply agreement after pool deletion or composition standardization.','')],
92:[('IFNG–Stem','IFN-γ–Stem'),('TNFA','TNF')],
93:[('the saved depleted NES','the depleted NES'),('The near-zero positive disease-up values for TA donor DNW9 in the latter representations are near-zero sign changes, not evidence of a large reverse response.','The positive disease-up values for TA donor DNW9 in the latter representations were close to zero.'),(' and are shown in their own units and are not converted into a percentage of recovery',' in their corresponding score units')],
98:[('the saved non-rank','the non-rank'),('; the corresponding mean-response scatterplots do not add a fitted model or significance test','.')],
110:[('ECDF, empirical cumulative distribution function; ',''),(' Source programme identifiers are retained in source tables.','')],
119:[('latest saved organoid sensitivity results','organoid sensitivity results')],
179:[('; arrows indicate analysis links, not causality',''),('These are definition-set descriptions, not independent validation.','These describe the programme-definition cohort.')],
181:[('Saved external','External')],
183:[('Saved joint','Joint')],
185:[('Saved non-rank','Non-rank'),('after IFNG or TNFA treatment','after IFN-γ (IFNG) or TNF (TNFA) treatment'),('The saved GSEA','The GSEA'),(" Displaying non-rank effects here does not redefine the study's primary inference or establish independent replication.",'Non-rank effects provide an expression-scale display alongside the rank-based sensitivity results.')],
187:[('All saved available','All available'),('These panels use saved mean-cell','These panels use mean-cell'),('; the previous ECDF remains in the project archive','')],
189:[('24 saved non-rank','24 non-rank'),('Every point is a saved-score three-donor mean.','Every point is a three-donor mean.'),('The annotated Spearman values come from the existing result table; no regression, P value or new test is added.','Spearman coefficients summarize the displayed ligand responses.'),('Saved depleted','Depleted'),('; the previous matrices remain in the project archive','')]
}
for i,rs in replacements.items():
 t=edits.get(i,original[i])
 for a,b in rs:
  assert a in t,(i,a);t=t.replace(a,b)
 edits[i]=t
for i,t in edits.items():settext(ps[i],t.replace('..','.').replace('tests.Non-rank','tests. Non-rank'))
# Add the requested study-to-resource comparison without inventing a mechanism.
comparison="Capeling et al. established an IBD-linked atlas of epithelial responses to secreted factors, including differentiation-associated BMP responses [14]. ICI-colitis studies characterized patient epithelial and immune changes [8, 11, 12]. Our analysis connects these resources through ICI-derived memberships held fixed across external patient assessments and the ligand atlas. Joint disease-up/depleted responses distinguish isolated depleted-programme increases from concordant opposing changes, while donor-level comparisons and common-support weighting place those responses in inflammatory and maturation context."
np=copy.deepcopy(ps[102]);settext(np,comparison);ps[101].addnext(np)
# Renumber only true bracketed numeric citations, preserving reference identities.
refs={};refnodes={}
for p in ps:
 m=re.match(r'^\[(\d+)\]\s',text(p))
 if m:refs[int(m.group(1))]=text(p);refnodes[int(m.group(1))]=p
body=[p for p in r.findall('.//w:p',NS) if p not in refnodes.values()]
pat=re.compile(r'\[(\d+(?:\s*[,–-]\s*\d+)*)\]')
def nums(s):
 out=[]
 for part in s.split(','):
  a=re.split('[–-]',part.strip());out.extend(range(int(a[0]),int(a[-1])+1))
 return out
order=[]
for p in body:
 for m in pat.finditer(text(p)):
  for n in nums(m.group(1)):
   assert n in refs,('Missing reference',n)
   if n not in order:order.append(n)
mapping={old:new for new,old in enumerate(order,1)}
for p in body:
 if not pat.search(text(p)):continue
 t=pat.sub(lambda m:'['+', '.join(str(mapping[n]) for n in nums(m.group(1)))+']',text(p))
 if t!=text(p):settext(p,t)
# Existing bibliography nodes are preserved except for their leading numbers and order.
parent=refnodes[1].getparent();insert=parent.index(refnodes[1])
for p in refnodes.values():parent.remove(p)
for j,n in enumerate(order):
 p=refnodes[n];old='['+str(n)+']';new='['+str(mapping[n])+']'
 ts=p.findall('.//w:t',NS);done=False
 for t in ts:
  if old in (t.text or ''):t.text=t.text.replace(old,new,1);done=True;break
 if not done:settext(p,re.sub(r'^\[\d+\]',new,text(p),count=1))
 parent.insert(insert+j,p)
# Replace six embedded current previews with the final source-redrawn renderings.
rels=E.fromstring(d['word/_rels/document.xml.rels'])
targets={x.get('Id'):x.get('Target') for x in rels}
blips=r.findall('.//a:blip',NS);assert len(blips)==6
for i,b in enumerate(blips,1):
 target=targets[b.get('{'+NS['r']+'}embed')]
 f=F/'main_figures'/f'Figure_{i}_600dpi.png'
 if f.exists():d['word/'+target]=f.read_bytes()
# Restore clickable repository URLs if prose edits flattened their original runs.
urls=re.compile(r'https://[^\s)]+')
ridmax=max([int(x.get('Id')[3:]) for x in rels if x.get('Id','').startswith('rId') and x.get('Id')[3:].isdigit()]+[0])
for p in r.findall('.//w:p',NS):
 if p is not ps[119]:continue
 t=text(p);matches=list(urls.finditer(t))
 if not matches:continue
 prop=p.find('w:pPr',NS);rp=p.find('.//w:rPr',NS);rp=copy.deepcopy(rp) if rp is not None else None
 for c in list(p):
  if c is not prop:p.remove(c)
 def run(parent,s):
  rr=E.SubElement(parent,W+'r')
  if rp is not None:rr.append(copy.deepcopy(rp))
  tt=E.SubElement(rr,W+'t');tt.text=s;tt.set('{http://www.w3.org/XML/1998/namespace}space','preserve')
 last=0
 for m in matches:
  run(p,t[last:m.start()]);ridmax+=1;rid='rId'+str(ridmax)
  h=E.SubElement(p,W+'hyperlink');h.set('{'+NS['r']+'}id',rid);run(h,m.group())
  E.SubElement(rels,'{http://schemas.openxmlformats.org/package/2006/relationships}Relationship',Id=rid,Type=NS['r']+'/hyperlink',Target=m.group(),TargetMode='External');last=m.end()
 run(p,t[last:])
d['word/_rels/document.xml.rels']=E.tostring(rels,xml_declaration=True,encoding='UTF-8',standalone=True)
save(d,r,O/'Manuscript.docx')
# Synchronize just the permitted wording and caption note in the standalone table.
td,tr,tp=load(I/'Table_1.docx')
for p in tp:
 t=text(p)
 if t==original[44]:settext(p,edits[44])
 if t==original[53]:settext(p,edits[53])
save(td,tr,O/'Table_1.docx')
protected=list(range(12))+list(range(111,118))+list(range(120,126))
assert all(text(ps[i])==original[i] for i in protected)
assert not any('saved' in text(p).lower() for p in r.findall('.//w:p',NS) if p not in refnodes.values())
changes=[{'baseline_paragraph':i,'before':original[i],'after':text(ps[i])} for i in range(len(ps)) if text(ps[i])!=original[i]]
qa={'protected_paragraphs_unchanged':True,'reference_count':len(order),'old_to_new_references':mapping,'removed_uncited_references':sorted(set(refs)-set(order)),'comparison_paragraph':comparison,'changes':changes}
(F/'qa/manuscript_edit_record.json').write_text(json.dumps(qa,ensure_ascii=False,indent=2),encoding='utf8')
(F/'qa/final_paragraphs.txt').write_text('\n\n'.join(text(p) for p in r.findall('.//w:p',NS)),encoding='utf8')
print('Manuscript narrative revised; references',len(order),'removed',qa['removed_uncited_references'])
