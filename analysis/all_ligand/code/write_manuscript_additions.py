from common import *
from analysis_functions import *
import shutil,zipfile,xml.etree.ElementTree as ET
cfg=jr(ROOT/'analysis_config.json');g=pd.read_csv(ROOT/'results/all_ligand_programme_responses.tsv',sep='\t');s=pd.read_csv(ROOT/'results/response_pattern_summary.tsv',sep='\t');bench=pd.read_csv(ROOT/'results/programme_benchmark.tsv',sep='\t');scope=pd.read_csv(ROOT/'results/FDR_scope_counts.tsv',sep='\t');cross=pd.read_csv(ROOT/'results/cross_family_pattern_summary.tsv',sep='\t');li=pd.read_csv(ROOT/'results/ligand_family_pattern_inventory.tsv',sep='\t')
eligible=s[(s.condition_class=='ligand')&(s.pattern!='NE')];both=eligible[eligible.pattern=='BOTH'];counts=eligible.pattern.value_counts();fam=cfg['families'];fl={'Stem':'stem','Transit_Amplifying':'TA','Colonocyte':'colonocyte','Goblet':'goblet'}
def put(name,text):(ROOT/'text'/name).write_text(text.strip()+'\n',encoding='utf8')
def table(df,cols=None,formats=None):
    if cols is not None:df=df[cols]
    formats=formats or {};head='| '+' | '.join(df.columns)+' |\n| '+' | '.join(['---']*len(df.columns))+' |\n';lines=[]
    for _,r in df.iterrows():
        lines.append('| '+' | '.join('NE' if pd.isna(v) else formats.get(c,lambda x:str(x))(v) for c,v in r.items())+' |')
    return head+'\n'.join(lines)
def fm(x):return f'{x:.3g}'
def stat(l,c,m,col):return g[(g.ligand==l)&(g.celltype==c)&(g.module==m)].iloc[0][col]
tn=s[s.ligand=='TNFSF12'].set_index('celltype');main=scope[scope.FDR_scope=='all_ligand_primary_948'].iloc[0]
names='; '.join(fl[f]+': '+', '.join(both[both.celltype==f].ligand) for f in fam if len(both[both.celltype==f]))
common=li[li.n_eligible_families>=2];nvar=int(common.different_patterns_in_observed_families.sum())
enriched_both=[];not_both=[]
for r in both.itertuples():
    ok=stat(r.ligand,r.celltype,UP,'NES')<0 and stat(r.ligand,r.celltype,DEP,'NES')>0 and stat(r.ligand,r.celltype,UP,'FDR_atlas')<.05 and stat(r.ligand,r.celltype,DEP,'FDR_atlas')<.05
    (enriched_both if ok else not_both).append(r)
discordance_text=f"Of the {len(both)} donor-concordant dual-direction comparisons, {len(enriched_both)} also had negative disease-up and positive disease-depleted enrichment with both atlas FDR values below 0.05. The remaining disease-up estimates were "+'; '.join(f"{r.ligand} in {fl[r.celltype]} cells (NES {stat(r.ligand,r.celltype,UP,'NES'):.3f}, FDR {fm(stat(r.ligand,r.celltype,UP,'FDR_atlas'))})" for r in not_both)+'. In particular, the positive disease-up NES for TGFB2 in stem cells contrasted with the three negative raw rank-score changes. Raw within-sample rank changes and size-factor-adjusted gene-set enrichment summarize different quantities; neither result was substituted for the other.'
new_both=[r for r in both.itertuples() if stat(r.ligand,r.celltype,UP,'value_origin')=='FIT_NEW']
origin_text=f"Of these {len(both)} comparisons, {len(both)-len(new_both)} reuse earlier fitted estimates and {len(new_both)} belong to the newly fitted extension ("+'; '.join(r.ligand+' in '+fl[r.celltype]+' cells' for r in new_both)+'). The combined directional classification is a new synthesis, not a claim that all eight underlying responses were newly discovered.'
paragraph=f'''We extended the approved donor-paired organoid workflow to all 79 released ligand conditions. Of 316 planned ligand-by-family comparisons, {len(eligible)} met the existing three-donor coverage criterion and 38 remained non-estimable (NE); 78 ligands had at least one estimable family. The eligible comparisons comprised {', '.join(str((eligible.celltype==f).sum())+' '+fl[f] for f in fam)} comparisons. The no-cytomix condition remained a separate assay reference. Primary enrichment results were adjusted together across all 948 planned ligand-by-family-by-programme tests, including 834 estimated and 114 NE tests. The IFN-gene-excluded residual programme was adjusted in a separate 316-test family (278 estimated, 38 NE), and the assay reference retained its separate 16-test correction. Published doses were retained where unambiguously matched; the released IL22 label could not be reconciled with the IL11 entry in the published dose table, so no dose was assigned to IL22. Responses therefore describe the tested conditions and do not rank pharmacological potency.

Using unrounded treatment-minus-matched-control rank-score changes, {counts['BOTH']}/{len(eligible)} eligible comparisons ({counts['BOTH']/len(eligible)*100:.1f}%) showed a decrease in disease-up and an increase in disease-depleted scores in each of the same three donors. These eight comparisons involved seven ligands: {names}. A further {counts['UP_ONLY']}/{len(eligible)} comparisons showed a disease-up decrease in all three donors without a disease-depleted increase in all three, whereas {counts['DEPLETED_ONLY']}/{len(eligible)} showed the converse pattern. The remaining {counts['OTHER_OBSERVED']}/{len(eligible)} comparisons had other observed directional patterns, not evidence of absent response. These categories describe the three observed donor pairs and are separate from enrichment significance.

Among the {len(common)} ligands with at least two eligible epithelial families, {nvar} exhibited different directional categories across their observed families. Pairwise family comparisons were restricted to their jointly eligible ligands (see the accompanying cross-family table). These observed differences do not constitute a formal treatment-by-family interaction test, and NE comparisons were not treated as evidence of family specificity.'''
tnp=f'''Across the three organoid donors, TNFSF12 decreased disease-up programme scores in TA, colonocyte and goblet cells. A concordant increase in the disease-depleted score was observed only in TA cells among the four families examined. In TA cells, mean paired score changes were {tn.loc['Transit_Amplifying','delta_up_mean']:.6f} for disease-up and {tn.loc['Transit_Amplifying','delta_depleted_mean']:+.6f} for disease-depleted. Depleted-programme enrichment was positive (NES {stat('TNFSF12','Transit_Amplifying',DEP,'NES'):.3f}, atlas FDR {fm(stat('TNFSF12','Transit_Amplifying',DEP,'FDR_atlas'))}), whereas disease-up enrichment was negative (NES {stat('TNFSF12','Transit_Amplifying',UP,'NES'):.3f}, atlas FDR {fm(stat('TNFSF12','Transit_Amplifying',UP,'FDR_atlas'))}).

Colonocyte and goblet cells showed negative disease-up enrichment (NES {stat('TNFSF12','Colonocyte',UP,'NES'):.3f} and {stat('TNFSF12','Goblet',UP,'NES'):.3f}; atlas FDR {fm(stat('TNFSF12','Colonocyte',UP,'FDR_atlas'))} and {fm(stat('TNFSF12','Goblet',UP,'FDR_atlas'))}, respectively). Both had decreased disease-up scores in all three donors, but only two of three donors had an increased disease-depleted score. Stem cells showed positive disease-up enrichment (NES {stat('TNFSF12','Stem',UP,'NES'):.3f}, atlas FDR {fm(stat('TNFSF12','Stem',UP,'FDR_atlas'))}) while disease-up score changes were not uniformly positive across donors. Thus, a donor-concordant reduction in disease-up scores did not ensure a donor-concordant increase in disease-depleted scores (Fig. 4).

IFN-reference enrichment was positive in stem, TA and colonocyte cells and negative in goblet cells. IFN responses remained a separate readout and were not incorporated into a composite reversal score. TNFSF12 was one of seven ligands with a qualifying dual-direction comparison, rather than a uniquely effective or best-ranked candidate. The underlying TNFSF12 estimates and donor scores were reused unchanged; only their all-ligand multiple-testing context and presentation were updated.'''
bb=bench[bench.stratum=='MEAN_OF_3'].copy()
def rho(di,re):return ', '.join(f'{bb[(bb.celltype==f)&(bb.disease_readout==di)&(bb.reference==re)].iloc[0].spearman_rho:.3f}' for f in fam)
bp=f'''Across ligand mean paired changes, the relationship with the fixed IFN reference differed by epithelial family and readout. In stem, TA, colonocyte and goblet cells, respectively, Spearman correlations with IFN were {rho('delta_up','delta_IFN')} for disease-up and {rho('delta_depleted','delta_IFN')} for disease-depleted, based on 76, 77, 77 and 48 common eligible ligands. Pearson coefficients and donor-stratified results are reported alongside these summaries; for example, stem disease-up versus IFN had Pearson r = {bb[(bb.celltype=='Stem')&(bb.disease_readout=='delta_up')&(bb.reference=='delta_IFN')].iloc[0].pearson_r:.3f}, illustrating that linear and rank associations need not agree.

We also scored the published IBD1 and IBD2 membership sets on the same matched pseudobulks using the same mean expression-rank rule. Exact-symbol coverage was 106/113 unique IBD1 members and 86/91 unique IBD2 members. The published IBD2 list contained CA1 twice; it was counted once as a set member. These membership-based comparisons do not reproduce the publisher's original module-activity scores or learned programme weights. Disease-depleted changes correlated positively with IBD2 membership responses across all four families (Spearman ρ = {rho('delta_depleted','delta_IBD2')}); corresponding correlations with IBD1 were {rho('delta_depleted','delta_IBD1')}. The substantial overlap, especially for disease-depleted versus IBD2, argues against describing that readout as wholly novel or independent of established intestinal programmes. Disease-up correlations with IBD1 ({rho('delta_up','delta_IBD1')}) and IBD2 ({rho('delta_up','delta_IBD2')}) were weaker or differed in sign.

All continuous responses and directional discordances are supplied, with examples chosen by the rule locked before analysis: the alphabetically first eligible ligand within each family and directional category, plus the alphabetically first within each family and disease/IFN discordance category. Correlations are descriptive; shared donors and controls preclude treating ligand-by-family cells as independent biological replicates. No independent-grid P values, new patient-level analyses, or new clinical-performance estimates were calculated.'''
put('Results_additions.md',f'''# Results additions and Figure 4 replacement text

**Status:** generated from this incremental run. This file is a proposed manuscript patch; it does not overwrite the approved Word manuscript. Figure 4 is the replacement figure; other new figures use descriptive names pending whole-manuscript numbering.

## New subsection: Ligand responses separate disease-up attenuation from disease-depleted score recovery

{paragraph}

{discordance_text}

{origin_text}

## Replacement core of section 3.4: TNFSF12-induced programme changes differ across epithelial families

{tnp}

The separate no-cytomix reference showed both score directions in all three donors in stem, TA and colonocyte cells. In goblet cells, only the disease-depleted readout met the three-donor directional criterion, despite a negative mean disease-up change. These reference results were retained separately from the ligand screen.

## New subsection: Comparisons with IFN and published IBD programme membership

{bp}

## Exact cross-family comparison denominators

{table(cross,['family_a','family_b','n_common_ligands','n_different_pattern'])}

## Scope and retained material

- The 257 newly fitted comparisons are new analyses. Twenty-one ligand models and four no-cytomix models, including all TNFSF12 estimates, were reused without refitting.
- Historical targeted FDR remains in `old_targeted_FDR` and `old_targeted_FDR_family`; the text and Figure 4 above use `FDR_atlas` (948 planned primary tests).
- Preserve the original selected-set transfer results and original Figure 4c in a supplementary location; do not delete their historical results paragraph. All seven original nominated ligands and the assay reference remain in the complete tables.
- Preserve existing patient replication, IFN-conditional, TCR, tumour/off-target, negative, discordant and NE results. No claims in those sections were recalculated here.
''')
put('Methods_additions.md',f'''# Methods patch for the expanded analysis

## Replace the targeted-only description of multiplicity where the new atlas is introduced

The expanded organoid analysis retained the frozen programme definitions, released epithelial family labels, donor/pool matching and count filters used in the approved targeted workflow. Treatment and pooled cytomix-only control cells were matched within donor, pool and cell family, retaining a pool only when each arm contained at least 20 released cells. Counts from eligible pools and technical control wells were summed within donor. Models required three paired donors and used DESeq2 1.44.0 with design ~ donor + treatment, condition-blind total count ≥10, no log-fold-change shrinkage, no automatic count replacement, and independentFiltering=FALSE. All finite Wald statistics with available gene-level P values formed the rank background. Enrichment used fgseaSimple 1.30.0, 10,000 permutations, minimum set size 10, maximum 500, scoreType='std', nproc=1 and seed 20260826. Saved estimates and scores for 25 already fitted comparisons were reused; only the 257 missing eligible comparisons were fitted.

The primary BH correction contained 948 planned tests (79 ligands × 4 families × 3 programmes: disease-up, disease-depleted and IFN), retaining 114 NE tests in the correction denominator. Residual sensitivity used 316 planned tests, including 38 NE, and no-cytomix used a separate 16-test family. Historical targeted adjustments were archived separately. The label FDR in the expanded Figure 4 refers to the 948-test family.

Programme scores were the unweighted mean of within-pseudobulk expression ranks divided by the full 36,912-gene background, with average ranks for ties. Donor changes were treatment minus matched control on the raw score scale; disease-depleted scores were not inverted. Directional categories used strict inequalities on unrounded paired deltas from the same three donors. Zero was counted as neither positive nor negative. These categories introduce no new significance test.

The IBD membership comparison used unique members from the publisher's supplementary workbook (`41467_2025_68247_MOESM5_ESM.xlsx`, sheet `B Genes in IBD1 and IBD2`) with the same rank rule. CA1's duplicated IBD2 entry was collapsed once; exact symbols were used without alias conversion and missing members were reported. No weighted cNMF or publisher activity scores were silently converted into unweighted scores; this is explicitly a membership-set benchmark. Within-family Pearson and Spearman coefficients used common eligible ligands, with three-donor means and donor-stratified summaries. No P values assuming independent ligands were computed.

## Local source caveats

The published dose table (`41467_2025_68247_MOESM16_ESM.docx`) lists IL11, whereas the released object and pinned processing code label IL22. IL22's dose is unresolved; its analysis label is retained as released. All other mapped concentrations, vendors and catalogue entries are preserved verbatim, not corrected from assumptions. These data do not support potency comparisons across doses.

## Runtime report

Full commands, times, exit codes, input hashes and session records accompany the run. The current R installation exhibited an abnormal process exit after successful completion output; refer to the execution records for the exact process status. Completed scientific outputs are checked separately from runtime exit status. No package was installed or upgraded.
''')
put('Figure_legends_additions.md','''# Figure legends for the actual generated figures

## Figure 4. Donor-resolved TNFSF12 responses of disease-up and disease-depleted epithelial programmes

a, Normalized enrichment scores (NES) for the disease-up, disease-depleted and IFN-reference programmes across stem, transit-amplifying (TA), colonocyte and goblet cells. Genes were ranked by treatment Wald statistics from donor-adjusted DESeq2 models with three paired donors. Cell colours encode NES; printed values give NES and atlas FDR (q). Outlined cells have FDR < 0.05 under BH correction across all 948 planned ligand-by-family-by-primary-programme tests, including NE tests in the denominator. These are not the historical 16-test targeted TNFSF12 FDR values. b,c, Treatment-minus-matched-control changes in disease-up (b) and disease-depleted (c) mean expression-rank scores. Each point is one donor pair; black horizontal ticks show family means. d, Joint donor-level score changes, with −Δ disease-up on the horizontal axis and Δ disease-depleted on the vertical axis. The upper-right quadrant indicates a lower disease-up and higher disease-depleted score in the same donor–family observation. Colours denote families consistently across b–d, and point shapes denote the same three donors. Score axes in b–d are displayed in units of 10⁻³. Directional concordance is descriptive and distinct from enrichment significance. The IFN-excluded residual and no-cytomix reference are retained in the supplementary figure and complete source tables. TNFSF12 was tested at the published 250 ng/mL concentration; no potency comparison is implied.

## New figure: Complete all-ligand response atlas

All 79 released ligand conditions are shown alphabetically over two pages, with four epithelial families and three primary readouts. Colour indicates the mean treatment-minus-control raw rank-score change across three paired donors; each programme has its own symmetric colour scale, fixed across the two pages. Numbers inside cells give the number of donors with disease-up Δ < 0, disease-depleted Δ > 0, or IFN Δ > 0, respectively. Black cell outlines indicate enrichment FDR < 0.05 under the single 948-test primary family. Thus colour, enrichment significance and donor direction counts encode distinct quantities. Gray NE cells lack three eligible paired donors and are neither zero responses nor negative findings. All 316 planned ligand-by-family comparisons remain visible; no-cytomix is a separate reference. Exact dose matches and the unresolved IL11/IL22 dose-table discrepancy are recorded in the source data. Raw rank-score magnitudes are not calibrated effect sizes across programmes or potency measures across ligands.

## New figure: Dual-readout donor patterns across epithelial families

Colours classify strict unrounded paired score changes. BOTH requires disease-up Δ < 0 and disease-depleted Δ > 0 in each of the same three donors. UP_ONLY or DEPLETED_ONLY indicates that only the respective readout meets that three-donor directional condition. OTHER_OBSERVED retains all other observed directions; it does not mean absent response. Numbers in cells show donors satisfying both directions simultaneously. NE indicates fewer than three pairs. Summary bars count the 316 planned comparisons, of which 278 are eligible. Shared donors and controls do not create independent replicates across matrix cells. Cross-family pattern comparisons use jointly eligible ligands only.

## New figure: Reference programme comparisons

Page 1 shows within-family Spearman correlations between ligand mean paired changes in each ICI readout and IFN, IBD1 membership or IBD2 membership responses. n is the number of common eligible ligands. Pages 2–3 show every corresponding ligand mean as a point; the open black ring identifies TNFSF12. Raw axes use 10⁻³ rank-score units and retain programme-specific scales. IBD comparisons use the publisher's unique member lists with the same mean expression-rank rule and exact-symbol coverage of 106/113 and 86/91 members, respectively. They do not reproduce the publisher's original activity or weighted cNMF scores. The source tables include Pearson and Spearman coefficients, donor-stratified coefficients and all continuous/sign-discordant pairs. No correlation P values assuming independent ligands are reported.

## Supplementary figure: Prespecified donor-level display examples

Each panel shows the three actual donor pairs for one ligand–family combination in the joint −Δ disease-up versus Δ disease-depleted plane. Panel headings also report donor counts for UP Δ < 0 (UP−), depleted Δ > 0 (Dep+) and IFN Δ > 0 (IFN+), each out of three. Examples were chosen using the rule fixed before results: alphabetical first within each family and joint-response category, plus alphabetical first within each family and disease/IFN direction-discordance category. Some panels fulfil more than one selection category; all category assignments are preserved in the selection source table. Donor shapes and family colours match Figure 4. These are descriptive illustrations, not selections of the most favourable effect or FDR. Every eligible and NE combination remains in the complete tables.

## Supplementary figure: IFN-excluded residual sensitivity and no-cytomix reference

Page 1 retains residual-programme enrichment in TNFSF12 and the no-cytomix reference. Filled points have FDR < 0.05. TNFSF12 residual FDR uses the separate 316-test all-ligand residual family; no-cytomix FDR uses the separate 16-test reference family. All-ligand residual estimates and coverage are provided in source tables. Page 2 displays donor-resolved raw paired changes in the three primary programmes for no-cytomix; black ticks denote family means. No-cytomix is an assay reference and is excluded from the 79-ligand total.

## Historical Figure 4c transfer

Retain the approved historical selected-set transfer panel and its original data as a supplementary panel, explicitly labelled as the earlier selected-set analysis. Its existing legend and selection restrictions must move with it. Do not relabel those historical values as new all-ligand estimates. The current Word manuscript and its embedded figure media are preserved unchanged in inputs; no historical text or image has been deleted.
''')
put('Conclusions_options.md',f'''# Conclusions options based on the completed incremental outputs

## Recommended complete replacement for Conclusions

Fixed epithelial programmes of ICI colitis retained their disease-associated directions in the previously analysed external patient datasets. Extending the approved organoid workflow across the 79 released ligands yielded 278 eligible ligand–family comparisons, of which eight showed both a disease-up score decrease and a disease-depleted score increase in each of three observed donors. Eighty-seven comparisons met only the disease-up directional criterion and seven met only the disease-depleted criterion. Responses differed across observed epithelial families, and TNFSF12 was one example rather than a uniquely effective perturbation. The disease-depleted readout substantially overlapped with a published IBD2 membership reference, while relationships with IFN varied by family. Together, these findings support reporting the two transcriptional directions and their reference-programme overlap separately. They do not establish functional rescue, clinical efficacy, or pharmacological safety.

## Shorter option

Across 278 eligible ligand–epithelial-family comparisons, donor-concordant attenuation of disease-up scores often occurred without a concordant increase in disease-depleted scores. Eight comparisons showed both directions in all three observed donors. These family-dependent transcriptional patterns, including TNFSF12 responses, support joint assessment of the two readouts alongside IFN and IBD references. The results describe tested organoid conditions and do not establish functional recovery or therapeutic benefit.

## Interpretation constraints retained

- Three shared healthy donors and reused controls; no grid-level independence or formal family-interaction claim.
- Thirty-eight NE comparisons remain missing for inference, including all four families for TNFSF10.
- Heterogeneous tested doses; IL22 dose remains unresolved because of the source IL11/IL22 discrepancy.
- IBD2 overlap is acknowledged, not presented as independent novelty.
- Keep existing tumour/off-target and negative guardrails. No new safety evidence was generated.
''')
cc=pd.read_csv(ROOT/'results/response_pattern_counts.tsv',sep='\t');small=cc.pivot(index='celltype',columns='pattern',values='n').reindex(index=['ALL']+fam,columns=['BOTH','UP_ONLY','DEPLETED_ONLY','OTHER_OBSERVED','NE']).reset_index();bt=[]
for r in both.itertuples():
    bt.append(dict(ligand=r.ligand,family=r.celltype,mean_UP=r.delta_up_mean,mean_depleted=r.delta_depleted_mean,UP_NES=stat(r.ligand,r.celltype,UP,'NES'),UP_atlas_FDR=stat(r.ligand,r.celltype,UP,'FDR_atlas'),depleted_NES=stat(r.ligand,r.celltype,DEP,'NES'),depleted_atlas_FDR=stat(r.ligand,r.celltype,DEP,'FDR_atlas'),value_origin=stat(r.ligand,r.celltype,UP,'value_origin')))
both_table=pd.DataFrame(bt);save(both_table,'results/both_direction_combinations.tsv')
put('ANALYSIS_REPORT_ZH.md',f'''# 全配体增量分析结果

## 执行范围

工作目录：`{ROOT}`。复用 25 个已完成模型，只新增 257 个合格比较的拟合。原始数据、旧结果及当前 Word 主稿均未修改；未执行 RUN_ALL，未新增患者、TCR、肿瘤或安全性分析。

## 1. 双读出分类

{table(small)}

分母为 278 个可估计组合，而不是把 316 个计划格子都当作独立重复。BOTH 为同一三个供体都同时 UP<0、depleted>0。UP_ONLY=87、DEPLETED_ONLY=7，说明单侧一致并不保证双侧一致。OTHER_OBSERVED 不是“无效”。

双读出一致的 8 个组合涉及 7 个配体：

{table(both_table,formats={c:(lambda x:f'{x:.6g}') for c in both_table.columns if c not in ['ligand','family','value_origin']})}

其中 {len(both)-len(new_both)} 个组合复用旧拟合估计，{len(new_both)} 个属于本次新增拟合（{'；'.join(r.ligand+'/'+r.celltype for r in new_both)}）。联合分类是本次新整理，不把所有已有反应都说成新发现。

方向分类与富集 FDR 是不同信息。{len(both)} 个双方向组合中，只有 {len(enriched_both)} 个同时具备 UP 负富集、depleted 正富集且两者 atlas FDR<0.05。INS/colonocyte、TGFB2/stem、TNFSF12/TA 的 UP 未过校正；其中 TGFB2/stem 的 UP NES 为正而三个供体的原始秩变化为负。这种估计量间的不一致完整保留，不据此宣称临床有效或药理效力排序。

## 2. 家族差异

至少有两个可估计家族的 {len(common)} 个配体中，{nvar} 个出现不同的观察方向类别。这是直接观察的类别差异，不是交互作用检验。

{table(cross,['family_a','family_b','n_common_ligands','n_different_pattern'])}

## 3. 参考程序比较

{table(bb,['celltype','disease_readout','reference','n_common_ligands','pearson_r','spearman_rho'],formats={'pearson_r':lambda x:f'{x:.3f}','spearman_rho':lambda x:f'{x:.3f}'})}

Disease-depleted 与 IBD2 成员集合在各家族有较强重叠，TA 最高；不能将其包装成完全独立的新信息。UP 与 IFN 的关系较弱或异向，而且 Pearson 与 Spearman 并不总同号。供体分层系数、连续散点和全部不一致组合均已输出；未基于独立格子假设计算显著性。

IBD1 为 113 个唯一成员、可用 106 个；IBD2 原表 92 行含 CA1 重复，去重后 91 个、可用 86 个。本次明确使用集合型秩评分，不等同于原作者 cNMF/模块活动分数。缺失符号逐一列在 `IBD_reference_coverage.tsv`，没有自行改名补齐。

## 4. 多重检验

{table(scope)}

旧 targeted FDR 与本次 atlas FDR 在表中分列。No-cytomix 的原 16-test FDR 保持一致。新版 Figure 4 用 948-test 主范围；residual 用单独 316-test 范围。

## 5. 局部待确认与技术状态

- IL22：对象和原发布代码为 IL22，但配体剂量表为 IL11；剂量缺项不影响按已发布标签完成的比较，不能写成已确认的 IL22 浓度。
- R 运行状态：样例及进程的计算完成标志与退出码分别记录；本机存在载入 fgsea 后退出异常的可复现现象。不能把非零进程退出伪报为 PASS。详见 `review/FINAL_STATUS.json` 和 `logs/execution_*.json`。
- 文字交付为基于实际结果生成的 Results、Methods、图注和 Conclusions 修改稿；未擅自覆盖当前主稿。历史 Figure 4c 的图和段落保留，提供移入补充的明确建议。
''')
put('MANUSCRIPT_CHANGE_MAP.md','''# Manuscript patch map

| Location | Deliverable | Action |
| --- | --- | --- |
| Current Methods 2.8 organoid scope | Methods_additions.md | Expand targeted scope to the exact atlas and separate FDR families; retain original model rules |
| New Results atlas subsection | Results_additions.md, first subsection | Insert counts, dual-readout categories and family comparison |
| Current Results 3.4 TNFSF12 core | Results_additions.md, replacement core | Replace with actual atlas FDR and donor directions |
| New Results reference comparison | Results_additions.md, third subsection | Add IFN/IBD membership comparisons with coverage and overlap |
| Figure 4 | Figure4_TNFSF12.pdf and legend | Replace a–d layout in a future manuscript copy |
| Historical Figure 4c | Approved original media in inputs | Move intact to a supplementary location with original legend and selection restrictions |
| Residual and no-cytomix | Supplementary_Residual_and_AssayReference.pdf | Retain sensitivity and assay reference |
| Conclusions | Conclusions_options.md | Choose recommended or short evidence-based version |

Other new figures carry descriptive names. Do not silently renumber unchanged patient/TCR/tumour figures or strip their results. A final full-manuscript placement step must assign those figure numbers and update references after the author chooses the main/supplementary placement; this delivery provides actual figures and text, not a manual analysis to-do list.
''')
print('MANUSCRIPT_ADDITIONS_COMPLETE',flush=True)
