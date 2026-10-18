from common import *
import shutil
cl=read(ROOT/'results/classification_stability.tsv');e=read(ROOT/'results/three_method_donor_effects.tsv');comp=read(ROOT/'results/composition_standardized_expression.tsv');lp=read(ROOT/'results/leave_one_pool_out_classification.tsv')
def combo(b):return '; '.join(r.ligand+'-'+r.family.replace('Transit_Amplifying','TA') for r in b.itertuples())
def mdtable(d,columns=None):
 if columns:d=d[columns]
 vals=[[str(v) for v in row] for row in d.itertuples(index=False,name=None)];return '| '+' | '.join(d.columns)+' |\n| '+' | '.join(['---']*len(d.columns))+' |\n'+'\n'.join('| '+' | '.join(r)+' |' for r in vals)+'\n'
def write(name,s):(ROOT/'text'/name).write_text(s,encoding='utf8')
# Add explicit composition and pool-sensitivity classifications without changing the three primary classifications.
extra=[]
for r in cl.itertuples():
 b=comp[(comp.ligand==r.ligand)&(comp.celltype==r.family)];row={}
 for col in ['unstandardized_full_delta','unstandardized_common_delta','standardized_common_delta']:
  w=b[b.programme.isin([UP,DEP])].pivot(index='donor',columns='programme',values=col).reindex(index=DONORS,columns=[UP,DEP]);row['composition_'+col+'_class']=pattern(w[UP],w[DEP])
 row['composition_estimable_donors']=b.loc[b.status=='ESTIMATED','donor'].nunique()
 if not b.empty and (b.annotation_resolution=='MULTI_SUBSTATE').all() and not (b.status=='ESTIMATED').any():cl.loc[(cl.ligand==r.ligand)&(cl.family==r.family),'composition_assessment_status']='NE_COMMON_SUPPORT_ALL_DONORS'
 l=lp[(lp.ligand==r.ligand)&(lp.celltype==r.family)];est=l[l.status=='ESTIMATED'];row['pool_deletion_scenarios']=len(l);row['pool_deletion_estimable_scenarios']=len(est)
 for me in ['raw','bg_relative','nonrank']:row['pool_deletion_'+me+'_Both_count']=int((est[me+'_class']=='BOTH').sum())
 row['all_three_methods_Both_in_every_planned_pool_deletion']=bool(len(l)>0 and len(est)==len(l) and (est[['raw_class','bg_relative_class','nonrank_class']]=='BOTH').all().all())
 extra.append(row)
for k in extra[0]:cl[k]=[r[k] for r in extra]
save(cl,'results/classification_stability.tsv');save(cl[cl.original_both],'results/original_eight_Both_tracking.tsv')
orig=cl[cl.original_both];stable=orig[orig.raw_bg_nonrank_both];changed=orig[~orig.raw_bg_nonrank_both];new=cl[cl.new_bg_both|cl.new_nonrank_both]
counts=read(ROOT/'results/classification_counts_by_method.tsv');counts=counts[counts.family=='ALL'];ct=counts.pivot(index='classification',columns='method',values='n_combinations').reindex(index=['BOTH','UP_ONLY','DEPLETED_ONLY','OTHER_OBSERVED','NE'])
ct=ct[['raw','bg_relative','nonrank']].reset_index();marg=read(ROOT/'results/raw_marginal_direction_counts.tsv');marg=marg[marg.family=='ALL'];joint=read(ROOT/'results/raw_joint_direction_table.tsv');joint=joint[joint.family=='ALL'];dd=read(ROOT/'results/raw_same_donor_joint_counts.tsv').query("family=='ALL' and donor=='ALL'").iloc[0]
rs=read(ROOT/'results/random_reference_summary.tsv');rs0=rs[rs.stratum=='module'].set_index('module');balf=read(ROOT/'results/matching_balance.tsv');balmain=balf[~balf.module.str.contains('__without__')];imb=balmain.loc[balmain.standardized_mean_difference.abs().idxmax()]
fam=read(ROOT/'results/family_response_summary.tsv');fb=fam[(fam.family=='ALL')&(fam.classification=='BOTH')]
po=read(ROOT/'results/pool_background_summary.tsv');po=po[po.ligand!='No_Cytomix'];reuse=read(ROOT/'results/control_reuse_summary.tsv');design=read(ROOT/'results/pool_family_design_summary.tsv')
genes,full,sets=definitions();ov=read(ROOT/'results/gene_set_overlap.tsv');corr=read(ROOT/'results/reference_correlations_by_method.tsv')
def meansrow(lig,fam,mod,method):return e[(e.ligand==lig)&(e.celltype==fam)&(e.module==mod)][method]
def signstr(v):return f'{int((v>0).sum())} positive / {int((v<0).sum())} negative / {int((v==0).sum())} zero'
reft=[]
for (l,f),b in e[e.ligand.isin(['IFNG','TNFA'])].groupby(['ligand','celltype']):
 for me in ['raw','bg_relative','nonrank']:reft.append(dict(ligand=l,family=f,method=me,up=signstr(b[b.module==UP][me]),depleted=signstr(b[b.module==DEP][me]),IFN=signstr(b[b.module==IFN][me])))
rt=pd.DataFrame(reft);save(rt,'results/reference_strict_direction_counts.tsv')
coverage=pd.DataFrame([dict(programme=SHORT[m],frozen_members=len(full[m]),measured_members=len(sets[m]),missing_members=';'.join(g for g in full[m] if g not in genes)) for m in MODS]);save(coverage,'results/fixed_programme_coverage.tsv')
gsea=read(BASE/'results/all_ligand_programme_responses.tsv');shutil.copy2(BASE/'results/all_ligand_programme_responses.tsv',ROOT/'results/original_GSEA_results.tsv')
gs=gsea.merge(orig[['ligand','family']],left_on=['ligand','celltype'],right_on=['ligand','family']);gs=gs[gs.module.isin([UP,DEP])][['ligand','celltype','module','NES','FDR_atlas','FDR_scope','FDR_planned_n']];save(gs,'results/original_Both_GSEA_context.tsv')
depibd=corr[(corr.donor=='MEAN_OF_3')&(corr.programme_a==DEP)&(corr.programme_b==IBD2)];dpiv=depibd.pivot(index=['celltype','membership'],columns='method',values='Spearman').round(3).reset_index()
bm=read(ROOT/'results/within_family_substate_proportions.tsv');bmp=bm[(bm.ligand=='BMP4')&(bm.celltype=='Colonocyte')&(bm.substate=='Colonocyte_Mature')]
bc=comp[(comp.ligand=='BMP4')&(comp.celltype=='Colonocyte')&(comp.programme==DEP)]
body=f'''# Proposed Results update

Draft for author integration; the manuscript has not been overwritten. This is a post hoc extension performed after viewing the all-ligand atlas, the original Both calls and the TGF-beta concentration. It does not constitute preregistration or independent validation.

## Directional asymmetry also occurred in expression/detection-matched gene backgrounds

The unchanged inventory contained {len(cl)} planned ligand-family combinations, of which {int(cl.eligible.sum())} met the existing three-donor criteria and {int((~cl.eligible).sum())} were non-estimable. These combinations were not independent biological replicates. For 0/1/2/3 concordant donors, disease-up decreases occurred in {'/'.join(map(str,marg[marg.axis=='up_down'].sort_values('donor_count').n_combinations))} combinations, whereas disease-depleted increases occurred in {'/'.join(map(str,marg[marg.axis=='depleted_up'].sort_values('donor_count').n_combinations))}. The observed joint distribution included {int(joint[(joint.up_down_donor_count==3)&(joint.depleted_up_donor_count==0)].n_combinations.iloc[0])} combinations with disease-up decreasing in all donors but disease-depleted increasing in none. Both raw programmes decreased in the same donor in {int(dd.same_donor_both_decline)}/{int(dd.n_matched_pairs)} matched ligand-family-donor rows, demonstrating that these two directional summaries cannot be inferred from separate margins.

We retained the frozen members and generated 1,000 joint expression/detection-quantile-matched nonmember sets per programme and epithelial family, with donor-equal baselines from deduplicated matched cytomix controls. The background median B was negative in every eligible ligand donor pair for both disease programmes ({int(rs0.loc[UP,'n'])} pairs each). Averaged across pairs, {rs0.loc[UP,'random_set_negative_fraction_mean']:.1%} of disease-up reference-set deltas and {rs0.loc[DEP,'random_set_negative_fraction_mean']:.1%} of disease-depleted reference-set deltas were negative. This shared asymmetry was present across all four families and all three donors; pool-level medians were also negative in all {len(po[po.module==UP])} observed eligible pool pairs for both programmes. Pool summaries include partial donor coverage and therefore use a different denominator from the three-donor atlas. These findings demonstrate a widespread directional background; they do not identify its technical or biological cause.

Raw disease-up deltas were negative in {rs0.loc[UP,'raw_negative_fraction']:.1%} of eligible donor rows, compared with {rs0.loc[UP,'background_relative_negative_fraction']:.1%} after subtracting B. Corresponding disease-depleted proportions were {rs0.loc[DEP,'raw_negative_fraction']:.1%} and {rs0.loc[DEP,'background_relative_negative_fraction']:.1%}. The three representations answer distinct questions: paired rank change relative to cytomix, deviation from an expression/detection-matched gene background, and mean fixed-member log2 normalized-expression change. A positive background-relative disease-depleted delta can arise when the programme declines less than its background; it is not necessarily an actual increase.

Matching preserved the available joint quantile strata without replacement, but did not eliminate continuous covariate imbalance within those strata. The largest absolute baseline standardized mean difference among the six main programmes was {abs(imb.standardized_mean_difference):.3f} ({SHORT[imb.module]}, {imb.celltype}, {imb.dimension}; standardization uses the family-wide gene-universe SD). IBD2 was particularly imperfectly matched. These are coarse-stratum descriptive reference comparisons, not fully balanced causal adjustments. No new response-based rematching, P value, or donor confidence interval was used; gene-set percentiles remain descriptive. All strata, continuous balance, random memberships and raw effects are supplied.

## A subset of the original Both patterns was consistent across the three pooled scoring methods

The same strict unrounded three-donor rule yielded the following classifications:

{mdtable(ct)}

Of the {len(orig)} original Both combinations, {len(stable)} retained Both in both additional scoring methods: {combo(stable)}. {combo(changed)} became depleted-only under both background-relative and nonrank scoring. Each original combination is shown below; the complete planned atlas retains every NE and exclusion reason.

{mdtable(orig,['ligand','family','raw_class','bg_relative_class','nonrank_class','composition_assessment_status'])}

Additional Both calls, which were not original raw Both, were also retained:

{mdtable(new,['ligand','family','raw_class','bg_relative_class','nonrank_class'])}

Only TNFSF11-Colonocyte was newly Both in both background-relative and nonrank scoring; it was not an original raw Both. These are directional classifications, not significance claims. The original GSEA NES and q values are preserved separately in original_GSEA_results.tsv, with the original planned correction families. They were not used as significance for the added scores.

## Pool dependence and epithelial composition qualify the pooled patterns

The eligible pool data reused {len(reuse)} distinct donor-family-pool cytomix controls, with {int(reuse.n_ligands.min())}-{int(reuse.n_ligands.max())} tested ligands per control. Shared controls induce dependence but are not, by themselves, evidence of biased mean treatment effects. Across all requested deletion scenarios, {int((lp.status=='ESTIMATED').sum())}/{len(lp)} retained three matched donors; remaining scenarios were NE, including originally incomplete combinations. The eligible comparison set therefore changed after deletion.

BMP4-Colonocyte retained Both under all three methods after deletion of either of its two original pools. BMP6-Goblet retained this pattern in its one estimable deletion scenario, while deletion of Pool 3 left only two donors. INHBA-TA, INS-Colonocyte and TGFB1-TA did not retain all-method Both across both pool deletions. TGFB2-TA retained raw Both in each individual pool but not background-relative Both in either. TNFSF12-TA also varied across the two pools and scoring methods. The complete pool and deletion tables retain these disagreements; the pool-relative median is a separate descriptive comparator and never replaces cytomix controls or adds to B subtraction.

Released mature/immature substate annotations supported composition summaries for Colonocyte and Goblet. Stem and TA had no finer final substate annotation, so their nominally unchanged standardized scores do not exclude within-family heterogeneity. Common-support standardization was often unavailable for Goblet, including all three BMP6-Goblet donors. No missing substate expression was imputed. All 20 existing publisher perturbation-derived cNMF usages were summarized with their published annotations; they were not re-fitted or treated as newly defined maturation or stemness signatures.

For BMP4-Colonocyte, the recovered mature-cell fraction increased by {bmp.proportion_delta.min():.3f}-{bmp.proportion_delta.max():.3f} across donors. On the separate cell LogNormalize scale, disease-depleted deltas decreased from {bc.unstandardized_full_delta.min():.3f}-{bc.unstandardized_full_delta.max():.3f} before weighting to {bc.standardized_common_delta.min():.3f}-{bc.standardized_common_delta.max():.3f} after fixed control-substate weighting, while disease-up remained negative and disease-depleted remained positive in each donor. This is consistent with a substantial differentiation/composition-associated component plus remaining within-supported-substate transcriptional changes. It is not evidence of functional protection. INS-Colonocyte did not show a three-donor Both pattern on this cell-level expression scale even before weighting; disease-depleted change became negative in two donors after weighting. Its sensitivity therefore cannot be attributed solely to composition. TGFB2-TA likewise had mixed disease-up direction on the cell-level scale, despite its pooled-pseudobulk three-method Both classification.

## TGF-beta superfamily and disease-direction references

All {len(TGF)} task-specified BMP/GDF/INHBA/TGFB ligands were retained as a structural ligand-superfamily group, without implying a common receptor or SMAD pathway. Family coverage and unique-ligand hit counts were:

{mdtable(fb,['ligand_group','method','n_combinations','n_eligible_combinations','n_ligands_with_category','n_ligands_with_eligible_family'])}

The descriptive concentration remained visible across scoring methods, but did not establish superfamily enrichment under an independent-combination model. Group members occurred across all original pools, with each ligand generally assigned to one linked pool pair. Thus the observation was not confined to one pool, but the design does not fully cross ligand identity with pool; intrinsic superfamily and pool effects cannot be uniquely separated. No combination- or ligand-level Fisher P value was reported. Pool-dependent effects and substate changes remain part of the interpretation.

IFNG (IFN-gamma) and TNFA (TNF-alpha) were retained as internal disease-direction reference perturbations in the existing inflammatory background. In Stem, both had strictly positive disease-up and IFN changes in all three donors under all methods. IFNG also retained strictly negative disease-depleted changes in all three donors under all methods. TNFA retained negative disease-depleted changes in all three donors under raw and nonrank scoring, but its background-relative disease-depleted direction was positive in two donors and negative in one. Thus TNFA's complete disease-direction pattern was not invariant to background referencing. In TA, raw disease-up decreased in one donor for each perturbation; background-relative disease-up increased in all three for both, while nonrank disease-up increased in two donors for IFNG and all three for TNFA. TNFA-Colonocyte remained mixed; IFNG-Colonocyte and both Goblet comparisons were NE. Every strict sign count is supplied in reference_strict_direction_counts.tsv. These are not new independent positive-validation experiments.

## Reference-programme relationships were method-dependent and not explained solely by overlap

The measurable disease-up set contained {len(sets[UP])} members and the disease-depleted set {len(sets[DEP])}. Measurable disease-up overlap was 4 genes with IFN, 2 with IBD1 and none with IBD2; disease-depleted overlap was 6, 12 and 17 genes, respectively. Full frozen-set overlap and unavailable symbols are reported separately. All membership-based comparisons retained the membership-score label for IBD1/IBD2.

Disease-depleted versus IBD2 Spearman correlations across eligible ligand means remained positive after symmetrically removing shared genes, although both method and member-set size affected their values:

{mdtable(dpiv)}

Disease-up relationships with IFN/IBD1 were generally weak or negative and changed with method; removing the four shared IFN members did not establish statistical independence or ICI specificity. The common ligand IDs, donor-stratified coefficients, Pearson coefficients and residual-set coverage are fully supplied. Background subtraction varied by ligand/donor/family rather than shifting each programme by a single constant. No correlation P values were used because donors and controls were shared.

## Interpretation boundary

The most defensible conclusion is method- and context-qualified anti-disease-direction transcriptional change for selected combinations. Cross-scoring concordance does not imply absence of pool, composition or normalization dependence. Background-relative positivity alone is a deviation from a gene reference, not restoration. Neither these observations nor the unchanged GSEA establishes epithelial functional rescue, ICI-specificity, causal mechanism, pharmacologic safety or a reason to advance systemic TNFSF12/Fn14 targeting.
'''
write('results_update.md',body)
write('conclusions_update.md',f'''# Proposed Conclusions update

The all-ligand atlas showed a widespread downward rank-score background that was also present in expression/detection-stratified reference gene sets. Retaining raw, matched-background-relative and normalized-expression readouts separated actual programme increases from relative preservation against that background. {len(stable)} of the {len(orig)} original donor-consistent Both combinations remained Both across all three pooled scoring methods; TGFB2-Stem and TNFSF12-TA did not. Additional combinations appeared under the alternative readouts and are reported without replacing the original classification.

The post hoc TGF-beta superfamily concentration remained descriptive and was qualified by pool dependence, incomplete substate support and differentiation-associated remodelling. BMP4-Colonocyte retained both directions after pool deletion and fixed mature/immature-state weighting, but its disease-depleted change was substantially attenuated by weighting. Other combinations, including INS-Colonocyte, depended on the expression aggregation scale or pool. IFNG/TNFA supported an internal disease-direction comparison in Stem, with TNFA's disease-depleted direction differing after background referencing. Reference-programme relationships persisted beyond direct member overlap but remained method- and family-dependent.

These results support a graded, context-dependent interpretation of apparent transcriptional-restoration signals. They do not demonstrate functional epithelial protection, ICI specificity, causal independence from differentiation, or therapeutic safety. Coarse gene-background matching and incomplete substate identifiability remain limitations; none of the three score representations is a uniquely corrected truth.
''')
write('figure_legends_update.md','''# Proposed figure legends

## Figure 4. Directional background and cross-scoring stability in the all-ligand organoid atlas

(A) Numbers of eligible ligand-family combinations with 0, 1, 2 or 3 paired donors showing disease-up decrease or disease-depleted increase. (B) Joint counts from the same combinations, not a product of marginal probabilities. (C) Distributions of observed donor-pair rank changes and the median B of 1,000 expression/detection-quantile-matched nonmember sets. Boxes span the 25th-75th percentiles, orange lines indicate medians, and whiskers extend to observations within 1.5 interquartile ranges; outliers are omitted from this summary display but retained in source tables and full atlases. (D) Classification of every planned ligand-family comparison under raw rank, matched-background-relative rank and normalized-expression scoring. (E) Every original Both combination, with B/U/D/O/NE denoting Both, UP-only, depleted-only, other observed directions and non-estimable. (F) Both numerators and eligible-combination denominators for the post hoc structural TGF-beta ligand-superfamily grouping and other ligands. No independent-combination or ligand-level Fisher test is implied.

Both requires disease-up delta <0 and disease-depleted delta >0 in each of the same three donors, using unrounded values; zero satisfies neither direction. Raw delta is treatment minus matched cytomix control. Background-relative delta is raw delta minus B and can be positive when the raw programme declines less than its matched gene background. Nonrank delta is the unweighted fixed-member mean of log2(normalized treatment count +1) minus log2(normalized control count +1), with the original DESeq2 ratio normalization formula and no model refit. Background libraries preserve joint quantile strata but exhibit residual continuous covariate imbalance; they do not match coexpression and are not no-treatment-effect nulls or donor confidence intervals. No original GSEA q value is used as significance for these added scores. No_Cytomix is excluded from the planned ligand denominator and shown separately in the supporting figures.

## Complete three-method atlas

Twelve pages show four epithelial families under each of the three scoring methods. All 79 planned ligands and the No_Cytomix reference are retained; non-estimable combinations are grey with an NE class. Cells show three-donor means with actual signs for all six programmes, including the residual and IBD1/IBD2 membership scores. Each programme/method uses its own symmetric colour range equal to the global all-family absolute mean-effect 95th percentile, fixed across family pages. Asterisks mark display saturation; printed values are rounded but not clipped, and full-precision means/extrema are in the source tables. Asterisks do not denote significance. B/U/D/O/NE use the Figure 4 definitions. IFN direction is not coded as universally favourable or adverse.

## Donor examples and joint directions

The fixed example rule retains all original Both combinations, IFNG and TNFA in every family, TNFSF12 in every family and No_Cytomix, including NE. The first two pages show individual paired donor effects with actual signs in each method's own units. Subsequent pages show family-specific joint effects: x = -delta disease-up, y = +delta disease-depleted. The upper-right quadrant therefore denotes anti-disease direction for these two programmes. Colours identify families, shapes identify donors, and numbers identify ligand means using the printed key. Shared donor/control structure is retained; dots are not independent perturbation experiments. Background-relative positivity is not interchangeable with raw or normalized-expression increase. All newly Both combinations are included in the full atlas and source tables.

## Shared-control and pool background sensitivity

For each family and programme, pool-level raw rank deltas, matched-background medians B and background-relative deltas are shown across eligible observed ligand-donor pairs. Boxes/whiskers follow the Figure 4 convention. Original pools remain distinct. Identical controls are counted once in the reuse table, and No_Cytomix is excluded from the ligand distributions. Pool-level scores are a separate sensitivity because pooling counts before ranking is nonlinear. Pool-relative comparisons against other ligands and leave-one-pool-out results are separate source tables, with changed coverage and NE explicitly reported; neither substitutes for matched cytomix controls or identifies an additive pool correction.

## Composition and reference-programme context

Page 1 shows treatment-minus-control recovered mature-substate proportions within Colonocyte and Goblet for all specified context ligands. These are recovered-cell proportions, not absolute survival or cell yield. Page 2 compares unstandardized and fixed-control-weighted common-support expression changes on the separate cell LogNormalize scale, ln(1 +10000 counts/cell library), averaging all available fixed members. Substates require >=20 cells in both arms and >=80% retained family-cell coverage in each arm; unsupported comparisons are not imputed. Weights are the same in treatment and control. Stem/TA lack finer released final substates. Pages 3-4 show descriptive family-specific Spearman correlations across eligible ligand means before and after symmetric removal of shared members, under all three methods. Exact common-ligand lists, each donor, remaining member counts and Pearson coefficients are supplied. Low correlation does not establish independence or specificity; overlap removal also changes set size. No new inference P values are displayed.
''')
claims=[
 ('Directional asymmetry extends to matched gene backgrounds','raw_direction_summary.tsv;random_reference_summary.tsv;pool_background_summary.tsv','All disease-programme pooled and pool-level B medians negative; raw double decline counted within donor','Descriptive background only; no cause of asymmetry identified'),
 ('Original Both is scoring-dependent','classification_stability.tsv;original_eight_Both_tracking.tsv',f'{len(stable)}/{len(orig)} original Both retained under three pooled methods','No claim of all-sensitivity stability or significance'),
 ('New Both combinations are disclosed','additional_Both_combinations.tsv','Includes all new background/nonrank Both, not just selected examples','Original raw classification retained'),
 ('TGF-beta concentration remains descriptive','family_response_summary.tsv;ligand_family_membership.tsv;pool_family_design_summary.tsv','Structural group counts with combination and unique-ligand denominators','Post hoc; shared donors/controls; no Fisher inference or common receptor assumption'),
 ('Pool deletion changes some conclusions','leave_one_pool_out_classification.tsv;leave_one_pool_out_donor_effects.tsv','All 610 scenarios retained, including NE','Changing donor coverage; pooled ranks not averages of pool ranks'),
 ('BMP4-Colonocyte has composition-associated and residual within-substate changes','within_family_substate_proportions.tsv;composition_standardized_expression.tsv','Mature fraction increases; weighted depleted change attenuates but remains positive','Separate cell normalization scale; no causal-direct-effect or functional protection claim'),
 ('INS-Colonocyte is not stable across all expression summaries','composition_standardized_expression.tsv','Cell-level unstandardized UP mixed; weighted depleted becomes negative in two donors','Cannot attribute all differences solely to composition'),
 ('IFNG/TNFA disease-direction references require qualification','reference_strict_direction_counts.tsv;reference_perturbation_effects.tsv','IFNG-Stem three directions persist; TNFA-Stem depleted direction mixed after B subtraction','Existing inflammatory background; not independent validation'),
 ('Reference relationships persist beyond direct overlap','gene_set_overlap.tsv;reference_correlations_by_method.tsv','All fixed/measured overlap and six symmetric-removal comparisons supplied','Membership size also changes; low rho is not independence'),
 ('Gene-background matching is coarse and imperfect','matching_balance.tsv;matching_strata.tsv;matched_set_membership',f'Max main-programme absolute global-SD standardized mean difference {abs(imb.standardized_mean_difference):.3f}','No data-derived retuning or claim of exact covariate balance'),
]
save(pd.DataFrame(claims,columns=['claim','evidence_files','observed_evidence','interpretation_limit']),'text/claim_evidence_matrix.tsv')
limitations=[
 dict(module='E3',status='BLOCKED_MISSING_FINER_ANNOTATIONS',item='Within-Stem/TA substate sensitivity',reason='Released final_annotations have one state per family. Degenerate equal-weight summaries cannot exclude hidden substate shifts.'),
 dict(module='E3',status='PARTIAL_DATA_DERIVED_NE',item='Common-support substate standardization',reason='Some colonocyte and most goblet contexts fail >=20 cells per arm/substate or >=80% retained coverage; per-donor NE retained. BMP6-Goblet all three donors fail.'),
 dict(module='E3',status='BLOCKED_MISSING_VALIDATED_DEFINITION',item='Dedicated stemness signature',reason='No independently validated frozen stemness signature in the scoped approved assets. Published annotations and all 20 saved cNMF programmes used without inventing a marker list.'),
 dict(module='B/E4',status='COMPLETED_WITH_RESIDUAL_MATCHING_IMBALANCE',item='Continuous expression/detection balance',reason='Joint quantile membership is matched but within-bin continuous differences remain, particularly IBD2. Detailed balance supplied; no interpretation as a fully balanced correction.'),
 dict(module='D/E1',status='NOT_IDENTIFIABLE_FROM_DESIGN',item='Pool-independent intrinsic ligand-superfamily effect',reason='Ligands are tied to linked pool pairs; family spans all pools but the design is not fully crossed. No deconfounding model or Fisher test fitted.'),
 dict(module='D',status='PARTIAL_DATA_DERIVED_NE',item='Pool deletion scenarios',reason=f'{int((lp.status!="ESTIMATED").sum())} of {len(lp)} deletion scenarios lack three remaining matched donors; this includes originally NE combinations.'),
 dict(module='Inherited metadata',status='UNCHANGED_UNRESOLVED_SOURCE_DISCREPANCY',item='IL22 dose metadata',reason='Prior approved run identified source IL22 versus publisher dose-list IL11 mismatch. No dose guessed or condition relabelled; no new scientific result depends on dose normalization.'),
]
save(pd.DataFrame(limitations),'review/module_status_and_limitations.tsv')
write('methods_update.md',(ROOT/'analysis_amendment.md').read_text(encoding='utf8'))
write('matching_interpretation_note.md',f'''# Matching interpretation

All random sets have exact target size, no duplicated genes, no target members, and the required counts in each occupied joint quantile stratum. No sparse-bin merges were necessary in this dataset. This does not guarantee equality of continuous baseline means: the broad upper quantile contains a range of highly expressed/detected genes. The largest main-programme imbalance was {abs(imb.standardized_mean_difference):.3f} gene-universe SD ({SHORT[imb.module]}, {imb.celltype}, {imb.dimension}). Detailed target/reference means, quantiles and differences are retained. No additional match-quality threshold was selected after observing response classifications. Results remain the pre-specified coarse-quantile reference analysis and are explicitly qualified, especially for IBD2; they are not an exactly matched or uniquely corrected biological effect.
''')
write('READ_ME_FIRST_zh.md',f'''# 本轮定向增补：先读结论与限制

- 唯一复用基线：`{BASE}`。本轮输出：`{ROOT}`；UTC 目录日期为 10 月 16 日，对应本地 10 月 17 日。
- 原始/背景相对/非秩 Both 数分别为 {int((cl.raw_class=='BOTH').sum())}/{int((cl.bg_relative_class=='BOTH').sum())}/{int((cl.nonrank_class=='BOTH').sum())}；原始 {len(orig)} 个中的 {len(stable)} 个在三套池合并评分下保持 Both。
- 保持的组合：{combo(stable)}。
- 改为 depleted-only：{combo(changed)}。新出现的所有 Both 见 `additional_Both_combinations.tsv`。
- 不能把上述 {len(stable)} 个称为“所有敏感性都通过”：池删一不一致、INS 的细胞级表达方向不一致、BMP6-Goblet 组成标准化缺乏共同支持均已保留。
- 背景匹配有真实连续协变量不平衡（尤其 IBD2），不是已消除偏差；B 的广泛负值也不证明技术原因。
- BMP4-Colonocyte 的成熟细胞回收比例增加；固定亚状态权重后 depleted 改变量明显减弱但三供体仍为正。这支持分化/组成相关转录重塑的限定解释，不证明功能保护。
- IFNG-Stem 的疾病方向在三套评分中一致；TNFA-Stem 的 depleted 在背景相对评分中变为 2 正 1 负，不能省略。
- `classification_stability.tsv` 包含全部 {len(cl)} 个计划组合，原始 NES/q 及检验族、供体连续值、控制与池、组成状态、池删一状态。
- 所有主分析保持原合格规则；没有重跑 DESeq2/fgsea 或 FDR，没有修改原手稿，也没有公开推送。
- 未能完成/不可识别的具体项目见 `review/module_status_and_limitations.tsv`；其余模块已执行。这里交付的是带明示限制的定向增补，不是无条件的投稿就绪认证。
''')
print('TEXT_AND_CONTEXT_COMPLETE',flush=True)
