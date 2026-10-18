from common import *
import itertools
x=read(ROOT/'results/composition_marker_summary.tsv');context=read(ROOT/'inputs/published_GP_context.tsv').set_index('programme');rows=[]
for ligand,fam,gp in itertools.product(SELECTED,FAMS,['GP11','GP14']):
 b=x[(x.ligand==ligand)&(x.celltype==fam)&(x.programme==gp)];v=b.unstandardized_full_delta
 rows.append(dict(ligand=ligand,celltype=fam,programme=gp,published_annotation=context.loc[gp,'top_published_GO_terms'],n_observed_donors=b.donor.nunique(),n_positive=int((v>0).sum()),n_negative=int((v<0).sum()),n_zero=int((v==0).sum()),mean_observed_delta=v.mean(),donor_deltas=json.dumps(dict(zip(b.donor,v))),status='THREE_DONORS_DESCRIPTIVE' if len(b)==3 else 'NE_THREE_DONOR_SUMMARY_PARTIAL_VALUES_RETAINED',units='saved_publisher_cNMF_usage'))
save(pd.DataFrame(rows),'results/published_proliferation_context.tsv')
p=ROOT/'text/results_update.md';text=p.read_text(encoding='utf8');anchor='## TGF-beta superfamily and disease-direction references'
paragraph='''Published proliferation-related programmes were not uniformly concordant. The source GO annotations identify GP11 with DNA replication and GP14 with mitotic spindle/cell-cycle processes. Among the original Both combinations, GP11 decreased in all three donors for BMP4-Colonocyte, BMP6-Goblet, INHBA-TA, INS-Colonocyte, TGFB1-TA and TGFB2-TA, whereas GP14 decreased in all three only for BMP6-Goblet and increased in all three for INHBA-TA. Other GP14 directions were mixed. These divergent phase-related usages do not establish a single proliferation or differentiation mechanism. All selected ligands, including missing/partial contexts, appear in published_proliferation_context.tsv; all 20 publisher programmes remain in composition_marker_summary.tsv.\n\n'''
assert anchor in text;text=text.replace(anchor,paragraph+anchor);p.write_text(text,encoding='utf8')
manifest=read(ROOT/'review/input_manifest.tsv')
def role(p):
 if p.endswith('Codex_Directional_Background_Sensitivity.md'):return ('User task file','Authorized scope and prespecified context groups')
 if 'frozen_module' in p:return ('GSE253720 and curated reference definitions; exact derivation stored in member CSV','Frozen disease/IFN/residual members')
 if 'IBD_membership' in p:return ('Published Capeling et al. supplementary member lists; DOI 10.1038/s41467-025-68247-6','IBD1/IBD2 fixed membership benchmarks')
 if 'article_supplementary' in p or 'Capeling_publisher' in p:return ('Capeling et al.; DOI 10.1038/s41467-025-68247-6','Publisher gene-programme annotations and interpretation provenance')
 if p.endswith(('.py','.R','.Rmd')):return ('Approved project or pinned publisher code','Existing scoring, normalization, parsing or annotation mapping reference')
 if '/model_counts/' in p.replace('\\','/') or 'pseudobulk_counts' in p:return ('GSE313368','Approved paired counts or their sample matching metadata')
 if 'counts_' in p or 'column_pointers' in p:return ('GSE313368','Saved sparse count array, descriptor or cell-column index')
 if 'cell_metadata' in p:return ('GSE313368','Released cell identities, donor/pool mapping, final substates and saved cNMF usages')
 if 'GSEA' in p or 'model_inventory' in p:return ('Approved all-ligand run on GSE313368','Preserved statistical results or original normalization-factor extrema')
 if 'donor_delta' in p or 'pattern_summary' in p:return ('Approved all-ligand run on GSE313368','Raw paired effects and original classification baseline')
 return ('GSE313368 / approved all-ligand run','Approved eligibility, matching, gene order or job inventory')
for i,r in manifest.iterrows():manifest.loc[i,['source','role']]=role(r.path)
save(manifest,'review/input_manifest.tsv')
print('CONTEXT_SUMMARIES_ADDED; no new scientific analysis or member selection')
