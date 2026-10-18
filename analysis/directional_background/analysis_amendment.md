# Directional/background sensitivity amendment

Locked before computing the new sensitivity results, 2026-10-17 Asia/Shanghai (run ID uses UTC).

This is a post hoc, targeted extension after inspection of the all-ligand atlas, the eight original Both combinations and the apparent TGF-beta concentration. It is not study-wide preregistration. The user's explicit reuse authorization supersedes the general clean-run prohibition on reading earlier results. No DESeq2 model, differential-expression test, fgsea or FDR is rerun. No original manuscript or result is overwritten.

## Inputs and unchanged primary definitions

Approved baseline: `E:/ICI_irAE_Codex_Workspace/runs/CSBJ_ALL_LIGANDS_20261016T152734Z_r3`. Reuse its inventory, saved paired scores, exact-symbol frozen/IBD membership, pool eligibility, pseudobulk counts and GSEA results. The saved sparse counts and released annotations are read only for gene-background detection, pool sensitivity and composition summaries. The existing rule is at least 20 cells in each arm within donor/family/pool, with all three matched donors for a classified combination. Missing comparisons remain NE. Strict signs are evaluated before rounding; zero is neither positive nor negative. No new inferential P values are generated. GSEA q values remain in their original 948/316/16 planned test families and do not test the new scores.

## A. Direction

Reuse saved donor deltas. Describe all six programmes by family/donor/matched control/pool coverage, including zeros, magnitudes, same-donor double decline and the 4 x 4 joint donor-count table. Keep assay-reference No_Cytomix separate from the planned ligand denominator.

## B. Matched background

One compatible library per epithelial family and actual member set: all samples use the same released gene universe, exact members, average-tie ranks divided by universe size. Baseline expression is log2(1 + CPM) from pooled unique matched cytomix cells within each donor/family; detection is the fraction of those cells with positive counts. Average the donor-specific expression and detection equally across donors. Deduplicate donor/pool/family controls before aggregation.

Use at most 10 empirical quantile bins for each baseline dimension, collapsing duplicate boundaries; assign identical values together. Start with occupied joint bins. For each member set exclude its members from candidate background. If a joint bin has fewer candidate genes than required target members, merge it with the nearest occupied bin/cluster by Euclidean distance between gene-count-weighted bin-coordinate centroids, tie-breaking by the lexicographically smallest constituent bin. Repeat until every cluster has enough candidates, or report NE if total capacity is insufficient. Record all merges, stratum target/candidate counts and continuous matching balance; merging is not presented as exact continuous matching.

Seed 20260826, deterministic SHA-256-derived stream by family and membership identity; 1,000 sets, sampling without replacement within a set. Store every gene index and universe order. Reuse each library for all donors, ligands and pools in that family. Compute B as the median of random-set paired deltas; background-relative delta = raw delta - B. Store 2.5/25/50/75/97.5 percentiles and empirical midrank percentile of the raw delta. These are gene-set reference distributions, not biological confidence intervals or P values. Do not generate random category-count nulls, independent sign flips or enrichment Fisher tests. Positive background-relative depleted change can mean less decline rather than an actual rise.

## C. Nonrank expression

The saved inventory contains size-factor minima/maxima but no per-sample factors. Recompute only the original DESeq2 default ratio size factors: filter genes by total count >=10 across the six original pseudobulks, take finite geometric means of genes positive in all samples, then exp(median(log(count)-log(geometric mean))) per sample. Use the same formula for all comparisons; verify against saved size-factor extrema. Do not fit a negative-binomial model. For every available fixed member compute log2(count_treatment/size_factor_treatment + 1) minus log2(count_control/size_factor_control + 1). Report unweighted member mean, median, quantiles, strict signs, and fraction with normalized count <1 in both arms. Retain low-count and zero members. This is normalized-expression change, not absolute RNA change.

## D. Pool sensitivity

Define unique controls by donor/family/original pool. Pool-level ranks are calculated separately using the same universe and matched libraries; they do not equal the mean of pooled-pseudobulk ranks. Reuse only eligible pool pairs. Report repeated ligand measurements and all leave-one-pool-out combinations: after deleting each pool present in a combination, require the same three donors and original cell threshold, otherwise NE. Recompute only pooled scores and normalization factors for remaining counts. Original and background-relative/nonrank results remain separate.

Optional pool-relative reference is raw delta minus median of other eligible ligand deltas in the same donor/family/pool, with at least three other ligands, excluding No_Cytomix from the reference. Save reference IDs. Never add this subtraction to the matched-background result. Show pool-family overlap and lack of full identifiability; fit no family-vs-pool deconfounding model.

## E. Biological context and overlap

Use the task-specified 15-member BMP/GDF/INHBA/TGFB structural ligand-superfamily grouping, retain all other tested ligands in the membership table, and do not equate grouping with shared receptors or SMAD signalling. Report combination and unique-ligand denominators. Show IFNG/TNFA in all eligible families, preserving source IDs.

Composition: use released final_annotations (mature/immature colonocyte and goblet labels; Stem and TA have no finer released final substate). Use all original pools with matched recovered cells for broad family proportions; within-family summaries use the original eligible pools. These are recovered proportions, not survival or cell yield. For each selected ligand (all 15 superfamily members, INS, TNFSF12, IFNG, TNFA and No_Cytomix), describe all 20 saved publisher perturbation-derived cNMF usages X1_DEG...X20_DEG, without re-fitting or selecting response-favourable programmes. Preserve published GO annotations and source mapping from Supplementary Data 9 / MOESM11 and pinned Figure 4 code. Maturation uses the published annotations. A dedicated validated stemness signature is not invented if absent.

Composition-standardized expression is a separate sensitivity in the publisher's cell LogNormalize scale: ln(1 + 10000 * gene_count / total_cell_counts), mean over the same fixed programme members, then mean over cells. This follows the saved publisher NormalizeData workflow, not the pseudobulk size-factor scale in C. Within each matched donor/family's eligible pooled cells, retain substates with >=20 recovered cells in both arms; require >=80% of each arm's family cells retained. Standardized treatment and control use identical weights, the matched control's proportions renormalized over supported substates. Save unstandardized full-family, unstandardized common-support and standardized common-support changes, weights and lost coverage. Stem/TA single-state comparisons cannot resolve finer composition. No missing-state imputation, fractional-count model, new clustering or causal-direct-effect claim.

Overlap: full frozen and measurable pairwise overlap for all six programmes. Main correlations compare each disease programme to IFN/IBD1/IBD2, within family across ligands, both donor means and each donor, with common IDs saved, no P values. Symmetrically remove shared genes for those six pairs, keeping >=10 measurable genes per residual set for this sensitivity only; main eligibility never changes. Re-score all three methods including separately matched libraries for residual member sets. Remaining size/coverage and NE are explicit. Changes may reflect both overlap removal and altered size.

## F. Outputs and displays

One row per planned ligand/family in classification_stability.tsv; No_Cytomix separately. Include all original Both, all newly Both, NE reasons, donor effects, controls, pool coverage, composition status and original GSEA fields. Main examples fixed now: all original eight Both, IFNG/TNFA all families, TNFSF12 all families and No_Cytomix. All other combinations appear in complete supplemental atlases and source tables.

Plots use actual deltas unless axes explicitly say -delta disease-up / +delta disease-depleted. IFN retains actual direction. No GSEA significance boxes on new-score plots. Any atlas saturation uses the global programme/method absolute-effect 95th percentile, same across pages, with saturated cells marked and extrema supplied. Joint plots label quadrants, family colours and donor shapes. No efficacy ratios are calculated. Deliver Results, legends, conclusions, claim-evidence table and provenance/check records; no automatic manuscript replacement or public release.
