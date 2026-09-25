#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: audit_organoid_offtarget_pathways.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

screen <- fread(file.path(DIRS$tables, "Table_S23_GSE313368_irAE_rescue_screen.tsv"))
candidates <- unique(c(screen[grepl("^Tier_[AB]", evidence_tier)]$perturbation, "NoCytomix"))
de <- fread(file.path(DIRS$results, "GSE313368_celltype_perturbation_DE_all.tsv.gz"))[
  perturbation %chin% candidates &
    celltype %chin% c("Transit_Amplifying", "Colonocyte") &
    supplement_qc_exclude_gene_level == FALSE &
    is.finite(log2FoldChange) & is.finite(lfcSE) & lfcSE > 0
]
de[, rank_z := log2FoldChange / lfcSE]

hallmark_names <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_IL2_STAT5_SIGNALING",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_HYPOXIA",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_E2F_TARGETS",
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_P53_PATHWAY"
)
msig <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
pathways <- split(msig[gs_name %chin% hallmark_names]$gene_symbol,
                  msig[gs_name %chin% hallmark_names]$gs_name)

groups <- unique(de[, .(perturbation, celltype)])
rows <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
  block <- de[perturbation == groups$perturbation[[i]] & celltype == groups$celltype[[i]]]
  block <- block[nzchar(gene)][order(-abs(rank_z))][!duplicated(gene)]
  ranks <- block$rank_z
  names(ranks) <- block$gene
  ranks <- sort(ranks, decreasing = TRUE)
  if (length(ranks) < 1000L) next
  ans <- as.data.table(fgseaMultilevel(
    pathways = pathways, stats = ranks, minSize = 15L, maxSize = 500L,
    eps = 1e-20, nproc = 1L, nPermSimple = 10000L
  ))
  ans[, `:=`(
    perturbation = groups$perturbation[[i]],
    celltype = groups$celltype[[i]],
    leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1))
  )]
  rows[[i]] <- ans
}
gsea <- rbindlist(rows, use.names = TRUE, fill = TRUE)
setcolorder(gsea, c(
  "perturbation", "celltype", "pathway", "NES", "pval", "padj",
  "size", "leadingEdge", "log2err", "ES"
))
gsea[, abs_NES := abs(NES)]
setorder(gsea, perturbation, celltype, padj, -abs_NES)
atomic_fwrite(gsea, file.path(DIRS$tables, "Table_S47_organoid_candidate_offtarget_Hallmark_GSEA.tsv"))

inflammatory <- c(
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING", "HALLMARK_IL2_STAT5_SIGNALING"
)
tumor_context <- c(
  "HALLMARK_HYPOXIA", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_ANGIOGENESIS", "HALLMARK_TGF_BETA_SIGNALING",
  "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_MYC_TARGETS_V1", "HALLMARK_MTORC1_SIGNALING"
)
summary <- gsea[, .(
  n_tested_pathway_celltypes = .N,
  n_positive_inflammatory_fdr05 = sum(pathway %chin% inflammatory & NES > 0 & padj <= 0.05),
  positive_inflammatory_pathways = paste(unique(pathway[
    pathway %chin% inflammatory & NES > 0 & padj <= 0.05
  ]), collapse = ";"),
  n_positive_tumor_context_fdr05 = sum(pathway %chin% tumor_context & NES > 0 & padj <= 0.05),
  positive_tumor_context_pathways = paste(unique(pathway[
    pathway %chin% tumor_context & NES > 0 & padj <= 0.05
  ]), collapse = ";"),
  n_negative_IFN_fdr05 = sum(
    pathway %chin% c("HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE") &
      NES < 0 & padj <= 0.05
  ),
  max_positive_inflammatory_NES = {
    x <- NES[pathway %chin% inflammatory & NES > 0 & padj <= 0.05]
    if (length(x)) max(x) else NA_real_
  },
  max_positive_tumor_context_NES = {
    x <- NES[pathway %chin% tumor_context & NES > 0 & padj <= 0.05]
    if (length(x)) max(x) else NA_real_
  }
), by = perturbation]
summary[screen, on = "perturbation", evidence_tier := i.evidence_tier]
summary[, interpretation := fcase(
  n_negative_IFN_fdr05 > 0,
  "deprioritize: suppresses one or more interferon Hallmark programs",
  n_positive_inflammatory_fdr05 > 0 | n_positive_tumor_context_fdr05 > 0,
  "mechanistic/local-only candidate: repair signal coexists with inflammatory or tumor-context programs",
  default = "no prespecified off-target Hallmark activation detected; absence is not evidence of safety"
)]
setorder(summary, perturbation)
atomic_fwrite(summary, file.path(DIRS$tables, "Table_S48_organoid_candidate_offtarget_summary.tsv"))

cat("Organoid off-target Hallmark audit complete.\n")
print(summary[perturbation == "TNFSF12"])
