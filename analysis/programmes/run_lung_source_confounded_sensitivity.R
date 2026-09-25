#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
  library(fgsea)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_lung_source_confounded_sensitivity.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA",
         compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

prefix <- file.path(DIRS$results, "GSE277136_pseudobulk")
meta <- fread(paste0(prefix, "_metadata.tsv"))
genes <- fread(paste0(prefix, "_genes.tsv"))
con <- gzfile(paste0(prefix, "_counts.mtx.gz"), "rt")
counts <- t(readMM(con))
close(con)
stopifnot(nrow(counts) == nrow(genes), ncol(counts) == nrow(meta))
rownames(counts) <- make.unique(genes$gene)
colnames(counts) <- meta$pseudobulk_id

eligible <- meta[eligible_min_20_cells == TRUE,
                 .N, by = .(celltype, condition)][, if (all(c("CIP", "control") %chin% condition)) {
                   .(n_CIP = N[condition == "CIP"], n_control = N[condition == "control"])
                 }, by = celltype][n_CIP >= 3L & n_control >= 3L]
stopifnot(nrow(eligible) >= 3L)

de_rows <- list()
for (ct in eligible$celltype) {
  m <- meta[celltype == ct & eligible_min_20_cells == TRUE & condition %chin% c("CIP", "control")]
  y <- DGEList(counts = counts[, m$pseudobulk_id, drop = FALSE])
  keep <- filterByExpr(y, group = factor(m$condition), min.count = 5L)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  condition <- relevel(factor(m$condition), ref = "control")
  design <- model.matrix(~ condition)
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "conditionCIP")
  tab <- as.data.table(topTags(qlf, n = Inf, sort.by = "none")$table, keep.rownames = "gene")
  tab[, `:=`(
    dataset = "GSE277136", organ = "lung_BALF", celltype = ct,
    contrast = "CIP_vs_external_healthy_BALF_source_confounded",
    case = "CIP", control = "external_healthy_BALF",
    n_case = m[condition == "CIP", .N], n_control = m[condition == "control", .N],
    design = "condition_only_source_perfectly_confounded",
    inference_role = "sensitivity_only_not_causal_or_receiver_tissue_evidence"
  )]
  setcolorder(tab, c(
    "dataset", "organ", "celltype", "contrast", "case", "control",
    "n_case", "n_control", "design", "inference_role", "gene",
    "logFC", "logCPM", "F", "PValue", "FDR"
  ))
  de_rows[[ct]] <- tab
}
de <- rbindlist(de_rows, use.names = TRUE, fill = TRUE)
atomic_fwrite(de, file.path(DIRS$results, "GSE277136_source_confounded_patient_level_DE.tsv.gz"))

module_def <- fread(file.path(DIRS$tables, "Table_S21_frozen_module_definitions.tsv"))
pathways <- split(module_def$gene, module_def$module)
gsea_rows <- list()
for (ct in eligible$celltype) {
  block <- de[celltype == ct & is.finite(F) & is.finite(logFC)]
  block <- block[nzchar(gene)][order(-abs(logFC))][!duplicated(gene)]
  ranks <- sign(block$logFC) * sqrt(pmax(block$F, 0))
  names(ranks) <- block$gene
  ranks <- sort(ranks, decreasing = TRUE)
  ans <- as.data.table(fgseaMultilevel(
    pathways = pathways, stats = ranks, minSize = 10L, maxSize = 1000L,
    eps = 1e-20, nproc = 1L, nPermSimple = 10000L
  ))
  ans[, `:=`(
    dataset = "GSE277136", organ = "lung_BALF", celltype = ct,
    contrast = "CIP_vs_external_healthy_BALF_source_confounded",
    inference_role = "sensitivity_only_source_confounded",
    leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1))
  )]
  gsea_rows[[ct]] <- ans
}
gsea <- rbindlist(gsea_rows, use.names = TRUE, fill = TRUE)
gsea[, expected_direction_in_irAE := module_def[
  match(pathway, module), expected_direction_in_irAE
]]
gsea[, observed_direction := fifelse(NES > 0, "up", "down")]
gsea[, direction_concordant := expected_direction_in_irAE == observed_direction]
setcolorder(gsea, c(
  "dataset", "organ", "celltype", "contrast", "inference_role", "pathway",
  "expected_direction_in_irAE", "observed_direction", "direction_concordant",
  "NES", "pval", "padj", "size", "leadingEdge", "log2err", "ES"
))
setorder(gsea, celltype, pathway)
atomic_fwrite(gsea, file.path(DIRS$tables, "Table_S57_lung_source_confounded_module_sensitivity.tsv"))

summary <- gsea[, .(
  n_eligible_celltypes = .N,
  n_concordant = sum(direction_concordant, na.rm = TRUE),
  n_concordant_fdr05 = sum(direction_concordant & padj <= 0.05, na.rm = TRUE),
  median_NES = median(NES, na.rm = TRUE)
), by = pathway]
summary[, limitation := paste(
  "all CIP samples came from the study cohort and all healthy BALF controls from an external atlas;",
  "source and disease are perfectly confounded; use only as a direction-of-effect sensitivity check"
)]
setorder(summary, pathway)
atomic_fwrite(summary, file.path(DIRS$tables, "Table_S58_lung_source_confounded_module_summary.tsv"))

cat("Lung source-confounded sensitivity analysis complete.\n")
print(summary)
