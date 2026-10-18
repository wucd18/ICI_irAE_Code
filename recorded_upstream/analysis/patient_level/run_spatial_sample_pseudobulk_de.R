#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_spatial_sample_pseudobulk_de.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA",
         compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

meta <- fread(file.path(DIRS$results, "GSE210037_spatial_pseudobulk_metadata.tsv"))
genes <- fread(file.path(DIRS$results, "GSE210037_spatial_pseudobulk_genes.tsv"))
con <- gzfile(file.path(DIRS$results, "GSE210037_spatial_pseudobulk_counts.mtx.gz"), "rt")
counts <- t(readMM(con))
close(con)
stopifnot(nrow(counts) == nrow(genes), ncol(counts) == nrow(meta))
rownames(counts) <- make.unique(genes$gene)
colnames(counts) <- meta$sample

# The dermatitis comparison has only two biological samples per group and is
# deliberately excluded from inferential DE. It remains available for module
# localization after modules are frozen in independent cohorts.
keep <- meta[tissue %chin% c("healthy_colon", "irAE_colitis")]
keep[, condition := factor(tissue, levels = c("healthy_colon", "irAE_colitis"))]
stopifnot(uniqueN(keep[tissue == "healthy_colon", sample]) == 8L)
stopifnot(uniqueN(keep[tissue == "irAE_colitis", sample]) == 7L)
design <- model.matrix(~ condition, data = keep)

y <- DGEList(counts = counts[, keep$sample, drop = FALSE])
keep_gene <- filterByExpr(y, design = design, min.count = 10, min.total.count = 15)
y <- y[keep_gene, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y, method = "TMM")
y <- estimateDisp(y, design, robust = TRUE)
fit <- glmQLFit(y, design, robust = TRUE)
test <- glmQLFTest(fit, coef = grep("^condition", colnames(design)))
de <- as.data.table(topTags(test, n = Inf, sort.by = "none")$table, keep.rownames = "gene")
de[, `:=`(
  dataset = "GSE210037", organ = "colon_spatial", feature_level = "whole_spot_sample",
  contrast = "irAE_colitis_vs_healthy_colon", n_case = 7L, n_control = 8L,
  design = paste(colnames(design), collapse = ";"),
  passes_effect_gate = FDR <= FDR_THRESHOLD & abs(logFC) >= MIN_ABS_LOG2FC
)]
setcolorder(de, c("dataset", "organ", "feature_level", "contrast", "n_case", "n_control",
                  "design", "gene", "logFC", "logCPM", "F", "PValue", "FDR",
                  "passes_effect_gate"))
atomic_fwrite(de, file.path(DIRS$results, "spatial_colon_patient_level_DE.tsv.gz"))

summary <- de[, .(
  n_genes_tested = .N, n_up_fdr_effect = sum(passes_effect_gate & logFC > 0),
  n_down_fdr_effect = sum(passes_effect_gate & logFC < 0),
  median_abs_logFC = median(abs(logFC)), min_FDR = min(FDR)
), by = .(dataset, organ, feature_level, contrast, n_case, n_control, design)]
atomic_fwrite(summary, file.path(DIRS$tables, "Table_S16_spatial_colon_DE_summary.tsv"))

top <- de[passes_effect_gate == TRUE]
top[, abs_logFC_order := abs(logFC)]
setorder(top, FDR, -abs_logFC_order)
top <- head(top, 100L)
top[, abs_logFC_order := NULL]
atomic_fwrite(top, file.path(DIRS$tables, "Table_S17_spatial_colon_top_DE_genes.tsv"))
cat("Completed sample-level spatial colon edgeR analysis.\n")
print(summary)
