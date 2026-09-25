#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: validate_TNFSF12_signature_reciprocity.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

feature <- fread(file.path(DIRS$tables, "Table_S25_ICB_guardrail_feature_definitions.tsv"))[
  feature_id == "PERTSIG_TNFSF12_TA"
]
stopifnot(nrow(feature) == 1L)
rescue_feature <- fread(file.path(DIRS$tables, "Table_S25_ICB_guardrail_feature_definitions.tsv"))[
  feature_id == "TNFSF12_ICI_LOSS_RESCUE_CORE"
]
stopifnot(nrow(rescue_feature) == 1L)
pathways <- list(
  TNFSF12_ORGANOID_UP = strsplit(feature$up_genes, ";", fixed = TRUE)[[1]],
  TNFSF12_ORGANOID_DOWN = strsplit(feature$down_genes, ";", fixed = TRUE)[[1]],
  TNFSF12_FUNCTION_RESCUE_CORE = strsplit(rescue_feature$up_genes, ";", fixed = TRUE)[[1]]
)

make_rank <- function(block, gene_col = "gene") {
  block <- copy(block[is.finite(F) & is.finite(logFC)])
  block[, symbol := sub("^.*\\|", "", get(gene_col))]
  block <- block[nzchar(symbol)][order(-abs(logFC))][!duplicated(symbol)]
  stats <- sign(block$logFC) * sqrt(pmax(block$F, 0))
  names(stats) <- block$symbol
  sort(stats, decreasing = TRUE)
}

run_one <- function(stats, dataset, feature_level, celltype, contrast, role) {
  ans <- as.data.table(fgseaMultilevel(
    pathways = pathways, stats = stats, minSize = 10L, maxSize = 200L,
    eps = 1e-20, nproc = 1L, nPermSimple = 10000L
  ))
  ans[, `:=`(
    dataset = dataset, feature_level = feature_level, celltype = celltype,
    contrast = contrast, cohort_role = role,
    leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1))
  )]
  ans
}

rows <- list()
discovery <- fread(file.path(DIRS$results, "patient_level_pseudobulk_DE_all.tsv.gz"))[
  dataset == "GSE253720" & celltype == "Epithelial" &
    contrast == "CPI_colitis_vs_HC"
]
rows[[length(rows) + 1L]] <- run_one(
  make_rank(discovery), "GSE253720", "broad_celltype", "Epithelial",
  "CPI_colitis_vs_HC", "derivation_cohort_reciprocity_check"
)

independent <- fread(file.path(DIRS$results, "independent_colon_patient_level_DE_all.tsv.gz"))[
  feature_level == "family" & contrast == "irColitis_vs_ICI_control_PD1_only" &
    celltype %chin% c(
      "Absorptive epithelial", "Immature epithelial",
      "Mature absorptive epithelial", "Secretory epithelial"
    )
]
for (ct in sort(unique(independent$celltype))) {
  block <- independent[celltype == ct]
  rows[[length(rows) + 1L]] <- run_one(
    make_rank(block), "GSE206300", "epithelial_family", ct,
    "irColitis_vs_ICI_control_PD1_only", "independent_primary"
  )
}

spatial <- fread(file.path(DIRS$results, "spatial_colon_patient_level_DE.tsv.gz"))
rows[[length(rows) + 1L]] <- run_one(
  make_rank(spatial), "GSE210037", "spatial_sample_pseudobulk", "whole_spot_sample",
  "irAE_colitis_vs_healthy_colon", "independent_spatial_platform"
)

result <- rbindlist(rows, use.names = TRUE, fill = TRUE)
result[, expected_direction_for_reciprocal_repair := fcase(
  pathway == "TNFSF12_ORGANOID_DOWN", "up_in_colitis",
  default = "down_in_colitis"
)]
result[, observed_direction_in_colitis := fifelse(NES > 0, "up_in_colitis", "down_in_colitis")]
result[, reciprocal_direction_concordant :=
  expected_direction_for_reciprocal_repair == observed_direction_in_colitis]
setcolorder(result, c(
  "dataset", "cohort_role", "feature_level", "celltype", "contrast",
  "pathway", "expected_direction_for_reciprocal_repair",
  "observed_direction_in_colitis", "reciprocal_direction_concordant",
  "NES", "pval", "padj", "size", "leadingEdge", "log2err", "ES"
))
setorder(result, dataset, celltype, pathway)
atomic_fwrite(
  result,
  file.path(DIRS$tables, "Table_S49_TNFSF12_organoid_signature_disease_reciprocity.tsv")
)

summary <- result[, .(
  n_tests = .N,
  n_reciprocal = sum(reciprocal_direction_concordant, na.rm = TRUE),
  n_reciprocal_fdr05 = sum(reciprocal_direction_concordant & padj <= 0.05, na.rm = TRUE),
  median_NES = median(NES, na.rm = TRUE)
), by = .(pathway, cohort_role)]
atomic_fwrite(
  summary,
  file.path(DIRS$tables, "Table_S50_TNFSF12_organoid_signature_reciprocity_summary.tsv")
)
cat("TNFSF12 signature reciprocity analysis complete.\n")
print(summary)
