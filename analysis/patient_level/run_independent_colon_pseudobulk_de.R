#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_independent_colon_pseudobulk_de.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA",
         compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

read_pseudobulk <- function(prefix) {
  meta <- fread(paste0(prefix, "_metadata.tsv"))
  genes <- fread(paste0(prefix, "_genes.tsv"))
  con <- gzfile(paste0(prefix, "_counts.mtx.gz"), "rt")
  on.exit(close(con), add = TRUE)
  counts <- t(readMM(con))
  stopifnot(nrow(counts) == nrow(genes), ncol(counts) == nrow(meta))
  rownames(counts) <- make.unique(genes$gene)
  colnames(counts) <- meta$pseudobulk_id
  list(counts = counts, meta = meta)
}

# Mapping is taken from the source study's public analysis/R/14-figures.R.
cluster_family <- c(
  "8" = "Immature epithelial", "3" = "Immature epithelial", "14" = "Immature epithelial",
  "1" = "Absorptive epithelial", "5" = "Absorptive epithelial", "6" = "Absorptive epithelial",
  "16" = "Absorptive epithelial", "18" = "Absorptive epithelial", "9" = "Absorptive epithelial",
  "12" = "Absorptive epithelial", "2" = "Mature absorptive epithelial",
  "11" = "Mature absorptive epithelial", "20" = "Mature absorptive epithelial",
  "4" = "Secretory epithelial", "10" = "Secretory epithelial", "17" = "Secretory epithelial",
  "15" = "Secretory epithelial", "19" = "Secretory epithelial",
  "7" = "Mesenchymal", "23" = "Mesenchymal", "21" = "Mesenchymal",
  "22" = "Mesenchymal", "13" = "Mesenchymal"
)

aggregate_to_family <- function(counts, meta) {
  meta <- copy(meta)
  meta[, family := unname(cluster_family[as.character(celltype)])]
  if (anyNA(meta$family)) stop("Unmapped GSE206300 clusters: ", paste(unique(meta[is.na(family), celltype]), collapse = ","))
  meta[, group_key := paste(patient_id, condition, family, sep = "||")]
  keys <- unique(meta$group_key)
  membership <- sparseMatrix(
    i = match(meta$group_key, keys), j = seq_len(nrow(meta)), x = 1,
    dims = c(length(keys), nrow(meta))
  )
  family_counts <- counts %*% t(membership)
  family_meta <- meta[, .(
    patient_id = first(patient_id), condition = first(condition), celltype = first(family),
    n_cells = sum(n_cells), drug = paste(sort(unique(drug)), collapse = ";"),
    chemistry = paste(sort(unique(chemistry)), collapse = ";"),
    class2 = paste(sort(unique(class2)), collapse = ";"),
    class_short = paste(sort(unique(class_short)), collapse = ";")
  ), by = group_key]
  family_meta[, pseudobulk_id := paste0("GSE206300_FAMILY_", sprintf("%04d", .I))]
  family_meta[, eligible_min_20_cells := n_cells >= MIN_CELLS_PER_SAMPLE_CELLTYPE]
  colnames(family_counts) <- family_meta$pseudobulk_id
  list(counts = family_counts, meta = family_meta)
}

fit_one <- function(counts, meta, celltype_name, feature_level, contrast_name,
                    pd1_only = FALSE) {
  keep <- meta[
    celltype == celltype_name & eligible_min_20_cells == TRUE &
      condition %chin% c("On ICI therapy", "irColitis")
  ]
  if (pd1_only) keep <- keep[drug == "PD-1"]
  n_case <- uniqueN(keep[condition == "irColitis", patient_id])
  n_control <- uniqueN(keep[condition == "On ICI therapy", patient_id])
  if (n_case < MIN_SAMPLES_PER_GROUP || n_control < MIN_SAMPLES_PER_GROUP) return(NULL)
  keep[, condition_factor := factor(condition, levels = c("On ICI therapy", "irColitis"))]
  design <- model.matrix(~ condition_factor, data = keep)
  y <- DGEList(counts = counts[, keep$pseudobulk_id, drop = FALSE])
  keep_gene <- filterByExpr(y, design = design, min.count = 10, min.total.count = 15)
  if (sum(keep_gene) < 200L) return(NULL)
  y <- y[keep_gene, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  test <- glmQLFTest(fit, coef = grep("^condition_factor", colnames(design)))
  ans <- as.data.table(topTags(test, n = Inf, sort.by = "none")$table, keep.rownames = "gene")
  ans[, `:=`(
    dataset = "GSE206300", organ = "colon", feature_level = feature_level,
    celltype = celltype_name, contrast = contrast_name,
    case = "irColitis", control = "On ICI therapy", n_case = n_case,
    n_control = n_control, regimen_restriction = if (pd1_only) "PD-1_only" else "all_regimens",
    design = paste(colnames(design), collapse = ";"),
    passes_effect_gate = FDR <= FDR_THRESHOLD & abs(logFC) >= MIN_ABS_LOG2FC
  )]
  setcolorder(ans, c(
    "dataset", "organ", "feature_level", "celltype", "contrast", "case", "control",
    "n_case", "n_control", "regimen_restriction", "design", "gene", "logFC", "logCPM",
    "F", "PValue", "FDR", "passes_effect_gate"
  ))
  ans
}

pb <- read_pseudobulk(file.path(DIRS$results, "GSE206300_pseudobulk"))
pb$meta[, celltype := as.character(celltype)]
family <- aggregate_to_family(pb$counts, pb$meta)

results <- list()
for (level in c("cluster", "family")) {
  obj <- if (level == "cluster") pb else family
  for (ct in sort(unique(obj$meta$celltype))) {
    primary <- fit_one(obj$counts, obj$meta, ct, level,
                       "irColitis_vs_ICI_control_PD1_only", pd1_only = TRUE)
    sensitivity <- fit_one(obj$counts, obj$meta, ct, level,
                           "irColitis_vs_ICI_control_all_regimens", pd1_only = FALSE)
    if (!is.null(primary)) results[[length(results) + 1L]] <- primary
    if (!is.null(sensitivity)) results[[length(results) + 1L]] <- sensitivity
  }
}
if (!length(results)) stop("No independent colon contrasts passed coverage gates")
de <- rbindlist(results, use.names = TRUE, fill = TRUE)
atomic_fwrite(de, file.path(DIRS$results, "independent_colon_patient_level_DE_all.tsv.gz"))

summary <- de[, .(
  n_genes_tested = .N,
  n_up_fdr_effect = sum(passes_effect_gate & logFC > 0),
  n_down_fdr_effect = sum(passes_effect_gate & logFC < 0),
  median_abs_logFC = median(abs(logFC)), min_FDR = min(FDR)
), by = .(dataset, organ, feature_level, celltype, contrast, n_case, n_control,
          regimen_restriction, design)]
atomic_fwrite(summary, file.path(DIRS$tables, "Table_S11_independent_colon_DE_summary.tsv"))

top <- de[passes_effect_gate == TRUE]
top[, abs_logFC_order := abs(logFC)]
setorder(top, contrast, feature_level, celltype, FDR, -abs_logFC_order)
top <- top[, head(.SD, 25L), by = .(contrast, feature_level, celltype)]
top[, abs_logFC_order := NULL]
atomic_fwrite(top, file.path(DIRS$tables, "Table_S12_independent_colon_top_DE_genes.tsv"))

mapping <- data.table(cluster = names(cluster_family), cell_family = unname(cluster_family))
mapping <- mapping[order(as.integer(cluster))]
atomic_fwrite(mapping, file.path(DIRS$tables, "Table_S13_GSE206300_cluster_family_mapping.tsv"))
cat("Completed independent GSE206300 patient-level edgeR analyses.\n")
print(summary)
