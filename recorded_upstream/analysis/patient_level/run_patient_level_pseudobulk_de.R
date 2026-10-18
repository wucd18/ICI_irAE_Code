#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(edgeR)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_patient_level_pseudobulk_de.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA", compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
  invisible(path)
}

read_pseudobulk <- function(prefix) {
  meta <- fread(paste0(prefix, "_metadata.tsv"))
  genes <- fread(paste0(prefix, "_genes.tsv"))
  con <- gzfile(paste0(prefix, "_counts.mtx.gz"), open = "rt")
  on.exit(close(con), add = TRUE)
  counts <- readMM(con)
  stopifnot(nrow(counts) == nrow(meta), ncol(counts) == nrow(genes))
  counts <- t(counts)
  rownames(counts) <- make.unique(genes$gene)
  colnames(counts) <- meta$pseudobulk_id
  list(counts = counts, meta = meta, genes = genes)
}

fit_one_celltype <- function(counts, meta, celltype_name, dataset, organ, contrast_name,
                             case_label, control_label, covariates = character()) {
  keep_meta <- meta[meta$celltype == celltype_name &
                      meta$eligible_min_20_cells == TRUE &
                      meta$condition %chin% c(control_label, case_label)]
  n_case <- uniqueN(keep_meta[condition == case_label, patient_id])
  n_control <- uniqueN(keep_meta[condition == control_label, patient_id])
  if (n_case < MIN_SAMPLES_PER_GROUP || n_control < MIN_SAMPLES_PER_GROUP) return(NULL)

  keep_meta[, condition_factor := factor(condition, levels = c(control_label, case_label))]
  valid_covariates <- covariates[covariates %chin% names(keep_meta)]
  valid_covariates <- valid_covariates[vapply(valid_covariates, function(x) uniqueN(keep_meta[[x]]) > 1L, logical(1))]
  formula_text <- paste("~", paste(c(valid_covariates, "condition_factor"), collapse = " + "))
  design <- model.matrix(as.formula(formula_text), data = keep_meta)
  if (qr(design)$rank < ncol(design)) {
    valid_covariates <- character()
    design <- model.matrix(~ condition_factor, data = keep_meta)
  }

  y <- DGEList(counts = counts[, keep_meta$pseudobulk_id, drop = FALSE])
  keep_gene <- filterByExpr(y, design = design, min.count = 10, min.total.count = 15)
  if (sum(keep_gene) < 200L) return(NULL)
  y <- y[keep_gene, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design = design, robust = TRUE)
  fit <- glmQLFit(y, design = design, robust = TRUE)
  coef_name <- grep("^condition_factor", colnames(design), value = TRUE)
  if (length(coef_name) != 1L) stop("Could not identify disease coefficient for ", dataset, "/", celltype_name)
  test <- glmQLFTest(fit, coef = coef_name)
  ans <- as.data.table(topTags(test, n = Inf, sort.by = "none")$table, keep.rownames = "gene")
  ans[, `:=`(
    dataset = dataset,
    organ = organ,
    celltype = celltype_name,
    contrast = contrast_name,
    case = case_label,
    control = control_label,
    n_case = n_case,
    n_control = n_control,
    design = paste(colnames(design), collapse = ";"),
    passes_effect_gate = FDR <= FDR_THRESHOLD & abs(logFC) >= MIN_ABS_LOG2FC
  )]
  setcolorder(ans, c("dataset", "organ", "celltype", "contrast", "case", "control",
                     "n_case", "n_control", "design", "gene", "logFC", "logCPM",
                     "F", "PValue", "FDR", "passes_effect_gate"))
  ans
}

heart <- read_pseudobulk(file.path(DIRS$results, "GSE228597_pseudobulk"))
heart_meta <- copy(heart$meta)
# Primary heart inference: all controls plus only pre-steroid myocarditis samples.
heart_meta <- heart_meta[
  condition == "control" |
    (condition == "myocarditis" & on_steroids == "False")
]
heart_meta[, institution := factor(institution)]

colon <- read_pseudobulk(file.path(DIRS$results, "GSE253720_pseudobulk"))

results <- list()
idx <- 0L
for (ct in sort(unique(heart_meta$celltype))) {
  x <- fit_one_celltype(
    heart$counts, heart_meta, ct, "GSE228597", "heart",
    "myocarditis_vs_control_pre_steroid", "myocarditis", "control",
    covariates = "institution"
  )
  if (!is.null(x)) { idx <- idx + 1L; results[[idx]] <- x }
}
for (ct in sort(unique(colon$meta$celltype))) {
  x <- fit_one_celltype(
    colon$counts, colon$meta, ct, "GSE253720", "colon",
    "CPI_colitis_vs_HC", "CPI_colitis", "HC"
  )
  if (!is.null(x)) { idx <- idx + 1L; results[[idx]] <- x }
}

if (!length(results)) stop("No cell-type contrast passed the patient/cell coverage gates")
de <- rbindlist(results, use.names = TRUE, fill = TRUE)
atomic_fwrite(de, file.path(DIRS$results, "patient_level_pseudobulk_DE_all.tsv.gz"))

summary <- de[, .(
  n_genes_tested = .N,
  n_up_fdr_effect = sum(passes_effect_gate & logFC > 0),
  n_down_fdr_effect = sum(passes_effect_gate & logFC < 0),
  median_abs_logFC = median(abs(logFC)),
  min_FDR = min(FDR)
), by = .(dataset, organ, celltype, contrast, n_case, n_control, design)]
atomic_fwrite(summary, file.path(DIRS$tables, "Table_S5_patient_level_DE_summary.tsv"))

top <- de[passes_effect_gate == TRUE]
top[, abs_logFC_order := abs(logFC)]
setorder(top, dataset, celltype, FDR, -abs_logFC_order)
top <- top[, head(.SD, 25L), by = .(dataset, organ, celltype, contrast)]
top[, abs_logFC_order := NULL]
atomic_fwrite(top, file.path(DIRS$tables, "Table_S6_top_patient_level_DE_genes.tsv"))

cat("Completed patient-level edgeR QL analyses.\n")
print(summary)
