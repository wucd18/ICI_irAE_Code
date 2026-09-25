#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: extract_GSE313368_module_effects.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA",
         compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

source_file <- file.path(DIRS$data, "41467_2025_68247_MOESM14_ESM.xls")
# libxls cannot open the Unicode target path on this Windows/R build. Copying to
# R's ASCII temporary directory preserves the source bytes and is reproducible.
temp_file <- tempfile(fileext = ".xls")
stopifnot(file.copy(source_file, temp_file, overwrite = TRUE))
on.exit(unlink(temp_file), add = TRUE)

sheets <- excel_sheets(temp_file)
out <- vector("list", length(sheets))
audit <- vector("list", length(sheets))

for (i in seq_along(sheets)) {
  sheet_name <- sheets[[i]]
  x <- as.data.table(read_excel(temp_file, sheet = sheet_name))
  setnames(x, names(x), make.names(names(x), unique = TRUE))

  if (sheet_name %in% c("A_IBD1_IBD2", "B_cNMF_Gene_Programs")) {
    y <- x[, .(
      dataset = "GSE313368",
      source_sheet = sheet_name,
      analysis_scope = "injury_vs_control",
      module_family = fifelse(sheet_name == "A_IBD1_IBD2", "IBD1_IBD2", "cNMF"),
      perturbation = as.character(perturbation),
      control = as.character(control),
      module = as.character(module),
      celltype = as.character(cell),
      pvalue = suppressWarnings(as.numeric(p_value)),
      padj = suppressWarnings(as.numeric(padj)),
      mean_pert = NA_real_,
      mean_control = NA_real_,
      mean_difference = suppressWarnings(as.numeric(mean_diff)),
      standard_deviation = suppressWarnings(as.numeric(stdev)),
      effect_size = suppressWarnings(as.numeric(effect_size))
    )]
  } else {
    y <- x[, .(
      dataset = "GSE313368",
      source_sheet = sheet_name,
      analysis_scope = "perturbation_vs_cytomix",
      module_family = fifelse(sheet_name == "C_S12D_IBD1_IBD2", "IBD1_IBD2", "cNMF"),
      perturbation = as.character(perturbation),
      control = "cytomix_only",
      module = as.character(Module),
      celltype = as.character(celltype),
      pvalue = suppressWarnings(as.numeric(pval)),
      padj = suppressWarnings(as.numeric(padj)),
      mean_pert = suppressWarnings(as.numeric(mean_pert)),
      mean_control = suppressWarnings(as.numeric(mean_control)),
      mean_difference = suppressWarnings(as.numeric(mean_pert) - as.numeric(mean_control)),
      standard_deviation = NA_real_,
      effect_size = suppressWarnings(as.numeric(cohen))
    )]
  }
  out[[i]] <- y
  audit[[i]] <- data.table(
    dataset = "GSE313368",
    source_file = basename(source_file),
    source_sheet = sheet_name,
    data_level = "donor_aggregated_module_effect",
    celltype = "multiple_epithelial_celltypes",
    n_rows = nrow(y),
    n_perturbations = uniqueN(y$perturbation),
    sha256_of_tab_separated_data_rows = NA_character_,
    qc_issue = NA_character_,
    analysis_use = fifelse(
      unique(y$analysis_scope) == "perturbation_vs_cytomix",
      "real_perturbation_module_validation",
      "injury_model_context_only"
    )
  )
}

effects <- rbindlist(out, use.names = TRUE, fill = TRUE)
setorder(effects, analysis_scope, module_family, module, perturbation, celltype)
atomic_fwrite(effects, file.path(DIRS$results, "GSE313368_celltype_module_effects.tsv.gz"))

audit_path <- file.path(DIRS$tables, "Table_S20_GSE313368_perturbation_data_audit.tsv")
existing_audit <- fread(audit_path)
combined_audit <- rbindlist(c(list(existing_audit), audit), use.names = TRUE, fill = TRUE)
setorder(combined_audit, source_file, source_sheet)
atomic_fwrite(combined_audit, audit_path)

cat("Extracted GSE313368 module effects.\n")
print(effects[, .(
  n_rows = .N,
  n_perturbations = uniqueN(perturbation),
  n_celltypes = uniqueN(celltype),
  n_modules = uniqueN(module)
), by = .(analysis_scope, module_family)])
