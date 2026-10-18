#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: derive_and_validate_irAE_modules.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA",
         compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

core_de <- fread(file.path(DIRS$results, "patient_level_pseudobulk_DE_all.tsv.gz"))
discovery_epi <- core_de[
  dataset == "GSE253720" & celltype == "Epithelial" & contrast == "CPI_colitis_vs_HC"
]
stopifnot(nrow(discovery_epi) > 1000L)

strict_up <- discovery_epi[passes_effect_gate == TRUE & logFC > 0, unique(gene)]
strict_loss <- discovery_epi[passes_effect_gate == TRUE & logFC < 0, unique(gene)]

hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
hallmark_sets <- split(hallmark$gene_symbol, hallmark$gs_name)
get_hallmark_union <- function(names) unique(unlist(hallmark_sets[names], use.names = FALSE))
ifn_visibility_genes <- get_hallmark_union(c(
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"
))

modules <- list(
  ICI_COLITIS_EPITHELIAL_UP_STRICT = strict_up,
  ICI_COLITIS_EPITHELIAL_UP_NON_IFN = setdiff(strict_up, ifn_visibility_genes),
  ICI_COLITIS_EPITHELIAL_LOSS_STRICT = strict_loss,
  CURATED_IFN_VISIBILITY = ifn_visibility_genes,
  CURATED_COLON_METABOLIC_FUNCTION = get_hallmark_union(c(
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
    "HALLMARK_FATTY_ACID_METABOLISM",
    "HALLMARK_BILE_ACID_METABOLISM",
    "HALLMARK_XENOBIOTIC_METABOLISM",
    "HALLMARK_PEROXISOME"
  )),
  CURATED_EPITHELIAL_BARRIER = get_hallmark_union(c(
    "HALLMARK_APICAL_JUNCTION",
    "HALLMARK_PROTEIN_SECRETION"
  ))
)
modules <- lapply(modules, function(x) sort(unique(x[nzchar(x)])))

definition_rules <- data.table(
  module = names(modules),
  derivation_source = c(
    "GSE253720 patient-level epithelial pseudobulk",
    "GSE253720 patient-level epithelial pseudobulk plus MSigDB Hallmark exclusion",
    "GSE253720 patient-level epithelial pseudobulk",
    "MSigDB Hallmark predefined gene sets",
    "MSigDB Hallmark predefined gene sets",
    "MSigDB Hallmark predefined gene sets"
  ),
  selection_rule = c(
    "FDR<=0.05, abs(log2FC)>=0.5, log2FC>0 in CPI colitis vs healthy control",
    "strict epithelial up module after removing Hallmark interferon alpha/gamma genes",
    "FDR<=0.05, abs(log2FC)>=0.5, log2FC<0 in CPI colitis vs healthy control",
    "union of Hallmark interferon alpha and gamma response",
    "union of Hallmark oxidative phosphorylation, fatty acid, bile acid, xenobiotic metabolism and peroxisome",
    "union of Hallmark apical junction and protein secretion"
  ),
  expected_direction_in_irAE = c("up", "up", "down", "up", "down", "down"),
  analysis_role = c(
    "data_derived_disease_up_module",
    "data_derived_non_IFN_disease_associated_residual",
    "data_derived_function_loss_module",
    "shared_attack_axis_context_not_primary_novelty",
    "predefined_function_axis",
    "predefined_function_axis"
  )
)

definitions <- rbindlist(lapply(names(modules), function(nm) {
  data.table(module = nm, gene = modules[[nm]])
}))
definitions <- merge(definitions, definition_rules, by = "module", all.x = TRUE, sort = FALSE)
setcolorder(definitions, c(
  "module", "gene", "derivation_source", "selection_rule",
  "expected_direction_in_irAE", "analysis_role"
))
atomic_fwrite(definitions, file.path(DIRS$tables, "Table_S21_frozen_module_definitions.tsv"))

make_stats <- function(block, gene_col = "gene") {
  block <- block[is.finite(F) & is.finite(logFC)]
  stats <- sign(block$logFC) * sqrt(pmax(block$F, 0))
  names(stats) <- block[[gene_col]]
  stats <- stats[nzchar(names(stats)) & !duplicated(names(stats))]
  sort(stats, decreasing = TRUE)
}

validation_blocks <- list()
validation_meta <- list()
add_block <- function(key, block, cohort_role, celltype, contrast, n_case, n_control, gene_col = "gene") {
  validation_blocks[[key]] <<- make_stats(block, gene_col)
  validation_meta[[key]] <<- data.table(
    validation_dataset = unique(block$dataset)[1],
    cohort_role = cohort_role,
    feature_level = if ("feature_level" %in% names(block)) unique(block$feature_level)[1] else "celltype",
    celltype = celltype,
    contrast = contrast,
    n_case = n_case,
    n_control = n_control
  )
}

add_block(
  "GSE253720_Epithelial", discovery_epi, "module_derivation_self_check", "Epithelial",
  "CPI_colitis_vs_HC", unique(discovery_epi$n_case), unique(discovery_epi$n_control)
)

independent_de <- fread(file.path(DIRS$results, "independent_colon_patient_level_DE_all.tsv.gz"))
independent_de[, gene_symbol := sub("^.*\\|", "", gene)]
for (contrast_name in c(
  "irColitis_vs_ICI_control_PD1_only",
  "irColitis_vs_ICI_control_all_regimens"
)) {
  for (ct in c(
    "Absorptive epithelial", "Immature epithelial",
    "Mature absorptive epithelial", "Secretory epithelial"
  )) {
    block <- independent_de[
      feature_level == "family" & contrast == contrast_name & celltype == ct
    ]
    if (!nrow(block)) next
    key <- paste("GSE206300", contrast_name, ct, sep = "::")
    add_block(
      key, block,
      ifelse(grepl("PD1_only", contrast_name), "independent_primary", "independent_sensitivity"),
      ct, contrast_name, unique(block$n_case), unique(block$n_control), "gene_symbol"
    )
  }
}

spatial_de <- fread(file.path(DIRS$results, "spatial_colon_patient_level_DE.tsv.gz"))
add_block(
  "GSE210037_spatial", spatial_de, "independent_spatial_validation", "whole_spot_sample",
  unique(spatial_de$contrast), unique(spatial_de$n_case), unique(spatial_de$n_control)
)

validation <- rbindlist(lapply(names(validation_blocks), function(key) {
  stats <- validation_blocks[[key]]
  fg <- as.data.table(fgseaMultilevel(
    pathways = modules, stats = stats, minSize = 10, maxSize = 500,
    eps = 0, scoreType = "std", nproc = 1
  ))
  fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1))]
  cbind(validation_meta[[key]], fg)
}), use.names = TRUE, fill = TRUE)
validation <- merge(
  validation, definition_rules[, .(module, expected_direction_in_irAE)],
  by.x = "pathway", by.y = "module", all.x = TRUE, sort = FALSE
)
validation[, observed_direction := fifelse(NES > 0, "up", "down")]
validation[, direction_concordant := observed_direction == expected_direction_in_irAE]
validation[, passes_fdr := padj <= FDR_THRESHOLD]
setnames(validation, "pathway", "module")
setcolorder(validation, c(
  "validation_dataset", "cohort_role", "feature_level", "celltype", "contrast",
  "n_case", "n_control", "module", "expected_direction_in_irAE",
  "observed_direction", "direction_concordant", "NES", "pval", "padj",
  "passes_fdr", "size", "leadingEdge"
))
setorder(validation, module, cohort_role, validation_dataset, celltype)
atomic_fwrite(validation, file.path(DIRS$tables, "Table_S22_frozen_module_validation.tsv"))
atomic_fwrite(validation, file.path(DIRS$results, "frozen_module_validation_all.tsv.gz"))

cat("Frozen modules and independent validation complete.\n")
print(definitions[, .(n_genes = .N), by = .(module, expected_direction_in_irAE)])
print(validation[, .(
  n_tests = .N,
  n_concordant = sum(direction_concordant),
  n_concordant_fdr = sum(direction_concordant & passes_fdr),
  median_NES = median(NES)
), by = .(module, cohort_role)])
