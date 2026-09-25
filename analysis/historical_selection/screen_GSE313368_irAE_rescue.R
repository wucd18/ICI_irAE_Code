#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: screen_GSE313368_irAE_rescue.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA",
         compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

module_def <- fread(file.path(DIRS$tables, "Table_S21_frozen_module_definitions.tsv"))
selected_modules <- c(
  "ICI_COLITIS_EPITHELIAL_UP_STRICT",
  "ICI_COLITIS_EPITHELIAL_UP_NON_IFN",
  "ICI_COLITIS_EPITHELIAL_LOSS_STRICT",
  "CURATED_IFN_VISIBILITY",
  "CURATED_COLON_METABOLIC_FUNCTION",
  "CURATED_EPITHELIAL_BARRIER"
)
pathways <- split(module_def[module %chin% selected_modules]$gene,
                  module_def[module %chin% selected_modules]$module)

de <- fread(file.path(DIRS$results, "GSE313368_celltype_perturbation_DE_all.tsv.gz"))
de <- de[supplement_qc_exclude_gene_level == FALSE]
numeric_cols <- c("baseMean", "log2FoldChange", "lfcSE", "pvalue", "padj")
de[, (numeric_cols) := lapply(.SD, as.numeric), .SDcols = numeric_cols]

groups <- de[, .(n_genes = uniqueN(gene)), by = .(celltype, perturbation)]
eligible <- groups[n_genes >= 1000L]
gsea_out <- vector("list", nrow(eligible))
for (i in seq_len(nrow(eligible))) {
  g <- eligible[i]
  block <- de[celltype == g$celltype & perturbation == g$perturbation]
  block <- block[is.finite(log2FoldChange) & is.finite(lfcSE) & lfcSE > 0]
  block[, rank_stat := log2FoldChange / lfcSE]
  block[, abs_rank_stat := abs(rank_stat)]
  setorder(block, gene, -abs_rank_stat)
  block <- block[!duplicated(gene)]
  stats <- block$rank_stat
  names(stats) <- block$gene
  stats <- sort(stats, decreasing = TRUE)
  set.seed(20260826L + i)
  fg <- as.data.table(fgseaSimple(
    pathways = pathways, stats = stats, minSize = 10, maxSize = 500,
    nperm = 10000, scoreType = "std", nproc = 1
  ))
  fg[, `:=`(
    dataset = "GSE313368",
    celltype = g$celltype,
    perturbation = g$perturbation,
    n_ranked_genes = length(stats),
    leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1))
  )]
  gsea_out[[i]] <- fg
  if (i %% 10L == 0L || i == nrow(eligible)) {
    cat(sprintf("%d/%d perturbation-celltype GSEA groups\n", i, nrow(eligible)))
  }
}

gsea <- rbindlist(gsea_out, use.names = TRUE, fill = TRUE)
gsea[, screen_padj := p.adjust(pval, method = "BH"), by = .(celltype, pathway)]
setnames(gsea, "pathway", "module")
setcolorder(gsea, c(
  "dataset", "celltype", "perturbation", "n_ranked_genes", "module",
  "ES", "NES", "pval", "padj", "screen_padj", "size", "leadingEdge"
))
gsea[, abs_NES_order := abs(NES)]
setorder(gsea, module, celltype, screen_padj, -abs_NES_order)
gsea[, abs_NES_order := NULL]
atomic_fwrite(gsea, file.path(DIRS$results, "GSE313368_irAE_module_GSEA_all.tsv.gz"))

wide <- dcast(
  gsea,
  perturbation + celltype + n_ranked_genes ~ module,
  value.var = "NES"
)

summarize_metric <- function(dt, value_col, prefix) {
  value <- dt[[value_col]]
  result <- data.table(
    perturbation = dt$perturbation,
    celltype = dt$celltype,
    value = value
  )[, {
    finite_value <- value[is.finite(value)]
    if (!length(finite_value)) {
      list(n_celltypes = 0L, n_positive = 0L, n_negative = 0L,
           median_value = NA_real_, min_value = NA_real_, max_value = NA_real_)
    } else {
      list(
        n_celltypes = length(finite_value),
        n_positive = sum(finite_value > 0),
        n_negative = sum(finite_value < 0),
        median_value = median(finite_value),
        min_value = min(finite_value),
        max_value = max(finite_value)
      )
    }
  }, by = perturbation]
  setnames(result, c(
    "perturbation", paste0(prefix, c(
      "_n_celltypes", "_n_positive", "_n_negative", "_median",
      "_min", "_max"
    ))
  ))
  result
}

metric_cols <- intersect(selected_modules, names(wide))
gene_summaries <- lapply(metric_cols, function(module_name) {
  prefix <- switch(module_name,
    ICI_COLITIS_EPITHELIAL_UP_STRICT = "ici_up",
    ICI_COLITIS_EPITHELIAL_UP_NON_IFN = "ici_up_non_ifn",
    ICI_COLITIS_EPITHELIAL_LOSS_STRICT = "ici_loss",
    CURATED_IFN_VISIBILITY = "ifn_visibility",
    CURATED_COLON_METABOLIC_FUNCTION = "metabolic_function",
    CURATED_EPITHELIAL_BARRIER = "epithelial_barrier"
  )
  summarize_metric(wide, module_name, prefix)
})
gene_summary <- Reduce(function(x, y) merge(x, y, by = "perturbation", all = TRUE), gene_summaries)

effects <- fread(file.path(DIRS$results, "GSE313368_celltype_module_effects.tsv.gz"))
effects <- effects[analysis_scope == "perturbation_vs_cytomix"]
effect_targets <- effects[module %chin% c("IBD1", "IBD2", "GP20")]
effect_summary <- effect_targets[, .(
  module_n_celltypes = sum(is.finite(effect_size)),
  module_n_positive = sum(effect_size > 0, na.rm = TRUE),
  module_n_negative = sum(effect_size < 0, na.rm = TRUE),
  module_median_effect = median(effect_size, na.rm = TRUE),
  module_min_effect = min(effect_size, na.rm = TRUE),
  module_max_effect = max(effect_size, na.rm = TRUE),
  module_min_padj = min(padj, na.rm = TRUE)
), by = .(perturbation, module)]
effect_wide <- dcast(
  effect_summary, perturbation ~ module,
  value.var = c(
    "module_n_celltypes", "module_n_positive", "module_n_negative",
    "module_median_effect", "module_min_effect", "module_max_effect", "module_min_padj"
  )
)
setnames(effect_wide, names(effect_wide), gsub("module_", "", names(effect_wide), fixed = TRUE))

all_perturbations <- data.table(perturbation = sort(unique(effects$perturbation)))
screen <- Reduce(
  function(x, y) merge(x, y, by = "perturbation", all = TRUE),
  list(all_perturbations, gene_summary, effect_wide)
)

# Pre-specified transparent gates. IBD2 is the published metabolism/differentiation
# program; the ICI-loss module is independently validated in two irAE cohorts.
screen[, ici_loss_consistent_rescue :=
  fcoalesce(ici_loss_n_celltypes >= 1 & ici_loss_n_positive == ici_loss_n_celltypes & ici_loss_median > 0, FALSE)]
screen[, ibd2_multicelltype_rescue :=
  fcoalesce(n_celltypes_IBD2 == 4 & n_positive_IBD2 >= 2 & median_effect_IBD2 > 0, FALSE)]
screen[, function_rescue_gate := ici_loss_consistent_rescue & ibd2_multicelltype_rescue]

# Strong visibility erasure is defined on two orthogonal representations: the
# curated IFN gene set (NES < -1) and cNMF GP20 (median Cohen's d < -1).
screen[, strong_ifn_gene_set_suppression :=
  fcoalesce(is.finite(ifn_visibility_median) & ifn_visibility_median < -1, FALSE)]
screen[, strong_gp20_suppression :=
  fcoalesce(is.finite(median_effect_GP20) & median_effect_GP20 < -1, FALSE)]
screen[, visibility_preserved_gate :=
  !strong_ifn_gene_set_suppression & !strong_gp20_suppression]
screen[, strong_non_IFN_disease_program_aggravation :=
  fcoalesce(is.finite(ici_up_non_ifn_median) & ici_up_non_ifn_median > 1, FALSE)]
screen[, strong_IBD1_aggravation :=
  fcoalesce(is.finite(median_effect_IBD1) & median_effect_IBD1 > 1, FALSE)]
screen[, disease_program_not_aggravated_gate :=
  !strong_non_IFN_disease_program_aggravation & !strong_IBD1_aggravation]
screen[, full_decoupling_gate :=
  function_rescue_gate & visibility_preserved_gate & disease_program_not_aggravated_gate]
screen[, is_control := perturbation == "NoCytomix"]
screen[, evidence_tier := fifelse(
  is_control, "positive_control",
  fifelse(
    full_decoupling_gate,
    "Tier_A_function_rescue_visibility_preserved_no_disease_program_aggravation",
    fifelse(
      function_rescue_gate & visibility_preserved_gate,
      "Tier_B_function_rescue_but_disease_program_aggravated",
      fifelse(
        function_rescue_gate,
        "Tier_C_function_rescue_visibility_suppressed",
        fifelse(
          ibd2_multicelltype_rescue,
          "Tier_D_IBD2_rescue_without_full_ICI_module_support",
          "not_prioritized"
        )
      )
    )
  )
)]

screen[, rank_tier := fcase(
  evidence_tier == "Tier_A_function_rescue_visibility_preserved_no_disease_program_aggravation", 1L,
  evidence_tier == "Tier_B_function_rescue_but_disease_program_aggravated", 2L,
  evidence_tier == "Tier_C_function_rescue_visibility_suppressed", 3L,
  evidence_tier == "Tier_D_IBD2_rescue_without_full_ICI_module_support", 4L,
  evidence_tier == "positive_control", 5L,
  default = 6L
)]
setorder(screen, rank_tier, -ici_loss_n_celltypes, -ici_loss_median, -median_effect_IBD2, perturbation)
screen[, rank_tier := NULL]

atomic_fwrite(screen, file.path(DIRS$tables, "Table_S23_GSE313368_irAE_rescue_screen.tsv"))
atomic_fwrite(screen, file.path(DIRS$results, "GSE313368_irAE_rescue_screen_all.tsv.gz"))

cat("GSE313368 irAE rescue screen complete.\n")
print(screen[, .N, by = evidence_tier][order(evidence_tier)])
print(screen[evidence_tier != "not_prioritized", .(
  perturbation, evidence_tier, ici_loss_n_celltypes, ici_loss_median,
  median_effect_IBD2, ici_up_non_ifn_median, ifn_visibility_median, median_effect_GP20,
  median_effect_IBD1
)])
