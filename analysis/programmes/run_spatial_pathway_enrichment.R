#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_spatial_pathway_enrichment.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA",
         compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

de <- fread(file.path(DIRS$results, "spatial_colon_patient_level_DE.tsv.gz"))
hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
reactome <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
gene_sets <- unique(rbindlist(list(
  hallmark[, .(collection = "HALLMARK", pathway = gs_name, gene = gene_symbol)],
  reactome[, .(collection = "REACTOME", pathway = gs_name, gene = gene_symbol)]
)))
pathways <- split(gene_sets$gene, paste(gene_sets$collection, gene_sets$pathway, sep = "::"))

stats <- sign(de$logFC) * sqrt(pmax(de$F, 0))
names(stats) <- de$gene
stats <- sort(stats[!duplicated(names(stats))], decreasing = TRUE)
fg <- as.data.table(fgseaMultilevel(
  pathways, stats, minSize = 10, maxSize = 500, eps = 0,
  scoreType = "std", nproc = 1
))
fg[, c("collection", "pathway") := tstrsplit(pathway, "::", fixed = TRUE)]
fg[, `:=`(
  dataset = "GSE210037", organ = "colon_spatial", feature_level = "whole_spot_sample",
  contrast = "irAE_colitis_vs_healthy_colon", n_case = 7L, n_control = 8L,
  direction = fifelse(NES > 0, "up_in_irAE", "down_in_irAE"),
  passes_fdr = padj <= FDR_THRESHOLD,
  leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1))
)]
setcolorder(fg, c(
  "dataset", "organ", "feature_level", "contrast", "n_case", "n_control",
  "collection", "pathway", "direction", "NES", "pval", "padj", "size",
  "passes_fdr", "leadingEdge"
))
fg[, abs_NES_order := abs(NES)]
setorder(fg, padj, -abs_NES_order)
fg[, abs_NES_order := NULL]
atomic_fwrite(fg, file.path(DIRS$results, "spatial_colon_pathway_enrichment_all.tsv.gz"))

summary <- fg[passes_fdr == TRUE, .(
  n_significant_pathways = .N, n_up = sum(NES > 0), n_down = sum(NES < 0),
  min_padj = min(padj)
), by = .(dataset, organ, feature_level, contrast, collection)]
atomic_fwrite(summary, file.path(DIRS$tables, "Table_S18_spatial_colon_pathway_summary.tsv"))
top <- fg[passes_fdr == TRUE]
top[, abs_NES_order := abs(NES)]
setorder(top, collection, padj, -abs_NES_order)
top <- top[, head(.SD, 25L), by = collection]
top[, abs_NES_order := NULL]
atomic_fwrite(top, file.path(DIRS$tables, "Table_S19_spatial_colon_top_pathways.tsv"))
cat("Completed spatial colon pathway enrichment.\n")
print(summary)
