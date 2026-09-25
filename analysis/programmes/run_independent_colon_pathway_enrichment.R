#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_independent_colon_pathway_enrichment.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA",
         compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

de <- fread(file.path(DIRS$results, "independent_colon_patient_level_DE_all.tsv.gz"))
de[, gene_symbol := sub("^.*\\|", "", gene)]

hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
reactome <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
gene_sets <- unique(rbindlist(list(
  hallmark[, .(collection = "HALLMARK", pathway = gs_name, gene = gene_symbol)],
  reactome[, .(collection = "REACTOME", pathway = gs_name, gene = gene_symbol)]
)))
pathways <- split(gene_sets$gene, paste(gene_sets$collection, gene_sets$pathway, sep = "::"))

groups <- unique(de[, .(dataset, organ, feature_level, celltype, contrast,
                        n_case, n_control, regimen_restriction)])
out <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
  g <- groups[i]
  block <- de[
    feature_level == g$feature_level & celltype == g$celltype & contrast == g$contrast
  ][is.finite(F) & is.finite(logFC)]
  stats <- sign(block$logFC) * sqrt(pmax(block$F, 0))
  names(stats) <- block$gene_symbol
  stats <- sort(stats[!duplicated(names(stats))], decreasing = TRUE)
  fg <- as.data.table(fgseaMultilevel(
    pathways, stats, minSize = 10, maxSize = 500, eps = 0,
    scoreType = "std", nproc = 1
  ))
  fg[, c("collection", "pathway") := tstrsplit(pathway, "::", fixed = TRUE)]
  fg[, `:=`(
    dataset = g$dataset, organ = g$organ, feature_level = g$feature_level,
    celltype = g$celltype, contrast = g$contrast, n_case = g$n_case,
    n_control = g$n_control, regimen_restriction = g$regimen_restriction,
    direction = fifelse(NES > 0, "up_in_irColitis", "down_in_irColitis"),
    passes_fdr = padj <= FDR_THRESHOLD,
    leadingEdge = vapply(leadingEdge, paste, collapse = ";", character(1))
  )]
  out[[i]] <- fg
  cat(sprintf("%d/%d %s %s %s\n", i, nrow(groups), g$feature_level, g$celltype, g$regimen_restriction))
}

enrichment <- rbindlist(out, use.names = TRUE, fill = TRUE)
setcolorder(enrichment, c(
  "dataset", "organ", "feature_level", "celltype", "contrast", "n_case", "n_control",
  "regimen_restriction", "collection", "pathway", "direction", "NES", "pval", "padj",
  "size", "passes_fdr", "leadingEdge"
))
enrichment[, abs_NES_order := abs(NES)]
setorder(enrichment, contrast, feature_level, celltype, padj, -abs_NES_order)
enrichment[, abs_NES_order := NULL]
atomic_fwrite(enrichment, file.path(DIRS$results, "independent_colon_pathway_enrichment_all.tsv.gz"))

summary <- enrichment[passes_fdr == TRUE, .(
  n_significant_pathways = .N, n_up = sum(NES > 0), n_down = sum(NES < 0),
  min_padj = min(padj)
), by = .(dataset, organ, feature_level, celltype, contrast, regimen_restriction, collection)]
atomic_fwrite(summary, file.path(DIRS$tables, "Table_S14_independent_colon_pathway_summary.tsv"))

top <- enrichment[passes_fdr == TRUE]
top[, abs_NES_order := abs(NES)]
setorder(top, contrast, feature_level, celltype, padj, -abs_NES_order)
top <- top[, head(.SD, 15L), by = .(contrast, feature_level, celltype, collection)]
top[, abs_NES_order := NULL]
atomic_fwrite(top, file.path(DIRS$tables, "Table_S15_independent_colon_top_pathways.tsv"))
cat("Completed independent colon pathway enrichment.\n")
print(summary)
