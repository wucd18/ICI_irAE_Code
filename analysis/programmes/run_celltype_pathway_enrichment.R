#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(msigdbr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_celltype_pathway_enrichment.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA", compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

de <- fread(file.path(DIRS$results, "patient_level_pseudobulk_DE_all.tsv.gz"))

hallmark <- as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))
reactome <- as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))
gene_sets <- rbindlist(list(
  hallmark[, .(collection = "HALLMARK", pathway = gs_name, gene = gene_symbol)],
  reactome[, .(collection = "REACTOME", pathway = gs_name, gene = gene_symbol)]
), use.names = TRUE)
gene_sets <- unique(gene_sets[nzchar(gene)])
pathways <- split(gene_sets$gene, paste(gene_sets$collection, gene_sets$pathway, sep = "::"))

groups <- unique(de[, .(dataset, organ, celltype, contrast, n_case, n_control)])
enrichments <- vector("list", nrow(groups))
for (i in seq_len(nrow(groups))) {
  g <- groups[i]
  block <- de[
    dataset == g$dataset & organ == g$organ & celltype == g$celltype & contrast == g$contrast
  ]
  block <- block[is.finite(F) & is.finite(logFC)]
  stats <- sign(block$logFC) * sqrt(pmax(block$F, 0))
  names(stats) <- block$gene
  stats <- sort(stats[!duplicated(names(stats))], decreasing = TRUE)
  fg <- as.data.table(fgseaMultilevel(
    pathways = pathways,
    stats = stats,
    minSize = 10,
    maxSize = 500,
    eps = 0,
    scoreType = "std",
    nproc = 1
  ))
  fg[, c("collection", "pathway") := tstrsplit(pathway, "::", fixed = TRUE, keep = c(1L, 2L))]
  fg[, `:=`(
    dataset = g$dataset,
    organ = g$organ,
    celltype = g$celltype,
    contrast = g$contrast,
    n_case = g$n_case,
    n_control = g$n_control,
    direction = fifelse(NES > 0, "up_in_irAE", "down_in_irAE"),
    passes_fdr = padj <= FDR_THRESHOLD
  )]
  fg[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", character(1))]
  enrichments[[i]] <- fg
  cat(sprintf("%d/%d %s %s\n", i, nrow(groups), g$dataset, g$celltype))
}

enrichment <- rbindlist(enrichments, use.names = TRUE, fill = TRUE)
setcolorder(enrichment, c(
  "dataset", "organ", "celltype", "contrast", "n_case", "n_control",
  "collection", "pathway", "direction", "NES", "pval", "padj", "size",
  "passes_fdr", "leadingEdge"
))
enrichment[, abs_NES_order := abs(NES)]
setorder(enrichment, dataset, celltype, padj, -abs_NES_order)
enrichment[, abs_NES_order := NULL]
atomic_fwrite(enrichment, file.path(DIRS$results, "celltype_pathway_enrichment_all.tsv.gz"))

summary <- enrichment[passes_fdr == TRUE, .(
  n_significant_pathways = .N,
  n_up = sum(NES > 0),
  n_down = sum(NES < 0),
  min_padj = min(padj)
), by = .(dataset, organ, celltype, contrast, collection)]
atomic_fwrite(summary, file.path(DIRS$tables, "Table_S8_pathway_enrichment_summary.tsv"))

top <- enrichment[passes_fdr == TRUE]
top[, abs_NES_order := abs(NES)]
setorder(top, dataset, celltype, padj, -abs_NES_order)
top <- top[, head(.SD, 15L), by = .(dataset, organ, celltype, contrast, collection)]
top[, abs_NES_order := NULL]
atomic_fwrite(top, file.path(DIRS$tables, "Table_S9_top_celltype_pathways.tsv"))

cat("Completed cell-type pathway enrichment.\n")
print(summary)
