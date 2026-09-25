cat(R.version.string, '\n')
for (p in c('Seurat','SeuratObject','Matrix','edgeR','limma','DESeq2','fgsea','ggplot2','cowplot','jsonlite')) {
  cat(p, if (requireNamespace(p,quietly=TRUE)) as.character(packageVersion(p)) else 'MISSING', '\n')
}
print(sessionInfo())
