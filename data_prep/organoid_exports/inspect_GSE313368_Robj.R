options(stringsAsFactors = FALSE, width = 200)

project_root <- Sys.getenv("WCD_PROJECT_ROOT", unset = normalizePath(".", winslash = "/", mustWork = TRUE))
input_file <- file.path(project_root, "01_data", "GSE313368_sunbio_3_donor_3_atlas.Robj")
output_file <- file.path(project_root, "02_results", "GSE313368_object_audit.txt")

audit <- character()
emit <- function(...) {
  line <- paste0(...)
  audit <<- c(audit, line)
  cat(line, "\n", sep = "")
}

e <- new.env(parent = emptyenv())
loaded <- load(input_file, envir = e)
emit("Loaded objects: ", paste(loaded, collapse = ", "))

for (nm in loaded) {
  obj <- e[[nm]]
  emit("\nOBJECT ", nm)
  emit("class: ", paste(class(obj), collapse = ", "))
  emit("size_bytes: ", format(as.numeric(object.size(obj)), scientific = FALSE))
  if (!is.null(dim(obj))) emit("dim: ", paste(dim(obj), collapse = " x "))

  if (inherits(obj, "Seurat")) {
    emit("assays: ", paste(names(obj@assays), collapse = ", "))
    emit("active_assay: ", Seurat::DefaultAssay(obj))
    emit("reductions: ", paste(names(obj@reductions), collapse = ", "))
    emit("meta_columns: ", paste(colnames(obj@meta.data), collapse = " | "))
    emit("meta_dim: ", paste(dim(obj@meta.data), collapse = " x "))
    for (col in colnames(obj@meta.data)) {
      vals <- obj@meta.data[[col]]
      n_unique <- length(unique(vals))
      if (n_unique <= 120L) {
        tab <- sort(table(vals, useNA = "ifany"), decreasing = TRUE)
        text <- paste(paste(names(tab), as.integer(tab), sep = "="), collapse = "; ")
        emit("META ", col, " [", n_unique, "]: ", text)
      } else {
        emit("META ", col, " [", n_unique, "]")
      }
    }
    for (assay_nm in names(obj@assays)) {
      assay_obj <- obj@assays[[assay_nm]]
      emit("ASSAY ", assay_nm, " class=", paste(class(assay_obj), collapse = ","))
      emit("ASSAY ", assay_nm, " features=", nrow(assay_obj), " cells=", ncol(assay_obj))
      if ("counts" %in% slotNames(assay_obj)) {
        counts <- methods::slot(assay_obj, "counts")
        emit("ASSAY ", assay_nm, " counts_dim=", paste(dim(counts), collapse = "x"),
             " nnzero=", Matrix::nnzero(counts))
      }
    }
  }
}

writeLines(audit, output_file, useBytes = TRUE)
