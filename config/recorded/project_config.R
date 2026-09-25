options(stringsAsFactors = FALSE)
required_root <- function(key) {
 value <- Sys.getenv(key, unset = "")
 if (!nzchar(value)) stop(paste("Missing explicit environment path:", key))
 normalizePath(value, winslash = "/", mustWork = TRUE)
}
PROJECT_ROOT <- required_root("WCD_PROJECT_ROOT")
INPUT_ROOT <- required_root("ICI_INPUT_ROOT")
WORK_ROOT <- required_root("ICI_WORK_ROOT")
if (!startsWith(paste0(PROJECT_ROOT, "/"), paste0(WORK_ROOT, "/"))) stop("Project must be within this fresh run")
if (startsWith(paste0(PROJECT_ROOT, "/"), paste0(INPUT_ROOT, "/"))) stop("Cannot write under inputs")
DIRS <- list(
 scripts = file.path(PROJECT_ROOT, "00_scripts"), data = INPUT_ROOT,
 results = file.path(PROJECT_ROOT, "02_results"), figures = file.path(PROJECT_ROOT, "figures_for_article"),
 tables = file.path(PROJECT_ROOT, "tables_for_article"), manuscript = file.path(PROJECT_ROOT, "manuscript"),
 logs = file.path(PROJECT_ROOT, "logs")
)
stopifnot(all(vapply(DIRS, dir.exists, logical(1))))

SEED <- 20260826L
set.seed(SEED)

MIN_CELLS_PER_SAMPLE_CELLTYPE <- 20L
MIN_SAMPLES_PER_GROUP <- 3L
FDR_THRESHOLD <- 0.05
MIN_ABS_LOG2FC <- log2(1.25)

write_tsv_atomic <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  data.table::fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
  invisible(path)
}

save_rds_atomic <- function(x, path, compress = "xz") {
  tmp <- paste0(path, ".tmp")
  saveRDS(x, tmp, compress = compress)
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
  invisible(path)
}
