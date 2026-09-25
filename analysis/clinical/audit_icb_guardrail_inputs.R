#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(data.table))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: audit_icb_guardrail_inputs.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

strip_geo_quotes <- function(x) {
  x <- sub('^"', "", x)
  sub('"$', "", x)
}

safe_field_name <- function(x) {
  x <- tolower(trimws(x))
  x <- gsub("[^a-z0-9]+", "_", x)
  gsub("^_+|_+$", "", x)
}

parse_series_metadata <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con), add = TRUE)
  lines <- readLines(con, warn = FALSE)
  sample_lines <- lines[grepl("^!Sample_", lines)]
  fields <- lapply(sample_lines, function(line) {
    strip_geo_quotes(strsplit(line, "\t", fixed = TRUE)[[1]])
  })
  geo_row <- fields[[which(vapply(fields, function(x) x[[1]] == "!Sample_geo_accession", logical(1)))[[1]]]]
  n <- length(geo_row) - 1L
  records <- lapply(seq_len(n), function(i) list(geo_accession = geo_row[[i + 1L]]))

  add_value <- function(record, key, value) {
    if (!nzchar(key) || !nzchar(value)) return(record)
    old <- record[[key]]
    if (is.null(old) || !nzchar(old)) {
      record[[key]] <- value
    } else if (!value %chin% strsplit(old, " | ", fixed = TRUE)[[1]]) {
      record[[key]] <- paste(old, value, sep = " | ")
    }
    record
  }

  for (field in fields) {
    key <- field[[1]]
    values <- field[-1L]
    if (length(values) < n) values <- c(values, rep("", n - length(values)))
    if (length(values) > n) values <- values[seq_len(n)]
    if (key == "!Sample_characteristics_ch1") {
      for (i in seq_len(n)) {
        value <- trimws(values[[i]])
        colon <- regexpr(":", value, fixed = TRUE)[[1]]
        if (colon > 0L) {
          label <- safe_field_name(substr(value, 1L, colon - 1L))
          item <- trimws(substr(value, colon + 1L, nchar(value)))
          records[[i]] <- add_value(records[[i]], label, item)
        }
      }
    } else if (key %chin% c("!Sample_title", "!Sample_source_name_ch1")) {
      label <- safe_field_name(sub("^!Sample_", "", key))
      for (i in seq_len(n)) records[[i]] <- add_value(records[[i]], label, trimws(values[[i]]))
    }
  }
  rbindlist(records, fill = TRUE)
}

out <- file.path(DIRS$results, "ICB_guardrail_input_audit.txt")
sink(out, split = TRUE)
on.exit(sink(), add = TRUE)

for (gse in c("GSE78220", "GSE91061", "GSE135222", "GSE126044")) {
  cat("\n### ", gse, " metadata\n", sep = "")
  meta <- parse_series_metadata(file.path(DIRS$data, paste0(gse, "_series_matrix.txt.gz")))
  cat("samples=", nrow(meta), " fields=", ncol(meta), "\n", sep = "")
  cat("field_names=", paste(names(meta), collapse = "|"), "\n", sep = "")
  response_fields <- grep("response|survival|progress|biopsy|treatment|sample|patient", names(meta), value = TRUE)
  for (field in response_fields) {
    values <- unique(meta[[field]])
    values <- values[!is.na(values) & nzchar(values)]
    cat(field, "=", paste(head(values, 30L), collapse = "|"), "\n", sep = "")
  }
}

cat("\n### Expression files\n")
for (file in c(
  "GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz",
  "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz",
  "GSE126044_counts.txt.gz"
)) {
  x <- fread(file.path(DIRS$data, file), nrows = 5L, check.names = FALSE)
  cat(file, ": first_id=", x[[1]][[1]], "; sample_columns=", ncol(x) - 1L,
      "; first_samples=", paste(head(names(x)[-1L], 8L), collapse = "|"), "\n", sep = "")
}

if (requireNamespace("readxl", quietly = TRUE)) {
  source_xlsx <- file.path(DIRS$data, "GSE78220_PatientFPKM.xlsx")
  local_xlsx <- tempfile(pattern = "GSE78220_", fileext = ".xlsx")
  copied <- file.copy(source_xlsx, local_xlsx, overwrite = TRUE)
  if (copied) {
    x <- readxl::read_excel(local_xlsx, n_max = 5L)
    unlink(local_xlsx)
    cat("GSE78220_PatientFPKM.xlsx: sample_columns=", ncol(x) - 1L,
        "; first_samples=", paste(head(names(x)[-1L], 8L), collapse = "|"), "\n", sep = "")
  } else {
    cat("GSE78220_PatientFPKM.xlsx: temporary ASCII-path copy failed\n")
  }
}

cat("\nAudit complete.\n")
