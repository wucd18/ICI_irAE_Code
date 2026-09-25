#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
options(error = function() {
  traceback(8)
  quit(save = "no", status = 1L)
})
suppressPackageStartupMessages({
  library(data.table)
  library(edgeR)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(Biobase)
  library(survival)
  library(metafor)
  library(pROC)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: run_icb_efficacy_guardrail.R <project_dir>")
source(file.path(args[[1]], "00_scripts", "project_config.R"))

atomic_fwrite <- function(x, path) {
  compressed <- grepl("\\.gz$", path, ignore.case = TRUE)
  tmp <- if (compressed) paste0(path, ".tmp.gz") else paste0(path, ".tmp")
  fwrite(x, tmp, sep = "\t", quote = FALSE, na = "NA",
         compress = if (compressed) "gzip" else "none")
  if (file.exists(path)) file.remove(path)
  stopifnot(file.rename(tmp, path))
}

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
  geo_index <- which(vapply(fields, function(x) x[[1]] == "!Sample_geo_accession", logical(1)))[[1]]
  geo_row <- fields[[geo_index]]
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

map_gene_ids <- function(ids, keytype) {
  ids <- as.character(ids)
  valid <- !is.na(ids) & nzchar(ids)
  out <- rep(NA_character_, length(ids))
  keys <- unique(ids[valid])
  mapped <- suppressMessages(mapIds(
    org.Hs.eg.db, keys = keys, column = "SYMBOL", keytype = keytype,
    multiVals = "first"
  ))
  out[valid] <- unname(mapped[ids[valid]])
  out
}

collapse_by_symbol <- function(mat, symbols, mode = c("mean", "sum")) {
  mode <- match.arg(mode)
  keep <- !is.na(symbols) & nzchar(symbols)
  mat <- mat[keep, , drop = FALSE]
  symbols <- symbols[keep]
  summed <- rowsum(mat, group = symbols, reorder = FALSE, na.rm = TRUE)
  if (mode == "mean") {
    counts <- as.numeric(table(symbols)[rownames(summed)])
    summed <- summed / counts
  }
  storage.mode(summed) <- "double"
  summed
}

logcpm_counts <- function(counts) {
  counts[!is.finite(counts)] <- 0
  counts[counts < 0] <- 0
  dge <- DGEList(counts = counts)
  dge <- calcNormFactors(dge, method = "TMM")
  cpm(dge, log = TRUE, prior.count = 1)
}

normalize_response <- function(x) {
  y <- tolower(trimws(as.character(x)))
  out <- rep(NA_integer_, length(y))
  nonresponse <- grepl("non[- _]?responder|progressive disease|stable disease|^pd$|^sd$|^sd/pd$", y)
  response <- grepl("partial response|complete response|^pr$|^cr$|^prcr$|^cr/pr$|^responder$", y)
  out[nonresponse] <- 0L
  out[response] <- 1L
  out
}

aggregate_patient_expression <- function(expr, patient_id, meta) {
  stopifnot(ncol(expr) == length(patient_id), nrow(meta) == length(patient_id))
  patient_id <- as.character(patient_id)
  patient_levels <- unique(patient_id)
  group <- factor(patient_id, levels = patient_levels)
  summed <- rowsum(t(expr), group = group, reorder = FALSE, na.rm = TRUE)
  averaged <- t(summed / as.numeric(table(group)[rownames(summed)]))
  colnames(averaged) <- patient_levels
  first <- match(patient_levels, patient_id)
  patient_meta <- copy(meta[first])
  patient_meta[, patient_id := patient_levels]
  list(expr = averaged, meta = patient_meta)
}

finalize_cohort <- function(name, cancer, therapy, expression_type, expr, meta,
                            patient_id, response = NULL, survival_time = NULL,
                            survival_event = NULL, survival_endpoint = NA_character_,
                            caveat = "") {
  meta <- copy(meta)
  response_value <- response
  survival_time_value <- survival_time
  survival_event_value <- survival_event
  meta[, patient_id := as.character(patient_id)]
  meta[, responder := if (is.null(response_value)) NA_integer_ else as.integer(response_value)]
  meta[, survival_time := if (is.null(survival_time_value)) NA_real_ else as.numeric(survival_time_value)]
  meta[, survival_event := if (is.null(survival_event_value)) NA_integer_ else as.integer(survival_event_value)]
  aggregated <- aggregate_patient_expression(expr, meta$patient_id, meta)
  stopifnot(identical(colnames(aggregated$expr), aggregated$meta$patient_id))
  list(
    name = name, cancer = cancer, therapy = therapy,
    expression_type = expression_type, expr = aggregated$expr,
    meta = aggregated$meta, survival_endpoint = survival_endpoint,
    caveat = caveat
  )
}

cat("Preparing public ICI efficacy cohorts...\n")

# GSE78220: pre-treatment melanoma; one on-treatment sample is excluded and
# the two Pt27 baseline biopsies are averaged at the patient level.
g782_meta <- parse_series_metadata(file.path(DIRS$data, "GSE78220_series_matrix.txt.gz"))
g782_xlsx <- file.path(DIRS$data, "GSE78220_PatientFPKM.xlsx")
g782_tmp <- tempfile(pattern = "GSE78220_", fileext = ".xlsx")
stopifnot(file.copy(g782_xlsx, g782_tmp, overwrite = TRUE))
g782_dt <- as.data.table(readxl::read_excel(g782_tmp))
unlink(g782_tmp)
g782_genes <- as.character(g782_dt[[1]])
g782_expr <- as.matrix(g782_dt[, -1L])
storage.mode(g782_expr) <- "double"
rownames(g782_expr) <- g782_genes
colnames(g782_expr) <- names(g782_dt)[-1L]
g782_keep <- grepl("\\.baseline$", colnames(g782_expr), ignore.case = TRUE)
g782_expr <- log2(g782_expr[, g782_keep, drop = FALSE] + 1)
g782_title <- sub("\\.baseline$", "", colnames(g782_expr), ignore.case = TRUE)
g782_match <- match(g782_title, g782_meta$title)
stopifnot(!anyNA(g782_match))
g782_m <- g782_meta[g782_match]
g782 <- finalize_cohort(
  "GSE78220", "melanoma", "pembrolizumab", "log2(FPKM+1)",
  g782_expr, g782_m, g782_m$patient_id,
  normalize_response(g782_m$anti_pd_1_response),
  suppressWarnings(as.numeric(g782_m$overall_survival_days)),
  fifelse(tolower(g782_m$vital_status) == "dead", 1L,
          fifelse(tolower(g782_m$vital_status) == "alive", 0L, NA_integer_)),
  "OS", "Pt27 duplicate baseline biopsies averaged; Pt16 on-treatment sample excluded"
)

# GSE91061: pre-treatment melanoma samples only; unknown response excluded in
# response models but retained in the expression audit.
g910_meta <- parse_series_metadata(file.path(DIRS$data, "GSE91061_series_matrix.txt.gz"))
g910_dt <- fread(file.path(DIRS$data, "GSE91061_BMS038109Sample.hg19KnownGene.rld.csv.gz"),
                 check.names = FALSE)
g910_ids <- as.character(g910_dt[[1]])
g910_mat <- as.matrix(g910_dt[, -1L])
storage.mode(g910_mat) <- "double"
g910_symbols <- map_gene_ids(g910_ids, "ENTREZID")
g910_expr <- collapse_by_symbol(g910_mat, g910_symbols, "mean")
colnames(g910_expr) <- names(g910_dt)[-1L]
g910_match <- match(colnames(g910_expr), g910_meta$title)
stopifnot(!anyNA(g910_match))
g910_m <- g910_meta[g910_match]
g910_keep <- tolower(g910_m$visit_pre_or_on_treatment) == "pre"
g910_expr <- g910_expr[, g910_keep, drop = FALSE]
g910_m <- g910_m[g910_keep]
g910_patient <- sub("_.*$", "", g910_m$title)
g910 <- finalize_cohort(
  "GSE91061", "melanoma", "nivolumab", "regularized-log expression",
  g910_expr, g910_m, g910_patient, normalize_response(g910_m$response),
  survival_endpoint = NA_character_,
  caveat = "pre-treatment subset; GEO response UNK excluded from response models"
)

# GSE135222: advanced NSCLC anti-PD-(L)1 cohort with PFS event/time metadata.
g135_meta <- parse_series_metadata(file.path(DIRS$data, "GSE135222_series_matrix.txt.gz"))
g135_dt <- fread(file.path(DIRS$data, "GSE135222_GEO_RNA-seq_omicslab_exp.tsv.gz"),
                 check.names = FALSE)
g135_ids <- sub("\\..*$", "", as.character(g135_dt[[1]]))
g135_mat <- as.matrix(g135_dt[, -1L])
storage.mode(g135_mat) <- "double"
g135_symbols <- map_gene_ids(g135_ids, "ENSEMBL")
g135_expr <- log2(collapse_by_symbol(g135_mat, g135_symbols, "mean") + 1)
colnames(g135_expr) <- names(g135_dt)[-1L]
g135_key <- gsub("[^A-Za-z0-9]", "", g135_meta$title)
g135_match <- match(colnames(g135_expr), g135_key)
stopifnot(!anyNA(g135_match))
g135_m <- g135_meta[g135_match]
g135 <- finalize_cohort(
  "GSE135222", "NSCLC", "anti-PD-1/PD-L1", "log2(processed expression+1)",
  g135_expr, g135_m, g135_key[g135_match], response = NULL,
  survival_time = suppressWarnings(as.numeric(g135_m$pfs_time)),
  survival_event = suppressWarnings(as.integer(g135_m$progression_free_survival_pfs)),
  survival_endpoint = "PFS",
  caveat = "GEO provides PFS but not RECIST response"
)

# GSE126044: pre-treatment NSCLC raw counts; sample material is retained for
# sensitivity analysis because the cohort mixes fresh and FFPE specimens.
g126_meta <- parse_series_metadata(file.path(DIRS$data, "GSE126044_series_matrix.txt.gz"))
g126_dt <- fread(file.path(DIRS$data, "GSE126044_counts.txt.gz"), check.names = FALSE)
g126_genes <- as.character(g126_dt[[1]])
g126_counts <- as.matrix(g126_dt[, -1L])
storage.mode(g126_counts) <- "double"
g126_counts <- collapse_by_symbol(g126_counts, g126_genes, "sum")
g126_expr <- logcpm_counts(g126_counts)
colnames(g126_expr) <- names(g126_dt)[-1L]
g126_key <- sub("^RNA-seq_", "", g126_meta$title, ignore.case = TRUE)
g126_match <- match(colnames(g126_expr), g126_key)
stopifnot(!anyNA(g126_match))
g126_m <- g126_meta[g126_match]
g126 <- finalize_cohort(
  "GSE126044", "NSCLC", "anti-PD-1", "TMM logCPM",
  g126_expr, g126_m, g126_key[g126_match], normalize_response(g126_m$patient_response),
  survival_endpoint = NA_character_,
  caveat = "11 fresh and 5 FFPE samples; material-adjusted sensitivity reported"
)

# IMvigor210: raw counts and trial annotations extracted without requiring the
# retired DESeq package by accessing the serialized ExpressionSet-compatible slots.
imv_env <- new.env(parent = emptyenv())
load(file.path(DIRS$data, "IMvigor210_cds.RData"), envir = imv_env)
imv_cds <- imv_env[["cds"]]
imv_counts <- get("counts", envir = attr(imv_cds, "assayData"))
imv_pd <- as.data.table(pData(attr(imv_cds, "phenoData")), keep.rownames = "sample_id")
imv_fd <- as.data.table(pData(attr(imv_cds, "featureData")), keep.rownames = "feature_id")
imv_counts <- collapse_by_symbol(imv_counts, imv_fd$Symbol, "sum")
imv_expr <- logcpm_counts(imv_counts)
imv_match <- match(colnames(imv_expr), imv_pd$sample_id)
stopifnot(!anyNA(imv_match))
imv_m <- imv_pd[imv_match]
imvigor <- finalize_cohort(
  "IMvigor210", "metastatic urothelial carcinoma", "atezolizumab", "TMM logCPM",
  imv_expr, imv_m, imv_m$ANONPT_ID, normalize_response(imv_m$binaryResponse),
  imv_m$os, imv_m$censOS, "OS",
  "non-evaluable RECIST cases excluded only from response models"
)

cohorts <- list(g782, g910, g135, g126, imvigor)

cohort_audit <- rbindlist(lapply(cohorts, function(x) {
  m <- x$meta
  data.table(
    cohort = x$name, cancer = x$cancer, therapy = x$therapy,
    expression_type = x$expression_type,
    n_patients_expression = nrow(m),
    n_response_evaluable = sum(!is.na(m$responder)),
    n_responders = sum(m$responder == 1L, na.rm = TRUE),
    n_nonresponders = sum(m$responder == 0L, na.rm = TRUE),
    survival_endpoint = x$survival_endpoint,
    n_survival_evaluable = sum(is.finite(m$survival_time) & !is.na(m$survival_event)),
    n_survival_events = sum(m$survival_event == 1L, na.rm = TRUE),
    n_genes = nrow(x$expr), caveat = x$caveat
  )
}))
atomic_fwrite(cohort_audit, file.path(DIRS$tables, "Table_S24_ICB_guardrail_cohort_audit.tsv"))

cat("Deriving candidate features from the independent organoid perturbation atlas...\n")
screen <- fread(file.path(DIRS$tables, "Table_S23_GSE313368_irAE_rescue_screen.tsv"))
candidates <- screen[grepl("^Tier_[AB]", evidence_tier)]$perturbation
candidates <- unique(candidates)
stopifnot("TNFSF12" %chin% candidates)

receptor_path <- Sys.getenv('ICI_RECEPTOR_MAP', unset='')
if (!nzchar(receptor_path) || !file.exists(receptor_path)) stop('BLOCKED_INPUT: ICI_RECEPTOR_MAP')
receptor_map <- jsonlite::fromJSON(receptor_path, simplifyVector=FALSE)
if (!length(receptor_map) || is.null(names(receptor_map)) || anyDuplicated(names(receptor_map))) stop('Invalid receptor map')
receptor_map <- lapply(receptor_map, function(x) {
 x <- unlist(x, use.names=FALSE)
 if (!length(x) || anyNA(x) || any(!nzchar(x)) || anyDuplicated(x)) stop('Invalid receptor members')
 x
})

feature_rows <- list()
add_feature <- function(feature_id, candidate, feature_type, celltype,
                        up_genes, down_genes = character(), derivation) {
  feature_rows[[length(feature_rows) + 1L]] <<- data.table(
    feature_id = feature_id, candidate = candidate,
    feature_type = feature_type, celltype = celltype,
    up_genes = list(sort(unique(up_genes))),
    down_genes = list(sort(unique(down_genes))),
    derivation = derivation
  )
}

for (candidate in candidates) {
  add_feature(
    paste0("GENE_", candidate), candidate, "direct_ligand_gene", NA_character_,
    candidate, derivation = "within-cohort percentile rank of the candidate ligand gene"
  )
  receptors <- receptor_map[[candidate]]
  if (length(receptors)) {
    add_feature(
      paste0("AXIS_", candidate), candidate, "ligand_receptor_axis", NA_character_,
      c(candidate, receptors),
      derivation = "mean within-cohort percentile rank of ligand and canonical receptor genes"
    )
  }
}
add_feature(
  "GENE_TNFRSF12A", "TNFSF12", "direct_receptor_gene", NA_character_,
  "TNFRSF12A", derivation = "within-cohort percentile rank of the TWEAK receptor gene"
)

org_de <- fread(file.path(DIRS$results, "GSE313368_celltype_perturbation_DE_all.tsv.gz"))
org_de <- org_de[
  perturbation %chin% candidates & supplement_qc_exclude_gene_level == FALSE &
    celltype %chin% c("Transit_Amplifying", "Colonocyte") &
    is.finite(log2FoldChange) & is.finite(lfcSE) & lfcSE > 0
]
org_de[, rank_z := log2FoldChange / lfcSE]
for (candidate in candidates) {
  for (ct in c("Transit_Amplifying", "Colonocyte")) {
    block <- copy(org_de[perturbation == candidate & celltype == ct])
    if (nrow(block) < 1000L) next
    block <- block[is.finite(padj) & padj <= 0.10]
    up <- block[rank_z > 0][order(-rank_z)]$gene
    down <- block[rank_z < 0][order(rank_z)]$gene
    up <- head(unique(up), 75L)
    down <- head(unique(down), 75L)
    if (length(up) < 15L || length(down) < 15L) next
    ct_short <- fifelse(ct == "Transit_Amplifying", "TA", "COLONOCYTE")
    add_feature(
      paste0("PERTSIG_", candidate, "_", ct_short), candidate,
      "bidirectional_organoid_perturbation_signature", ct, up, down,
      "top positive and negative GSE313368 genes at publisher padj<=0.10; ranked by log2FC/SE; capped at 75 per direction"
    )
  }
}

gsea <- fread(file.path(DIRS$results, "GSE313368_irAE_module_GSEA_all.tsv.gz"))
tweak_rescue <- gsea[
  perturbation == "TNFSF12" & celltype == "Transit_Amplifying" &
    module == "ICI_COLITIS_EPITHELIAL_LOSS_STRICT"
]
if (nrow(tweak_rescue) == 1L) {
  rescue_genes <- strsplit(tweak_rescue$leadingEdge[[1]], ";", fixed = TRUE)[[1]]
  add_feature(
    "TNFSF12_ICI_LOSS_RESCUE_CORE", "TNFSF12", "function_rescue_leading_edge",
    "Transit_Amplifying", rescue_genes,
    derivation = "leading edge of positive enrichment of the frozen ICI-colitis epithelial-loss module after TNFSF12 perturbation"
  )
}

features <- rbindlist(feature_rows, use.names = TRUE, fill = TRUE)
feature_table <- copy(features)
feature_table[, `:=`(
  n_up_genes = lengths(up_genes), n_down_genes = lengths(down_genes),
  up_genes = vapply(up_genes, paste, collapse = ";", character(1)),
  down_genes = vapply(down_genes, paste, collapse = ";", character(1))
)]
setcolorder(feature_table, c(
  "feature_id", "candidate", "feature_type", "celltype", "n_up_genes",
  "n_down_genes", "up_genes", "down_genes", "derivation"
))
atomic_fwrite(feature_table, file.path(DIRS$tables, "Table_S25_ICB_guardrail_feature_definitions.tsv"))

score_cohort <- function(cohort, features) {
  expr <- cohort$expr
  keep <- rowSums(is.finite(expr)) == ncol(expr)
  expr <- expr[keep, , drop = FALSE]
  rank_mat <- apply(expr, 2L, rank, ties.method = "average") / nrow(expr)
  if (is.null(dim(rank_mat))) rank_mat <- matrix(rank_mat, ncol = 1L)
  rownames(rank_mat) <- rownames(expr)
  colnames(rank_mat) <- colnames(expr)
  out <- vector("list", nrow(features))
  for (i in seq_len(nrow(features))) {
    up <- intersect(features$up_genes[[i]], rownames(rank_mat))
    down <- intersect(features$down_genes[[i]], rownames(rank_mat))
    min_up <- if (features$feature_type[[i]] %chin% c(
      "direct_ligand_gene", "direct_receptor_gene", "downstream_rescue_gene"
    )) 1L else 2L
    if (length(up) < min_up || (length(features$down_genes[[i]]) > 0L && length(down) < 5L)) next
    score <- colMeans(rank_mat[up, , drop = FALSE])
    if (length(down)) score <- score - colMeans(rank_mat[down, , drop = FALSE])
    score_sd <- sd(score, na.rm = TRUE)
    score_z <- if (is.finite(score_sd) && score_sd > 0) as.numeric(scale(score)) else rep(NA_real_, length(score))
    out[[i]] <- data.table(
      cohort = cohort$name, cancer = cohort$cancer,
      patient_id = colnames(rank_mat), feature_id = features$feature_id[[i]],
      candidate = features$candidate[[i]], feature_type = features$feature_type[[i]],
      score = score, score_z = score_z,
      n_up_available = length(up), n_down_available = length(down)
    )
  }
  scores <- rbindlist(out, use.names = TRUE, fill = TRUE)
  cohort_meta <- cohort[["meta"]]
  scores[cohort_meta, on = "patient_id", `:=`(
    responder = i.responder,
    survival_time = i.survival_time,
    survival_event = i.survival_event
  )]
  endpoint_value <- cohort[["survival_endpoint"]]
  scores[, survival_endpoint := endpoint_value]
  scores
}

cat("Scoring candidate features in each ICI cohort...\n")
scores <- rbindlist(lapply(cohorts, score_cohort, features = features), use.names = TRUE, fill = TRUE)
atomic_fwrite(scores, file.path(DIRS$results, "ICB_guardrail_patient_scores.tsv.gz"))

response_effect_one <- function(dt) {
  x <- dt[!is.na(responder) & is.finite(score_z)]
  n1 <- sum(x$responder == 1L)
  n0 <- sum(x$responder == 0L)
  if (n1 < 3L || n0 < 3L) return(NULL)
  s1 <- x[responder == 1L]$score_z
  s0 <- x[responder == 0L]$score_z
  esc <- escalc(
    measure = "SMD", m1i = mean(s1), sd1i = sd(s1), n1i = n1,
    m2i = mean(s0), sd2i = sd(s0), n2i = n0, vtype = "UB"
  )
  roc_obj <- suppressMessages(roc(x$responder, x$score_z, levels = c(0, 1), direction = "<", quiet = TRUE))
  glm_fit <- suppressWarnings(glm(responder ~ score_z, family = binomial(), data = x))
  glm_coef <- summary(glm_fit)$coefficients
  data.table(
    n = nrow(x), n_responders = n1, n_nonresponders = n0,
    mean_score_z_responders = mean(s1), mean_score_z_nonresponders = mean(s0),
    hedges_g = as.numeric(esc$yi), var_hedges_g = as.numeric(esc$vi),
    se_hedges_g = sqrt(as.numeric(esc$vi)),
    ci_low = as.numeric(esc$yi - 1.96 * sqrt(esc$vi)),
    ci_high = as.numeric(esc$yi + 1.96 * sqrt(esc$vi)),
    p_smd = 2 * pnorm(-abs(as.numeric(esc$yi) / sqrt(esc$vi))),
    p_wilcoxon = suppressWarnings(wilcox.test(s1, s0, exact = FALSE)$p.value),
    auc_higher_predicts_response = as.numeric(auc(roc_obj)),
    log_odds_ratio_per_sd = glm_coef["score_z", "Estimate"],
    se_log_odds_ratio = glm_coef["score_z", "Std. Error"],
    p_logistic = glm_coef["score_z", "Pr(>|z|)"]
  )
}

response_results <- scores[, response_effect_one(.SD),
                           by = .(cohort, cancer, feature_id, candidate, feature_type)]
response_results[, fdr_smd_within_cohort := p.adjust(p_smd, method = "BH"), by = cohort]
setorder(response_results, candidate, feature_id, cohort)
atomic_fwrite(response_results, file.path(DIRS$tables, "Table_S26_ICB_guardrail_response_associations.tsv"))

meta_response <- response_results[, {
  if (.N < 2L) NULL else {
    fit <- do.call(rma.uni, list(
      yi = as.numeric(.SD[["hedges_g"]]),
      vi = as.numeric(.SD[["var_hedges_g"]]), method = "REML"
    ))
    data.table(
      n_cohorts = .N, total_n = sum(n),
      estimate = as.numeric(fit$b), se = fit$se,
      ci_low = fit$ci.lb, ci_high = fit$ci.ub, p = fit$pval,
      tau2 = fit$tau2, I2 = fit$I2, Q_p = fit$QEp,
      cohorts = paste(cohort, collapse = ";")
    )
  }
}, by = .(feature_id, candidate, feature_type)]
meta_response[, fdr := p.adjust(p, method = "BH")]
meta_response[, direction := fifelse(estimate > 0, "higher_in_responders", "higher_in_nonresponders")]
setorder(meta_response, candidate, feature_type, feature_id)
atomic_fwrite(meta_response, file.path(DIRS$tables, "Table_S27_ICB_guardrail_response_meta_analysis.tsv"))

# Material-adjusted response sensitivity for GSE126044.
g126_scores <- scores[cohort == "GSE126044" & !is.na(responder) & is.finite(score_z)]
g126_material <- g126$meta[, .(patient_id, sample_material = factor(sample))]
g126_scores[g126_material, on = "patient_id", sample_material := i.sample_material]
g126_adjusted <- g126_scores[, {
  if (uniqueN(sample_material) < 2L) NULL else {
    fit <- suppressWarnings(glm(responder ~ score_z + sample_material, family = binomial(), data = .SD))
    co <- summary(fit)$coefficients
    data.table(
      n = .N, log_odds_ratio_per_sd = co["score_z", "Estimate"],
      se_log_odds_ratio = co["score_z", "Std. Error"],
      p = co["score_z", "Pr(>|z|)"]
    )
  }
}, by = .(feature_id, candidate, feature_type)]
g126_adjusted[, fdr := p.adjust(p, method = "BH")]
atomic_fwrite(g126_adjusted, file.path(DIRS$results, "GSE126044_material_adjusted_response_sensitivity.tsv.gz"))

survival_effect_one <- function(dt) {
  x <- dt[is.finite(survival_time) & survival_time > 0 &
            survival_event %in% c(0L, 1L) & is.finite(score_z)]
  if (nrow(x) < 10L || sum(x$survival_event == 1L) < 5L) return(NULL)
  fit <- coxph(Surv(survival_time, survival_event) ~ score_z, data = x, ties = "efron")
  sm <- summary(fit)
  data.table(
    n = nrow(x), events = sum(x$survival_event == 1L),
    log_hazard_ratio_per_sd = sm$coefficients["score_z", "coef"],
    se_log_hazard_ratio = sm$coefficients["score_z", "se(coef)"],
    hazard_ratio = sm$coefficients["score_z", "exp(coef)"],
    ci_low_hr = sm$conf.int["score_z", "lower .95"],
    ci_high_hr = sm$conf.int["score_z", "upper .95"],
    p = sm$coefficients["score_z", "Pr(>|z|)"],
    concordance = sm$concordance[[1]]
  )
}

survival_results <- scores[!is.na(survival_endpoint), survival_effect_one(.SD),
                           by = .(cohort, cancer, survival_endpoint, feature_id, candidate, feature_type)]
survival_results[, fdr_within_cohort_endpoint := p.adjust(p, method = "BH"),
                 by = .(cohort, survival_endpoint)]
setorder(survival_results, candidate, feature_id, cohort)
atomic_fwrite(survival_results, file.path(DIRS$tables, "Table_S28_ICB_guardrail_survival_associations.tsv"))

meta_survival_scope <- function(dt, scope_name) {
  dt[, {
    if (.N < 2L) NULL else {
      fit <- do.call(rma.uni, list(
        yi = as.numeric(.SD[["log_hazard_ratio_per_sd"]]),
        sei = as.numeric(.SD[["se_log_hazard_ratio"]]), method = "REML"
      ))
      data.table(
        meta_scope = scope_name, n_cohorts = .N, total_n = sum(n), total_events = sum(events),
        log_hazard_ratio = as.numeric(fit$b), se = fit$se,
        hazard_ratio = exp(as.numeric(fit$b)), ci_low_hr = exp(fit$ci.lb),
        ci_high_hr = exp(fit$ci.ub), p = fit$pval, tau2 = fit$tau2,
        I2 = fit$I2, Q_p = fit$QEp, cohorts = paste(cohort, collapse = ";")
      )
    }
  }, by = .(feature_id, candidate, feature_type)]
}

meta_survival <- rbindlist(list(
  meta_survival_scope(survival_results[survival_endpoint == "OS"], "OS_only"),
  meta_survival_scope(survival_results, "OS_PFS_sensitivity")
), use.names = TRUE, fill = TRUE)
meta_survival[, fdr := p.adjust(p, method = "BH"), by = meta_scope]
meta_survival[, direction := fifelse(hazard_ratio > 1, "higher_is_worse", "higher_is_better")]
setorder(meta_survival, candidate, feature_type, feature_id, meta_scope)
atomic_fwrite(meta_survival, file.path(DIRS$tables, "Table_S29_ICB_guardrail_survival_meta_analysis.tsv"))

candidate_summary <- data.table(candidate = candidates)
feature_counts <- features[, .(n_tested_features = uniqueN(feature_id)), by = candidate]
candidate_summary[feature_counts, on = "candidate", n_tested_features := i.n_tested_features]
candidate_summary[, `:=`(
  response_significant_harm = vapply(candidate, function(z) {
    any(meta_response[candidate == z, estimate < 0 & ci_high < 0 & fdr <= 0.10])
  }, logical(1)),
  os_significant_harm = vapply(candidate, function(z) {
    any(meta_survival[candidate == z & meta_scope == "OS_only",
                      hazard_ratio > 1 & ci_low_hr > 1 & fdr <= 0.10])
  }, logical(1)),
  pfs_significant_harm = vapply(candidate, function(z) {
    any(survival_results[candidate == z & survival_endpoint == "PFS",
                         hazard_ratio > 1 & ci_low_hr > 1 & fdr_within_cohort_endpoint <= 0.10])
  }, logical(1)),
  response_worst_meta_g = vapply(candidate, function(z) {
    x <- meta_response[candidate == z]$estimate
    if (length(x)) min(x, na.rm = TRUE) else NA_real_
  }, numeric(1)),
  os_worst_meta_hr = vapply(candidate, function(z) {
    x <- meta_survival[candidate == z & meta_scope == "OS_only"]$hazard_ratio
    if (length(x)) max(x, na.rm = TRUE) else NA_real_
  }, numeric(1)),
  max_nominal_harmful_response_cohorts = vapply(candidate, function(z) {
    x <- response_results[candidate == z & hedges_g < 0 & p_smd < 0.10, .N, by = feature_id]$N
    if (length(x)) max(x) else 0L
  }, integer(1))
)]
candidate_summary[, guardrail_status := fcase(
  response_significant_harm | os_significant_harm | pfs_significant_harm,
  "FAIL_detected_tumor_efficacy_harm",
  response_worst_meta_g < -0.20 | os_worst_meta_hr > 1.20 |
    max_nominal_harmful_response_cohorts >= 2L,
  "CONCERN_systemic_use_not_supported",
  default = "NO_DETECTED_HARM_not_equivalent_to_safety"
)]
candidate_summary[, translation_constraint := fcase(
  guardrail_status == "FAIL_detected_tumor_efficacy_harm",
  "reject systemic translation; retain only as mechanistic clue if biologically useful",
  guardrail_status == "CONCERN_systemic_use_not_supported",
  "consider only epithelial-local or downstream tissue-restricted strategies pending stronger validation",
  default = "requires orthogonal tumor perturbation and safety evidence before any systemic claim"
)]
setorder(candidate_summary, guardrail_status, candidate)
atomic_fwrite(candidate_summary, file.path(DIRS$tables, "Table_S30_candidate_tumor_efficacy_guardrail_summary.tsv"))

# Resolve the systemic TWEAK/Fn14 concern by testing the 66-gene epithelial
# rescue leading edge one gene at a time. This is explicitly a downstream
# prioritization analysis, not evidence that absence of an association proves
# oncologic safety.
cat("Ranking TNFSF12 downstream epithelial rescue genes...\n")
stopifnot(exists("rescue_genes"), length(rescue_genes) >= 20L)
downstream_features <- data.table(
  feature_id = paste0("DOWNSTREAM_GENE_", rescue_genes),
  candidate = "TNFSF12_downstream",
  feature_type = "downstream_rescue_gene",
  celltype = "Transit_Amplifying",
  up_genes = lapply(rescue_genes, function(z) z),
  down_genes = lapply(rescue_genes, function(z) character()),
  derivation = "single gene from the positive leading edge of the frozen epithelial function-loss module after TNFSF12 perturbation"
)
downstream_definition <- copy(downstream_features)
downstream_definition[, `:=`(
  gene = rescue_genes,
  n_up_genes = 1L,
  n_down_genes = 0L,
  up_genes = vapply(up_genes, paste, collapse = ";", character(1)),
  down_genes = ""
)]
setcolorder(downstream_definition, c(
  "feature_id", "gene", "candidate", "feature_type", "celltype",
  "n_up_genes", "n_down_genes", "up_genes", "down_genes", "derivation"
))
atomic_fwrite(
  downstream_definition,
  file.path(DIRS$tables, "Table_S38_TNFSF12_downstream_feature_definitions.tsv")
)

downstream_scores <- rbindlist(
  lapply(cohorts, score_cohort, features = downstream_features),
  use.names = TRUE, fill = TRUE
)
atomic_fwrite(
  downstream_scores,
  file.path(DIRS$results, "TNFSF12_downstream_ICB_patient_scores.tsv.gz")
)

downstream_response <- downstream_scores[, response_effect_one(.SD),
  by = .(cohort, cancer, feature_id, candidate, feature_type)]
downstream_response[, gene := sub("^DOWNSTREAM_GENE_", "", feature_id)]
downstream_response[, fdr_smd_within_cohort := p.adjust(p_smd, method = "BH"), by = cohort]
setorder(downstream_response, gene, cohort)
atomic_fwrite(
  downstream_response,
  file.path(DIRS$tables, "Table_S39_TNFSF12_downstream_ICB_response_associations.tsv")
)

downstream_response_meta <- downstream_response[, {
  if (.N < 2L) NULL else {
    fit <- do.call(rma.uni, list(
      yi = as.numeric(.SD[["hedges_g"]]),
      vi = as.numeric(.SD[["var_hedges_g"]]), method = "REML"
    ))
    data.table(
      n_cohorts = .N, total_n = sum(n), estimate = as.numeric(fit$b),
      se = fit$se, ci_low = fit$ci.lb, ci_high = fit$ci.ub, p = fit$pval,
      tau2 = fit$tau2, I2 = fit$I2, Q_p = fit$QEp,
      cohorts = paste(cohort, collapse = ";")
    )
  }
}, by = .(feature_id, candidate, feature_type, gene)]
downstream_response_meta[, fdr := p.adjust(p, method = "BH")]
downstream_response_meta[, direction := fifelse(
  estimate > 0, "higher_in_responders", "higher_in_nonresponders"
)]
setorder(downstream_response_meta, p)
atomic_fwrite(
  downstream_response_meta,
  file.path(DIRS$tables, "Table_S40_TNFSF12_downstream_ICB_response_meta_analysis.tsv")
)

downstream_survival <- downstream_scores[!is.na(survival_endpoint), survival_effect_one(.SD),
  by = .(cohort, cancer, survival_endpoint, feature_id, candidate, feature_type)]
downstream_survival[, gene := sub("^DOWNSTREAM_GENE_", "", feature_id)]
downstream_survival[, fdr_within_cohort_endpoint := p.adjust(p, method = "BH"),
                    by = .(cohort, survival_endpoint)]
setorder(downstream_survival, gene, cohort)
atomic_fwrite(
  downstream_survival,
  file.path(DIRS$tables, "Table_S41_TNFSF12_downstream_ICB_survival_associations.tsv")
)

downstream_survival_meta <- rbindlist(list(
  meta_survival_scope(downstream_survival[survival_endpoint == "OS"], "OS_only"),
  meta_survival_scope(downstream_survival, "OS_PFS_sensitivity")
), use.names = TRUE, fill = TRUE)
downstream_survival_meta[, gene := sub("^DOWNSTREAM_GENE_", "", feature_id)]
downstream_survival_meta[, fdr := p.adjust(p, method = "BH"), by = meta_scope]
downstream_survival_meta[, direction := fifelse(
  hazard_ratio > 1, "higher_is_worse", "higher_is_better"
)]
setorder(downstream_survival_meta, meta_scope, p)
atomic_fwrite(
  downstream_survival_meta,
  file.path(DIRS$tables, "Table_S42_TNFSF12_downstream_ICB_survival_meta_analysis.tsv")
)

organoid_gene <- org_de[
  perturbation == "TNFSF12" & celltype == "Transit_Amplifying" & gene %chin% rescue_genes,
  .(gene, organoid_log2FC = log2FoldChange, organoid_lfcSE = lfcSE,
    organoid_p = pvalue, organoid_fdr = padj, organoid_rank_z = rank_z)
]
independent_gene <- fread(file.path(DIRS$results, "independent_colon_patient_level_DE_all.tsv.gz"))
independent_gene[, symbol := sub("^.*\\|", "", gene)]
independent_gene <- independent_gene[
  feature_level == "family" & contrast == "irColitis_vs_ICI_control_PD1_only" &
    celltype %chin% c(
      "Absorptive epithelial", "Immature epithelial",
      "Mature absorptive epithelial", "Secretory epithelial"
    ) & symbol %chin% rescue_genes,
  .(
    independent_n_families = .N,
    independent_n_down = sum(logFC < 0, na.rm = TRUE),
    independent_median_logFC = median(logFC, na.rm = TRUE),
    independent_min_p = min(PValue, na.rm = TRUE),
    independent_min_fdr = min(FDR, na.rm = TRUE),
    independent_family_logFC = paste(celltype, sprintf("%.3f", logFC), sep = ":", collapse = ";")
  ), by = .(gene = symbol)
]
spatial_gene <- fread(file.path(DIRS$results, "spatial_colon_patient_level_DE.tsv.gz"))[
  gene %chin% rescue_genes,
  .(gene, spatial_logFC = logFC, spatial_p = PValue, spatial_fdr = FDR)
]
response_gene <- downstream_response_meta[, .(
  gene, response_meta_g = estimate, response_ci_low = ci_low,
  response_ci_high = ci_high, response_p = p, response_fdr = fdr
)]
os_gene <- downstream_survival_meta[meta_scope == "OS_only", .(
  gene, os_meta_hr = hazard_ratio, os_ci_low = ci_low_hr,
  os_ci_high = ci_high_hr, os_p = p, os_fdr = fdr
)]

downstream_rank <- Reduce(
  function(x, y) merge(x, y, by = "gene", all.x = TRUE, sort = FALSE),
  list(data.table(gene = rescue_genes), organoid_gene, independent_gene,
       spatial_gene, response_gene, os_gene)
)
downstream_rank[, detected_response_harm :=
  is.finite(response_ci_high) & response_ci_high < 0]
downstream_rank[, detected_os_harm :=
  is.finite(os_ci_low) & os_ci_low > 1]
downstream_rank[, no_detected_icb_harm := !(detected_response_harm | detected_os_harm)]
downstream_rank[, evidence_tier := fcase(
  organoid_log2FC > 0 & organoid_fdr <= 0.10 &
    independent_n_down >= 3L & independent_median_logFC < 0 &
    spatial_logFC < 0 & spatial_fdr <= 0.10 & no_detected_icb_harm,
  "Tier_A_convergent_downstream_candidate",
  organoid_log2FC > 0 & organoid_p <= 0.05 &
    independent_n_down >= 3L & independent_median_logFC < 0 &
    spatial_logFC < 0 & spatial_fdr <= 0.10 & no_detected_icb_harm,
  "Tier_B_supportive_downstream_candidate",
  default = "not_prioritized"
)]
downstream_rank[, ranking_score :=
  3 * as.integer(is.finite(organoid_fdr) & organoid_fdr <= 0.10) +
  as.integer(is.finite(organoid_p) & organoid_p <= 0.05) +
  pmin(fifelse(is.na(independent_n_down), 0, independent_n_down), 4L) / 2 +
  2 * as.integer(is.finite(spatial_fdr) & spatial_logFC < 0 & spatial_fdr <= 0.10) +
  as.integer(is.finite(response_meta_g) & response_meta_g >= 0) +
  as.integer(is.finite(os_meta_hr) & os_meta_hr <= 1) -
  4 * as.integer(detected_response_harm | detected_os_harm)
]
downstream_rank[, interpretation_constraint := fifelse(
  no_detected_icb_harm,
  "no detected efficacy harm in available bulk cohorts is not proof of safety; requires causal and tumor-context validation",
  "deprioritized because a 95% confidence interval indicates possible tumor-efficacy harm"
)]
downstream_rank[, evidence_tier_order := fcase(
  evidence_tier == "Tier_A_convergent_downstream_candidate", 1L,
  evidence_tier == "Tier_B_supportive_downstream_candidate", 2L,
  default = 3L
)]
setorder(downstream_rank, evidence_tier_order, -ranking_score, organoid_fdr)
atomic_fwrite(
  downstream_rank,
  file.path(DIRS$tables, "Table_S43_TNFSF12_downstream_multicohort_ranking.tsv")
)

cat("ICI efficacy guardrail complete.\n")
print(cohort_audit)
print(candidate_summary)
print(downstream_rank[evidence_tier != "not_prioritized"])
