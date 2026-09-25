#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(cowplot)
  library(ggrepel)
  library(scales)
})


# Recorded display annotations are external scientific metadata, never code constants.
archived_label <- function(key) {
 p<-Sys.getenv('ICI_HISTORICAL_DISPLAY_LABELS',unset='')
 if(!nzchar(p)||!file.exists(p))stop('BLOCKED_INPUT: ICI_HISTORICAL_DISPLAY_LABELS')
 labels<-jsonlite::fromJSON(p)
 if(is.null(labels[[key]])||length(labels[[key]])!=1)stop(paste('Missing archived label',key))
 as.character(labels[[key]])
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) stop("Usage: make_manuscript_figures.R <project_dir>")
# Keep the ASCII junction path supplied by the runner. Resolving the junction to
# the Chinese target path breaks some Windows graphics devices under R 4.2.
root <- gsub("\\\\", "/", args[[1]])
if (!dir.exists(root)) stop("Project directory does not exist: ", root)
tabdir <- file.path(root, "tables_for_article")
figdir <- file.path(root, "figures_for_article")

read_tab <- function(name) read.delim(file.path(tabdir, name), check.names = FALSE, stringsAsFactors = FALSE)
clean_pathway <- function(x) gsub("_", " ", gsub("^HALLMARK_", "", x))
short_module <- function(x) recode(x,
  CURATED_IFN_VISIBILITY = "IFN visibility",
  ICI_COLITIS_EPITHELIAL_UP_STRICT = "Disease-up",
  ICI_COLITIS_EPITHELIAL_UP_NON_IFN = "Non-IFN residual",
  ICI_COLITIS_EPITHELIAL_LOSS_STRICT = "Function-loss",
  CURATED_COLON_METABOLIC_FUNCTION = "Metabolic function",
  CURATED_EPITHELIAL_BARRIER = "Generic barrier",
  TNFSF12_FUNCTION_RESCUE_CORE = "TNFSF12 rescue core",
  TNFSF12_ORGANOID_DOWN = "TNFSF12-down genes",
  TNFSF12_ORGANOID_UP = "TNFSF12-up genes",
  .default = x
)

ink <- "#17212B"
blue <- "#2474A6"
red <- "#C6453D"
gold <- "#D99922"
green <- "#2C8C69"
purple <- "#7253A6"
grey <- "#87939D"

theme_pub <- theme_classic(base_size = 9, base_family = "sans") +
  theme(
    text = element_text(colour = ink),
    axis.text = element_text(colour = ink),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 10, hjust = 0),
    plot.subtitle = element_text(size = 8.5, colour = "#4B5964"),
    strip.background = element_rect(fill = "#EEF2F4", colour = NA),
    strip.text = element_text(face = "bold", colour = ink),
    legend.title = element_text(face = "bold"),
    legend.key.height = unit(3.5, "mm"),
    plot.margin = margin(5.5, 7, 5.5, 7)
  )

save_figure <- function(plot, stem, width, height) {
  ggsave(file.path(figdir, paste0(stem, ".pdf")), plot, width = width, height = height,
         units = "in", device = cairo_pdf, bg = "white")
  ggsave(file.path(figdir, paste0(stem, ".tiff")), plot, width = width, height = height,
         units = "in", dpi = 600, compression = "lzw", bg = "white", device = grDevices::tiff)
  ggsave(file.path(figdir, paste0(stem, "_preview.png")), plot, width = width, height = height,
         units = "in", dpi = 180, bg = "white", device = grDevices::png)
}

panel_label <- function(...) plot_grid(..., labels = LETTERS, label_fontfamily = "sans",
                                        label_fontface = "bold", label_size = 12, hjust = -0.2)

# Figure 1: evidence architecture --------------------------------------------------
nodes <- data.frame(
  x = c(1, 3, 5, 7, 9), y = 1,
  label = c(
    "Human irAE discovery\nHeart + colon\npatient pseudobulk",
    "Independent transfer\nColon sc/snRNA + spatial\nsample-level validation",
    archived_label('organoid_overview'),
    "Preservation guardrails\n5 tumour-ICI cohorts\nresponse + survival",
    "Evidence-gated decision\nReject systemic ligands\nretain local repair core"
  ),
  fill = c("#DDECF4", "#E7F1EC", "#F8ECD5", "#F2E8F4", "#F3E2E0")
)
fig1 <- ggplot() +
  geom_segment(data = data.frame(x = c(1.82, 3.82, 5.82, 7.82), xend = c(2.18, 4.18, 6.18, 8.18)),
               aes(x = x, xend = xend, y = 1, yend = 1), linewidth = 1, colour = "#65737E",
               arrow = arrow(length = unit(2.5, "mm"))) +
  geom_tile(data = nodes, aes(x, y, fill = fill), width = 1.55, height = 0.72,
            colour = "white", linewidth = 1) +
  geom_text(data = nodes, aes(x, y, label = label), size = 2.65, lineheight = 1.05,
            fontface = "bold", colour = ink) +
  scale_fill_identity() +
  annotate("text", x = 1, y = 0.35, label = "TCR repertoire\nexploratory support", size = 3, colour = blue) +
  annotate("segment", x = 1, xend = 1, y = 0.62, yend = 0.47, colour = blue,
           arrow = arrow(length = unit(2, "mm"))) +
  annotate("text", x = 3, y = 0.35, label = "Frozen modules\n(no refitting)", size = 3, colour = green) +
  annotate("segment", x = 3, xend = 3, y = 0.62, yend = 0.47, colour = green,
           arrow = arrow(length = unit(2, "mm"))) +
  annotate("text", x = 7, y = 0.35, label = "Off-target Hallmarks\n+ efficacy preservation", size = 3, colour = purple) +
  annotate("segment", x = 7, xend = 7, y = 0.62, yend = 0.47, colour = purple,
           arrow = arrow(length = unit(2, "mm"))) +
  annotate("label", x = 5, y = 1.63,
           label = "Question: can tissue repair be separated from shared immune attack without weakening antitumour immunity?",
           size = 3.5, fontface = "bold", fill = "white", colour = ink, linewidth = 0) +
  annotate("text", x = 5, y = -0.02,
           label = "Predefined boundaries: patient-level inference | independent validation | real perturbation\nexplicit rejection rules and preserved uncertainty",
           size = 2.9, colour = "#4B5964") +
  coord_cartesian(xlim = c(0.25, 9.75), ylim = c(-0.2, 1.85), clip = "off") +
  theme_void(base_family = "sans") +
  theme(plot.margin = margin(12, 12, 12, 12))
save_figure(fig1, "Figure_1_study_design_and_evidence_gates", 11.0, 3.4)

# Figure 2: shared attack, divergent receivers -----------------------------------
s52 <- read_tab("Table_S52_cross_organ_receiver_Hallmark_matrix.tsv")
s53 <- read_tab("Table_S53_cross_organ_receiver_program_summary.tsv")
s58 <- read_tab("Table_S58_lung_source_confounded_module_summary.tsv")
focus_paths <- s53 %>% filter(cross_organ_class %in% c("shared_up_receiver_program", "organ_divergent_receiver_program")) %>% pull(pathway)
heat <- s52 %>%
  filter(pathway %in% focus_paths) %>%
  mutate(compartment = paste(ifelse(receiver_organ == "heart", "Heart", "Colon"), celltype, sep = " | "),
         pathway_label = clean_pathway(pathway), sig = ifelse(padj < 0.05, "*", ""))
path_order <- s53 %>% filter(pathway %in% focus_paths) %>% arrange(class_order, desc(heart_median_NES + colon_median_NES)) %>% pull(pathway) %>% clean_pathway()
heat$pathway_label <- factor(heat$pathway_label, levels = rev(unique(path_order)))
comp_order <- c("Heart | Endothelial cells", "Heart | Fibroblasts", "Heart | Mural cells",
                "Colon | Endothelial", "Colon | Epithelial", "Colon | Mesenchymal Stromal")
heat$compartment <- factor(heat$compartment, levels = comp_order)
p2a <- ggplot(heat, aes(compartment, pathway_label, fill = NES)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = sig), size = 3.2, colour = ink) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0,
                       limits = c(-3.3, 3.3), oob = squish, name = "NES") +
  labs(title = "Shared attack with divergent receiver remodeling",
       x = NULL, y = NULL, subtitle = "* = FDR < 0.05 within the tested compartment") +
  theme_pub + theme(axis.text.x = element_text(angle = 42, hjust = 1), legend.position = "right",
                    plot.margin = margin(5.5, 7, 5.5, 24))

scatter <- s53 %>% mutate(
  class = recode(cross_organ_class,
    shared_up_receiver_program = "Shared up",
    organ_divergent_receiver_program = "Organ-divergent",
    .default = "Other"),
  label = ifelse(class != "Other", clean_pathway(pathway), "")
)
p2b <- ggplot(scatter, aes(heart_median_NES, colon_median_NES, colour = class)) +
  geom_hline(yintercept = 0, colour = "#C8D0D5", linewidth = 0.4) +
  geom_vline(xintercept = 0, colour = "#C8D0D5", linewidth = 0.4) +
  geom_point(alpha = 0.8, size = 2.0) +
  geom_text_repel(data = subset(scatter, class != "Other"), aes(label = label), size = 2.3,
                  max.overlaps = Inf, min.segment.length = 0, box.padding = 0.25, show.legend = FALSE) +
  scale_colour_manual(values = c("Shared up" = red, "Organ-divergent" = blue, "Other" = grey)) +
  labs(title = "Whole-Hallmark concordance is absent", subtitle = archived_label('receiver_correlation'),
       x = "Heart median NES", y = "Colon median NES", colour = NULL) +
  theme_pub + theme(legend.position = "top")

lung <- s58 %>% filter(pathway %in% c("CURATED_IFN_VISIBILITY", "ICI_COLITIS_EPITHELIAL_UP_STRICT",
                                     "ICI_COLITIS_EPITHELIAL_UP_NON_IFN", "ICI_COLITIS_EPITHELIAL_LOSS_STRICT")) %>%
  mutate(module = short_module(pathway), module = factor(module, levels = rev(c("IFN visibility", "Disease-up", "Non-IFN residual", "Function-loss"))))
p2c <- ggplot(lung, aes(median_NES, module, fill = n_concordant_fdr05)) +
  geom_col(width = 0.62) +
  geom_text(aes(label = paste0(n_concordant_fdr05, "/", n_eligible_celltypes)), hjust = -0.15, size = 2.8) +
  scale_fill_gradient(low = "#D8E6EE", high = blue, limits = c(0, 4), name = "FDR-concordant\ncell types") +
  coord_cartesian(xlim = c(-0.1, 2.35)) +
  labs(title = "Lung reproduces shared attack only", subtitle = "Sensitivity analysis: disease and source are confounded",
       x = "Median NES (CIP vs external healthy)", y = NULL) + theme_pub

fig2 <- plot_grid(p2a, plot_grid(p2b, p2c, ncol = 1, rel_heights = c(1.15, 0.85)),
                  labels = c("A", "B"), label_fontface = "bold", label_size = 12,
                  rel_widths = c(1.65, 1))
save_figure(fig2, "Figure_2_shared_attack_and_organ_divergent_receiver_programs", 11.5, 7.6)

# Figure 3: frozen colitis module validation -------------------------------------
s21 <- read_tab("Table_S21_frozen_module_definitions.tsv")
s22 <- read_tab("Table_S22_frozen_module_validation.tsv")
s56 <- read_tab("Table_S56_patient_level_frozen_module_classification_performance.tsv")
key_modules <- c("CURATED_IFN_VISIBILITY", "ICI_COLITIS_EPITHELIAL_UP_STRICT",
                 "ICI_COLITIS_EPITHELIAL_UP_NON_IFN", "ICI_COLITIS_EPITHELIAL_LOSS_STRICT")
val <- s22 %>% filter(module %in% key_modules, cohort_role %in% c("independent_primary", "independent_spatial_validation")) %>%
  mutate(oriented_NES = ifelse(expected_direction_in_irAE == "down", -NES, NES),
         validation = ifelse(validation_dataset == "GSE210037", "Spatial sample pseudobulk", celltype),
         module_label = short_module(module), sig = ifelse(padj < 0.05 & direction_concordant, "*", ""))
val$module_label <- factor(val$module_label, levels = c("IFN visibility", "Disease-up", "Non-IFN residual", "Function-loss"))
val$validation <- factor(val$validation, levels = rev(c("Absorptive epithelial", "Immature epithelial",
  "Mature absorptive epithelial", "Secretory epithelial", "Spatial sample pseudobulk")))
p3a <- ggplot(val, aes(module_label, validation, fill = oriented_NES)) +
  geom_tile(colour = "white", linewidth = 0.35) + geom_text(aes(label = sig), size = 4) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, limits = c(-3.6, 3.6), oob = squish,
                       name = "Direction-oriented\nNES") +
  labs(title = "Frozen modules transfer without refitting",
       subtitle = "Positive = prespecified irAE direction; * = concordant FDR < 0.05", x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 25, hjust = 1))

auc <- s56 %>% filter(module %in% key_modules, dataset %in% c("GSE206300", "GSE210037")) %>%
  mutate(module_label = short_module(module), validation = ifelse(dataset == "GSE210037", "Spatial sample pseudobulk", celltype))
p3b <- ggplot(auc, aes(auc_disease_higher, validation, colour = module_label)) +
  geom_vline(xintercept = 0.5, linetype = 2, colour = grey) +
  geom_errorbarh(aes(xmin = bootstrap_auc_ci_low, xmax = bootstrap_auc_ci_high), height = 0.15, alpha = 0.7) +
  geom_point(size = 2.2) +
  facet_wrap(~module_label, ncol = 2) +
  scale_colour_manual(values = c("IFN visibility" = red, "Disease-up" = gold,
                                 "Non-IFN residual" = purple, "Function-loss" = blue), guide = "none") +
  coord_cartesian(xlim = c(0.25, 1.02)) +
  labs(title = "Patient-level separation is reproducible",
       subtitle = "Repeatability estimate, not a diagnostic test",
       x = "AUC (irAE higher after direction orientation)", y = NULL) + theme_pub +
  theme(axis.text.y = element_text(size = 7.2))

sizes <- s21 %>% filter(module %in% key_modules) %>% distinct(module, gene) %>% count(module, name = "n_genes") %>%
  mutate(module_label = short_module(module), module_label = factor(module_label,
    levels = rev(c("IFN visibility", "Disease-up", "Non-IFN residual", "Function-loss"))))
p3c <- ggplot(sizes, aes(n_genes, module_label, fill = module_label)) +
  geom_col(width = 0.65) + geom_text(aes(label = n_genes), hjust = -0.15, size = 3) +
  scale_fill_manual(values = c("IFN visibility" = red, "Disease-up" = gold,
                               "Non-IFN residual" = purple, "Function-loss" = blue), guide = "none") +
  coord_cartesian(xlim = c(0, 315)) + labs(title = "Frozen module sizes", x = "Genes", y = NULL) + theme_pub

fig3 <- plot_grid(p3a, plot_grid(p3b, p3c, ncol = 1, rel_heights = c(1.55, 0.7)),
                  labels = c("A", "B"), label_fontface = "bold", label_size = 12,
                  rel_widths = c(1.1, 1.35))
save_figure(fig3, "Figure_3_independent_validation_of_frozen_colitis_modules", 11.5, 7.2)

# Figure 4: human organoid perturbation and de-risking ----------------------------
s23 <- read_tab("Table_S23_GSE313368_irAE_rescue_screen.tsv")
s47 <- read_tab("Table_S47_organoid_candidate_offtarget_Hallmark_GSEA.tsv")
s50 <- read_tab("Table_S50_TNFSF12_organoid_signature_reciprocity_summary.tsv")
candidate_names <- c("TNFSF12", "INHBA", "BMP2", "BMP4", "TGFB1", "TGFB2", "TGFB3")
screen_plot <- s23 %>% filter(!is.na(ici_loss_median), !is.na(ici_up_non_ifn_median)) %>%
  mutate(label = ifelse(perturbation %in% candidate_names, perturbation, ""),
         class = case_when(perturbation == "TNFSF12" ~ "TNFSF12", perturbation %in% candidate_names ~ "Other function-rescue ligands", TRUE ~ "Other perturbations"))
p4a <- ggplot(screen_plot, aes(ici_loss_median, -ici_up_non_ifn_median, colour = class)) +
  geom_hline(yintercept = 0, colour = "#C8D0D5") + geom_vline(xintercept = 0, colour = "#C8D0D5") +
  geom_point(size = 2.1, alpha = 0.8) +
  geom_text_repel(aes(label = label), size = 2.7, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = c("TNFSF12" = red, "Other function-rescue ligands" = gold, "Other perturbations" = grey)) +
  labs(title = "Organoid screen nominates apparent repair signals",
       subtitle = "Upper-right: function restoration with suppression of the non-IFN disease residual",
       x = "Function-loss module rescue (median NES)", y = "Disease-residual suppression (-median NES)", colour = NULL) +
  theme_pub + theme(legend.position = "top")

tn_off <- s47 %>% filter(perturbation == "TNFSF12", padj < 0.05) %>%
  arrange(desc(abs_NES)) %>% slice_head(n = 14) %>%
  mutate(pathway_label = clean_pathway(pathway), pathway_label = factor(pathway_label, levels = rev(unique(pathway_label))))
p4b <- ggplot(tn_off, aes(NES, pathway_label, fill = NES > 0)) +
  geom_col(width = 0.7) + geom_vline(xintercept = 0, colour = ink, linewidth = 0.35) +
  scale_fill_manual(values = c(`TRUE` = red, `FALSE` = blue), labels = c("Suppressed", "Activated"), name = NULL) +
  labs(title = "TNFSF12 fails the expanded audit",
       subtitle = "TNFSF12, significant pathways (FDR < 0.05)", x = "NES", y = NULL) +
  theme_pub + theme(legend.position = "top", axis.text.y = element_text(size = 7.3))

recip <- s50 %>% mutate(signature = short_module(pathway),
  cohort = recode(cohort_role, independent_primary = "Independent sc/snRNA",
                  independent_spatial_platform = "Spatial sample pseudobulk",
                  derivation_cohort_reciprocity_check = "Discovery reciprocity"),
  mark = ifelse(n_reciprocal_fdr05 > 0, paste0(n_reciprocal_fdr05, "/", n_tests), ""))
p4c <- ggplot(recip, aes(cohort, signature, fill = median_NES)) +
  geom_tile(colour = "white", linewidth = 0.35) + geom_text(aes(label = mark), size = 2.8) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, name = "Median NES") +
  labs(title = archived_label('reciprocity_title'),
       subtitle = "Text: FDR-concordant tests / tests", x = NULL, y = NULL) + theme_pub +
  theme(axis.text.x = element_text(angle = 28, hjust = 1))

fig4 <- plot_grid(p4a, p4b, p4c, labels = c("A", "B", "C"), label_fontface = "bold",
                  label_size = 12, ncol = 3, rel_widths = c(1.15, 1.1, 0.95))
save_figure(fig4, "Figure_4_human_organoid_screen_and_TNFSF12_derisking", 13.5, 5.8)

# Figure 5: TCR evidence ----------------------------------------------------------
s31 <- read_tab("Table_S31_colon_TCR_patient_metrics.tsv")
s35 <- read_tab("Table_S35_TCR_statistical_summary.tsv")
s36 <- read_tab("Table_S36_colon_TCR_epithelial_module_patient_coupling.tsv")
tcr_metrics <- c("expanded_cell_fraction_ge2", "normalized_clonality", "cross_cd4_cd8_cell_fraction")
tcr_labels <- c(expanded_cell_fraction_ge2 = "Expanded-cell fraction",
                normalized_clonality = "Normalized clonality",
                cross_cd4_cd8_cell_fraction = "CD4-CD8 shared-clone cell fraction")
tcr_long <- s31 %>% filter(disease %in% c("CPI_colitis", "HC")) %>%
  select(disease, patient_id, all_of(tcr_metrics)) %>% pivot_longer(all_of(tcr_metrics), names_to = "metric", values_to = "value") %>%
  mutate(metric_label = tcr_labels[metric], disease = recode(disease, CPI_colitis = "ICI colitis", HC = "Healthy"))
tcr_stats <- s35 %>% filter(analysis == "colon_CPI_colitis_vs_HC_patient_level", metric %in% tcr_metrics) %>%
  mutate(metric_label = tcr_labels[metric], label = paste0("P=", signif(p, 2), "; q=", signif(fdr_within_analysis, 2)))
p5a <- ggplot(tcr_long, aes(disease, value, fill = disease)) +
  geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.45) +
  geom_point(position = position_jitter(width = 0.08, height = 0), size = 2, shape = 21, colour = ink) +
  facet_wrap(~metric_label, scales = "free_y", nrow = 1) +
  geom_text(data = tcr_stats, aes(x = 1.5, y = Inf, label = label), inherit.aes = FALSE,
            vjust = 1.5, size = 2.7) +
  scale_fill_manual(values = c("ICI colitis" = red, "Healthy" = "#B9C2C8"), guide = "none") +
  labs(title = "Patient-level TCR expansion and cross-lineage sharing accompany ICI colitis",
       subtitle = archived_label('tcr_sample_coverage'), x = NULL, y = NULL) + theme_pub +
  theme(axis.text.x = element_text(angle = 18, hjust = 1))

couple <- s36 %>% filter(disease == "CPI_colitis")
p5b <- ggplot(couple, aes(cross_cd4_cd8_cell_fraction, ICI_COLITIS_EPITHELIAL_LOSS_STRICT)) +
  geom_smooth(method = "lm", se = TRUE, colour = grey, fill = "#DDE4E8", linewidth = 0.7) +
  geom_point(size = 2.6, colour = blue) + geom_text_repel(aes(label = patient_id), size = 2.4, max.overlaps = Inf) +
  labs(title = "TCR-epithelial coupling is hypothesis-generating",
       subtitle = archived_label('tcr_coupling_statistics'),
       x = "CD4-CD8 shared-clone cell fraction", y = "Epithelial function-loss module score") + theme_pub

state_stats <- s35 %>% filter(analysis == "colon_CPI_colitis_vs_HC_patient_level",
                              metric %in% c("rarefied_clonality_200", "rarefied_top10_fraction_200")) %>%
  mutate(metric_label = recode(metric, rarefied_clonality_200 = "Rarefied clonality",
                               rarefied_top10_fraction_200 = "Rarefied top-10 fraction"))
p5c <- ggplot(state_stats, aes(median_difference, metric_label)) +
  geom_vline(xintercept = 0, colour = grey, linetype = 2) +
  geom_segment(aes(x = 0, xend = median_difference, yend = metric_label), linewidth = 1, colour = gold) +
  geom_point(size = 3, colour = gold) +
  geom_text(aes(label = paste0("P=", signif(p, 2), "; q=", signif(fdr_within_analysis, 2))), hjust = -0.1, size = 2.8) +
  coord_cartesian(xlim = c(-0.005, max(state_stats$median_difference) * 1.75)) +
  labs(title = "Depth-controlled metrics are weaker", x = "Median difference: ICI colitis - healthy", y = NULL) + theme_pub

fig5 <- plot_grid(p5a, plot_grid(p5b, p5c, ncol = 1, rel_heights = c(1.2, 0.7)),
                  labels = c("A", "B"), label_fontface = "bold", label_size = 12,
                  rel_widths = c(1.65, 1))
save_figure(fig5, "Figure_5_patient_level_TCR_evidence", 11.8, 6.8)

# Figure 6: preservation guardrails and final decisions ---------------------------
s30 <- read_tab("Table_S30_candidate_tumor_efficacy_guardrail_summary.tsv")
s43 <- read_tab("Table_S43_TNFSF12_downstream_multicohort_ranking.tsv")
s59 <- read_tab("Table_S59_final_candidate_decision_matrix.tsv")
guard_plot <- s30 %>% mutate(status = case_when(grepl("FAIL", guardrail_status) ~ "Detected efficacy harm",
  grepl("CONCERN", guardrail_status) ~ "Concern", TRUE ~ "No detected harm is not safety"))
p6a <- ggplot(guard_plot, aes(response_worst_meta_g, log(os_worst_meta_hr), colour = status, label = candidate)) +
  geom_hline(yintercept = 0, colour = grey, linetype = 2) + geom_vline(xintercept = 0, colour = grey, linetype = 2) +
  geom_point(size = 3) + geom_text_repel(size = 2.7, max.overlaps = Inf) +
  scale_colour_manual(values = c("Detected efficacy harm" = red, "Concern" = gold, "No detected harm is not safety" = blue)) +
  labs(title = "Tumour-ICI guardrails constrain rescue ligands",
       x = "Worst response meta-effect (negative is adverse)", y = "log worst OS hazard ratio", colour = NULL) +
  theme_pub + theme(legend.position = "top")

ligand_dec <- s59 %>% filter(candidate_level == "systemic ligand perturbation") %>%
  transmute(candidate,
    `Function rescue` = ifelse(organoid_function_rescue == "yes", "pass", "fail"),
    `Legacy module gate` = ifelse(legacy_module_decoupling_gate == "yes", "pass", "fail"),
    `Expanded off-target gate` = expanded_offtarget_gate,
    `Tumour efficacy gate` = case_when(grepl("FAIL", tumor_efficacy_guardrail) ~ "fail",
                                      grepl("CONCERN", tumor_efficacy_guardrail) ~ "concern", TRUE ~ "uncertain")) %>%
  pivot_longer(-candidate, names_to = "gate", values_to = "status") %>%
  mutate(gate = factor(gate, levels = c("Function rescue", "Legacy module gate",
                                        "Expanded off-target gate", "Tumour efficacy gate")))
p6b <- ggplot(ligand_dec, aes(gate, candidate, fill = status)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_text(aes(label = case_when(status == "pass" ~ "PASS", status == "fail" ~ "FAIL", status == "concern" ~ "!", TRUE ~ "?")),
            size = 2.7, fontface = "bold") +
  scale_fill_manual(values = c(pass = "#B9DDCF", fail = "#E7B4AF", concern = "#F1D59C", uncertain = "#CBDCE8"), name = NULL) +
  labs(title = "No systemic ligand survives the expanded decision rules", x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 28, hjust = 1), legend.position = "top")

sorl1 <- s43 %>% filter(gene == "SORL1")
sorl1_tissue <- data.frame(
  evidence = c("Organoid induction", "Independent epithelial depletion", "Spatial sample depletion"),
  estimate = c(sorl1$organoid_log2FC, -sorl1$independent_median_logFC, -sorl1$spatial_logFC),
  note = c(paste0("FDR ", signif(sorl1$organoid_fdr, 2)),
           paste0(sorl1$independent_n_down, "/", sorl1$independent_n_families, " families"),
           paste0("FDR ", signif(sorl1$spatial_fdr, 2)))
)
p6c <- ggplot(sorl1_tissue, aes(estimate, reorder(evidence, estimate))) +
  geom_segment(aes(x = 0, xend = estimate, yend = reorder(evidence, estimate)), colour = purple, linewidth = 1) +
  geom_point(size = 3, colour = purple) + geom_text(aes(label = note), hjust = -0.12, size = 2.8) +
  coord_cartesian(xlim = c(0, max(sorl1_tissue$estimate) * 1.42)) +
  labs(title = "SORL1 is convergent but exploratory",
       subtitle = paste0("Response g ", sprintf("%.2f", sorl1$response_meta_g), " [", sprintf("%.2f", sorl1$response_ci_low), ", ", sprintf("%.2f", sorl1$response_ci_high),
                         "]\nOS HR ", sprintf("%.2f", sorl1$os_meta_hr), " [", sprintf("%.2f", sorl1$os_ci_low), ", ", sprintf("%.2f", sorl1$os_ci_high), "]"),
       x = "Direction-oriented effect", y = NULL) + theme_pub

fig6 <- plot_grid(p6a, p6b, p6c, labels = c("A", "B", "C"), label_fontface = "bold", label_size = 12,
                  ncol = 3, rel_widths = c(1.15, 1.05, 1.05))
save_figure(fig6, "Figure_6_tumor_efficacy_guardrails_and_final_decisions", 13.5, 5.7)

# Supplementary Figure 1: complete Hallmark receiver matrix -----------------------
full_heat <- s52 %>% mutate(compartment = paste(ifelse(receiver_organ == "heart", "Heart", "Colon"), celltype, sep = " | "),
                            pathway_label = clean_pathway(pathway), sig = ifelse(padj < 0.05, "*", ""))
full_heat$compartment <- factor(full_heat$compartment, levels = comp_order)
path_levels <- s53 %>% arrange(class_order, desc(abs(heart_median_NES) + abs(colon_median_NES))) %>% pull(pathway) %>% clean_pathway()
full_heat$pathway_label <- factor(full_heat$pathway_label, levels = rev(unique(path_levels)))
sf1 <- ggplot(full_heat, aes(compartment, pathway_label, fill = NES)) + geom_tile(colour = "white", linewidth = 0.15) +
  geom_text(aes(label = sig), size = 1.8) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, limits = c(-3.3, 3.3), oob = squish, name = "NES") +
  labs(title = "Supplementary Figure 1 | Complete receiver-compartment Hallmark landscape",
       subtitle = "* = FDR < 0.05", x = NULL, y = NULL) + theme_pub +
  theme(axis.text.x = element_text(angle = 40, hjust = 1), axis.text.y = element_text(size = 5.6))
save_figure(sf1, "Supplementary_Figure_1_complete_receiver_Hallmark_matrix", 8.5, 12.0)

# Supplementary Figure 2: all module validation -----------------------------------
all_val <- s22 %>% mutate(oriented_NES = ifelse(expected_direction_in_irAE == "down", -NES, NES),
                          module_label = short_module(module),
                          validation = paste(validation_dataset, celltype, sep = " | "), sig = ifelse(padj < 0.05 & direction_concordant, "*", ""))
sf2 <- ggplot(all_val, aes(module_label, validation, fill = oriented_NES)) + geom_tile(colour = "white", linewidth = 0.2) +
  geom_text(aes(label = sig), size = 2) + scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, name = "Oriented NES") +
  labs(title = "Supplementary Figure 2 | Complete frozen-module transfer results", subtitle = "* = concordant FDR < 0.05", x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6))
save_figure(sf2, "Supplementary_Figure_2_complete_frozen_module_validation", 9.5, 10.5)

# Supplementary Figure 3: source-confounded lung sensitivity ----------------------
sf3 <- ggplot(s58 %>% mutate(module = short_module(pathway)), aes(median_NES, reorder(module, median_NES), fill = n_concordant_fdr05)) +
  geom_col(width = 0.65) + geom_text(aes(label = paste0(n_concordant_fdr05, "/", n_eligible_celltypes)), hjust = -0.1, size = 3) +
  scale_fill_gradient(low = "#D8E6EE", high = blue, name = "FDR-concordant") +
  coord_cartesian(xlim = c(min(0, min(s58$median_NES)), max(s58$median_NES) * 1.25)) +
  labs(title = "Supplementary Figure 3 | Lung BALF module sensitivity",
       subtitle = "All CIP cases and healthy controls originate from different sources; no causal inference",
       x = "Median NES", y = NULL) + theme_pub
save_figure(sf3, "Supplementary_Figure_3_lung_source_confounded_sensitivity", 7.2, 4.8)

# Supplementary Figure 4: candidate off-target matrix -----------------------------
off_sig <- s47 %>% filter(padj < 0.05, perturbation %in% candidate_names) %>% group_by(perturbation, pathway) %>%
  summarise(NES = NES[which.max(abs(NES))], .groups = "drop") %>%
  mutate(pathway_label = clean_pathway(pathway))
sf4 <- ggplot(off_sig, aes(perturbation, pathway_label, fill = NES)) + geom_tile(colour = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = blue, mid = "white", high = red, midpoint = 0, name = "NES") +
  labs(title = "Supplementary Figure 4 | Significant organoid off-target pathways across rescue candidates",
       subtitle = "Maximum absolute NES across retained epithelial cell types; FDR < 0.05", x = NULL, y = NULL) +
  theme_pub + theme(axis.text.x = element_text(angle = 35, hjust = 1), axis.text.y = element_text(size = 6.5))
save_figure(sf4, "Supplementary_Figure_4_organoid_candidate_offtarget_matrix", 8.5, 9.5)

# Supplementary Figure 5: paired heart-tumour repertoire sensitivity --------------
s34 <- read_tab("Table_S34_heart_tumor_TCR_overlap.tsv")
sf5 <- s34 %>% filter(clone_definition == "nucleotide_exact") %>%
  ggplot(aes(compartment_b, morisita_horn_abundance, group = donor, colour = donor)) +
  geom_line(alpha = 0.55) + geom_point(size = 2.5) +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Supplementary Figure 5 | Paired heart-tumour/adjacent repertoire overlap",
       subtitle = "Four donors; no significant compartment contrast. Descriptive sensitivity only.",
       x = "Comparator compartment", y = "Morisita-Horn abundance overlap", colour = "Donor") + theme_pub
save_figure(sf5, "Supplementary_Figure_5_heart_tumor_TCR_overlap_sensitivity", 7.2, 4.8)

# Supplementary Figure 6: downstream repair-core ranking --------------------------
rankplot <- s43 %>%
  filter(!is.na(independent_median_logFC), !is.na(spatial_logFC)) %>%
  mutate(highlight = evidence_tier %in% c("Tier_A_convergent_downstream_candidate", "Tier_B_supportive_downstream_candidate"),
                           label = ifelse(highlight, gene, ""))
sf6 <- ggplot(rankplot, aes(-independent_median_logFC, -spatial_logFC, colour = highlight, size = pmax(organoid_log2FC, 0))) +
  geom_hline(yintercept = 0, colour = grey) + geom_vline(xintercept = 0, colour = grey) +
  geom_point(alpha = 0.75) + geom_text_repel(aes(label = label), size = 2.7, max.overlaps = Inf, show.legend = FALSE) +
  scale_colour_manual(values = c(`TRUE` = purple, `FALSE` = "#A8B1B7"), guide = "none") +
  scale_size_continuous(name = "Organoid induction\n(log2FC)", range = c(1.5, 5)) +
  labs(title = "Supplementary Figure 6 | TNFSF12 repair-core downstream ranking",
       subtitle = "Genes with both independent epithelial and spatial estimates; positive axes indicate depletion",
       x = "Independent epithelial depletion (-median logFC)", y = "Spatial sample depletion (-logFC)") + theme_pub
save_figure(sf6, "Supplementary_Figure_6_TNFSF12_downstream_ranking", 7.2, 5.4)

message("Created 6 main and 6 supplementary figures in ", figdir)
