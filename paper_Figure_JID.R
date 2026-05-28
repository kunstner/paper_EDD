# Credentials -------------------------------------------------------------

#
# Author: Axel Künstner
# Project: Ichtyosis
# Data: Microbiome
# Description: Restructured figures for manuscript submission

# Libraries ---------------------------------------------------------------

library(sf)
library(tidyverse)
library(phyloseq)
library(patchwork)
library(ggpubr)
library(flextable)
library(vegan)
library(breakaway)
library(microbiome)
library(Maaslin2)
library(WriteXLS)

source("code/helper.R")

# Data and variables ------------------------------------------------------

SEED <- 1253

ps <- readRDS(file = "data/phyloseq.OTU.RDS")

sample_data(ps) <- sample_data(ps) |>
    data.frame() |>
    dplyr::mutate(Control = case_when(Control == TRUE ~ 'Control',
                                      Control == FALSE ~ 'Patient'))

sample_data(ps) <- sample_data(ps) |>
    data.frame() |>
    dplyr::mutate(Age_group = case_when(Age_group == 'old' ~ 'older',
                                        .default = Age_group)) |>
    dplyr::mutate(Age_group = factor(Age_group,
                                     levels = c('young', 'middle', 'older')))

cova <- ps |> sample_data() |> data.frame()

# Plotting themes ---------------------------------------------------------

theme_alpha <- theme_classic() +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())

theme_def <- theme_classic() +
    theme(
        axis.line = element_line(color = "black"),
        legend.position = "none",
        axis.ticks = element_line(),
        plot.title = element_text(size = 18, face = "bold"),
        axis.title.x = element_text(size = 14),
        axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        strip.text = element_text(face = 'italic') )

# =============================================================================
# NEW FIGURE 1 PANELS
# =============================================================================

# Fig1a (original Fig1d) - plot_phylum -------------------------------------

plot_phylum <- ps |>
    tax_glom(taxrank = "Phylum") |>
    transform_sample_counts(function(x) {x/sum(x)} ) |>
    psmelt() |>
    arrange(Phylum) |>
    dplyr::select( Phylum, Location, Abundance, Age_group, Control) |>
    dplyr::group_by(Phylum, Location, Age_group, Control) |>
    dplyr::reframe(mean = mean(Abundance), groups = c(Phylum)) |>
    unique() |>
    dplyr::mutate(Location = gsub(pattern = "_", replacement = " ", x = Location)) |>
    dplyr::mutate(Location = factor(x = Location, levels = c('Arm left', 'Arm right', 'Abdomen'))) |>
    dplyr::mutate(Phylum = gsub(pattern = "_", replacement = " ", x = Phylum)) |>
    ggplot(aes(x = Age_group, y = mean, fill = Phylum)) +
    geom_bar(stat = "identity") +
    ggh4x::facet_nested_wrap(facets = ~Control + Location,
                             ncol = 4) +
    scale_fill_manual(values = cols) +
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance \n") +
    ggtitle("") +
    theme_classic() +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          axis.text = element_text(size=12),
          strip.text = element_text(size=12)) +
    ylim(0,1.01) + xlab('') +
    guides(fill = guide_legend(title = "Phylum", ncol = 1))
plot_phylum

# Fig1b (original Fig1b) - plot_genus -------------------------------------

# Genera
means_genus <- ps |>
    tax_glom(taxrank = "Genus") |> # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} ) |>
    psmelt() |> # Melt to long format
    dplyr::arrange(Genus) |>
    dplyr::select( Genus, Location, Abundance, Age_group, Control) |>
    dplyr::group_by(Genus, Location, Age_group, Control) |>
    dplyr::reframe(mean = mean(Abundance), groups = c(Genus)) |>
    unique() |>
    dplyr::mutate(Genus = gsub(pattern = "_", replacement = " ", x = Genus))
# filter(mean >= 0.02) # Filter out low abundance taxa

# summarise lowly abundant genera
means_genus$Genus[means_genus$mean <= 0.015] <- "Other"
means_genus <- means_genus |>
    dplyr::group_by(Genus, Location, Age_group, Control) |>
    reframe(mean = sum(mean), groups = c(Genus)) |>
    unique()

plot_genus <- means_genus |>
    dplyr::mutate(Location = gsub(pattern = "_", replacement = " ", x = Location)) |>
    dplyr::mutate(Location = factor(x = Location, levels = c('Arm left', 'Arm right', 'Abdomen'))) |>
    #dplyr::mutate(Genus = gsub(pattern = "f__", replacement = "", x = Genus)) |>
    #dplyr::mutate(Genus = gsub(pattern = "_unclassified", replacement = " uncl.", x = Genus)) |>
    ggplot(aes(x = Age_group, y = mean, fill = Genus)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = cols) +
    # facet_wrap(~Control + Location, scales = "free_x", ncol = 4) +
    ggh4x::facet_nested_wrap(facets = ~Control + Location,
                             ncol = 4) +
    # Remove x axis title
    theme(axis.title.x = element_blank()) +
    guides(fill = guide_legend(reverse = FALSE, keywidth = 1, keyheight = 1)) +
    ylab("Relative Abundance") +
    ggtitle("") +
    theme_classic() +
    theme(axis.line = element_line(colour = "black"),
          legend.text = element_text(face = "italic"),
          legend.justification = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, size=12),
          #axis.text.x = element_blank(),
          axis.text = element_text(size=12),
          # axis.ticks.x = element_blank(),
          strip.text = element_text(size=12)) +
    ylim(0,1.01) + xlab('') +
    guides(fill = guide_legend(title = "Genus", ncol = 1))
plot_genus

# Fig1c (original Fig2a) - p_alpha1 ----------------------------------------

# Subset the samples to a genus level for Fig2 panels
ps_genus_fig2 <- ps |>
    microbiome::aggregate_taxa(level = "Genus", verbose = T)  |>
    subset_samples(Location == "Arm_left" & Age_group != 'older')

ps_genus_fig2@sam_data$Ichthyosis_type <- ifelse(
    is.na(ps_genus_fig2@sam_data$Ichthyosis_type),
    "Control", as.character(ps_genus_fig2@sam_data$Ichthyosis_type))

count_data_fig2 <-  ps_genus_fig2 %>%
    phylosmith::taxa_filter(frequency = 0.05) |>
    otu_table() |>
    data.frame()

metadata_fig2 <- ps_genus_fig2 |> sample_data() |> data.frame()
rownames(metadata_fig2) <- gsub(pattern = "-", replacement = ".", x = rownames(metadata_fig2))

p_alpha1 <- metadata_fig2 |>
    # dplyr::mutate(Control = case_when( Control == TRUE ~ 'Control',
    #                                    Control == FALSE ~ 'Patient')) |>
    ggplot(data = , aes(x = Control, y = dn_est, fill = Control)) +
    geom_violin(trim = TRUE, scale = "count", width = 0.5, alpha = 0.7) +
    geom_boxplot(width = 0.25, fill = "white", color = "black", outlier.shape = NA, linewidth = 1) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5) +
    scale_fill_manual(values = colv3) +
    labs(x = "", y = "Shannon Diversity") +
    facet_wrap(~Age_group) +
    ggtitle("") +
    ylim(0,4)  +
    theme_def +
    theme(strip.text = element_text(face = 'plain'))
p_alpha1

# Fig1d (original Fig2c) - p_maaslin1 -------------------------------------

Maaslin2_results_fig2 <- Maaslin2(input_data = count_data_fig2,
                                  input_metadata = metadata_fig2 |>
                                      dplyr::mutate(Control = factor(x = Control,
                                                                     levels = c('Control', 'Patient'))),
                                  output = "maaslin2_output_tmp/",
                                  min_abundance = 0.0,
                                  min_prevalence = 0.2,
                                  min_variance = 0.0,
                                  normalization = "TSS",
                                  transform = "AST",
                                  analysis_method = "LM",
                                  max_significance = 1,
                                  random_effects = NULL,
                                  fixed_effects = c("Control", "Age"),
                                  correction = "BH",
                                  standardize = TRUE,
                                  cores = 10,
                                  plot_heatmap = F,
                                  plot_scatter = F,
                                  heatmap_first_n = 50,
                                  reference = "FALSE")

Maaslin2_results_tb_fig2 <- Maaslin2_results_fig2[["results"]] |> dplyr::filter(metadata != "Age")

maaslin_filter_fig2 <- Maaslin2_results_tb_fig2 |>
    dplyr::filter(metadata == 'Control' & pval < 0.05) |>
    dplyr::select(significant_taxa = feature, effect_size = coef, standard_error = stderr, pval, qval)

p_maaslin1 <- maaslin_filter_fig2 |>
    ggplot(aes(x = effect_size, y = factor(significant_taxa, levels = significant_taxa))) +
    geom_point(size = 3, color = "black") +
    geom_errorbar(aes(xmin = effect_size - standard_error, xmax = effect_size + standard_error), height = 0.1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    theme_def +
    theme(axis.text.y = element_text(face = "italic"),
          plot.title = element_text(size = 14)) +
    labs(x = "Effect Size", y = "", title = "") +
    xlim(-1,1)
p_maaslin1

# Fig1e (original Fig3a) - p_alpha1_fig3a ----------------------------------

ps_genus_fig3a <- ps |>
    microbiome::aggregate_taxa(level = "Genus", verbose = T) |>
    subset_samples(Location == "Arm_left" &
                       (Ichthyosis_type == "ARCI" | Control == 'Control'))

ps_genus_fig3a@sam_data$Ichthyosis_type <- ifelse(
    is.na(ps_genus_fig3a@sam_data$Ichthyosis_type),
    "Control", as.character(ps_genus_fig3a@sam_data$Ichthyosis_type))

metadata_fig3a <- ps_genus_fig3a |> sample_data() |> data.frame()
rownames(metadata_fig3a) <- gsub(pattern = "-", replacement = ".", x = rownames(metadata_fig3a))

p_alpha1_fig3a <- ggplot(data = metadata_fig3a, aes(x = Ichthyosis_type, y = dn_est, fill = Ichthyosis_type)) +
    geom_violin(trim = TRUE, scale = "count", width = 0.3, alpha = 0.7) +
    geom_boxplot(width = 0.07, fill = "white", color = "black", outlier.shape = NA, linewidth = 1) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5) +
    scale_fill_manual(values = colv3) +
    labs(x = "", y = "Shannon Diversity") +
    ggtitle("") +
    ylim(0,4)  +
    theme_def
p_alpha1_fig3a

# Fig1f (original Fig3b) - p_beta1_fig3b -----------------------------------

ps_genus_fig3b <- ps |>
    microbiome::aggregate_taxa(level = "Genus", verbose = T) |>
    subset_samples(Location == "Arm_left" &
                       (Ichthyosis_type == "IV" | Control == 'Control'))

ps_genus_fig3b@sam_data$Ichthyosis_type <- ifelse(
    is.na(ps_genus_fig3b@sam_data$Ichthyosis_type),
    "Control", as.character(ps_genus_fig3b@sam_data$Ichthyosis_type))

metadata_fig3b <- ps_genus_fig3b |> sample_data() |> data.frame()

p_alpha2_fig3b <- metadata_fig3b |>
    dplyr::mutate(Ichthyosis_type = factor(Ichthyosis_type, levels = c('IV', 'Control'))) |>
    ggplot(aes(x = Ichthyosis_type, y = dn_est, fill = Ichthyosis_type)) +
    geom_violin(trim = TRUE, scale = "count", width = 0.3, alpha = 0.7) +
    geom_boxplot(width = 0.07, fill = "white", color = "black", outlier.shape = NA, linewidth = 1) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5) +
    scale_fill_manual(values = colv3) +
    labs(x = "", y = "Shannon Diversity") +
    ggtitle("") +
    ylim(0,4)  +
    theme_def
p_alpha2_fig3b

# Fig1g (original Fig3e) - p_maaslin1_fig3e --------------------------------

ps_genus_fig3e <- ps |>
    microbiome::aggregate_taxa(level = "Genus", verbose = T) |>
    subset_samples(Location == "Arm_left" &
                       (Ichthyosis_type == "ARCI" | Control == 'Control'))

ps_genus_fig3e@sam_data$Ichthyosis_type <- ifelse(
    is.na(ps_genus_fig3e@sam_data$Ichthyosis_type),
    "Control", as.character(ps_genus_fig3e@sam_data$Ichthyosis_type))

count_data_fig3e <-  ps_genus_fig3e |>
    phylosmith::taxa_filter(frequency = 0.05) |>
    otu_table() |>
    data.frame()

metadata_fig3e <- ps_genus_fig3e |> sample_data() |> data.frame()
rownames(metadata_fig3e) <- gsub(pattern = "-", replacement = ".", x = rownames(metadata_fig3e))

Maaslin2_results_fig3e <- Maaslin2(input_data = count_data_fig3e,
                                   input_metadata = metadata_fig3e |>
                                       dplyr::mutate(Ichthyosis_type = factor(x = Ichthyosis_type,
                                                                              levels = c('Control', 'ARCI'))),
                                   output = "maaslin2_output_tmp_fig3e/",
                                   min_abundance = 0.0,
                                   min_prevalence = 0.2,
                                   min_variance = 0.0,
                                   normalization = "TSS",
                                   transform = "AST",
                                   analysis_method = "LM",
                                   max_significance = 1,
                                   random_effects = NULL,
                                   fixed_effects = c("Ichthyosis_type", "Age"),
                                   correction = "BH",
                                   standardize = TRUE,
                                   cores = 10,
                                   plot_heatmap = F,
                                   plot_scatter = F,
                                   heatmap_first_n = 50,
                                   reference = "Control")

Maaslin2_results_tb_fig3e <- Maaslin2_results_fig3e[["results"]] |> dplyr::filter(metadata != "Age")

maaslin_filter_fig3e <- Maaslin2_results_tb_fig3e |>
    dplyr::filter(metadata == 'Ichthyosis_type' & pval < 0.05) |>
    dplyr::select(significant_taxa = feature, effect_size = coef, standard_error = stderr, pval, qval)

p_maaslin1_fig3e <- maaslin_filter_fig3e |>
    ggplot(aes(x = effect_size, y = factor(significant_taxa, levels = rev(sort(significant_taxa)) )) ) +
    geom_point(size = 3, color = "black") +
    geom_errorbar(aes(xmin = effect_size - standard_error, xmax = effect_size + standard_error), height = 0.1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    theme_def +
    theme(axis.text.y = element_text(face = "italic"),
          plot.title = element_text(size = 14)) +
    labs(x = "Effect Size", y = "", title = "") +
    xlim(-1,1)
p_maaslin1_fig3e

# Save NEW Figure 1 --------------------------------------------------------

layout_fig1 <- "
AAABBB
AAABBB
CCCCDD
EEFFGG
"

wrap_elements(full = plot_phylum) +
    wrap_elements(full = plot_genus) +
    wrap_elements(full = p_alpha1) +
    wrap_elements(full = p_maaslin1) +
    wrap_elements(full = p_alpha1_fig3a) +
    wrap_elements(full = p_alpha2_fig3b) +
    wrap_elements(full = p_maaslin1_fig3e) +
    plot_layout(design = layout_fig1) +
    plot_annotation(tag_levels = 'a')
ggsave(filename = "plots_JID/Figure_1.pdf", height = 16, width = 15)
