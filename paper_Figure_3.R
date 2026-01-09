# Credentials -------------------------------------------------------------

#
# Author: Ahmed Abelhamid
# Changes: Axel Künstner
# Project: Ichtyosis
# Data: Microbiome

# Libraires ---------------------------------------------------------------

library(sf)
library(tidyverse)
library(phyloseq)
library(vegan)
library(patchwork)
library(microbiome)
library(breakaway)
library(ggpubr)
library(Maaslin2)
library(WriteXLS)

source(file = "code/helper.R")

# Get data ARCI vs Controls -----------------------------------------------

ps <- readRDS(file = "data/phyloseq.OTU.RDS")

sample_data(ps) <- sample_data(ps) %>%
    data.frame() %>%
    dplyr::mutate(Age_group = case_when(Age_group == 'old' ~ 'older',
                                        .default = Age_group)) %>%
    dplyr::mutate(Age_group = factor(Age_group,
                                     levels = c('young', 'middle', 'older')))

# Plotting theme ----------------------------------------------------------

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

# Subset the samples to a genus level -------------------------------------

# Remove rows corresponding to abdomen & right_arm & old age from sample data
# Subset the samples to Ichthyosis_type of ARCI

ps_genus <- ps %>%
    microbiome::aggregate_taxa(x = ., level = "Genus", verbose = T) %>%
    subset_samples(Location == "Arm_left" &
                       (Ichthyosis_type == "ARCI" | Control == TRUE))

# check if each value in the "Ichthyosis_type" column is NA using is.na().
#If it is NA, it replaces it with the string "TRUE".
#If it is not NA, it converts the value to a character using as.character() to ensure it remains as the original value.
ps_genus@sam_data$Ichthyosis_type <- ifelse(
    is.na(ps_genus@sam_data$Ichthyosis_type),
    "Control", as.character(ps_genus@sam_data$Ichthyosis_type))

table(ps_genus@sam_data$Ichthyosis_type)

# Extract count data from phyloseq object for Differential abundance analysis
count_data <-  ps_genus %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 0.05) %>%
    otu_table() %>%
    data.frame()

# Extract metadata from phyloseq object for Differential abundance analysis
metadata <- ps_genus %>% sample_data() %>% data.frame()
rownames(metadata) <- gsub(pattern = "-", replacement = ".", x = rownames(metadata))

# Alpha diversity ---------------------------------------------------------

p_alpha1 <- ggplot(data = metadata, aes(x = Ichthyosis_type, y = dn_est, fill = Ichthyosis_type)) +
    geom_violin(trim = TRUE, scale = "count", width = 0.3, alpha = 0.7) +
    geom_boxplot(width = 0.07, fill = "white", color = "black", outlier.shape = NA, linewidth = 1) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5) +
    scale_fill_manual(values = colv3) +
    labs(x = "", y = "Shannon Diversity") +
    ggtitle("") +
    ylim(0,4)  +
    theme_def
p_alpha1

# betta
betta_Value <- betta(formula = dn_est ~ Ichthyosis_type + Age,
                     ses = dn_error, p.digits = 8,
                     data = metadata %>%
                         dplyr::mutate(Ichthyosis_type = factor(x = Ichthyosis_type,
                                                                levels = c('Control', 'ARCI'))))
betta_Value$table

# Beta diversity ----------------------------------------------------------

ps_genus_clr <- ps_genus %>% microbiome::transform(transform = 'clr')
otu.table_clr <- otu_table(ps_genus_clr) %>% t()
ps_genus_clr_dist <- dist(otu.table_clr, method="euclidean") # for PERMANOVA
ps_genus_clr_ord <- phyloseq::ordinate(ps_genus_clr, "RDA", distance = "euclidean")
p_beta1 <- plot_ordination(
    physeq = ps_genus_clr,
    ordination = ps_genus_clr_ord, color='Ichthyosis_type')  +
    scale_color_manual('Cohort', values = colv3) +
    geom_point(size = 3.6) +
    ggtitle("") +
    theme_def +
    theme(legend.position = c(0.1, 0.8))
p_beta1

# PERMANOVA
adonis2(ps_genus_clr_dist ~ Ichthyosis_type + Age,
        data = metadata, permutations = 9999)

# Differential abundance analysis by MaAsLin2 -----------------------------

Maaslin2_results <- Maaslin2(input_data = count_data,
                             input_metadata = metadata %>%
                                 dplyr::mutate(Ichthyosis_type = factor(x = Ichthyosis_type,
                                                                        levels = c('Control', 'ARCI'))),
                             output = "maaslin2_output_tmp/",
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
                             plot_heatmap = TRUE,
                             plot_scatter = TRUE,
                             heatmap_first_n = 50,
                             reference = "Control")
Maaslin2_results_tb <- Maaslin2_results[["results"]] %>% dplyr::filter(metadata != "Age")
Maaslin2_results_tb
# filter maaslin results
maaslin_filter <- Maaslin2_results_tb %>%
    dplyr::filter(metadata == 'Ichthyosis_type' & pval < 0.05) %>%
    dplyr::select(significant_taxa = feature, effect_size = coef, standard_error = stderr, pval, qval)

maaslin_filter

mylist <- list()
mylist[["ARCI vs Controls"]] <- maaslin_filter
WriteXLS::WriteXLS(x = mylist, ExcelFileName = 'maaslin_ARCI.xlsx', BoldHeaderRow = TRUE, row.names = TRUE, FreezeRow = 1)

# Create the error bar plot
p_maaslin1 <- maaslin_filter %>%
    ggplot(aes(x = effect_size, y = factor(significant_taxa, levels = rev(sort(significant_taxa)) )) ) +
    geom_point(size = 3, color = "black") +
    geom_errorbarh(aes(xmin = effect_size - standard_error, xmax = effect_size + standard_error), height = 0.1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # coord_cartesian(xlim = c(-0.5, 0.5)) +
    theme_def +
    theme(axis.text.y = element_text(face = "italic"),  # Make taxa names italic
          plot.title = element_text(size = 14)) +  # Increase title size
    labs(x = "Effect Size", y = "", title = "") +
    xlim(-1,1)
p_maaslin1

# plot individuals clr values
plot_data3 <- ps_genus %>%
    microbiome::abundances(x = ., transform = "clr") %>%
    otu_table(taxa_are_rows = T) %>% t()%>% data.frame()

plot_data3$Group <- ps_genus %>%
    sample_data() %>% data.frame() %>% .$Ichthyosis_type
colnames(plot_data3) <- gsub(pattern = "unclassified", replacement = "uncl.", x = colnames(plot_data3))

p_da1 <- plot_data3 %>%
    rownames_to_column('Sample') %>%
    pivot_longer(cols = -c("Sample", "Group"), names_to = "Genus", values_to = "CLR") %>%
    mutate(Genus = gsub(pattern = "\\.", replacement = "-", x = Genus)) %>%
    dplyr::filter(Genus %in% maaslin_filter$significant_taxa) %>%
    # left_join(maaslin_filter, by = c("Genus" = "significant_taxa")) %>%
    ggplot(data = ., aes(x = Group, y = CLR, fill = Group)) +
    geom_violin(width = 0.5) +
    geom_boxplot(outlier.shape = NA, width = 0.15, fill = "grey95") +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
    scale_fill_manual(values = colv3) +
    facet_wrap(~ Genus, ncol = 4) +
    theme_def +
    xlab('') +
    ylab('CLR transformed abundances')
p_da1

# Get data: IV vs Controls ------------------------------------------------

# IV vs Controls

ps <- readRDS(file = "data/phyloseq.OTU.RDS")

# Subset the samples to a genus level -------------------------------------

# Remove rows corresponding to abdomen & right_arm & old age from sample data
# Subset the samples to Ichthyosis_type of IV

ps_genus <- ps %>%
    microbiome::aggregate_taxa(x = ., level = "Genus", verbose = T) %>%
    subset_samples(Location == "Arm_left" &
                       (Ichthyosis_type == "IV" | Control == TRUE))

# check if each value in the "Ichthyosis_type" column is NA using is.na().
#If it is NA, it replaces it with the string "TRUE".
#If it is not NA, it converts the value to a character using as.character() to ensure it remains as the original value.
ps_genus@sam_data$Ichthyosis_type <- ifelse(
    is.na(ps_genus@sam_data$Ichthyosis_type),
    "Control", as.character(ps_genus@sam_data$Ichthyosis_type))

table(ps_genus@sam_data$Ichthyosis_type)

# Extract count data from phyloseq object for Differential abundance analysis
count_data <-  ps_genus %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 0.05) %>%
    otu_table() %>%
    data.frame()

# Extract metadata from phyloseq object for Differential abundance analysis
metadata <- ps_genus %>% sample_data() %>% data.frame()
rownames(metadata) <- gsub(pattern = "-", replacement = ".", x = rownames(metadata))

# Alpha diversity ---------------------------------------------------------

p_alpha2 <- metadata %>%
    dplyr::mutate(Ichthyosis_type = factor(Ichthyosis_type, levels = c("IV", "Control"))) %>%
    ggplot(data = , aes(x = Ichthyosis_type, y = dn_est, fill = Ichthyosis_type)) +
    geom_violin(trim = TRUE, scale = "count", width = 0.3, alpha = 0.7) +
    geom_boxplot(width = 0.07, fill = "white", color = "black", outlier.shape = NA, linewidth = 1) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5) +
    scale_fill_manual(values = colv3) +
    labs(x = "", y = "Shannon Diversity") +
    ggtitle("") +
    ylim(0,4)  +
    theme_def
p_alpha2

# betta
betta_Value <- betta(formula = dn_est ~ Ichthyosis_type + Age,
                     ses = dn_error, p.digits = 8,
                     data = metadata %>%
                         dplyr::mutate(Ichthyosis_type = factor(x = Ichthyosis_type,
                                                                levels = c('Control', 'IV'))))
betta_Value$table

# Beta diversity ----------------------------------------------------------

ps_genus_clr <- ps_genus %>% microbiome::transform(transform = 'clr')
otu.table_clr <- otu_table(ps_genus_clr) %>% t()
ps_genus_clr_dist <- dist(otu.table_clr, method="euclidean") # for PERMANOVA
ps_genus_clr_ord <- phyloseq::ordinate(ps_genus_clr, "RDA", distance = "euclidean")
p_beta2 <- plot_ordination(
    physeq = ps_genus_clr,
    ordination = ps_genus_clr_ord, color='Ichthyosis_type')  +
    scale_color_manual('Cohort', values = colv3) +
    geom_point(size = 3.6) +
    ggtitle("") +
    theme_def +
    theme(legend.position = c(0.8, 0.8))
p_beta2

# PERMANOVA
adonis2(ps_genus_clr_dist ~ Ichthyosis_type + Age,
        data = metadata, permutations = 9999)

# Differential abundance analysis by MaAsLin2 -----------------------------

Maaslin2_results <- Maaslin2(input_data = count_data,
                             input_metadata = metadata %>%
                                 dplyr::mutate(Ichthyosis_type = factor(x = Ichthyosis_type,
                                                                        levels = c('Control', 'IV'))),
                             output = "maaslin2_output_tmp/",
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
                             plot_heatmap = TRUE,
                             plot_scatter = TRUE,
                             heatmap_first_n = 50,
                             reference = "Control")
Maaslin2_results_tb <- Maaslin2_results[["results"]] %>% dplyr::filter(metadata != "Age")
Maaslin2_results_tb
# filter maaslin results
maaslin_filter <- Maaslin2_results_tb %>%
    dplyr::filter(metadata == 'Ichthyosis_type' & pval < 0.05) %>%
    dplyr::select(significant_taxa = feature, effect_size = coef, standard_error = stderr, pval, qval)

maaslin_filter

mylist <- list()
mylist[["IV vs Controls"]] <- maaslin_filter
WriteXLS::WriteXLS(x = mylist, ExcelFileName = 'maaslin_IV.xlsx', BoldHeaderRow = TRUE, row.names = TRUE, FreezeRow = 1)

# Get data: IV vs ARCI-----------------------------------------------------

ps <- readRDS(file = "data/phyloseq.OTU.RDS")

# Subset the samples to a genus level -------------------------------------

# Remove rows corresponding to abdomen & right_arm & old age from sample data
# Subset the samples to Ichthyosis_type of ARCI and IV

ps_genus <- ps %>%
    microbiome::aggregate_taxa(x = ., level = "Genus", verbose = T) %>%
    subset_samples(Location == "Arm_left" &
                       (Ichthyosis_type == "ARCI" | Ichthyosis_type == "IV"))

# check if each value in the "Ichthyosis_type" column is NA using is.na().
#If it is NA, it replaces it with the string "TRUE".
#If it is not NA, it converts the value to a character using as.character() to ensure it remains as the original value.
# ps_genus@sam_data$Ichthyosis_type <- ifelse(
#     is.na(ps_genus@sam_data$Ichthyosis_type),
#     "Control", as.character(ps_genus@sam_data$Ichthyosis_type))

table(ps_genus@sam_data$Ichthyosis_type)

# Extract count data from phyloseq object for Differential abundance analysis
count_data <-  ps_genus %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 0.05) %>%
    otu_table() %>%
    data.frame()

# Extract metadata from phyloseq object for Differential abundance analysis
metadata <- ps_genus %>% sample_data() %>% data.frame()
rownames(metadata) <- gsub(pattern = "-", replacement = ".", x = rownames(metadata))

# Alpha diversity ---------------------------------------------------------

p_alpha3 <- ggplot(data = metadata, aes(x = Ichthyosis_type, y = dn_est, fill = Ichthyosis_type)) +
    geom_violin(trim = TRUE, scale = "count", width = 0.3, alpha = 0.7) +
    geom_boxplot(width = 0.07, fill = "white", color = "black", outlier.shape = NA, linewidth = 1) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.5) +
    scale_fill_manual(values = colv3) +
    labs(x = "", y = "Shannon Diversity") +
    ggtitle("") +
    ylim(0,4)  +
    theme_def
p_alpha3

# betta
betta_Value <- betta(formula = dn_est ~ Ichthyosis_type + Age,
                     ses = dn_error, p.digits = 8,
                     data = metadata %>%
                         dplyr::mutate(Ichthyosis_type = factor(x = Ichthyosis_type,
                                                                levels = c('ARCI', 'IV'))))
betta_Value$table

# Beta diversity ----------------------------------------------------------

ps_genus_clr <- ps_genus %>% microbiome::transform(transform = 'clr')
otu.table_clr <- otu_table(ps_genus_clr) %>% t()
ps_genus_clr_dist <- dist(otu.table_clr, method="euclidean") # for PERMANOVA
ps_genus_clr_ord <- phyloseq::ordinate(ps_genus_clr, "RDA", distance = "euclidean")
p_beta3 <- plot_ordination(
    physeq = ps_genus_clr,
    ordination = ps_genus_clr_ord, color='Ichthyosis_type')  +
    scale_color_manual('Cohort', values = colv3) +
    geom_point(size = 3.6) +
    ggtitle("") +
    theme_def +
    theme(legend.position = c(0.8, 0.8))
p_beta3

# PERMANOVA
adonis2(ps_genus_clr_dist ~ Ichthyosis_type + Age,
        data = metadata, permutations = 9999)

# Differential abundance analysis by MaAsLin2 -----------------------------

Maaslin2_results <- Maaslin2(input_data = count_data,
                             input_metadata = metadata %>%
                                 dplyr::mutate(Ichthyosis_type = factor(x = Ichthyosis_type,
                                                                        levels = c('ARCI', 'IV'))),
                             output = "maaslin2_output_tmp/",
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
                             plot_heatmap = TRUE,
                             plot_scatter = TRUE,
                             heatmap_first_n = 50,
                             reference = "ARCI")
Maaslin2_results_tb <- Maaslin2_results[["results"]] %>% dplyr::filter(metadata != "Age")
Maaslin2_results_tb
# filter maaslin results
maaslin_filter <- Maaslin2_results_tb %>%
    dplyr::filter(metadata == 'Ichthyosis_type' & pval < 0.05) %>%
    dplyr::select(significant_taxa = feature, effect_size = coef, standard_error = stderr, pval, qval)

maaslin_filter

mylist <- list()
mylist[["ARCI vs IV"]] <- maaslin_filter
WriteXLS::WriteXLS(x = mylist, ExcelFileName = 'maaslin_ARCI_IV.xlsx', BoldHeaderRow = TRUE, row.names = TRUE, FreezeRow = 1)

# Create the error bar plot
p_maaslin3 <- maaslin_filter %>%
    ggplot(aes(x = effect_size, y = factor(significant_taxa, levels = rev(sort(significant_taxa)) )) ) +
    geom_point(size = 3, color = "black") +
    geom_errorbarh(aes(xmin = effect_size - standard_error, xmax = effect_size + standard_error), height = 0.1) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    # coord_cartesian(xlim = c(-0.5, 0.5)) +
    theme_def +
    theme(axis.text.y = element_text(face = "italic"),  # Make taxa names italic
          plot.title = element_text(size = 14)) +  # Increase title size
    labs(x = "Effect Size", y = "", title = "") +
    xlim(-1,1)
p_maaslin3

# plot individuals clr values
plot_data3 <- ps_genus %>%
    microbiome::abundances(x = ., transform = "clr") %>%
    otu_table(taxa_are_rows = T) %>% t()%>% data.frame()

plot_data3$Group <- ps_genus %>%
    sample_data() %>% data.frame() %>% .$Ichthyosis_type
colnames(plot_data3) <- gsub(pattern = "unclassified", replacement = "uncl.", x = colnames(plot_data3))

p_da3 <- plot_data3 %>%
    rownames_to_column('Sample') %>%
    pivot_longer(cols = -c("Sample", "Group"), names_to = "Genus", values_to = "CLR") %>%
    mutate(Genus = gsub(pattern = "\\.", replacement = "-", x = Genus)) %>%
    dplyr::filter(Genus %in% maaslin_filter$significant_taxa) %>%
    # left_join(maaslin_filter, by = c("Genus" = "significant_taxa")) %>%
    ggplot(data = ., aes(x = Group, y = CLR, fill = Group)) +
    geom_violin(width = 0.5) +
    geom_boxplot(outlier.shape = NA, width = 0.15, fill = "grey95") +
    geom_jitter(width = 0.1, alpha = 0.5, size = 2) +
    scale_fill_manual(values = colv3) +
    facet_wrap(~ Genus, ncol = 5) +
    theme_def +
    xlab('') +
    ylab('CLR transformed abundances')
p_da3


# Final plot --------------------------------------------------------------

# layout <- "
# ABBC
# DEE#
# FGGH
# "
#
# p_alpha1 + p_beta1 + theme(legend.position = 'none') +
#     p_maaslin1 + #p_da1 +
#     p_alpha2 + p_beta2 + theme(legend.position = 'none') +
#     p_alpha3 + p_beta3 + theme(legend.position = 'none') +
#     p_maaslin3 +
#     plot_layout(design = layout, guides = 'collect') +
#     plot_annotation(tag_levels = 'A')
#
# ggsave(filename = "plots/Figure_3.pdf", height = 13.5, width = 15)

# layout <- "
# AABBCC
# DDEEFF
# GGGHHH
# "
# p_alpha1 + p_alpha2 + p_alpha3 +
#     p_beta1 + theme(legend.position = 'none') +
#     p_beta2 + theme(legend.position = 'none') +
#     p_beta3 + theme(legend.position = 'none') +
#     wrap_elements(full = p_maaslin1) +
#     wrap_elements(full = p_maaslin3) +
#     plot_layout(design = layout, guides = 'collect') +
#     plot_annotation(tag_levels = 'a')
# ggsave(filename = "plots/Figure_3.pdf", height = 13.5, width = 15)

layout <- "
AABBCCDD
EEFFGGHH
"
p_alpha1 + p_alpha2 +
    p_beta1 + theme(legend.position = 'none') +
    p_beta2 + theme(legend.position = 'none') +
    wrap_elements(full = p_maaslin1) +

    p_alpha3 +
    p_beta3 + theme(legend.position = 'none') +
    wrap_elements(full = p_maaslin3) +
    plot_layout(design = layout, guides = 'collect') +
    plot_annotation(tag_levels = 'a')
ggsave(filename = "plots/Figure_3.pdf", height = 7.5, width = 15)

layout <- "
AAAA#
AAAA#
BBBBB
"

p_da1 + p_da3 +
    plot_layout(design = layout, guides = 'collect') +
    plot_annotation(tag_levels = 'a')
ggsave(filename = "plots/Figure_S1.pdf", height = 9, width = 15)

# Session info ------------------------------------------------------------

sessioninfo::session_info()
