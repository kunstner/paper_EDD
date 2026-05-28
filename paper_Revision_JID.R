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

# Subset the samples to a genus level -------------------------------------

# Remove rows corresponding to abdomen & right_arm & old age from sample data
# Subset the samples to Ichthyosis_type of ARCI

ps_genus <- ps %>%
    microbiome::aggregate_taxa(x = ., level = "Genus", verbose = T)  %>%
    subset_samples(Location == "Arm_left" & Age_group != 'older')

# check if each value in the "Ichthyosis_type" column is NA using is.na().
#If it is NA, it replaces it with the string "TRUE".
#If it is not NA, it converts the value to a character using as.character() to ensure it remains as the original value.
ps_genus@sam_data$Ichthyosis_type <- ifelse(
    is.na(ps_genus@sam_data$Ichthyosis_type),
    "Control", as.character(ps_genus@sam_data$Ichthyosis_type))

table(ps_genus@sam_data$Control)
table(ps_genus@sam_data$Control, ps_genus@sam_data$Age_group)

# Extract count data from phyloseq object for Differential abundance analysis
count_data <-  ps_genus %>%
    phylosmith::taxa_filter(phyloseq_obj = ., frequency = 0.05) %>%
    otu_table() %>%
    data.frame()

# Extract metadata from phyloseq object for Differential abundance analysis
metadata <- ps_genus %>% sample_data() %>% data.frame()
rownames(metadata) <- gsub(pattern = "-", replacement = ".", x = rownames(metadata))

# Sex visualization in beta diversity (Reviewer 1, 4) ---------------------

ps_genus_clr <- ps_genus %>% microbiome::transform(transform = 'clr')
otu.table_clr <- otu_table(ps_genus_clr) %>% t()
ps_genus_clr_dist <- dist(otu.table_clr, method="euclidean") # for PERMANOVA
ps_genus_clr_ord <- phyloseq::ordinate(ps_genus_clr, "RDA", distance = "euclidean")

sample_data(ps_genus_clr)$Age_Sex <- paste(
    sample_data(ps_genus_clr)$Age_group,
    sample_data(ps_genus_clr)$Gender,
    sep = " ")

sample_data(ps_genus_clr)$Age_Sex <- gsub(' f', ' female', sample_data(ps_genus_clr)$Age_Sex)
sample_data(ps_genus_clr)$Age_Sex <- gsub(' m', ' male', sample_data(ps_genus_clr)$Age_Sex)

p_beta1 <- plot_ordination(
    physeq = ps_genus_clr,
    ordination = ps_genus_clr_ord,
    color = 'Control',
    shape = 'Age_Sex') +
    scale_color_manual('Cohort', values = colv3) +
    scale_shape_manual('Age / Sex',
                       values = c('young female'  = 16,  # filled circle
                                  'young male'  = 1,   # open circle
                                  'middle female' = 17,  # filled triangle
                                  'middle male' = 2)) +  # open triangle
    geom_point(size = 3.6) +
    ggtitle("") +
    theme_def +
    theme(legend.position = 'right')
p_beta1 +
    plot_annotation(tag_levels = 'a')
ggsave(filename = 'plots_JID/Figure_S2a.pdf', height = 5, width = 6)

# PERMANOVA
adonis2(ps_genus_clr_dist ~ Control,
        data = metadata, permutations = 9999, by = "margin")

perm_result <- adonis2(ps_genus_clr_dist ~ Control + Age_group + Gender,
        data = metadata, permutations = 9999, by = "margin")
perm_result

# 2. Systemic therapy sensitivity analysis (Reviewers 1, 4) ---------------

# Sensitivity analysis: exclude patients on systemic therapy
systemic_therapy_ids <- ps |>
    sample_data() |>
    data.frame() |>
    dplyr::filter(is.na(Therapy) | Therapy == 'emolients') |>
    pull(SampleID)

ps_sensitivity <- ps_genus |>
    subset_samples( SampleID %in% systemic_therapy_ids )

ps_sens_genus <- ps_sensitivity |>
    microbiome::aggregate_taxa(level = "Genus", verbose = TRUE) |>
    subset_samples(Location == "Arm_left" & Age_group != 'older')

ps_sens_clr <- ps_sens_genus |> microbiome::transform(transform = 'clr')
otu_sens_clr <- otu_table(ps_sens_clr) |> t()
dist_sens_clr <- dist(otu_sens_clr, method = "euclidean")
meta_sens <- ps_sens_genus |> sample_data() |> data.frame()

adonis2(dist_sens_clr ~ Control,
        data = meta_sens, permutations = 9999)
adonis2(dist_sens_clr ~ Control,
        data = meta_sens, permutations = 9999, by = 'margin')


