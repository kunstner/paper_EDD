# Credentials -------------------------------------------------------------

#
# Author: Axel Künstner
# Project: Ichtyosis
# Data: Microbiome

# Libraries ---------------------------------------------------------------

library(tidyverse)
library(phyloseq)

library(patchwork)
library(ggpubr)
library(flextable)

library(vegan)
library(breakaway)

source("code/helper.R")

# Data and variables ------------------------------------------------------

SEED <- 1253

ps   <- readRDS(file = "data/phyloseq.OTU.RDS")

sample_data(ps) <- sample_data(ps) %>%
    data.frame() %>%
    dplyr::mutate(Control = case_when(Control == TRUE ~ 'Control',
                               Control == FALSE ~ 'Patient'))

sample_data(ps) <- sample_data(ps) %>%
    data.frame() %>%
    dplyr::mutate(Age_group = case_when(Age_group == 'old' ~ 'older',
                                 .default = Age_group)) %>%
    dplyr::mutate(Age_group = factor(Age_group,
                                     levels = c('young', 'middle', 'older')))

cova <- ps %>% sample_data %>% data.frame

# Alpha diversity ---------------------------------------------------------

# p_shannon_1 <- ps %>%
#     sample_data() %>%
#     data.frame() %>%
#     dplyr::mutate(Location = gsub(pattern = "_", replacement = " ", x = Location)) %>%
#     dplyr::filter( Control != "Control") %>%
#     ggpubr::ggviolin(data = .,
#                      x = "Location", y = "dn_est",
#                      trim = TRUE, width = 0.5,
#                      palette = colv3, fill = "Location",
#                      alpha = 0.75,
#                      add = c("median_q1q3", "jitter"),
#                      add.params = list(alpha = 0.5)) +
#     # ggpubr::stat_compare_means(method = 'kruskal.test') +
#     ylim(0,5) +
#     xlab("") + ylab("DivNet estimate of Shannon") +
#     theme_alpha
# p_shannon_1

p_shannon_1 <- ps %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::mutate(Location = gsub(pattern = "_", replacement = " ", x = Location)) %>%
    dplyr::filter( Control != "Control") %>%
    ggplot(data = , aes(x = Location, y = dn_est, fill = Location)) +
    geom_violin(trim = TRUE, scale = "count", width = 0.5, alpha = 0.7) +
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, linewidth = 1) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.3) +
    scale_fill_manual(values = colv3) +
    labs(x = "", y = "Shannon Diversity") +
    # facet_wrap(~Age_group) +
    ggtitle("") +
    ylim(0,4)  +
    theme_alpha +
    theme(strip.text = element_text(face = 'plain'))
p_shannon_1

p_shannon_2 <- ps %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::mutate(Location = gsub(pattern = "_", replacement = " ", x = Location)) %>%
    dplyr::filter( Control != "Control") %>%
    ggpubr::ggviolin(data = .,
                     x = "Location", y = "dn_est",
                     trim = TRUE, width = 0.5,
                     palette = colv3, fill = "Location",
                     alpha = 0.75,
                     add = c("median_q1q3", "jitter"),
                     add.params = list(alpha = 0.5)) +
    facet_wrap(Age_group~., scales = "free_x") +
    # ggpubr::stat_compare_means() +
    ylim(0,5) +
    xlab("") + ylab("DivNet estimate of Shannon") +
    theme_alpha
p_shannon_2

# p_shannon_3 <- ps %>%
#     sample_data() %>%
#     data.frame() %>%
#     dplyr::mutate(Location = gsub(pattern = "_", replacement = " ", x = Location)) %>%
#     dplyr::filter( Control != "Control") %>%
#     ggpubr::ggviolin(data = .,
#                      x = "Age_group", y = "dn_est",
#                      trim = TRUE, width = 0.5,
#                      palette = colv3, fill = "Age_group",
#                      alpha = 0.25,
#                      add = c("median_q1q3", "jitter"),
#                      add.params = list(alpha = 0.5)) +
#     # ggpubr::stat_compare_means() +
#     ylim(0,5) +
#     xlab("") + ylab("DivNet estimate of Shannon") +
#     theme_alpha
# p_shannon_3

p_shannon_3 <- ps %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::mutate(Location = gsub(pattern = "_", replacement = " ", x = Location)) %>%
    dplyr::filter( Control != "Control") %>%
    ggplot(data = , aes(x = Age_group, y = dn_est, fill = Age_group)) +
    geom_violin(trim = TRUE, scale = "count", width = 0.5, alpha = 0.50) +
    geom_boxplot(width = 0.15, fill = "white", color = "black", outlier.shape = NA, linewidth = 1) +
    geom_jitter(position = position_jitter(width = 0.1), size = 2, alpha = 0.3) +
    scale_fill_manual(values = colv3) +
    labs(x = "", y = "Shannon Diversity") +
    # facet_wrap(~Age_group) +
    ggtitle("") +
    ylim(0,4)  +
    theme_alpha +
    theme(strip.text = element_text(face = 'plain'))
p_shannon_3

p_shannon_4 <- ps %>%
    sample_data() %>%
    data.frame() %>%
    dplyr::mutate(Location = gsub(pattern = "_", replacement = " ", x = Location)) %>%
    dplyr::filter( Control != "Control") %>%
    ggpubr::ggviolin(data = .,
                     x = "Age_group", y = "dn_est",
                     trim = TRUE, width = 0.5,
                     palette = colv3, fill = "Age_group",
                     alpha = 0.25,
                     add = c("median_q1q3", "jitter"),
                     add.params = list(alpha = 0.5)) +
    facet_wrap(Location~., scales = "free_x") +
    # ggpubr::stat_compare_means() +
    ylim(0,5) +
    xlab("") + ylab("DivNet estimate of Shannon") +
    theme_alpha
p_shannon_4

betta_1 <- betta_random(formula = dn_est ~ (Location) | Age_group,
                        ses = dn_error, p.digits = 6,
                        data = cova %>%
                            dplyr::filter(Control != 'Control') %>%
                            dplyr::mutate(Location = factor(x = Location,
                                                            levels = c("Arm_left", "Abdomen", "Arm_right")))
)
betta_1$table

betta_2 <- betta_random(formula = dn_est ~ ( Age_group ) | Location,
                        ses = dn_error, p.digits = 6,
                        data = cova %>%
                            dplyr::filter(Control != 'Control') )
betta_2$table

betta_random(formula = dn_est ~ ( Age_group ) | Location,
             ses = dn_error, p.digits = 6,
             data = cova %>%
                 dplyr::filter(Control != 'Control' & Age_group != 'young') %>%
                 dplyr::mutate(Age_group = factor(Age_group)))$table

# Beta diversity ----------------------------------------------------------

ps_clr <- microbiome::transform(ps %>%
                                    subset_samples(Control != 'Control'), "clr")
otu.table_clr <- otu_table(ps_clr) %>% t()
ps_clr_dist <- dist(otu.table_clr, method="euclidean")
ps_clr_ord <- phyloseq::ordinate(ps_clr, "RDA", distance = "euclidean")
p_beta_all <- plot_ordination(
    physeq = ps_clr,
    ordination = ps_clr_ord, color='Location', shape = "Age_group")  +
    scale_color_manual(values = colv3) +
    geom_point(size=3) +
    geom_line(mapping = aes(group = Patient), color = "grey45", linetype = "dashed", linewidth = 0.25) +
    #geom_text(mapping = aes(label = Volunteer), size = 3, nudge_x = 0.25) +
    ggtitle("RDA of Aitchison distance") +
    theme(legend.key = element_rect(fill = "white"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black"))
p_beta_all

p_beta_all +
    geom_text(mapping = aes(label = Patient), color = 'black', size = 3, check_overlap = T, nudge_x = 0.05)
ggsave(filename = "plots/diversity_beta_IDs.pdf", height = 12, width = 12)

beta_df_all <- vegan::adonis2(formula = ps_clr_dist ~ Location * Age_group,
                              data = data.frame(ps_clr %>% sample_data),
                              permutations = 9999) %>% data.frame() %>%
    round(4) %>%
    rownames_to_column("Variable")
beta_df_all

perm_ait_pairw <- RVAideMemoire::pairwise.perm.manova(resp = ps_clr_dist,
                                                      fact = ps_clr %>% sample_data() %>%
                                                          data.frame() %>%
                                                          .$Age_group %>% as.character,
                                                      test = "Pillai", nperm = 9999, progress = TRUE,
                                                      p.method = "fdr", F = TRUE, R2 = TRUE)
pairwiseAdonis::pairwise.adonis2(x = ps_clr_dist ~ Age_group, nperm = 9999,
                                 data = ps_clr %>% sample_data() %>% data.frame())

# Taxonomics --------------------------------------------------------------

# Plot Phylum
plot_phylum <- ps %>%
    tax_glom(taxrank = "Phylum") %>%                     # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} ) %>%
    psmelt() %>%                                         # Melt to long format
    #filter(Abundance > 0.05) %>%                         # Filter out low abundance taxa
    arrange(Phylum) %>%
    dplyr::select( Phylum, Location, Abundance, Age_group, Control) %>%
    dplyr::group_by(Phylum, Location, Age_group, Control) %>%
    dplyr::reframe(mean = mean(Abundance), groups = c(Phylum)) %>%
    unique() %>%
    dplyr::mutate(Location = gsub(pattern = "_", replacement = " ", x = Location)) %>%
    dplyr::mutate(Location = factor(x = Location, levels = c('Arm left', 'Arm right', 'Abdomen'))) %>%
    dplyr::mutate(Phylum = gsub(pattern = "_", replacement = " ", x = Phylum)) %>%
    ggplot(., aes(x = Age_group, y = mean, fill = Phylum)) +
    geom_bar(stat = "identity") +
    # facet_wrap(~Control + Location, scales = "free_x", ncol = 4) +
    ggh4x::facet_nested_wrap(facets = ~Control + Location,
                             ncol = 4) +
    scale_fill_manual(values = cols) +
    # Remove x axis title
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
          #axis.text.x = element_blank(),
          axis.text = element_text(size=12),
          # axis.ticks.x = element_blank(),
          strip.text = element_text(size=12)) +
    ylim(0,1.01) + xlab('') +
    guides(fill = guide_legend(title = "Phylum", ncol = 1))
plot_phylum

# Genera
means_genus <- ps %>%
    tax_glom(taxrank = "Genus") %>% # agglomerate at phylum level
    transform_sample_counts(function(x) {x/sum(x)} ) %>%
    psmelt() %>% # Melt to long format
    dplyr::arrange(Genus) %>%
    dplyr::select( Genus, Location, Abundance, Age_group, Control) %>%
    dplyr::group_by(Genus, Location, Age_group, Control) %>%
    dplyr::reframe(mean = mean(Abundance), groups = c(Genus)) %>%
    unique() %>%
    dplyr::mutate(Genus = gsub(pattern = "_", replacement = " ", x = Genus))
# filter(mean >= 0.02) # Filter out low abundance taxa

# summarise lowly abundant genera
means_genus$Genus[means_genus$mean <= 0.015] <- "Other"
means_genus <- means_genus %>%
    dplyr::group_by(Genus, Location, Age_group, Control) %>%
    reframe(mean = sum(mean), groups = c(Genus)) %>%
    unique()

plot_genus <- means_genus %>%
    dplyr::mutate(Location = gsub(pattern = "_", replacement = " ", x = Location)) %>%
    dplyr::mutate(Location = factor(x = Location, levels = c('Arm left', 'Arm right', 'Abdomen'))) %>%
    #dplyr::mutate(Genus = gsub(pattern = "f__", replacement = "", x = Genus)) %>%
    #dplyr::mutate(Genus = gsub(pattern = "_unclassified", replacement = " uncl.", x = Genus)) %>%
    ggplot(data = ., aes(x = Age_group, y = mean, fill = Genus)) +
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

# Pretty plotting ---------------------------------------------------------

layout <- "
ABCC
ABCC
ABCC
ABCC
DDEE
DDEE
DDEE
DDEE
DDEE
DDEE
DDEE
"
p_shannon_1 +
    p_shannon_3 +
    p_beta_all + ggtitle('') +
    plot_phylum +
    plot_genus +
    plot_layout(design = layout) +
    plot_annotation(tag_levels = 'a')
ggsave(filename = "plots/Figure_1.pdf", height = 12, width = 15)
