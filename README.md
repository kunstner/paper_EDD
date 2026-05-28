# 2026_EDD/Ichthyosis_paper

Repository contains code for the following paper:

# Cutaneous dysbiosis in inherited ichthyoses / epidermal differentiation disorders

> **Cutaneous dysbiosis in inherited ichthyoses/epidermal differentiation disorders: a prospective case-control study**

This repository contains all R code used to produce the figures and statistical analyses for the manuscript published in the *Journal of Investigative Dermatology* (JID).

---

## Abstract

Inherited ichthyoses, recently renamed as epidermal differentiation disorders (EDD), are characterised by cutaneous dryness, scaling, hyperkeratosis, and varying degrees of inflammation. Patients suffer from hypohidrosis, recurrent cutaneous infections, eczema, atopy, and intense pruritus. We investigated the cutaneous microbiome in 41 patients with different types of ichthyosis/EDD using a prospective case-control design. Overall, the patient group shows a significant reduction in beta diversity and an increased abundance of *Corynebacterium* and *Staphylococcus* on the skin. The major finding is the early age at which microbiome imbalances are detectable — in childhood, before age 6. Autosomal recessive congenital ichthyosis (ARCI) was associated with more pronounced alterations than ichthyosis vulgaris (IV). Clinical severity correlated with dysbiosis: higher scaling scores correlated with *Corynebacteria* presence, higher disease scores with *Streptococcus* presence, and eczema/atopy with reduced beta diversity.

---

## Repository structure

```
.
├── code/
│   └── helper.R                   # Shared colour palettes, themes, and utility functions
├── data/
│   └── phyloseq.OTU.RDS           # Processed phyloseq object (OTU table + taxonomy + metadata)
├── plots/                         # Output directory for main and supplementary figures
├── paper_Figure_1.R               # Alpha & beta diversity overview; phylum/genus composition
├── paper_Figure_2.R               # Patients vs. controls: alpha/beta diversity + MaAsLin2 DA
├── paper_Figure_3.R               # ARCI vs. IV vs. controls: diversity + MaAsLin2 DA
├── paper_Figure_JID.R             # Restructured composite figure for JID submission
├── paper_Revision_JID.R           # Additional analyses added during peer review
└── R_Ichthyosis_with_ctrl.Rproj   # RStudio project file
```

---

## Script overview

### `code/helper.R`
Defines shared colour palettes (`cols`, `colv3`), ggplot2 themes (`theme_alpha`, `theme_def`), and any utility functions sourced by the figure scripts.

### `paper_Figure_1.R`
Produces the overview diversity and taxonomic composition figure.

- **Alpha diversity** — DivNet-estimated Shannon diversity by body site (left arm, right arm, abdomen) and age group (young, middle, older), using `breakaway::betta_random()` for inference
- **Beta diversity** — CLR-transformed Aitchison distance, ordinated with RDA (`phyloseq::ordinate`); PERMANOVA via `vegan::adonis2` with pairwise tests
- **Phylum- and genus-level composition** — stacked bar charts faceted by patient/control status, body site, and age group; low-abundance genera (< 1.5%) collapsed into "Other"

Output: `plots/Figure_1.pdf`

### `paper_Figure_2.R`
Compares all patients vs. healthy controls, restricted to the left arm in young and middle age groups.

- **Alpha diversity** — violin/boxplot of DivNet Shannon estimates by group and age; `breakaway::betta()` for statistical testing
- **Beta diversity** — CLR + Aitchison RDA ordination; PERMANOVA adjusted for ichthyosis type, age, and sex
- **Differential abundance** — `Maaslin2` (TSS normalisation, AST transform, LM, BH correction); results exported to `maaslin_Patients_Controls.xlsx`
- **Individual CLR abundances** — violin/jitter plots for significant genera

Output: `plots/Figure_2.pdf`

### `paper_Figure_3.R`
Subtype-specific comparisons, all on the left arm:

- **ARCI vs. controls** — alpha/beta diversity + MaAsLin2; exports `maaslin_ARCI.xlsx`
- **IV vs. controls** — alpha/beta diversity + MaAsLin2; exports `maaslin_IV.xlsx`
- **ARCI vs. IV** — alpha/beta diversity + MaAsLin2; exports `maaslin_ARCI_IV.xlsx`
- **Supplementary Figure S1** — individual CLR violin plots for ARCI and ARCI vs. IV significant genera

Output: `plots/Figure_3.pdf`, `plots/Figure_S1.pdf`

### `paper_Figure_JID.R`
Restructured composite figure prepared for the final JID manuscript submission. Recombines panels from Figures 1–3 into a single layout:

- Fig 1a/1b — phylum and genus composition (all patients vs. controls)
- Fig 1c/1d — alpha diversity by group/age; MaAsLin2 effect sizes (patients vs. controls)
- Fig 1e/1f — alpha diversity for ARCI vs. controls and IV vs. controls
- Fig 1g — MaAsLin2 effect sizes (ARCI vs. controls)

Output: `plots_JID/Figure_1.pdf`

### `paper_Revision_JID.R`
Supplementary analyses added in response to reviewer comments during peer review (e.g., sensitivity analyses for beta diversity, adjusted models).

---

## Data

The main input is `data/phyloseq.OTU.RDS`, a [`phyloseq`](https://joey711.github.io/phyloseq/) object containing:

- OTU table (16S amplicon sequencing)
- Taxonomy table
- Sample metadata (patient/control status, ichthyosis type, body-site location, age group, DivNet alpha diversity estimates)

Raw sequencing data are deposited at [ENA study PRJEB112279](https://www.ebi.ac.uk/ena/browser/view/PRJEB112279).

---

## Dependencies

All analyses were performed in R. Key packages:

| Package | Purpose |
|---|---|
| `phyloseq` | Microbiome data handling |
| `microbiome` | CLR transformation, taxa aggregation |
| `vegan` | PERMANOVA (`adonis2`), Aitchison distances |
| `breakaway` | DivNet alpha diversity estimation (`betta`, `betta_random`) |
| `Maaslin2` | Differential abundance analysis |
| `tidyverse` | Data wrangling and visualisation |
| `ggpubr` | Violin/box plot helpers |
| `patchwork` | Multi-panel figure assembly |
| `ggh4x` | Nested facets |
| `phylosmith` | Prevalence filtering |

---

## Reproducibility

A fixed random seed (1253) is used where applicable. Scripts are self-contained except for the shared `code/helper.R` and the `data/phyloseq.OTU.RDS` input.

---

## Citation

If you use this code or data, please cite the original publication (citation will be updated upon journal assignment of a DOI).

---

## Authors

- Ahmed Abdelhamid (primary analysis, Figures 2 & 3)
- Axel Künstner (analysis, revisions, figure restructuring)
