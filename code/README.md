# FoodSeq-FLOW Wastewater Dietary DNA Analysis

## Project Overview

This repository contains a complete reproducible analysis pipeline for the manuscript:

**"Dietary DNA in municipal wastewater reveals signatures of wealth, immigration, and coastal proximity"**

*Dong et al., PNAS 2025, 122(5), e2313516122*

### Key Statistics

- **Total Samples**: 183 wastewater samples from 19 municipalities across North Carolina
- **Population Monitored**: ~2.1 million people
- **Cost per Person**: <$0.01 for sequencing and analysis
- **Taxonomic Resolution**:
  - 184 unique plant ASVs (amplicon sequence variants)
  - 123 unique animal ASVs
  - 98.9% of animal reads mapped to food taxa

### Key Findings

1. **Wealth Signature**: Craft beer ingredients (hops, barley) indicate higher income areas; Southern staples (okra, pecans) indicate lower income
2. **Immigration Signature**: Tropical fruits and Asian legumes (chickpeas, pigeon peas) prevalent in immigrant communities
3. **Geographic Patterns**: Coastal communities consume more wild-caught Atlantic fish; inland urban areas show farmed seafood preference (Spearman ρ = -0.65, p = 0.0024)
4. **Dietary Diversity**: Education, immigration, and income are strong predictors of dietary diversity (VIP scores 1.68-1.95)

## Repository Structure

```
Wastewater/
├── README.md                          # This file
├── ANALYSIS_REPRODUCTION_GUIDE.md     # Detailed analysis workflow documentation
├── run_full_analysis.R                # Main analysis script (run this first!)
├── Rcode.Rmd                          # Original RMarkdown analysis file
├── Rcode_reproduced.R                 # Modular analysis script version
├── data/
│   ├── NCWW_allsamples_animal.rds    # Master phyloseq animal DNA object
│   ├── NCWW_allsamples_trnL.rds      # Master phyloseq plant DNA object
│   ├── metadata_animal.csv            # Sample metadata (37 variables)
│   ├── metadata_plant.csv             # Plant sample metadata
│   ├── animal_taxa.csv                # Animal taxonomy with food annotations
│   ├── plant_taxa.csv                 # Plant taxonomy with food annotations
│   └── [additional supplementary files]
├── figures/
│   ├── Figure_1A_locations.png        # Geographic distribution map
│   ├── Figure_1B_food_composition.png # Food vs non-food composition
│   ├── Figure_2_seasonal_patterns.png # Temporal dietary shifts
│   ├── Figure_3A_PAR_by_location.png  # Plant-to-Animal ratio
│   ├── Figure_3B_Shannon_diversity.png # Shannon diversity index
│   ├── Figure_3C_taxa_richness.png    # Taxa richness
│   ├── Figure_4A_plant_biplot.png     # Plant PCA biplot
│   ├── Figure_5A_fish_PCA.png         # Fish species patterns
│   └── Figure_5B_distance_correlation.png # Distance to coast correlation
├── results/
│   ├── summary_statistics.csv         # Key statistics table
│   ├── diversity_by_location.csv      # Diversity metrics by location
│   └── [statistical results tables]
└── docs/
    ├── ANALYSIS_REPRODUCTION_GUIDE.md # Detailed methodology
    └── QUICK_START.md                 # Quick start guide
```

## Quick Start Guide

### Prerequisites

- **R 4.0+** (download from https://cran.r-project.org/)
- **RStudio** (recommended, download from https://posit.co/download/rstudio-desktop/)
- **Git** (for version control)

### Installation Steps

1. **Clone the repository**:
   ```bash
   git clone https://github.com/YOUR_USERNAME/Wastewater.git
   cd Wastewater
   ```

2. **Open R or RStudio** and install required packages:
   ```R
   # Run this in R console or RStudio
   packages <- c(
     "phyloseq", "vegan", "ggplot2", "ggpubr", "RColorBrewer",
     "tidyr", "dplyr", "ape", "microbiome", "stats", "PathoStat",
     "here", "ggthemes", "tibble", "textshape", "envalysis",
     "vegan3d", "scatterplot3d", "readxl", "pls", "pheatmap",
     "factoextra", "grid", "gridExtra"
   )

   install.packages(packages, dependencies = TRUE)
   ```

3. **Run the analysis**:
   ```bash
   # From command line
   Rscript run_full_analysis.R
   ```

   OR in RStudio:
   ```R
   source("run_full_analysis.R")
   ```

4. **View results**:
   - Figures saved to: `figures/`
   - Summary tables saved to: `results/`

## File Descriptions

### Main Analysis Scripts

**`run_full_analysis.R`** (Recommended)
- Complete analysis pipeline in a single executable script
- Loads data, performs all analyses, generates all figures and tables
- Output: 8 publication-quality figures + 2 CSV summary tables
- Runtime: ~5-10 minutes on standard machine

**`Rcode.Rmd`**
- Original RMarkdown file from manuscript authors
- Includes detailed comments and intermediate outputs
- Can be rendered to HTML report: `rmarkdown::render("Rcode.Rmd")`

**`Rcode_reproduced.R`**
- Alternative modular version organized by analysis type
- Useful for running specific analyses independently

### Data Files

All data files located in `data/` directory:

**Phyloseq Objects (R format, .rds files)**:
- `NCWW_allsamples_animal.rds`: Combined animal DNA data (123 unique taxa, 183 samples)
- `NCWW_allsamples_trnL.rds`: Combined plant DNA data (184 unique taxa, 183 samples)
- `NCWW_allsamples_animal_norm.rds`: Normalized animal reads
- `NCWW_allsamples_trnL_norm.rds`: Normalized plant reads

**Taxonomy Files**:
- `animal_taxa.csv`: Animal taxonomy with IsFood flag (Y/N)
- `plant_taxa.csv`: Plant taxonomy with Streptophyta classification
- `animal_taxa_update.csv`: Updated animal taxonomy
- `plant_taxa_update.csv`: Updated plant taxonomy

**Metadata (37 variables)**:
- `metadata_animal.csv`: Complete sample metadata
  - Sample info: SampleID, Location, SampleType, Date, Month
  - Sequencing: plant_reads, animal_reads, read counts
  - Geographic: Coast_Inland, City_Town, County, lon, lat, DistancetoCoast
  - Demographic: Per_capita_income_k, Population_density, Foreign_born_percent, etc.
  - Health: Obesity_20yr_percent, Diabetes_20yr_percent, FoodInsecure_percent, etc.
- `metadata_plant.csv`: Plant sample metadata

## Analysis Workflow

### Stage 1: Data Quality Assessment
- Assess food vs non-food reads (expected: 98.9% animal food, 76.5% plant food)
- Read count statistics
- Taxonomic resolution validation

### Stage 2: Temporal Analysis (Longitudinal Samples)
- 147 samples from June 2020 - June 2021
- PERMANOVA: County effect (R² = 0.090, p = 0.001), Time effect (R² = 0.019, p = 0.001)
- Seasonal food patterns: summer fruits/vegetables vs. holiday turkey

### Stage 3: Spatial Analysis (2021 Statewide Samples)
- 39 simultaneous samples across 19 municipalities
- PERMANOVA: Coast/Inland effect (R² = 0.064, p = 0.003), City/Town effect (R² = 0.043, p = 0.031)
- Geographic food sourcing patterns

### Stage 4: Diversity Analysis
- Alpha diversity metrics: Shannon index, taxa richness
- Plant-to-Animal Ratio (PAR): higher in urban centers
- PLSR predictors: income, education, density, food insecurity

### Stage 5: Demographic Associations
- Income predictors: Hops (VIP=2.31), Okra (VIP=2.29)
- Cultural signatures: Immigration (tropical fruits), Southern (pecans, black-eyed peas)
- Chi-square test: χ² = 11.3, p = 0.0008

### Stage 6: Geographic Food Sourcing
- Coastal: Atlantic fish (Menhaden, croaker, striped mullet)
- Inland urban: Farmed salmon and tilapia
- Correlation: ρ = -0.65, p = 0.0024

### Stage 7: Validation Against Stool Samples
- Durham wastewater vs. individual stool samples
- Spearman correlation: ρ = 0.64, p < 0.0001
- Validates wastewater as population-level surveillance tool

## Generated Outputs

### Figures (Publication Quality, 300 DPI)

1. **Figure 1A**: Geographic distribution of 19 WWTP sampling sites
   - Shows population served and municipalities
   - Demonstrates nationwide representativeness

2. **Figure 1B**: Food vs Non-Food DNA composition
   - Pie charts showing 98.9% animal food, 76.5% plant food
   - Validates dietary signal in wastewater

3. **Figure 2**: Temporal & Geographic Patterns
   - PCA of seasonal food composition
   - County and month clustering
   - Seasonal dietary shifts

4. **Figure 3A-C**: Diversity Metrics
   - Plant-to-Animal Ratio by location
   - Shannon diversity index
   - Taxa richness
   - Higher in urban, affluent areas

5. **Figure 4A**: Plant Taxa Biplot
   - PCA of plant DNA abundance
   - Shows local vs. specialty foods gradient
   - Cultural dietary patterns

6. **Figure 5A-B**: Fish Consumption Patterns
   - Fish species PCA (coastal vs. farmed)
   - Distance to coast correlation (ρ = -0.65)
   - Demonstrates local food sourcing

### Tables

1. **summary_statistics.csv**
   - Key statistics: samples, municipalities, population, taxa counts

2. **diversity_by_location.csv**
   - Shannon diversity, taxa richness, PAR by location

3. **Statistical results** (generated by analysis):
   - PERMANOVA tables
   - Correlation matrices
   - PLSR variable importance scores

## Statistical Methods

### Data Processing
- **CLR Transformation**: Centered-log ratio for compositional data
- **Normalization**: Read normalization for pooling factors

### Ordination
- **PCA**: Principal Component Analysis for exploratory analysis
- **PERMANOVA**: Permutational MANOVA with Euclidean distance (999 permutations)

### Correlation & Regression
- **Spearman Correlation**: Non-parametric associations with distance to coast
- **PLSR**: Partial Least Squares Regression with 10-fold CV
- **VIP Scores**: Variable Importance in Projection (VIP > 0.8 = important)

### Hypothesis Tests
- **Shapiro-Wilk**: Test normality
- **Chi-square**: Contingency table tests for categorical associations
- **Kruskal-Wallis**: Non-parametric ANOVA for multiple group comparisons

## Key Findings Summary

### Income-Diet Relationships
**Positive (Wealthy areas)**:
- Hops (VIP=2.31) - craft beer ingredient
- Kratom (VIP=2.30) - specialty supplement
- Barley (VIP=2.15) - beer production

**Negative (Lower income)**:
- Okra (VIP=2.29) - Southern staple
- Pecan (VIP=2.25) - regional specialty
- Potato (VIP=2.11) - affordable carbohydrate

### Geographic-Diet Relationships
**Coastal Communities**:
- Wild-caught Atlantic species: Menhaden, mackerels, croaker
- Striped mullet, weakfish
- Higher frequency of local species

**Inland Urban Areas**:
- Farmed Atlantic salmon
- Farmed tilapia
- Imported/commercial species
- Spearman correlation: ρ = -0.65 (p = 0.0024)

### Cultural-Diet Relationships
**Immigration Signature**:
- Tropical fruits: Mango, coconut/palm
- Legumes: Vigna beans, chickpeas, pigeon peas
- Spices indicating Asian cuisine (VIP=1.68)

**Southern/Rural Signature**:
- Vegetables: Celery, carrot, onion
- Nuts: Pecans
- Legumes: Black-eyed peas
- Higher in lower-income areas

## Validation

**Stool vs. Wastewater Comparison (Durham)**:
- 96 overlapping plant taxa across food categories
- Spearman correlation: ρ = 0.64, p < 0.0001
- **Conclusion**: Wastewater DNA accurately reflects community dietary consumption

## Citation

If you use this analysis or reproduce these results, please cite:

```bibtex
@article{Dong2025,
  author = {Dong, M. and others},
  title = {Dietary DNA in municipal wastewater reveals signatures of wealth, immigration, and coastal proximity},
  journal = {Proceedings of the National Academy of Sciences},
  year = {2025},
  volume = {122},
  number = {5},
  pages = {e2313516122}
}
```

## Troubleshooting

### Issue: Package installation fails
**Solution**: Install dependencies separately:
```R
install.packages("phyloseq", dependencies = TRUE)
# For phyloseq, may need BiocManager
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")
```

### Issue: Data files not found
**Solution**: Ensure data files are in `data/` subdirectory and run from repository root

### Issue: Memory error with large phyloseq objects
**Solution**: Reduce sample set by running individual analyses instead of all at once

### Issue: Different output than expected
**Solution**: Check R version (need 4.0+) and that all packages loaded successfully

## Contributing

To contribute improvements or bug fixes:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/improvement`)
3. Commit changes (`git commit -am 'Add improvement'`)
4. Push to branch (`git push origin feature/improvement`)
5. Create Pull Request

## License

This analysis is provided for research and educational purposes. Please refer to the original manuscript for licensing information.

## Questions & Support

For questions about:
- **Analysis methodology**: See ANALYSIS_REPRODUCTION_GUIDE.md
- **Data structure**: See data dictionary in metadata CSV files
- **Code implementation**: Check inline comments in R scripts
- **Manuscript**: Refer to Dong et al. PNAS 2025

---

**Analysis prepared**: November 2025

**Data location**: `data/` directory in this repository

**Code location**: R scripts in repository root

**Figures location**: `figures/` directory (generated by run_full_analysis.R)

**Results location**: `results/` directory (generated by run_full_analysis.R)
