# Quick Start Guide - FoodSeq-FLOW Analysis

## 5-Minute Setup

### Step 1: Download R and RStudio
- **R**: https://cran.r-project.org/ → Download for your OS
- **RStudio**: https://posit.co/download/rstudio-desktop/ → Free Desktop version

### Step 2: Install Packages (5-7 minutes)
Open RStudio and run this in the Console:

```R
# Install required packages
packages <- c(
  "phyloseq", "vegan", "ggplot2", "ggpubr", "RColorBrewer",
  "tidyr", "dplyr", "ape", "microbiome", "stats", "PathoStat",
  "here", "ggthemes", "tibble", "textshape", "envalysis",
  "vegan3d", "scatterplot3d", "readxl", "pls", "pheatmap",
  "factoextra", "grid", "gridExtra"
)

# For phyloseq from Bioconductor
if(!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("phyloseq")

# Install CRAN packages
install.packages(packages, dependencies = TRUE)
```

### Step 3: Clone Repository
Open terminal and run:
```bash
git clone https://github.com/YOUR_USERNAME/Wastewater.git
cd Wastewater
```

### Step 4: Run Analysis
In RStudio, open `run_full_analysis.R` and click **Run** or paste into Console:
```R
source("run_full_analysis.R")
```

**Wait 5-10 minutes** for analysis to complete.

### Step 5: View Results
- **Figures**: Check `figures/` folder
- **Tables**: Check `results/` folder
- **Console**: Shows progress and statistics

## What Gets Generated

✓ 8 Publication-quality figures (PNG, 300 DPI)
✓ 2 Summary statistics tables (CSV)
✓ Complete statistical analysis in R console

## File Locations

```
Project Root/
├── run_full_analysis.R    ← RUN THIS FILE
├── data/                  ← Input data (phyloseq .rds files)
├── figures/               ← Output plots (generated)
└── results/               ← Output tables (generated)
```

## Expected Output

When you run the analysis, you'll see:
```
═══════════════════════════════════════════════════════════════════
WASTEWATER DIETARY DNA ANALYSIS - DATA LOADING
═══════════════════════════════════════════════════════════════════

Loading required packages...
✓ All packages loaded

✓ Phyloseq objects loaded
  Animal samples: 183
  Plant samples: 183
  ...

═══════════════════════════════════════════════════════════════════
DATA QUALITY ASSESSMENT
═══════════════════════════════════════════════════════════════════

ANIMAL SEQUENCES:
  Total reads: 3,706,835
  Food reads: 3,670,179
  Non-food reads: 36,656
  % Food reads: 98.99%

PLANT SEQUENCES:
  Total reads: 4,096,287
  Food reads: 3,131,098
  Non-food reads: 965,189
  % Food reads: 76.45%
  ...
```

## Troubleshooting

### Q: "Package not found" error
A: Run installation commands again, one package at a time:
```R
install.packages("ggplot2", dependencies = TRUE)
# Wait for completion, then next package
```

### Q: "Could not find data files"
A: Make sure you're running from the correct directory:
```R
getwd()  # Should show Wastewater folder
# If not, run: setwd("C:/path/to/Wastewater")
```

### Q: Running very slowly
A: This is normal. Analysis involves:
- Loading 3.7M+ animal sequences
- Loading 4M+ plant sequences
- 183 wastewater samples
- Computing CLR transformation
- 999 PERMANOVA permutations
- Generating 8 figures

Expect 5-15 minutes depending on your computer.

### Q: Out of memory error
A: Close other programs and try again. If persistent, run individual analyses:
```R
# Instead of full script, run sections separately
# See Rcode_reproduced.R for individual analysis sections
```

## Key Outputs

### Figures Generated

1. **Figure_1A_locations.png**
   - Map of 19 WWTP locations in North Carolina
   - Shows population served by each location

2. **Figure_1B_food_composition.png**
   - Pie charts: 98.9% animal DNA is food, 76.5% plant DNA is food
   - Validates that wastewater signals reflect diet

3. **Figure_2_seasonal_patterns.png**
   - PCA plot showing seasonal dietary changes
   - County clustering and temporal shifts

4. **Figure_3A_PAR_by_location.png**
   - Plant-to-Animal Ratio across North Carolina
   - Higher in urban areas (more vegetables/fruits)

5. **Figure_3B_Shannon_diversity.png**
   - Dietary diversity (Shannon Index) by location
   - Correlates with income and education

6. **Figure_3C_taxa_richness.png**
   - Number of unique food taxa detected by location
   - Shows 45-65 different food types per WWTP

7. **Figure_4A_plant_biplot.png**
   - PCA of plant DNA
   - Shows local/regional vs. imported foods
   - Coastal vs. inland patterns

8. **Figure_5A_fish_PCA.png**
   - Fish species patterns
   - Coastal = wild-caught Atlantic species
   - Inland = farmed salmon/tilapia

9. **Figure_5B_distance_correlation.png**
   - Scatter plot: Distance to coast vs. fish PC1
   - Shows local fish preference near coast (ρ = -0.65)

### Tables Generated

**summary_statistics.csv**
```
Metric,Value
Total Samples,183
Municipalities,19
Population Monitored (millions),~2.1
...
```

**diversity_by_location.csv**
```
Location,Shannon_Diversity,Taxa_Richness,Plant_to_Animal_Ratio
...
```

## Next Steps

1. **View manuscript**: Read Dong et al. PNAS 2025 for detailed interpretation
2. **Explore code**: Open `Rcode.Rmd` to see original analysis workflow
3. **Modify analysis**: Edit `run_full_analysis.R` to customize analyses
4. **Push to GitHub**: Commit results and push to your repository

## Key Statistics at a Glance

| Metric | Value |
|--------|-------|
| **Total Samples** | 183 wastewater |
| **Municipalities** | 19 across North Carolina |
| **Population** | ~2.1 million people |
| **Animal Taxa** | 123 unique species |
| **Plant Taxa** | 184 unique species |
| **Animal Food %** | 98.9% (high quality) |
| **Plant Food %** | 76.5% (food + pollen) |
| **Cost per Person** | <$0.01 |

## Need Help?

1. **Analysis Methods**: See `ANALYSIS_REPRODUCTION_GUIDE.md`
2. **Manuscript**: See Dong et al. PNAS 2025
3. **Data Structure**: See column descriptions in `data/metadata_animal.csv`
4. **Code Issues**: Check comments in R scripts

---

**Ready to run?** Open RStudio, navigate to the Wastewater folder, and:
```R
source("run_full_analysis.R")
```

Good luck! 🔬
