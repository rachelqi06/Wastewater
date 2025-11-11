###############################################################################
# Wastewater Dietary DNA Analysis - Reproduced Analysis Script
# Based on: "Dietary DNA in municipal wastewater reveals signatures of wealth,
#           immigration, and coastal proximity"
# Analysis Date: 2025-11-04
###############################################################################

# Set working directory
setwd("C:/Users/rache/Box/project_davidlab/LAD_LAB_Personnel/Rachel_Q/Code and Data")

# Set seed for reproducibility
set.seed(001)

###############################################################################
# 1. LOAD REQUIRED PACKAGES
###############################################################################

packages <- c("phyloseq", "vegan", "ggplot2", "ggpubr", "RColorBrewer",
              "tidyr", "dplyr", "ape", "microbiome", "stats", "PathoStat",
              "here", "ggthemes", "tibble", "textshape", "envalysis",
              "vegan3d", "scatterplot3d", "readxl", "pls", "pheatmap",
              "factoextra")

# Install missing packages
new_packages <- packages[!(packages %in% installed.packages()[,"Package"])]
if(length(new_packages) > 0) {
  install.packages(new_packages, dependencies = TRUE)
}

# Load all packages
invisible(lapply(packages, library, character.only = TRUE))

###############################################################################
# 2. LOAD PHYLOSEQ OBJECTS
###############################################################################

cat("Loading phyloseq objects...\n")

# Load animal and plant data
data_path <- "Data for code optimization_Do not submit/"
NCWW_allsamples_animal <- readRDS(paste0(data_path, "NCWW_allsamples_animal.rds"))
NCWW_allsamples_trnL <- readRDS(paste0(data_path, "NCWW_allsamples_trnL.rds"))

cat("✓ Phyloseq objects loaded\n")
cat("  Animal samples:", nsamples(NCWW_allsamples_animal), "\n")
cat("  Plant samples:", nsamples(NCWW_allsamples_trnL), "\n")

###############################################################################
# 3. DATA QUALITY ASSESSMENT (Fig 1C, Table S3)
###############################################################################

cat("\n=== DATA QUALITY ASSESSMENT ===\n")

# Animal reads - food vs non-food
food_animal <- subset_taxa(NCWW_allsamples_animal, IsFood == "Y")
total_animal_reads <- sum(taxa_sums(NCWW_allsamples_animal))
food_animal_reads <- sum(taxa_sums(food_animal))
nonfood_animal_reads <- total_animal_reads - food_animal_reads
animal_food_percent <- (food_animal_reads / total_animal_reads) * 100

cat("\nANIMAL SEQUENCES:\n")
cat("  Total animal reads:", total_animal_reads, "\n")
cat("  Food animal reads:", food_animal_reads, "\n")
cat("  Non-food animal reads:", nonfood_animal_reads, "\n")
cat("  % Food reads:", round(animal_food_percent, 2), "%\n")

# Plant reads - food vs non-food
food_plant <- subset_taxa(NCWW_allsamples_trnL, phylum == "Streptophyta")
total_plant_reads <- sum(taxa_sums(NCWW_allsamples_trnL))
food_plant_reads <- sum(taxa_sums(food_plant))
nonfood_plant_reads <- total_plant_reads - food_plant_reads
plant_food_percent <- (food_plant_reads / total_plant_reads) * 100

cat("\nPLANT SEQUENCES:\n")
cat("  Total plant reads:", total_plant_reads, "\n")
cat("  Food plant reads:", food_plant_reads, "\n")
cat("  Non-food plant reads:", nonfood_plant_reads, "\n")
cat("  % Food reads:", round(plant_food_percent, 2), "%\n")

###############################################################################
# 4. GENERATE SUBSETS FOR ANALYSIS
###############################################################################

cat("\nGenerating analysis subsets...\n")

# Remove low abundance samples
samples_to_remove <- c("NCWW121720_1","NCWW121720_2","NCWW121720_10",
                        "NCWW121720_56","NCWW121720_52","NCWW121720_18",
                        "NCWW121720_21","NCWW022621_16","NCWW121720_27",
                        "NCWW121720_6")

# Create food subsets
ps_animal <- subset_taxa(NCWW_allsamples_animal, IsFood == "Y" & CommonName != "Emu")
ps_animal_clr <- microbiome::transform(ps_animal, "clr")

ps_plant <- subset_taxa(NCWW_allsamples_trnL, phylum == "Streptophyta")
ps_plant_clr <- microbiome::transform(ps_plant, "clr")

# Longitudinal sample set (seasonal)
NCWW_Seasonal_plant <- subset_samples(ps_plant, County %in% c("carteret","mecklenburg","pitt"))
NCWW_Seasonal_plant <- prune_samples(!(sample_names(NCWW_Seasonal_plant) %in% samples_to_remove),
                                      NCWW_Seasonal_plant)
NCWW_Seasonal_plant <- prune_taxa(taxa_sums(NCWW_Seasonal_plant) > 0, NCWW_Seasonal_plant)
NCWW_Seasonal_plant_clr <- microbiome::transform(NCWW_Seasonal_plant, "clr")

NCWW_Seasonal_animal <- subset_samples(ps_animal, County %in% c("carteret","mecklenburg","pitt"))
NCWW_Seasonal_animal <- prune_samples(!(sample_names(NCWW_Seasonal_animal) %in% samples_to_remove),
                                       NCWW_Seasonal_animal)
NCWW_Seasonal_animal <- prune_taxa(taxa_sums(NCWW_Seasonal_animal) > 0, NCWW_Seasonal_animal)
NCWW_Seasonal_animal_clr <- microbiome::transform(NCWW_Seasonal_animal, "clr")

NCWW_Seasonal_fish <- subset_taxa(NCWW_Seasonal_animal, class == "Actinopteri")
NCWW_Seasonal_fish_clr <- subset_taxa(NCWW_Seasonal_animal_clr, class == "Actinopteri")

# Merge seasonal plant and animal
ps_notree_1 <- phyloseq(NCWW_Seasonal_animal@otu_table, NCWW_Seasonal_animal@tax_table,
                         NCWW_Seasonal_animal@sam_data)
ps_notree_2 <- phyloseq(NCWW_Seasonal_plant@otu_table, NCWW_Seasonal_plant@tax_table,
                         NCWW_Seasonal_plant@sam_data)
NCWW_Seasonal_food <- merge_phyloseq(ps_notree_1, ps_notree_2)

ps_notree_3 <- phyloseq(NCWW_Seasonal_animal_clr@otu_table, NCWW_Seasonal_animal_clr@tax_table,
                         NCWW_Seasonal_animal_clr@sam_data)
ps_notree_4 <- phyloseq(NCWW_Seasonal_plant_clr@otu_table, NCWW_Seasonal_plant_clr@tax_table,
                         NCWW_Seasonal_plant_clr@sam_data)
NCWW_Seasonal_food_clr <- merge_phyloseq(ps_notree_3, ps_notree_4)

# Spatial sample set (2021 samples)
NCWW_2021_animal <- subset_samples(ps_animal, Plate == "NCWW2021")
NCWW_2021_animal <- prune_taxa(taxa_sums(NCWW_2021_animal) > 0, NCWW_2021_animal)
NCWW_2021_animal_clr <- microbiome::transform(NCWW_2021_animal, "clr")
NCWW_2021_fish <- subset_taxa(NCWW_2021_animal, class == "Actinopteri")
NCWW_2021_fish_clr <- subset_taxa(NCWW_2021_animal_clr, class == "Actinopteri")

NCWW_2021_plant <- subset_samples(ps_plant, Plate == "NCWW2021")
NCWW_2021_plant <- prune_taxa(taxa_sums(NCWW_2021_plant) > 0, NCWW_2021_plant)
NCWW_2021_plant_clr <- microbiome::transform(NCWW_2021_plant, "clr")

# Merge 2021 spatial data
ps_notree_5 <- phyloseq(NCWW_2021_animal@otu_table, NCWW_2021_animal@tax_table,
                         NCWW_2021_animal@sam_data)
ps_notree_6 <- phyloseq(NCWW_2021_plant@otu_table, NCWW_2021_plant@tax_table,
                         NCWW_2021_plant@sam_data)
NCWW_2021 <- merge_phyloseq(ps_notree_5, ps_notree_6)

ps_notree_7 <- phyloseq(NCWW_2021_animal_clr@otu_table, NCWW_2021_animal_clr@tax_table,
                         NCWW_2021_animal_clr@sam_data)
ps_notree_8 <- phyloseq(NCWW_2021_plant_clr@otu_table, NCWW_2021_plant_clr@tax_table,
                         NCWW_2021_plant_clr@sam_data)
NCWW_2021_clr <- merge_phyloseq(ps_notree_7, ps_notree_8)

cat("✓ Data subsets created\n")
cat("  Longitudinal samples (seasonal):", nsamples(NCWW_Seasonal_food), "\n")
cat("  Spatial samples (2021):", nsamples(NCWW_2021), "\n")

###############################################################################
# 5. STATISTICS: PERMANOVA (Table S6, S13)
###############################################################################

cat("\n=== PERMANOVA ANALYSIS ===\n")

# Temporal and spatial pattern of food
ps <- NCWW_Seasonal_food_clr
df <- data.frame(ps@otu_table)
metadata <- data.frame(ps@sam_data)
metadata$Month_Date <- as.Date(metadata$Date, format = "%m/%d")
metadata$Month_Date <- as.numeric(metadata$Month_Date)

permanova_seasonal <- adonis2(df ~ County + Month_Date, data = metadata,
                               method = "euclidean", by = "terms", permutations = 999)
cat("\nSeasonal Food Composition (County + Month):\n")
print(permanova_seasonal)

# Spatial fish consumption
ps <- NCWW_2021_fish_clr
df <- data.frame(ps@otu_table)
metadata <- data.frame(ps@sam_data)

permanova_fish <- adonis2(df ~ Coast_Inland + City_Town, data = metadata,
                           method = "euclidean", by = "terms", permutations = 999)
cat("\nSpatial Fish Composition (Coast/Inland + City/Town):\n")
print(permanova_fish)

###############################################################################
# 6. DIVERSITY METRICS (Fig 3B,D; Supplementary Fig 6A)
###############################################################################

cat("\n=== DIVERSITY ANALYSIS ===\n")

NCWW_2021_int <- NCWW_2021
otu_table(NCWW_2021_int) <- round(otu_table(NCWW_2021_int))
alpha_diversity <- estimate_richness(NCWW_2021_int, measures=c("Observed", "Shannon"))

# Plant-animal ratio
ps_phylum <- tax_glom(NCWW_2021, taxrank = "phylum")
ps_phylum <- microbiome::transform(ps_phylum, "compositional")
otu_df <- as.data.frame(otu_table(ps_phylum))
tax_df <- as.data.frame(tax_table(ps_phylum))
colnames(otu_df) <- tax_df$phylum
otu_df$PAR <- otu_df$Streptophyta / otu_df$Chordata

df_diversity <- data.frame(otu_df, alpha_diversity)
df_diversity$Location <- ps_phylum@sam_data$Location

cat("\nDiversity Metrics Summary:\n")
cat("  Number of locations:", length(unique(df_diversity$Location)), "\n")
cat("  Mean Shannon diversity:", round(mean(df_diversity$Shannon, na.rm=T), 3), "\n")
cat("  Mean taxa richness:", round(mean(df_diversity$Observed, na.rm=T), 3), "\n")
cat("  Mean Plant-to-Animal Ratio:", round(mean(df_diversity$PAR, na.rm=T), 3), "\n")

# Statistical tests
cat("\nShapiro-Wilk normality test:\n")
cat("  Shannon diversity p-value:", round(shapiro.test(df_diversity$Shannon)$p.value, 4), "\n")
cat("  Observed (richness) p-value:", round(shapiro.test(df_diversity$Observed)$p.value, 4), "\n")
cat("  PAR p-value:", round(shapiro.test(df_diversity$PAR)$p.value, 4), "\n")

###############################################################################
# 7. SAVE RESULTS FOR VISUALIZATION
###############################################################################

cat("\n=== SAVING ANALYSIS RESULTS ===\n")

# Save diversity dataframe
saveRDS(df_diversity, paste0(data_path, "df_diversity.rds"))

# Save metadata and sample information
sample_summary <- data.frame(
  Total_Samples = nsamples(NCWW_allsamples_animal),
  Animal_Samples = nsamples(ps_animal),
  Plant_Samples = nsamples(ps_plant),
  Longitudinal_Samples = nsamples(NCWW_Seasonal_food),
  Spatial_2021_Samples = nsamples(NCWW_2021)
)

print(sample_summary)

cat("\n✓ Analysis complete! Results ready for visualization\n")

###############################################################################
# 8. SUMMARY STATISTICS TABLE
###############################################################################

cat("\n=== SUMMARY STATISTICS ===\n\n")

# Table 1: Sample Distribution
cat("TABLE 1: SAMPLE DISTRIBUTION\n")
cat("─────────────────────────────────\n")
cat(sprintf("Total wastewater samples analyzed: %d\n", nsamples(NCWW_allsamples_animal)))
cat(sprintf("Total municipalities: 19\n"))
cat(sprintf("Population monitored: ~2.1 million\n"))
cat(sprintf("Longitudinal samples (temporal): %d\n", nsamples(NCWW_Seasonal_food)))
cat(sprintf("Spatial samples (2021): %d\n\n", nsamples(NCWW_2021)))

# Table 2: Read counts and quality
cat("TABLE 2: SEQUENCE READ STATISTICS\n")
cat("─────────────────────────────────\n")
cat(sprintf("Mean plant reads per sample: %.0f\n", mean(sample_sums(NCWW_2021_plant))))
cat(sprintf("Mean animal reads per sample: %.0f\n", mean(sample_sums(NCWW_2021_animal))))
cat(sprintf("Food animal reads: %.1f%%\n", animal_food_percent))
cat(sprintf("Food plant reads: %.1f%%\n\n", plant_food_percent))

# Table 3: Taxonomic resolution
cat("TABLE 3: TAXONOMIC RESOLUTION\n")
cat("─────────────────────────────────\n")
cat(sprintf("Unique plant ASVs (trnL): %d\n", ntaxa(ps_plant)))
cat(sprintf("Unique animal ASVs (12SV5): %d\n", ntaxa(ps_animal)))
cat(sprintf("Unique fish taxa: %d\n\n", ntaxa(NCWW_2021_fish)))

cat("Analysis completed at:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
