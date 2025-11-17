#!/usr/bin/env Rscript
###############################################################################
# Figure 2A: PCA of CLR-transformed Plant and Animal DNA
# Samples from three locations: Beaufort, Charlotte, and Greenville
# PC3 vs PC4 highlighting temporal (monthly) separation
###############################################################################

set.seed(001)

# Load required packages
cat("Loading packages...\n")
packages <- c("phyloseq", "ggplot2", "dplyr", "lubridate", "microbiome", "ggrepel")
invisible(lapply(packages, library, character.only = TRUE))

# Set paths
data_path <- "C:/Users/rache/Box/project_davidlab/LAD_LAB_Personnel/Rachel_Q/Code and Data/Data for code optimization_Do not submit/"
output_path <- "C:/Users/rache/Box/project_davidlab/LAD_LAB_Personnel/Rachel_Q/Code and Data/code/figures/"

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

###############################################################################
# Load and Prepare Data
###############################################################################

cat("Loading data...\n")
NCWW_allsamples_trnL <- readRDS(paste0(data_path, "NCWW_allsamples_trnL.rds"))
NCWW_allsamples_animal <- readRDS(paste0(data_path, "NCWW_allsamples_animal.rds"))

# Create plant subset: food plants only, from three target locations
ps_plant <- subset_taxa(NCWW_allsamples_trnL, phylum == "Streptophyta" & IsFood == "Y")
ps_plant <- prune_taxa(taxa_sums(ps_plant) > 0, ps_plant)

# Create animal subset: food animals only, from three target locations
ps_animal <- subset_taxa(NCWW_allsamples_animal, IsFood == "Y" & CommonName != "Emu")
ps_animal <- prune_taxa(taxa_sums(ps_animal) > 0, ps_animal)

# Filter to specific LOCATIONS: Beaufort, Charlotte 1-4, Greenville
# (BEST configuration determined by PCA variance testing)
target_locations <- c("Beaufort", "Charlotte 1", "Charlotte 2", "Charlotte 3", "Charlotte 4", "Greenville")

ps_plant_filtered <- subset_samples(ps_plant, Location %in% target_locations)
ps_animal_filtered <- subset_samples(ps_animal, Location %in% target_locations)

# Apply sample removal list (problematic samples to exclude)
samples_to_remove <- c("NCWW121720_1","NCWW121720_2","NCWW121720_10",
                       "NCWW121720_56","NCWW121720_52","NCWW121720_18",
                       "NCWW121720_21","NCWW022621_16","NCWW121720_27",
                       "NCWW121720_6")

ps_plant_filtered <- prune_samples(!(sample_names(ps_plant_filtered) %in% samples_to_remove), ps_plant_filtered)
ps_animal_filtered <- prune_samples(!(sample_names(ps_animal_filtered) %in% samples_to_remove), ps_animal_filtered)

# Remove taxa with zero sums after filtering
ps_plant_filtered <- prune_taxa(taxa_sums(ps_plant_filtered) > 0, ps_plant_filtered)
ps_animal_filtered <- prune_taxa(taxa_sums(ps_animal_filtered) > 0, ps_animal_filtered)

cat("Plant samples after filtering:", nsamples(ps_plant_filtered), "\n")
cat("Animal samples after filtering:", nsamples(ps_animal_filtered), "\n")
cat("Plant taxa after filtering:", ntaxa(ps_plant_filtered), "\n")
cat("Animal taxa after filtering:", ntaxa(ps_animal_filtered), "\n\n")

# Merge plant and animal data
# Create phyloseq objects without trees for merging
ps_plant_notree <- phyloseq(otu_table(ps_plant_filtered), tax_table(ps_plant_filtered), sample_data(ps_plant_filtered))
ps_animal_notree <- phyloseq(otu_table(ps_animal_filtered), tax_table(ps_animal_filtered), sample_data(ps_animal_filtered))

ps_combined <- merge_phyloseq(ps_plant_notree, ps_animal_notree)
cat("Combined samples (plant + animal):", nsamples(ps_combined), "\n")
cat("Combined taxa:", ntaxa(ps_combined), "\n\n")

###############################################################################
# CLR Transformation
###############################################################################

cat("Applying CLR transformation...\n")

# Apply CLR transformation (microbiome package handles zero-inflation automatically)
ps_combined_clr <- microbiome::transform(ps_combined, "clr")

cat("CLR transformation completed\n\n")

###############################################################################
# Prepare Metadata
###############################################################################

metadata <- data.frame(sample_data(ps_combined_clr))
metadata$SamplingDate <- as.Date(metadata$Date, format = "%m/%d/%y")
metadata$Month <- month(metadata$SamplingDate)
metadata$MonthName <- month.abb[metadata$Month]

cat("Sampling date range:", min(metadata$SamplingDate), "to", max(metadata$SamplingDate), "\n")
cat("Months represented:", paste(unique(sort(metadata$Month)), collapse=", "), "\n")
cat("Month names:", paste(unique(metadata$MonthName[order(metadata$Month)]), collapse=", "), "\n")
cat("Locations:", paste(unique(metadata$Location), collapse=", "), "\n\n")

###############################################################################
# PCA Analysis
###############################################################################

cat("Performing PCA on CLR-transformed combined data...\n")

otu_data <- as.data.frame(otu_table(ps_combined_clr))
pca_result <- prcomp(otu_data, scale. = FALSE)

# Calculate variance explained
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2) * 100

cat("Variance explained by first 5 components:\n")
for(i in 1:5) {
  cat(sprintf("  PC%d: %.2f%%\n", i, var_explained[i]))
}

cat("\nPC3 variance explained: ", round(var_explained[3], 2), "%\n", sep="")
cat("PC4 variance explained: ", round(var_explained[4], 2), "%\n\n", sep="")

# Extract PC scores
pca_scores <- data.frame(
  PC1 = pca_result$x[,1],
  PC2 = pca_result$x[,2],
  PC3 = pca_result$x[,3],
  PC4 = pca_result$x[,4]
)

# Combine with metadata
plot_data <- cbind(metadata, pca_scores)

###############################################################################
# Create Figure 2A: PC3 vs PC4 colored by month, shaped by location
###############################################################################

cat("Creating Figure 2A visualization...\n")

# Define colors for months (rainbow palette)
month_colors <- c(
  "5" = "#FF6B6B",    # May - Red
  "6" = "#FFA500",    # June - Orange
  "7" = "#FFD700",    # July - Gold
  "8" = "#90EE90",    # August - Light green
  "9" = "#00FF00",    # September - Green
  "10" = "#00CED1",   # October - Cyan
  "11" = "#4169E1",   # November - Blue
  "12" = "#9370DB"    # December - Purple
)

# Define shapes for locations
shape_values <- c(
  "Beaufort" = 16,
  "Charlotte 1" = 17,
  "Charlotte 2" = 17,
  "Charlotte 3" = 17,
  "Charlotte 4" = 17,
  "Greenville" = 18
)

# Create plot
fig2a <- ggplot(plot_data, aes(x = PC3, y = PC4, color = as.factor(Month), shape = Location)) +
  geom_point(size = 4, alpha = 0.7, stroke = 1) +
  scale_color_manual(
    name = "Month",
    values = month_colors,
    labels = c("5" = "May", "6" = "June", "7" = "July", "8" = "August",
               "9" = "September", "10" = "October", "11" = "November", "12" = "December")
  ) +
  scale_shape_manual(
    name = "Location",
    values = shape_values
  ) +
  labs(
    title = "Figure 2A: Temporal Patterns in Food DNA (Plant & Animal Combined)",
    subtitle = "PCA of CLR-transformed food DNA from Beaufort, Charlotte, and Greenville",
    x = paste0("PC3 (", round(var_explained[3], 1), "% variance)"),
    y = paste0("PC4 (", round(var_explained[4], 1), "% variance)")
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11, face = "italic", hjust = 0),
    axis.title = element_text(size = 11, face = "bold"),
    axis.text = element_text(size = 10),
    legend.position = "right",
    panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

ggsave(paste0(output_path, "Figure_2A_temporal_PCA.png"), fig2a, width = 12, height = 8, dpi = 300)
cat("✓ Figure 2A saved\n\n")

###############################################################################
# Summary Statistics
###############################################################################

cat("═══════════════════════════════════════════════════════════════════\n")
cat("Figure 2A: Temporal PCA Analysis Summary\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("Data Summary:\n")
cat("  Combined samples (plant + animal):", nsamples(ps_combined_clr), "\n")
cat("  Combined taxa:", ntaxa(ps_combined_clr), "\n")
cat("  Locations: Beaufort, Charlotte (4 variants), Greenville\n")
cat("  Sampling period: May - December\n\n")

cat("Sample counts by location:\n")
print(table(plot_data$Location))

cat("\n\nSample counts by month:\n")
print(table(plot_data$Month))

cat("\n\nPCA Analysis Results:\n")
cat("  Transformation: CLR (centered log-ratio)\n")
cat("  Scaling: No scaling (data already CLR-transformed)\n")
cat("  Locations: Beaufort, Charlotte 1-4, Greenville\n")
cat("  Sample removal: Yes (10 problematic samples excluded)\n")
cat("  Configuration: BEST (determined by variance optimization testing)\n")
cat("  Number of components: 4\n\n")

cat("Variance Explained:\n")
for(i in 1:4) {
  cat(sprintf("  PC%d: %6.2f%%\n", i, var_explained[i]))
}

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("Script completed successfully!\n")
cat("═══════════════════════════════════════════════════════════════════\n")
