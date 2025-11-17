#!/usr/bin/env Rscript
###############################################################################
# Figure 2B: Top 20 Plant Taxa Most Strongly Associated with Sampling Time
# Using PLS Regression (first component) on CLR-transformed plant abundances
# Colors indicate each taxon's typical growing/harvest season
###############################################################################

set.seed(001)

# Load required packages
cat("Loading packages...\n")
packages <- c("phyloseq", "pls", "ggplot2", "dplyr", "lubridate", "microbiome")
invisible(lapply(packages, library, character.only = TRUE))

# Set paths
data_path <- "C:/Users/rache/Box/project_davidlab/LAD_LAB_Personnel/Rachel_Q/Code and Data/Data for code optimization_Do not submit/"
output_path <- "C:/Users/rache/Box/project_davidlab/LAD_LAB_Personnel/Rachel_Q/Code and Data/code/figures/"

if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

###############################################################################
# Load Data
###############################################################################

cat("Loading data...\n")
NCWW_allsamples_trnL <- readRDS(paste0(data_path, "NCWW_allsamples_trnL.rds"))
NCWW_allsamples_animal <- readRDS(paste0(data_path, "NCWW_allsamples_animal.rds"))

# Create seasonal subset
ps_plant <- subset_taxa(NCWW_allsamples_trnL, phylum == "Streptophyta")
ps_plant_clr <- microbiome::transform(ps_plant, "clr")

samples_to_remove <- c("NCWW121720_1","NCWW121720_2","NCWW121720_10",
                       "NCWW121720_56","NCWW121720_52","NCWW121720_18",
                       "NCWW121720_21","NCWW022621_16","NCWW121720_27",
                       "NCWW121720_6")

NCWW_Seasonal_plant <- subset_samples(ps_plant, County %in% c("carteret","mecklenburg","pitt"))
NCWW_Seasonal_plant <- prune_samples(!(sample_names(NCWW_Seasonal_plant) %in% samples_to_remove),
                                      NCWW_Seasonal_plant)
NCWW_Seasonal_plant <- prune_taxa(taxa_sums(NCWW_Seasonal_plant) > 0, NCWW_Seasonal_plant)
NCWW_Seasonal_plant_clr <- microbiome::transform(NCWW_Seasonal_plant, "clr")

###############################################################################
# Prepare data for PLS regression
###############################################################################

cat("Preparing data for PLS regression...\n")

# Filter to ONLY plant taxa
NCWW_Seasonal_plant_only <- subset_taxa(NCWW_Seasonal_plant_clr, phylum == "Streptophyta")

# Get metadata with Month
fig2b_metadata <- data.frame(sam_data(NCWW_Seasonal_plant_only))
fig2b_metadata$Month <- month(as.Date(fig2b_metadata$Date, format = "%m/%d/%y"))

# Get OTU table
otu_df <- as.data.frame(otu_table(NCWW_Seasonal_plant_only))
otu_ids <- colnames(otu_df)

# Get common names mapping
tax_tab <- tax_table(NCWW_Seasonal_plant_only)
common_names_map <- as.character(tax_tab[, "CommonName"])
names(common_names_map) <- rownames(tax_tab)

# Remove zero-abundance taxa
otu_df <- otu_df[, colSums(otu_df) > 0]
otu_ids <- colnames(otu_df)
otu_df[is.na(otu_df)] <- 0

cat("  Samples:", nrow(otu_df), "\n")
cat("  Taxa retained:", ncol(otu_df), "\n")
cat("  Sampling months:", min(fig2b_metadata$Month, na.rm=TRUE), "-", max(fig2b_metadata$Month, na.rm=TRUE), "\n\n")

# Create PLSR data frame
month_numeric <- as.numeric(fig2b_metadata$Month)
plsr_data <- cbind(data.frame(Month = month_numeric), otu_df)

###############################################################################
# Perform PLS Regression on first component
###############################################################################

cat("Performing PLS regression...\n")

# Fit PLSR model with first component
plsr_model <- plsr(Month ~ ., data = plsr_data, scale = TRUE, ncomp = 1, na.action = na.omit)

# Get loadings from first component
plsr_loadings <- plsr_model$loadings[, 1]

# Map to common names
all_common_names <- common_names_map[names(plsr_loadings)]
all_common_names[is.na(all_common_names)] <- names(plsr_loadings)[is.na(all_common_names)]

###############################################################################
# Select Top 20 Taxa by Absolute Loading Magnitude
###############################################################################

# Calculate absolute loadings
abs_loadings <- abs(plsr_loadings)

# Find top 20 by magnitude
sorted_idx <- order(abs_loadings, decreasing = TRUE)
top_20_idx <- sorted_idx[1:20]

# Get names and loadings for top 20
top_taxa_otu_ids <- names(plsr_loadings[top_20_idx])
top_taxa_common_names <- all_common_names[top_taxa_otu_ids]
top_taxa_loadings <- plsr_loadings[top_taxa_otu_ids]

# Sort by signed loading (positive to negative) for visualization
viz_order <- order(top_taxa_loadings, decreasing = TRUE)
top_taxa_common_names <- top_taxa_common_names[viz_order]
top_taxa_loadings <- top_taxa_loadings[viz_order]

###############################################################################
# Define Growing/Harvest Seasons
###############################################################################

season_map <- c(
  "Coriander" = "Year-round", "Pecan" = "Fall/Winter", "Black pepper" = "Year-round",
  "Canola" = "Fall/Winter", "Rice" = "Fall", "Cabbage" = "Fall/Winter",
  "Melon" = "Summer", "Olive" = "Fall/Winter", "Camellia" = "Fall/Winter",
  "Grape" = "Summer", "Asparagus" = "Spring/Summer", "Kiwifruit" = "Fall/Winter",
  "Barley" = "Summer", "Pistachio" = "Fall", "Pineapple" = "Year-round",
  "Mango" = "Summer", "Walnut" = "Fall", "Okra" = "Summer",
  "Blueberry" = "Summer", "Cocoa" = "Year-round", "Corn" = "Fall",
  "Wheat" = "Summer", "Lettuce" = "Year-round", "Tomato" = "Summer",
  "Citrus" = "Fall/Winter", "Apple" = "Fall", "Carrot" = "Fall/Winter"
)

# Assign seasons to top 20 taxa
taxa_seasons <- sapply(top_taxa_common_names, function(x) {
  for (season_name in names(season_map)) {
    if (grepl(season_name, x, ignore.case = TRUE)) {
      return(season_map[[season_name]])
    }
  }
  # Default classification
  if (grepl("berry|melon|grape|mango", x, ignore.case = TRUE)) {
    return("Summer")
  } else if (grepl("pecan|citrus|cabbage|carrot|onion|celery|pea|walnut|apple", x, ignore.case = TRUE)) {
    return("Fall/Winter")
  } else if (grepl("asparagus", x, ignore.case = TRUE)) {
    return("Spring/Summer")
  } else {
    return("Year-round")
  }
})

###############################################################################
# Create Figure 2B
###############################################################################

cat("Creating Figure 2B bar plot...\n")

# Create data frame for plotting
fig2b_data <- data.frame(
  CommonName = as.character(top_taxa_common_names),
  Loading = as.numeric(top_taxa_loadings),
  Season = as.character(taxa_seasons),
  stringsAsFactors = FALSE
)

# Order by loading magnitude (already done above)
fig2b_data$CommonName <- factor(fig2b_data$CommonName,
                                 levels = rev(fig2b_data$CommonName))

# Create bar chart
fig2b <- ggplot(fig2b_data, aes(x = Loading, y = CommonName, fill = Season)) +
  geom_col(color = "black", linewidth = 0.3, alpha = 0.85) +
  scale_fill_manual(
    values = c(
      "Summer" = "#F39C12",           # Orange
      "Fall" = "#E67E22",              # Dark orange
      "Fall/Winter" = "#E74C3C",       # Red
      "Spring/Summer" = "#27AE60",     # Green
      "Year-round" = "#3498DB"         # Blue
    ),
    name = "Growing/Harvest Season"
  ) +
  labs(
    title = "Figure 2B: Top 20 Plant Taxa Most Strongly Associated with Sampling Time",
    subtitle = "PLSR Loading on Component 1 (CLR-transformed abundances)",
    x = "PLSR Loading (Component 1)",
    y = "Plant Taxon"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold", hjust = 0),
    plot.subtitle = element_text(size = 11, face = "italic", hjust = 0),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 11, face = "bold"),
    legend.position = "right",
    panel.grid.major.x = element_line(color = "gray90", linewidth = 0.3)
  )

# Save figure
ggsave(paste0(output_path, "Figure_2B_plant_taxa_by_season.png"),
       fig2b, width = 12, height = 8, dpi = 300)

cat("✓ Figure 2B saved to:", paste0(output_path, "Figure_2B_plant_taxa_by_season.png\n\n"))

###############################################################################
# Output Top 20 Plant Taxa List
###############################################################################

cat("═══════════════════════════════════════════════════════════════════\n")
cat("TOP 20 PLANT TAXA MOST STRONGLY ASSOCIATED WITH SAMPLING TIME\n")
cat("Identified by PLS Regression (Component 1, CLR-transformed)\n")
cat("═══════════════════════════════════════════════════════════════════\n\n")

cat("Rank | Plant Taxon                  | PLSR Loading | Season\n")
cat("─────┼──────────────────────────────┼──────────────┼──────────────\n")

for (i in 1:nrow(fig2b_data)) {
  loading_val <- fig2b_data$Loading[i]
  taxa_name <- fig2b_data$CommonName[i]
  season <- fig2b_data$Season[i]

  cat(sprintf("%3d  | %-28s | %+10.4f   | %s\n",
              i, taxa_name, loading_val, season))
}

cat("\n═══════════════════════════════════════════════════════════════════\n\n")

# Summary statistics
cat("Summary Statistics:\n")
cat("  Model: PLS regression (Month ~ Plant taxa abundances)\n")
cat("  Data transformation: CLR (centered log-ratio)\n")
cat("  Component analyzed: 1st component (loadings extracted)\n")
cat("  Number of samples:", nrow(plsr_data), "\n")
cat("  Number of plant taxa analyzed:", ncol(otu_df), "\n")
cat("  Loading range (all taxa): [", round(min(plsr_loadings), 4), ", ",
      round(max(plsr_loadings), 4), "]\n", sep="")
cat("  Loading range (top 20): [", round(min(fig2b_data$Loading), 4), ", ",
      round(max(fig2b_data$Loading), 4), "]\n\n", sep="")

# Print the 20 taxa names in a simple list
cat("Plant Taxa Names (in order of association strength):\n\n")
for (i in 1:nrow(fig2b_data)) {
  cat(i, ". ", as.character(fig2b_data$CommonName[i]), "\n", sep="")
}

cat("\n═══════════════════════════════════════════════════════════════════\n")
cat("Script completed successfully!\n")
cat("═══════════════════════════════════════════════════════════════════\n")
