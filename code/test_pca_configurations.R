#!/usr/bin/env Rscript
###############################################################################
# Test Multiple PCA Configurations for Figure 2A
# Goal: Find settings that give PC3 ≈ 4.9% and PC4 ≈ 4.1%
###############################################################################

set.seed(001)

library(phyloseq)
library(microbiome)
library(lubridate)

# Set paths
data_path <- "C:/Users/rache/Box/project_davidlab/LAD_LAB_Personnel/Rachel_Q/Code and Data/Data for code optimization_Do not submit/"

# Load data
NCWW_allsamples_trnL <- readRDS(paste0(data_path, "NCWW_allsamples_trnL.rds"))
NCWW_allsamples_animal <- readRDS(paste0(data_path, "NCWW_allsamples_animal.rds"))

# Create plant and animal subsets
ps_plant <- subset_taxa(NCWW_allsamples_trnL, phylum == "Streptophyta" & IsFood == "Y")
ps_plant <- prune_taxa(taxa_sums(ps_plant) > 0, ps_plant)

ps_animal <- subset_taxa(NCWW_allsamples_animal, IsFood == "Y" & CommonName != "Emu")
ps_animal <- prune_taxa(taxa_sums(ps_animal) > 0, ps_animal)

target_counties <- c("carteret", "mecklenburg", "pitt")

samples_to_remove <- c("NCWW121720_1","NCWW121720_2","NCWW121720_10",
                       "NCWW121720_56","NCWW121720_52","NCWW121720_18",
                       "NCWW121720_21","NCWW022621_16","NCWW121720_27",
                       "NCWW121720_6")

cat("╔════════════════════════════════════════════════════════════════════════════╗\n")
cat("║ Testing PCA Configurations for Figure 2A                                   ║\n")
cat("║ Target: PC3 ≈ 4.9%, PC4 ≈ 4.1%                                           ║\n")
cat("╚════════════════════════════════════════════════════════════════════════════╝\n\n")

results <- data.frame(
  Config = character(),
  Description = character(),
  PC3 = numeric(),
  PC4 = numeric(),
  PC3_Error = numeric(),
  PC4_Error = numeric(),
  stringsAsFactors = FALSE
)

config_num <- 0

###############################################################################
# Configuration 1: County filtering + remove samples + scale = FALSE
###############################################################################
config_num <- config_num + 1

ps_p1 <- subset_samples(ps_plant, County %in% target_counties)
ps_a1 <- subset_samples(ps_animal, County %in% target_counties)

ps_p1 <- prune_samples(!(sample_names(ps_p1) %in% samples_to_remove), ps_p1)
ps_a1 <- prune_samples(!(sample_names(ps_a1) %in% samples_to_remove), ps_a1)

ps_p1 <- prune_taxa(taxa_sums(ps_p1) > 0, ps_p1)
ps_a1 <- prune_taxa(taxa_sums(ps_a1) > 0, ps_a1)

ps_p1_clr <- microbiome::transform(ps_p1, "clr")
ps_a1_clr <- microbiome::transform(ps_a1, "clr")

ps_notree_1 <- phyloseq(otu_table(ps_p1_clr), tax_table(ps_p1_clr), sample_data(ps_p1_clr))
ps_notree_2 <- phyloseq(otu_table(ps_a1_clr), tax_table(ps_a1_clr), sample_data(ps_a1_clr))
ps_combined_1 <- merge_phyloseq(ps_notree_1, ps_notree_2)

otu_data_1 <- as.data.frame(otu_table(ps_combined_1))
pca_1 <- prcomp(otu_data_1, scale. = FALSE)
var_1 <- pca_1$sdev^2 / sum(pca_1$sdev^2) * 100

results <- rbind(results, data.frame(
  Config = paste0("Config ", config_num),
  Description = "County filter + remove samples + scale=FALSE",
  PC3 = round(var_1[3], 2),
  PC4 = round(var_1[4], 2),
  PC3_Error = round(abs(var_1[3] - 4.9), 2),
  PC4_Error = round(abs(var_1[4] - 4.1), 2)
))

cat(sprintf("Config %d: County filter + remove samples + scale=FALSE\n", config_num))
cat(sprintf("  PC3: %.2f%% (error: %.2f%%)\n", var_1[3], abs(var_1[3] - 4.9)))
cat(sprintf("  PC4: %.2f%% (error: %.2f%%)\n\n", var_1[4], abs(var_1[4] - 4.1)))

###############################################################################
# Configuration 2: County filtering + NO sample removal + scale = FALSE
###############################################################################
config_num <- config_num + 1

ps_p2 <- subset_samples(ps_plant, County %in% target_counties)
ps_a2 <- subset_samples(ps_animal, County %in% target_counties)

ps_p2 <- prune_taxa(taxa_sums(ps_p2) > 0, ps_p2)
ps_a2 <- prune_taxa(taxa_sums(ps_a2) > 0, ps_a2)

ps_p2_clr <- microbiome::transform(ps_p2, "clr")
ps_a2_clr <- microbiome::transform(ps_a2, "clr")

ps_notree_3 <- phyloseq(otu_table(ps_p2_clr), tax_table(ps_p2_clr), sample_data(ps_p2_clr))
ps_notree_4 <- phyloseq(otu_table(ps_a2_clr), tax_table(ps_a2_clr), sample_data(ps_a2_clr))
ps_combined_2 <- merge_phyloseq(ps_notree_3, ps_notree_4)

otu_data_2 <- as.data.frame(otu_table(ps_combined_2))
pca_2 <- prcomp(otu_data_2, scale. = FALSE)
var_2 <- pca_2$sdev^2 / sum(pca_2$sdev^2) * 100

results <- rbind(results, data.frame(
  Config = paste0("Config ", config_num),
  Description = "County filter + NO sample removal + scale=FALSE",
  PC3 = round(var_2[3], 2),
  PC4 = round(var_2[4], 2),
  PC3_Error = round(abs(var_2[3] - 4.9), 2),
  PC4_Error = round(abs(var_2[4] - 4.1), 2)
))

cat(sprintf("Config %d: County filter + NO sample removal + scale=FALSE\n", config_num))
cat(sprintf("  PC3: %.2f%% (error: %.2f%%)\n", var_2[3], abs(var_2[3] - 4.9)))
cat(sprintf("  PC4: %.2f%% (error: %.2f%%)\n\n", var_2[4], abs(var_2[4] - 4.1)))

###############################################################################
# Configuration 3: County filtering + remove samples + scale = TRUE
###############################################################################
config_num <- config_num + 1

ps_p3 <- subset_samples(ps_plant, County %in% target_counties)
ps_a3 <- subset_samples(ps_animal, County %in% target_counties)

ps_p3 <- prune_samples(!(sample_names(ps_p3) %in% samples_to_remove), ps_p3)
ps_a3 <- prune_samples(!(sample_names(ps_a3) %in% samples_to_remove), ps_a3)

ps_p3 <- prune_taxa(taxa_sums(ps_p3) > 0, ps_p3)
ps_a3 <- prune_taxa(taxa_sums(ps_a3) > 0, ps_a3)

ps_p3_clr <- microbiome::transform(ps_p3, "clr")
ps_a3_clr <- microbiome::transform(ps_a3, "clr")

ps_notree_5 <- phyloseq(otu_table(ps_p3_clr), tax_table(ps_p3_clr), sample_data(ps_p3_clr))
ps_notree_6 <- phyloseq(otu_table(ps_a3_clr), tax_table(ps_a3_clr), sample_data(ps_a3_clr))
ps_combined_3 <- merge_phyloseq(ps_notree_5, ps_notree_6)

otu_data_3 <- as.data.frame(otu_table(ps_combined_3))
pca_3 <- prcomp(otu_data_3, scale. = TRUE)
var_3 <- pca_3$sdev^2 / sum(pca_3$sdev^2) * 100

results <- rbind(results, data.frame(
  Config = paste0("Config ", config_num),
  Description = "County filter + remove samples + scale=TRUE",
  PC3 = round(var_3[3], 2),
  PC4 = round(var_3[4], 2),
  PC3_Error = round(abs(var_3[3] - 4.9), 2),
  PC4_Error = round(abs(var_3[4] - 4.1), 2)
))

cat(sprintf("Config %d: County filter + remove samples + scale=TRUE\n", config_num))
cat(sprintf("  PC3: %.2f%% (error: %.2f%%)\n", var_3[3], abs(var_3[3] - 4.9)))
cat(sprintf("  PC4: %.2f%% (error: %.2f%%)\n\n", var_3[4], abs(var_3[4] - 4.1)))

###############################################################################
# Configuration 4: County filtering + NO sample removal + scale = TRUE
###############################################################################
config_num <- config_num + 1

ps_p4 <- subset_samples(ps_plant, County %in% target_counties)
ps_a4 <- subset_samples(ps_animal, County %in% target_counties)

ps_p4 <- prune_taxa(taxa_sums(ps_p4) > 0, ps_p4)
ps_a4 <- prune_taxa(taxa_sums(ps_a4) > 0, ps_a4)

ps_p4_clr <- microbiome::transform(ps_p4, "clr")
ps_a4_clr <- microbiome::transform(ps_a4, "clr")

ps_notree_7 <- phyloseq(otu_table(ps_p4_clr), tax_table(ps_p4_clr), sample_data(ps_p4_clr))
ps_notree_8 <- phyloseq(otu_table(ps_a4_clr), tax_table(ps_a4_clr), sample_data(ps_a4_clr))
ps_combined_4 <- merge_phyloseq(ps_notree_7, ps_notree_8)

otu_data_4 <- as.data.frame(otu_table(ps_combined_4))
pca_4 <- prcomp(otu_data_4, scale. = TRUE)
var_4 <- pca_4$sdev^2 / sum(pca_4$sdev^2) * 100

results <- rbind(results, data.frame(
  Config = paste0("Config ", config_num),
  Description = "County filter + NO sample removal + scale=TRUE",
  PC3 = round(var_4[3], 2),
  PC4 = round(var_4[4], 2),
  PC3_Error = round(abs(var_4[3] - 4.9), 2),
  PC4_Error = round(abs(var_4[4] - 4.1), 2)
))

cat(sprintf("Config %d: County filter + NO sample removal + scale=TRUE\n", config_num))
cat(sprintf("  PC3: %.2f%% (error: %.2f%%)\n", var_4[3], abs(var_4[3] - 4.9)))
cat(sprintf("  PC4: %.2f%% (error: %.2f%%)\n\n", var_4[4], abs(var_4[4] - 4.1)))

###############################################################################
# Configuration 5: Location filtering (Beaufort, Charlotte 1-4, Greenville) + remove samples + scale=FALSE
###############################################################################
config_num <- config_num + 1

target_locations <- c("Beaufort", "Charlotte 1", "Charlotte 2", "Charlotte 3", "Charlotte 4", "Greenville")

ps_p5 <- subset_samples(ps_plant, Location %in% target_locations)
ps_a5 <- subset_samples(ps_animal, Location %in% target_locations)

ps_p5 <- prune_samples(!(sample_names(ps_p5) %in% samples_to_remove), ps_p5)
ps_a5 <- prune_samples(!(sample_names(ps_a5) %in% samples_to_remove), ps_a5)

ps_p5 <- prune_taxa(taxa_sums(ps_p5) > 0, ps_p5)
ps_a5 <- prune_taxa(taxa_sums(ps_a5) > 0, ps_a5)

ps_p5_clr <- microbiome::transform(ps_p5, "clr")
ps_a5_clr <- microbiome::transform(ps_a5, "clr")

ps_notree_9 <- phyloseq(otu_table(ps_p5_clr), tax_table(ps_p5_clr), sample_data(ps_p5_clr))
ps_notree_10 <- phyloseq(otu_table(ps_a5_clr), tax_table(ps_a5_clr), sample_data(ps_a5_clr))
ps_combined_5 <- merge_phyloseq(ps_notree_9, ps_notree_10)

otu_data_5 <- as.data.frame(otu_table(ps_combined_5))
pca_5 <- prcomp(otu_data_5, scale. = FALSE)
var_5 <- pca_5$sdev^2 / sum(pca_5$sdev^2) * 100

results <- rbind(results, data.frame(
  Config = paste0("Config ", config_num),
  Description = "Location filter + remove samples + scale=FALSE",
  PC3 = round(var_5[3], 2),
  PC4 = round(var_5[4], 2),
  PC3_Error = round(abs(var_5[3] - 4.9), 2),
  PC4_Error = round(abs(var_5[4] - 4.1), 2)
))

cat(sprintf("Config %d: Location filter + remove samples + scale=FALSE\n", config_num))
cat(sprintf("  PC3: %.2f%% (error: %.2f%%)\n", var_5[3], abs(var_5[3] - 4.9)))
cat(sprintf("  PC4: %.2f%% (error: %.2f%%)\n\n", var_5[4], abs(var_5[4] - 4.1)))

###############################################################################
# Configuration 6: Location filtering + NO sample removal + scale=FALSE
###############################################################################
config_num <- config_num + 1

ps_p6 <- subset_samples(ps_plant, Location %in% target_locations)
ps_a6 <- subset_samples(ps_animal, Location %in% target_locations)

ps_p6 <- prune_taxa(taxa_sums(ps_p6) > 0, ps_p6)
ps_a6 <- prune_taxa(taxa_sums(ps_a6) > 0, ps_a6)

ps_p6_clr <- microbiome::transform(ps_p6, "clr")
ps_a6_clr <- microbiome::transform(ps_a6, "clr")

ps_notree_11 <- phyloseq(otu_table(ps_p6_clr), tax_table(ps_p6_clr), sample_data(ps_p6_clr))
ps_notree_12 <- phyloseq(otu_table(ps_a6_clr), tax_table(ps_a6_clr), sample_data(ps_a6_clr))
ps_combined_6 <- merge_phyloseq(ps_notree_11, ps_notree_12)

otu_data_6 <- as.data.frame(otu_table(ps_combined_6))
pca_6 <- prcomp(otu_data_6, scale. = FALSE)
var_6 <- pca_6$sdev^2 / sum(pca_6$sdev^2) * 100

results <- rbind(results, data.frame(
  Config = paste0("Config ", config_num),
  Description = "Location filter + NO sample removal + scale=FALSE",
  PC3 = round(var_6[3], 2),
  PC4 = round(var_6[4], 2),
  PC3_Error = round(abs(var_6[3] - 4.9), 2),
  PC4_Error = round(abs(var_6[4] - 4.1), 2)
))

cat(sprintf("Config %d: Location filter + NO sample removal + scale=FALSE\n", config_num))
cat(sprintf("  PC3: %.2f%% (error: %.2f%%)\n", var_6[3], abs(var_6[3] - 4.9)))
cat(sprintf("  PC4: %.2f%% (error: %.2f%%)\n\n", var_6[4], abs(var_6[4] - 4.1)))

###############################################################################
# Summary Table
###############################################################################

cat("\n╔════════════════════════════════════════════════════════════════════════════╗\n")
cat("║ SUMMARY: All Configurations Ranked by Total Error                          ║\n")
cat("╚════════════════════════════════════════════════════════════════════════════╝\n\n")

results$Total_Error <- results$PC3_Error + results$PC4_Error
results_sorted <- results[order(results$Total_Error), ]

print(results_sorted[, c("Config", "Description", "PC3", "PC4", "PC3_Error", "PC4_Error", "Total_Error")])

cat("\n╔════════════════════════════════════════════════════════════════════════════╗\n")
cat("║ BEST CONFIGURATION:                                                        ║\n")
cat("╚════════════════════════════════════════════════════════════════════════════╝\n\n")

best_idx <- which.min(results_sorted$Total_Error)
cat(sprintf("Config: %s\n", results_sorted$Config[best_idx]))
cat(sprintf("Description: %s\n", results_sorted$Description[best_idx]))
cat(sprintf("PC3: %.2f%% (target 4.9%%, error %.2f%%)\n", results_sorted$PC3[best_idx], results_sorted$PC3_Error[best_idx]))
cat(sprintf("PC4: %.2f%% (target 4.1%%, error %.2f%%)\n", results_sorted$PC4[best_idx], results_sorted$PC4_Error[best_idx]))
cat(sprintf("Total Error: %.2f%%\n\n", results_sorted$Total_Error[best_idx]))
