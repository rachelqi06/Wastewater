#!/usr/bin/env Rscript
# scripts/run_all.R
# Driver: sets working dir, loads RDS objects, runs scripts, writes figures.

# Abort on error
options(error = function() { quit(status = 1) })

# Set working directory to repo root (when run in Actions, cwd is repo root)
if (requireNamespace("here", quietly = TRUE)) {
  setwd(here::here())
} else {
  # fallback: assume script is in scripts/
  setwd(normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")))
}

cat("Working dir:", getwd(), "\n")

# Create results dir
dir.create("results/figures", recursive = TRUE, showWarnings = FALSE)

# Load phyloseq objects from data/
ps_files <- list.files("data", pattern = "\\.rds$|\\.RDS$", full.names = TRUE)
if (length(ps_files) == 0) {
  message("Warning: no .rds files found in data/ — continuing anyway.")
} else {
  ps_objects <- lapply(ps_files, readRDS)
  names(ps_objects) <- basename(ps_files)
  cat("Loaded", length(ps_objects), "RDS files from data/\n")
}

# Run analysis sections in order (create these scripts next)
scripts_to_run <- c(
  "scripts/01_data_prep.R",
  "scripts/02_alpha_beta.R",
  "scripts/03_diff_abundance.R",
  "scripts/04_make_figures.R"
)

for (f in scripts_to_run) {
  if (file.exists(f)) {
    cat("Sourcing", f, "...\n")
    source(f, local = new.env())
  } else {
    cat("Skipping missing file:", f, "\n")
  }
}

cat("Driver finished. Figures should be in results/figures\n")
