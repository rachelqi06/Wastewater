# scripts/install_packages_from_rmd.R
rmd <- "Rcode.Rmd"
if (!file.exists(rmd)) {
  message("Rcode.Rmd not found; skipping package install.")
} else {
  txt <- readLines(rmd, warn = FALSE)
  # crude extraction of library(...) and require(...)
  libs <- unique(gsub(".*\\b(library|require)\\((['\"]?)([[:alnum:].]+)\\2\\).*", "\\3",
                      grep("\\b(library|require)\\(", txt, value = TRUE)))
  # also look for pacman::p_load, :: calls, etc. (simple)
  libs <- libs[libs != ""]
  if (length(libs) > 0) {
    inst <- libs[!(libs %in% installed.packages()[, "Package"])]
    if (length(inst) > 0) {
      message("Installing packages: ", paste(inst, collapse = ", "))
      install.packages(inst, repos = "https://cloud.r-project.org")
    } else {
      message("All packages already installed.")
    }
  } else {
    message("No library()/require() calls found in Rcode.Rmd.")
  }
}
