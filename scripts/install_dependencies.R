# Install the R packages used by the attribution analysis.
# Copyright (c) 2026 Michiel van Boven. BSD-3-Clause; see LICENSE.
cran_packages <- c("readxl", "readr", "dplyr", "tidyr", "tibble", "ggplot2",
                   "posterior", "loo", "cowplot", "withr")
missing <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) install.packages(missing, repos = "https://cloud.r-project.org")
if (!requireNamespace("cmdstanr", quietly = TRUE)) {
  install.packages("cmdstanr", repos = c("https://stan-dev.r-universe.dev", "https://cloud.r-project.org"))
}
message("R packages installed. For model fitting, install CmdStan once with cmdstanr::install_cmdstan().")
