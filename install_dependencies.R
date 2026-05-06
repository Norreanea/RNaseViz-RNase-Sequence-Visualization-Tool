# Run this once locally or during deployment setup. Do NOT source this from ui.R/server.R.
cran_packages <- c("shiny", "shinyjs", "bslib", "htmltools", "DT", "ggplot2", "plotly", "ape")
bioc_packages <- c("Biostrings", "ggtree")

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

if (!requireNamespace("seqvisr", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("vragh/seqvisr")
}
