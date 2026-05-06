# RNaseViz package loading -----------------------------------------------------
# Do not install packages at app runtime. Install once before deployment using
# install_dependencies.R or renv, then keep runtime limited to library() calls.

required_packages <- c(
  "shiny", "shinyjs", "bslib", "htmltools", "DT",
  "ggplot2", "plotly", "ape", "ggtree", "Biostrings", "seqvisr", "ggiraph", "htmlwidgets",
  "ggimage", "ggrepel"
)

load_rnaseviz_packages <- function() {
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) {
    stop(
      "Missing packages: ", paste(missing, collapse = ", "),
      ". Install them before launching the app; do not install packages inside server.R.",
      call. = FALSE
    )
  }

  suppressPackageStartupMessages({
    library(shiny)
    library(shinyjs)
    library(bslib)
    library(htmltools)
    library(DT)
    library(ggplot2)
    library(plotly)
    library(ape)
    library(ggtree)
    library(Biostrings)
    library(seqvisr)
    library(ggiraph)
    library(htmlwidgets)
    library(ggimage)
    library(ggrepel)
  })
}
