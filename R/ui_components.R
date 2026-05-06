# RNaseViz UI -----------------------------------------------------------------

app_css <- "
.rnaseviz-title { margin-bottom: 0.2rem; }
.rnaseviz-subtitle { color: #586069; margin-top: 0; }
.rnaseviz-article-note { color: #4b5563; margin-bottom: 1rem; }
.card-lite { border: 1px solid #e5e7eb; border-radius: 12px; padding: 14px; background: #ffffff; box-shadow: 0 1px 4px rgba(0,0,0,0.04); margin-bottom: 12px; }
.metric { font-size: 1.35rem; font-weight: 700; }
.small-note { color: #6b7280; font-size: 0.9rem; }
.warning-note { color: #8a4b00; background: #fff8e6; border-left: 4px solid #f4b400; padding: 0.75rem; border-radius: 8px; }
.sidebar .form-group { margin-bottom: 1rem; }
.rnase-info-html a { text-decoration: none; }
.rnase-image-zoomable { cursor: zoom-in; transition: transform 0.15s ease, box-shadow 0.15s ease; border-radius: 8px; }
.rnase-image-zoomable:hover { transform: scale(1.015); box-shadow: 0 4px 18px rgba(0,0,0,0.16); }
.rnase-img-modal { display: none; position: fixed; z-index: 9999; padding-top: 40px; left: 0; top: 0; width: 100%; height: 100%; overflow: auto; background-color: rgba(0,0,0,0.88); }
.rnase-img-modal-content { margin: auto; display: block; max-width: 94%; max-height: 88vh; background: #ffffff; border-radius: 10px; }
.rnase-img-modal-close { position: absolute; top: 16px; right: 30px; color: #ffffff; font-size: 42px; font-weight: bold; cursor: pointer; }
.click-note { color: #4b5563; font-size: 0.95rem; margin-bottom: 8px; }
.rnaseviz-layout { display: flex; gap: 0; align-items: flex-start; }
.rnaseviz-sidebar { flex: 0 0 320px; width: 320px; min-width: 280px; max-width: 520px; overflow-x: hidden; overflow-y: auto; }
.rnaseviz-splitter { width: 14px; min-height: 100vh; cursor: col-resize; position: relative; flex: 0 0 14px; }
.rnaseviz-splitter::before { content: ''; position: absolute; left: 6px; top: 0; bottom: 0; width: 2px; background: #d1d5db; border-radius: 999px; }
.rnaseviz-main { flex: 1 1 auto; min-width: 0; }
.overview-context-card { min-height: 100%; }
.overview-tree-card { min-height: 100%; }
.overview-block { margin-bottom: 12px; }
.overview-organisms-card { margin-top: 12px; }
.figure-legend-card { margin-top: 12px; }
@media (max-width: 991px) {
  .rnaseviz-layout { display: block; }
  .rnaseviz-sidebar { width: 100%; max-width: none; margin-bottom: 16px; }
  .rnaseviz-splitter { display: none; }
}
"

rnaseviz_ui <- function() {
  shiny::fluidPage(
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly"),
    shinyjs::useShinyjs(),
    tags$head(tags$style(HTML(app_css)), tags$script(splitter_js)),

    tags$h1("RNaseViz", class = "rnaseviz-title"),
    tags$p(
      "Comparative RNase sequence visualization, disease-variant and residue-level conservation inspection.",
      class = "rnaseviz-subtitle"
    ),
    tags$p(
      class = "rnaseviz-article-note",
      HTML("This tool was developed as part of the article <a href='https://www.sciencedirect.com/science/article/pii/S2352304225001023' target='_blank'><em>Ribonucleases in Mendelian disease: Characterization and insight from model organisms</em></a>.")
    ),

    tags$div(
      id = "rnase-img-modal",
      class = "rnase-img-modal",
      onclick = "this.style.display='none';",
      tags$span(class = "rnase-img-modal-close", HTML("&times;")),
      tags$img(id = "rnase-img-modal-content", class = "rnase-img-modal-content")
    ),

    tags$div(
      class = "rnaseviz-layout",
      tags$div(
        id = "rnaseviz-sidebar",
        class = "sidebar rnaseviz-sidebar",
        tags$div(
          class = "card-lite",
          shiny::selectInput("rnase", "RNase family member", choices = rnase_order, selected = "ELAC2"),
          tags$hr(),
          shiny::selectInput("referenceOrganism", "Reference organism", choices = organism_lookup$organism, selected = "Homo sapiens"),
          shiny::textInput("mutationQuery", "Reference position or mutation label", value = "520", placeholder = "Example: 520 or T520I"),
          shiny::actionButton("highlightMutation", "Highlight / inspect residue", class = "btn-secondary w-100"),
          tags$br(), tags$br(),
          shiny::radioButtons(
            "downloadFormat",
            "Alignment download format",
            choices = c(
              "FASTA alignment (.fa)" = "fasta",
              "Alignment figure with marked pathogenic variants (PDF)" = "pdf"
            ),
            selected = "fasta",
            inline = FALSE
          ),
          shiny::downloadButton("downloadAlignment", "Download alignment", class = "w-100"),
          tags$hr(),
          shiny::uiOutput("rnaseImageContainer"),
          tags$div(
            class = "small-note",
            "Tip: drag the right edge of this panel to widen it. Select an RNase, then click the alignment to select the nearest aligned residue, or enter a reference-sequence position and press Highlight / inspect residue. Click the RNase image to zoom."
          )
        )
      ),
      tags$div(id = "rnaseviz-splitter", class = "rnaseviz-splitter"),
      tags$div(
        class = "rnaseviz-main",
        shiny::tabsetPanel(
          id = "mainTabs",
          shiny::tabPanel(
            "Overview",
            tags$br(),
            shiny::uiOutput("rnaseSummary"),
            shiny::h4("Curated variants / marked positions"),
            DT::DTOutput("mutationTable")
          ),
          shiny::tabPanel(
            "Alignment",
            tags$br(),
            tags$div(class = "click-note", "Click anywhere on the global alignment or the residue-window plot to update the selected residue."),
            tags$div(
              class = "row",
              tags$div(
                class = "col-lg-7",
                shiny::uiOutput("alignedPosInfo"),
                shiny::uiOutput("alignmentDomainLegendContainer")
              ),
              tags$div(class = "col-lg-5", shiny::uiOutput("alignmentDomainContainer"))
            ),
            shiny::plotOutput("seqPlot", height = "680px", click = "seqPlot_click"),
            tags$br(),
            shiny::plotOutput("windowPlot", height = "420px", click = "windowPlot_click")
          )
        )
      )
    )
  )
}

splitter_js <- HTML("
document.addEventListener('DOMContentLoaded', function () {
  var splitter = document.getElementById('rnaseviz-splitter');
  var sidebar = document.getElementById('rnaseviz-sidebar');
  if (!splitter || !sidebar) return;

  var dragging = false;

  splitter.addEventListener('mousedown', function (e) {
    dragging = true;
    document.body.style.userSelect = 'none';
    e.preventDefault();
  });

  document.addEventListener('mousemove', function (e) {
    if (!dragging) return;
    var minWidth = 280;
    var maxWidth = 520;
    var nextWidth = Math.min(maxWidth, Math.max(minWidth, e.clientX - sidebar.getBoundingClientRect().left));
    sidebar.style.width = nextWidth + 'px';
    sidebar.style.flexBasis = nextWidth + 'px';
  });

  document.addEventListener('mouseup', function () {
    dragging = false;
    document.body.style.userSelect = '';
  });
});
")
