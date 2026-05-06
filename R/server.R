# RNaseViz server --------------------------------------------------------------

rnaseviz_server <- function(input, output, session) {
  shiny::addResourcePath("rnase-domains", normalizePath(file.path("data", "domains"), winslash = "/", mustWork = TRUE))
  shiny::addResourcePath("rnase-trees", normalizePath(file.path("data", "trees"), winslash = "/", mustWork = TRUE))

  rv <- shiny::reactiveValues(
    selected_residue = NULL
  )

  selected_rnase <- shiny::reactive({
    input$rnase
  })

  shiny::observeEvent(input$rnase, {
    rv$selected_residue <- NULL
    orgs <- available_organisms_for_rnase(input$rnase)
    selected_org <- if ("Homo sapiens" %in% orgs) "Homo sapiens" else orgs[[1]]
    shiny::updateSelectInput(session, "referenceOrganism", choices = orgs, selected = selected_org)
  }, ignoreInit = FALSE)

  mapping_reactive <- shiny::reactive({
    shiny::req(selected_rnase(), input$referenceOrganism)
    position_mapping(selected_rnase(), input$referenceOrganism)
  })

  build_selected_residue <- function(row, source = "manual") {
    row <- row[1, , drop = FALSE]
    report <- residue_report(selected_rnase(), row$aligned_position)
    cons <- conservation_summary(report, row$original_aa)
    nearby <- nearby_catalog_features(selected_rnase(), row$aligned_position)

    list(
      rnase = selected_rnase(),
      organism = input$referenceOrganism,
      organism_code = get_organism_code(input$referenceOrganism),
      original_position = row$original_position,
      original_aa = row$original_aa,
      aligned_position = row$aligned_position,
      aligned_aa = row$aligned_aa,
      report = report,
      conservation = cons,
      nearby = nearby,
      source = source
    )
  }

  set_selected_from_original <- function(pos, source = "manual") {
    mapping <- mapping_reactive()
    row <- mapping[mapping$original_position == pos, , drop = FALSE]
    if (nrow(row) == 0) {
      shiny::showNotification(paste0("Position ", pos, " is outside the selected reference sequence."), type = "error")
      return(invisible(NULL))
    }

    rv$selected_residue <- build_selected_residue(row, source = source)
    shiny::updateTextInput(session, "mutationQuery", value = as.character(rv$selected_residue$original_position))
    shiny::updateTabsetPanel(session, "mainTabs", selected = "Alignment")
    invisible(rv$selected_residue)
  }

  set_selected_from_aligned <- function(aligned_pos, source = "alignment click") {
    if (is.null(aligned_pos) || is.na(aligned_pos) || !is.finite(aligned_pos)) return(invisible(NULL))

    aligned_pos <- as.integer(round(aligned_pos))
    mapping <- mapping_reactive()
    row <- mapping[mapping$aligned_position == aligned_pos, , drop = FALSE]

    if (nrow(row) == 0) {
      nearest_i <- which.min(abs(mapping$aligned_position - aligned_pos))
      row <- mapping[nearest_i, , drop = FALSE]
      shiny::showNotification(
        paste0("Clicked alignment position ", aligned_pos, "; selected nearest reference residue ", row$original_aa, row$original_position, "."),
        type = "message"
      )
    }

    rv$selected_residue <- build_selected_residue(row, source = source)
    shiny::updateTextInput(session, "mutationQuery", value = as.character(rv$selected_residue$original_position))
    shiny::updateTabsetPanel(session, "mainTabs", selected = "Alignment")
    invisible(rv$selected_residue)
  }

  shiny::observeEvent(input$highlightMutation, {
    pos <- parse_position_query(input$mutationQuery)
    if (is.na(pos) || pos <= 0) {
      shiny::showNotification("Enter a valid residue position, e.g. 520 or T520I.", type = "error")
      return(invisible(NULL))
    }
    set_selected_from_original(pos, source = "manual")
  })

  shiny::observeEvent(input$seqPlot_click, {
    set_selected_from_aligned(input$seqPlot_click$x, source = "global alignment click")
  })

  shiny::observeEvent(input$windowPlot_click, {
    set_selected_from_aligned(input$windowPlot_click$x, source = "window alignment click")
  })

  output$rnaseSummary <- shiny::renderUI({
    shiny::req(selected_rnase())
    cfg <- rnase_config(selected_rnase())
    n_variants <- nrow(cfg$features)
    n_org <- length(available_organisms_for_rnase(selected_rnase()))
    organism_key <- paste(paste0(organism_lookup$code, " = ", organism_lookup$organism), collapse = "; ")

    htmltools::tagList(
      tags$div(
        class = "row",
        tags$div(class = "col-md-4", tags$div(class = "card-lite", tags$div(class = "metric", cfg$title), tags$div("Selected RNase"))),
        tags$div(class = "col-md-4", tags$div(class = "card-lite", tags$div(class = "metric", n_variants), tags$div("Marked variants / regions"))),
        tags$div(class = "col-md-4", tags$div(class = "card-lite", tags$div(class = "metric", n_org), tags$div("Available orthologs")))
      ),
      tags$div(
        class = "row overview-block",
        tags$div(
          class = "col-lg-7",
          tags$div(
            class = "card-lite rnase-info-html overview-context-card",
            tags$h4("Scientific context"),
            HTML(rnase_info_html(selected_rnase()))
          )
        ),
        tags$div(
          class = "col-lg-5",
          tags$div(
            class = "card-lite overview-tree-card",
            tags$h4(paste(cfg$title, "ortholog tree")),
            tags$img(
              src = cfg$tree_image_url,
              alt = paste(selected_rnase(), "tree"),
              style = "width:100%; max-width:420px; height:auto; display:block; margin:0 auto;"
            ),
            tags$div(
              class = "small-note",
              "The length of the branches indicates the evolutionary distance between sequences (longer branches indicate greater distance). The positioning of the nodes shows the branching patterns and suggests common ancestors for groups of sequences."
            )
          )
        )
      ),
      tags$div(
        class = "card-lite overview-organisms-card",
        tags$h4("Organism abbreviations"),
        tags$p(class = "small-note", organism_key)
      )
    )
  })

  output$rnaseImageContainer <- shiny::renderUI({
    cfg <- rnase_config(input$rnase)
    legacy_exists <- file.exists(cfg$legacy_image)
    if (!legacy_exists) return(NULL)

    tags$div(
      class = "card-lite",
      tags$img(
        src = paste0("images/", input$rnase, ".png"),
        class = "rnase-image-zoomable",
        style = "max-width:100%; height:auto;",
        onclick = "document.getElementById('rnase-img-modal-content').src=this.src; document.getElementById('rnase-img-modal').style.display='block';"
      ),
      tags$div(class = "small-note", "RNase figure created with BioRender.com. Click to zoom.")
    )
  })

  output$alignmentDomainContainer <- shiny::renderUI({
    cfg <- rnase_config(input$rnase)
    if (!file.exists(cfg$image)) return(NULL)
    tags$div(
      class = "card-lite",
      tags$h4("Domain organization"),
      tags$img(
        src = cfg$image_url,
        class = "rnase-image-zoomable",
        style = "max-width:100%; height:auto; display:block; margin:0 auto;",
        onclick = "document.getElementById('rnase-img-modal-content').src=this.src; document.getElementById('rnase-img-modal').style.display='block';"
      ),
      tags$div(class = "small-note", "Domain organization figure from the article. Click to zoom.")
    )
  })

  output$alignmentDomainLegendContainer <- shiny::renderUI({
    cfg <- rnase_config(input$rnase)
    if (!file.exists(cfg$image)) return(NULL)
    tags$div(
      class = "card-lite figure-legend-card",
      tags$h4("Legend for the domain organization figure"),
      tags$p(
        class = "small-note",
        "Amino acid sequence alignments highlight the mutation sites (color background and position number), flanked by adjacent residues or deletion of a whole exon. Aligned protein sequences are depicted as gray rectangles with black borders, with particular domains shown as color rectangles according to the legend, and positions indicated on the axis. Mutation sites are indicated by red triangles."
      )
    )
  })

  output$mutationTable <- DT::renderDT({
    shiny::req(selected_rnase())
    cfg <- rnase_config(selected_rnase())
    out <- cfg$features[, c("organism_code", "span", "label", "variant_type")]
    names(out) <- c("Organism", "Aligned position", "Variant / feature", "Class")
    DT::datatable(
      out,
      rownames = FALSE,
      filter = "top",
      selection = "single",
      options = list(pageLength = 10, autoWidth = TRUE)
    )
  })

  shiny::observeEvent(input$mutationTable_rows_selected, {
    shiny::req(selected_rnase())
    idx <- input$mutationTable_rows_selected
    if (length(idx) != 1) return(invisible(NULL))
    cfg <- rnase_config(selected_rnase())
    feat <- cfg$features[idx, , drop = FALSE]
    if (feat$organism_code != get_organism_code(input$referenceOrganism)) {
      org <- organism_lookup$organism[match(feat$organism_code, organism_lookup$code)]
      if (!is.na(org) && org %in% available_organisms_for_rnase(selected_rnase())) {
        shiny::updateSelectInput(session, "referenceOrganism", selected = org)
      }
    }
    set_selected_from_aligned(feat$aligned_start, source = "variant table click")
  })

  output$seqPlot <- shiny::renderPlot({
    shiny::req(selected_rnase())
    selected <- NULL
    if (!is.null(rv$selected_residue)) {
      selected <- c(
        rv$selected_residue$organism_code,
        rv$selected_residue$aligned_position,
        paste0("Selected ", rv$selected_residue$original_aa, rv$selected_residue$original_position)
      )
    }
    plot_global_alignment(selected_rnase(), selected = selected)
  }, res = 120)

  output$alignedPosInfo <- shiny::renderUI({
    sel <- rv$selected_residue
    if (is.null(sel)) {
      return(tags$div(class = "warning-note", "Click the alignment or enter a residue to see mapped alignment coordinates and conservation."))
    }

    nearby_txt <- if (nrow(sel$nearby) == 0) {
      "No curated feature is within +/-20 aligned positions."
    } else {
      paste(sel$nearby$label, collapse = "; ")
    }

    tags$div(
      class = "card-lite",
      tags$p(tags$strong("Selection source: "), sel$source),
      tags$p(tags$strong("Reference residue: "), paste0(sel$original_aa, sel$original_position, " in ", sel$organism)),
      tags$p(tags$strong("Alignment coordinate: "), sel$aligned_position),
      tags$p(tags$strong("Conservation: "), paste0(sel$conservation$n_same, "/", sel$conservation$n_valid, " non-gap orthologs match the reference residue (", sel$conservation$percent_same, "%).")),
      tags$p(tags$strong("Nearby curated features: "), nearby_txt)
    )
  })

  output$windowPlot <- shiny::renderPlot({
    sel <- rv$selected_residue
    shiny::req(sel)
    plot_alignment_window(alignment_window(sel$rnase, sel$aligned_position, flank = 12L))
  }, res = 120)

  output$downloadAlignment <- shiny::downloadHandler(
    filename = function() {
      shiny::req(selected_rnase())
      if (identical(as.character(input$downloadFormat), "pdf")) {
        paste0(selected_rnase(), "_final.pdf")
      } else {
        paste0("reordered_for_ggmsa_", selected_rnase(), "Aln.fa")
      }
    },
    content = function(file) {
      shiny::req(selected_rnase())
      cfg <- rnase_config(selected_rnase())

      if (identical(as.character(input$downloadFormat), "pdf")) {
        src <- file.path("data/alignments", paste0(selected_rnase(), "_final.pdf"))
        validate_file(src)
        ok <- file.copy(src, file, overwrite = TRUE)
        if (!isTRUE(ok)) {
          stop("Failed to copy alignment PDF to the download path.", call. = FALSE)
        }
      } else {
        validate_file(cfg$msa_ggmsa)
        ok <- file.copy(cfg$msa_ggmsa, file, overwrite = TRUE)
        if (!isTRUE(ok)) {
          stop("Failed to copy alignment FASTA to the download path.", call. = FALSE)
        }
      }
    }
  )
}
