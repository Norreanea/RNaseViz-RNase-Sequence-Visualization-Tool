# Alignment visualization ------------------------------------------------------

seqvisr_feature_list <- function(rnase, extra_feature = NULL) {
  cfg <- rnase_config(rnase)
  features <- cfg$features_seqvisr
  if (!is.null(extra_feature)) {
    selected_org <- extra_feature[[1]]
    selected_pos <- as.integer(extra_feature[[2]])
    features <- Filter(function(feature) {
      feature_org <- feature[[1]]
      feature_positions <- suppressWarnings(as.integer(feature[2:(length(feature) - 1)]))
      !(identical(feature_org, selected_org) && selected_pos %in% feature_positions)
    }, features)
    features <- c(features, list(extra_feature))
  }
  features
}

plot_global_alignment <- function(rnase, selected = NULL) {
  cfg <- rnase_config(rnase)
  validate_file(cfg$msa_seqvisr)

  features <- cfg$features_seqvisr
  colors <- cfg$colors

  p <- seqvisr::msavisr(
    mymsa = cfg$msa_seqvisr,
    myref = "HS",
    myroi = features,
    roicolors = colors,
    basecolors = c("Snow1", "Snow2", "Snow3"),
    wnon = 1.0,
    wroi = 4.0,
    hroi = 1,
    hmat = 0.9,
    hnon = 0.9,
    cbfcols = FALSE
  ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    ) +
    ggplot2::labs(title = paste0(rnase, " multiple sequence alignment"))

  if (!is.null(selected)) {
    selected_pos <- as.integer(selected[[2]])
    p <- p +
      ggplot2::geom_rect(
        data = data.frame(xmin = selected_pos - 0.5, xmax = selected_pos + 0.5),
        ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
        inherit.aes = FALSE,
        fill = "#ffff00",
        alpha = 0.35,
        colour = NA
      )
  }

  p
}

plot_alignment_window <- function(window_df) {
  window_df$residue_label <- ifelse(window_df$residue == "-", ".", window_df$residue)
  window_df$status <- ifelse(window_df$is_query_position, "Selected position", "Context")

  ggplot2::ggplot(
    window_df,
    ggplot2::aes(x = aligned_position, y = organism_code, fill = status)
  ) +
    ggplot2::geom_tile(color = "white", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(label = residue_label), size = 3.2) +
    ggplot2::scale_fill_manual(values = c("Context" = "#f2f2f2", "Selected position" = "#ffff00")) +
    ggplot2::scale_x_continuous(breaks = sort(unique(window_df$aligned_position))) +
    ggplot2::labs(
      x = "Alignment position",
      y = "Organism",
      fill = NULL,
      title = "Residue-level alignment context"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

extract_tree_candidate <- function(tree) {
  if (inherits(tree, c("girafe", "htmlwidget", "gg", "ggplot", "ggtree", "phylo", "multiPhylo", "treedata", "hclust"))) {
    return(tree)
  }

  if (is.character(tree) && length(tree) == 1) {
    return(tree)
  }

  if (is.list(tree)) {
    for (candidate_name in c("tree", "phylo", "ggtree", "ggobj", "plot", "widget", "object", "x")) {
      if (!is.null(tree[[candidate_name]])) {
        candidate <- extract_tree_candidate(tree[[candidate_name]])
        if (!is.null(candidate)) return(candidate)
      }
    }
  }

  NULL
}

tree_class_report <- function(tree) {
  candidate <- extract_tree_candidate(tree)
  paste0(
    "RDS class: ", paste(class(tree), collapse = ", "),
    if (!is.null(candidate)) paste0("; extracted class: ", paste(class(candidate), collapse = ", ")) else "; no supported candidate extracted"
  )
}

coerce_tree_for_plot <- function(tree) {
  candidate <- extract_tree_candidate(tree)

  if (is.null(candidate)) {
    stop("Unsupported tree object stored in RDS. ", tree_class_report(tree), call. = FALSE)
  }

  if (inherits(candidate, c("gg", "ggplot", "ggtree"))) {
    return(candidate)
  }

  if (inherits(candidate, "multiPhylo")) {
    return(candidate[[1]])
  }

  if (inherits(candidate, c("phylo", "treedata"))) {
    return(candidate)
  }

  if (inherits(candidate, "hclust")) {
    return(ape::as.phylo(candidate))
  }

  if (is.character(candidate) && length(candidate) == 1) {
    if (file.exists(candidate)) {
      return(ape::read.tree(candidate))
    }
    return(ape::read.tree(text = candidate))
  }

  stop("Unsupported tree object for static plotting. ", tree_class_report(tree), call. = FALSE)
}

plot_tree <- function(rnase) {
  candidate <- extract_tree_candidate(read_tree(rnase))
  if (is.null(candidate) || is.null(candidate$data)) {
    stop("Tree object does not contain plottable node data. ", tree_class_report(candidate), call. = FALSE)
  }

  tree_df <- as.data.frame(candidate$data, stringsAsFactors = FALSE)
  parent_df <- tree_df[, c("node", "x", "y")]
  names(parent_df) <- c("parent", "parent_x", "parent_y")

  horiz <- merge(
    tree_df[!is.na(tree_df$parent), c("parent", "x", "y")],
    parent_df,
    by = "parent",
    all.x = TRUE,
    sort = FALSE
  )

  vertical <- do.call(rbind, lapply(split(horiz, horiz$parent), function(df) {
    data.frame(
      x = df$parent_x[[1]],
      xend = df$parent_x[[1]],
      y = min(df$y, na.rm = TRUE),
      yend = max(df$y, na.rm = TRUE)
    )
  }))

  tips <- tree_df[tree_df$isTip %in% TRUE, c("x", "y", "label"), drop = FALSE]
  tips$uid <- tree_df$uid[tree_df$isTip %in% TRUE]
  tips$Phylum <- tree_df$Phylum[tree_df$isTip %in% TRUE]
  tips$size <- tree_df$size[tree_df$isTip %in% TRUE]
  tips$image_url <- paste0("https://images.phylopic.org/images/", tips$uid, "/vector.svg")
  tips$image_size <- tips$size * 0.55

  branch_labels <- tree_df[!is.na(tree_df$branch.length) & tree_df$branch.length > 0, c("branch", "y", "branch.length"), drop = FALSE]
  branch_labels$label <- round(branch_labels$branch.length, 2)

  phylum_cols <- c(
    "Ascomycota" = "darkblue",
    "Arthropoda" = "#8a2be2",
    "Chordata" = "#0b7d20",
    "Nematoda" = "#8b0000"
  )
  tree_title <- if (!is.null(candidate$labels$title)) candidate$labels$title else paste0(rnase, " ortholog tree")

  label_offset <- max(tree_df$x, na.rm = TRUE) * 0.03
  if (!is.finite(label_offset) || label_offset <= 0) {
    label_offset <- 0.03
  }

  ggplot2::ggplot() +
    ggplot2::geom_segment(data = vertical, ggplot2::aes(x = x, y = y, xend = xend, yend = yend), linewidth = 0.5, colour = "black") +
    ggplot2::geom_segment(data = horiz, ggplot2::aes(x = parent_x, y = y, xend = x, yend = y), linewidth = 0.5, colour = "black") +
    ggimage::geom_image(
      data = tips,
      ggplot2::aes(x = x + label_offset, y = y, image = image_url, colour = Phylum),
      size = tips$image_size,
      by = "height"
    ) +
    ggrepel::geom_text_repel(
      data = branch_labels,
      ggplot2::aes(x = branch, y = y, label = label),
      size = 4,
      segment.colour = NA,
      hjust = 1.3,
      vjust = -1.5
    ) +
    ggplot2::scale_colour_manual(values = phylum_cols, drop = FALSE) +
    ggplot2::coord_cartesian(xlim = c(-0.01, max(tree_df$x, na.rm = TRUE) + label_offset * 10), clip = "off") +
    ggplot2::labs(title = tree_title, x = NULL, y = NULL, colour = "Phylum") +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      text = ggplot2::element_text(colour = "black"),
      legend.position = "right",
      plot.margin = ggplot2::margin(10, 120, 10, 10)
    )
}
