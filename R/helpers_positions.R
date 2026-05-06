# Position mapping and residue reports ----------------------------------------

map_amino_acid_positions <- function(original_seq, aligned_seq) {
  original_aas <- strsplit(original_seq, "", fixed = TRUE)[[1]]
  aligned_aas <- strsplit(aligned_seq, "", fixed = TRUE)[[1]]

  original_pos <- integer(0)
  aligned_pos <- integer(0)
  original_aa <- character(0)
  aligned_aa <- character(0)
  j <- 1L

  for (i in seq_along(aligned_aas)) {
    if (!identical(aligned_aas[[i]], "-")) {
      if (j <= length(original_aas)) {
        original_pos <- c(original_pos, j)
        aligned_pos <- c(aligned_pos, i)
        original_aa <- c(original_aa, original_aas[[j]])
        aligned_aa <- c(aligned_aa, aligned_aas[[i]])
        j <- j + 1L
      }
    }
  }

  data.frame(
    original_position = original_pos,
    original_aa = original_aa,
    aligned_position = aligned_pos,
    aligned_aa = aligned_aa,
    stringsAsFactors = FALSE
  )
}

position_mapping <- function(rnase, organism) {
  map_amino_acid_positions(
    original_seq = get_original_sequence(rnase, organism),
    aligned_seq = get_aligned_sequence(rnase, organism)
  )
}

parse_position_query <- function(query) {
  query <- trimws(query %||% "")
  if (!nzchar(query)) return(NA_integer_)
  hit <- regmatches(query, regexpr("[0-9]+", query))
  if (!length(hit) || identical(hit, character(0))) return(NA_integer_)
  as.integer(hit)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

residue_report <- function(rnase, aligned_position) {
  aln <- read_alignment(rnase)
  rows <- lapply(names(aln), function(code) {
    seq_chars <- strsplit(as.character(aln[[code]]), "", fixed = TRUE)[[1]]
    organism <- organism_lookup$organism[match(code, organism_lookup$code)]
    data.frame(
      organism_code = code,
      organism = ifelse(is.na(organism), code, organism),
      aligned_position = aligned_position,
      residue = if (aligned_position <= length(seq_chars)) seq_chars[[aligned_position]] else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

conservation_summary <- function(report, reference_residue) {
  valid <- report[!is.na(report$residue) & report$residue != "-", , drop = FALSE]
  if (nrow(valid) == 0 || is.na(reference_residue)) {
    return(list(n_valid = 0L, n_same = NA_integer_, percent_same = NA_real_))
  }
  n_same <- sum(valid$residue == reference_residue)
  list(
    n_valid = nrow(valid),
    n_same = n_same,
    percent_same = round(100 * n_same / nrow(valid), 1)
  )
}

alignment_window <- function(rnase, aligned_position, flank = 12L) {
  aln <- read_alignment(rnase)
  max_len <- max(width(aln))
  start <- max(1L, aligned_position - flank)
  end <- min(max_len, aligned_position + flank)

  rows <- list()
  k <- 1L
  for (code in names(aln)) {
    chars <- strsplit(as.character(aln[[code]]), "", fixed = TRUE)[[1]]
    pos <- start:end
    organism <- organism_lookup$organism[match(code, organism_lookup$code)]
    rows[[k]] <- data.frame(
      organism_code = code,
      organism = ifelse(is.na(organism), code, organism),
      aligned_position = pos,
      residue = chars[pos],
      is_query_position = pos == aligned_position,
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
  do.call(rbind, rows)
}

nearby_catalog_features <- function(rnase, aligned_position, window = 20L) {
  feats <- rnase_config(rnase)$features
  feats[abs(feats$aligned_start - aligned_position) <= window |
          abs(feats$aligned_end - aligned_position) <= window, , drop = FALSE]
}
