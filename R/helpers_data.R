# Cached data loading ----------------------------------------------------------

.rnaseviz_cache <- new.env(parent = emptyenv())

cache_get <- function(key, expr) {
  if (!exists(key, envir = .rnaseviz_cache, inherits = FALSE)) {
    assign(key, force(expr), envir = .rnaseviz_cache)
  }
  get(key, envir = .rnaseviz_cache, inherits = FALSE)
}

read_raw_sequences <- function(rnase) {
  cfg <- rnase_config(rnase)
  cache_get(paste0("raw_", rnase), {
    validate_file(cfg$raw_seq)
    readRDS(cfg$raw_seq)
  })
}

read_alignment <- function(rnase) {
  cfg <- rnase_config(rnase)
  cache_get(paste0("alignment_", rnase), {
    validate_file(cfg$msa_ggmsa)
    Biostrings::readAAStringSet(cfg$msa_ggmsa)
  })
}

read_tree <- function(rnase) {
  cfg <- rnase_config(rnase)
  cache_get(paste0("tree_", rnase), {
    validate_file(cfg$tree)
    readRDS(cfg$tree)
  })
}

validate_file <- function(path) {
  if (!file.exists(path)) {
    stop("Required file is missing: ", path, call. = FALSE)
  }
  invisible(TRUE)
}

get_original_sequence <- function(rnase, organism) {
  raw <- read_raw_sequences(rnase)
  nm <- get_organism_rds_name(organism)
  if (is.na(nm) || !nm %in% names(raw)) {
    stop("No raw sequence found for ", organism, " in ", rnase, call. = FALSE)
  }
  as.character(raw[[nm]])
}

get_aligned_sequence <- function(rnase, organism) {
  aln <- read_alignment(rnase)
  code <- get_organism_code(organism)
  if (is.na(code) || !code %in% names(aln)) {
    stop("No aligned sequence found for ", organism, " in ", rnase, call. = FALSE)
  }
  as.character(aln[[code]])
}

available_organisms_for_rnase <- function(rnase) {
  cfg <- rnase_config(rnase)
  if (file.exists(cfg$msa_ggmsa)) {
    aln <- read_alignment(rnase)
    organism_lookup$organism[organism_lookup$code %in% names(aln)]
  } else {
    cfg$organisms
  }
}

file_status_table <- function(rnase) {
  cfg <- rnase_config(rnase)
  data.frame(
    data_type = c("MSA for seqvisr", "MSA for residue window", "Raw sequences", "Tree image", "Domain image"),
    path = c(cfg$msa_seqvisr, cfg$msa_ggmsa, cfg$raw_seq, cfg$tree_image, cfg$image),
    status = c(
      file.exists(cfg$msa_seqvisr), file.exists(cfg$msa_ggmsa), file.exists(cfg$raw_seq),
      file.exists(cfg$tree_image), file.exists(cfg$image)
    ),
    stringsAsFactors = FALSE
  )
}
