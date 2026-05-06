# RNaseViz data configuration --------------------------------------------------

organism_lookup <- data.frame(
  organism = c(
    "Homo sapiens", "Mus musculus", "Rattus norvegicus", "Danio rerio",
    "Xenopus tropicalis", "Caenorhabditis elegans", "Drosophila melanogaster",
    "Saccharomyces cerevisiae"
  ),
  code = c("HS", "MM", "RN", "DR", "XT", "CE", "DM", "SC"),
  rds_name = c(
    "Homo_sapiens", "Mus_musculus", "Rattus_norvegicus", "Danio_rerio",
    "Xenopus_tropicalis", "Caenorhabditis_elegans", "Drosophila_melanogaster",
    "Saccharomyces_cerevisiae"
  ),
  stringsAsFactors = FALSE
)

rnase_order <- c(
  "RNASEH2A", "RNASEH2B", "RNASEH2C", "AGO2", "RNASET2",
  "DICER1", "DIS3L2", "ELAC2", "PRORP", "PARN", "RNASEH1"
)

rnase_palette <- c(
  "#1f77b4", "#b22222", "#556b2f", "#8b008b", "#8b4513",
  "#0f766e", "#1d4ed8", "#7c3aed", "#be123c", "#14532d",
  "#334155", "#9f1239", "#0369a1", "#4b5563", "#6b7280"
)

rnase_asset_map <- list(
  RNASEH1 = list(domain = "RNASEH1.png", tree = "RNASEH1.jpg"),
  RNASEH2A = list(domain = "RNaseH2A.JPG", tree = "RNaseH2A.jpg"),
  RNASEH2B = list(domain = "RNaseH2B.JPG", tree = "RNaseH2B.jpg"),
  RNASEH2C = list(domain = "RNASEH2C.png", tree = "RNaseH2C.jpg"),
  AGO2 = list(domain = "Ago.JPG", tree = "AGO2.jpg"),
  RNASET2 = list(domain = "RNASET2.png", tree = "RNaseT2.jpg"),
  DICER1 = list(domain = "Dicer.JPG", tree = "DICER1.jpg"),
  DIS3L2 = list(domain = "Dis3L.JPG", tree = "DIS3L2.jpg"),
  ELAC2 = list(domain = "Elac2.JPG", tree = "ELAC2.jpg"),
  PRORP = list(domain = "PRORP.png", tree = "PRORP.jpg"),
  PARN = list(domain = "PARN.png", tree = "PARN.jpg")
)

feature_df <- function(rnase, organism_code, start, label, end = start) {
  data.frame(
    rnase = rnase,
    organism_code = organism_code,
    aligned_start = as.integer(start),
    aligned_end = as.integer(end),
    label = label,
    stringsAsFactors = FALSE
  )
}

features_catalog <- do.call(rbind, list(
  feature_df("RNASEH1", "HS", 227, "V142I"),
  feature_df("RNASEH1", "HS", 242, "R157X"),
  feature_df("RNASEH1", "HS", 271, "A185V"),

  feature_df("RNASEH2A", "HS", 93, "G37S"),
  feature_df("RNASEH2A", "HS", 78, "R25R (R[CGC]>R[CGT])"),
  feature_df("RNASEH2A", "HS", 324, "R235Q"),
  feature_df("RNASEH2A", "HS", 73, "V23V (V[GTG]>V[GTA])"),
  feature_df("RNASEH2A", "HS", 262, "R186W"),
  feature_df("RNASEH2A", "HS", 288, "N212I"),
  feature_df("RNASEH2A", "SC", 93, "G42S"),
  feature_df("RNASEH2A", "MM", 93, "G37S"),

  feature_df("RNASEH2B", "HS", 203, "A177T"),
  feature_df("RNASEH2B", "HS", 211, "V185G"),
  feature_df("RNASEH2B", "MM", 203, "A174T"),
  feature_df("RNASEH2B", "MM", 234, "E202X"),

  feature_df("RNASEH2C", "HS", 106, "R69W"),
  feature_df("RNASEH2C", "HS", 182, "K143I"),

  feature_df("AGO2", "HS", 369, "L192P"),
  feature_df("AGO2", "HS", 535, "T357M"),
  feature_df("AGO2", "HS", 542, "M364T"),
  feature_df("AGO2", "HS", 945, "C751Y"),
  feature_df("AGO2", "HS", 927, "G733R"),
  feature_df("AGO2", "HS", 330, "F182del"),
  feature_df("AGO2", "CE", 330, "F180del"),

  feature_df("DICER1", "HS", 2248, "D1713V"),
  feature_df("DICER1", "HS", 2244, "D1709Y"),
  feature_df("DICER1", "HS", 316, "E292fs"),
  feature_df("DICER1", "HS", 1069, "I813_Y819del PAZ", 1075),
  feature_df("DICER1", "HS", 1095, "S839F"),
  feature_df("DICER1", "HS", 2115, "L1583R"),
  feature_df("DICER1", "HS", 599, "E503X"),
  feature_df("DICER1", "HS", 1206, "R944X"),
  feature_df("DICER1", "HS", 1054, "T798fs"),
  feature_df("DICER1", "HS", 640, "R544X"),
  feature_df("DICER1", "HS", 1835, "L1303VfsX4"),
  feature_df("DICER1", "HS", 1632, "Y1204LfsX29"),
  feature_df("DICER1", "HS", 1143, "EX18del", 1193),
  feature_df("DICER1", "DM", 1205, "Q1147X"),
  feature_df("DICER1", "HS", 2365, "S1826X"),

  feature_df("ELAC2", "HS", 264, "R211X"),
  feature_df("ELAC2", "HS", 678, "T520I"),
  feature_df("ELAC2", "HS", 198, "F154L"),
  feature_df("ELAC2", "HS", 532, "L423F"),
  feature_df("ELAC2", "HS", 271, "S217L"),
  feature_df("ELAC2", "HS", 700, "A541T"),
  feature_df("ELAC2", "HS", 707, "H548Afs"),
  feature_df("ELAC2", "HS", 987, "R781H"),
  feature_df("ELAC2", "HS", 805, "E622V"),
  feature_df("ELAC2", "MM", 700, "A537T"),
  feature_df("ELAC2", "DM", 198, "F155L"),
  feature_df("ELAC2", "DM", 678, "T494I"),
  feature_df("ELAC2", "SC", 678, "T513I"),

  feature_df("DIS3L2", "HS", 470, "82.8-KB del EX6del", 580),
  feature_df("DIS3L2", "HS", 707, "22-KB del EX9del RNB", 767),
  feature_df("DIS3L2", "HS", 891, "C489Y"),
  feature_df("DIS3L2", "HS", 1167, "EX19del", 1205),
  feature_df("DIS3L2", "DM", 7, "V7GfsX10"),
  feature_df("DIS3L2", "MM", 707, "22-KB del EX10del", 767),

  feature_df("RNASET2", "HS", 99, "2.5-KB del", 113),
  feature_df("RNASET2", "HS", 36, "15-BP del", 40),
  feature_df("RNASET2", "HS", 372, "C184R"),
  feature_df("RNASET2", "HS", 377, "Q189Q Q[CAG]>Q[CAA]"),
  feature_df("RNASET2", "CE", 2219, "G119E"),
  feature_df("RNASET2", "CE", 136, "P55X"),

  feature_df("PARN", "HS", 403, "A383V"),
  feature_df("PARN", "HS", 295, "PARTIAL ND2 del", 320),
  feature_df("PARN", "HS", 295, "G281TfsX4"),
  feature_df("PARN", "HS", 302, "N288KfsX23"),
  feature_df("PARN", "HS", 222, "PARTIAL R3H del", 234),
  feature_df("PARN", "HS", 363, "R349W"),
  feature_df("PARN", "HS", 321, "D307VfsX22"),
  feature_df("PARN", "HS", 321, "EX14_18del", 437),
  feature_df("PARN", "HS", 190, "Q177X"),
  feature_df("PARN", "HS", 201, "I188IfsX7 CAF1"),
  feature_df("PARN", "HS", 445, "K421R"),

  feature_df("PRORP", "HS", 549, "A485V"),
  feature_df("PRORP", "HS", 470, "N412S"),
  feature_df("PRORP", "HS", 498, "A434D"),
  feature_df("PRORP", "HS", 509, "R445Q"),
  feature_df("PRORP", "HS", 458, "S400IfsX6"),
  feature_df("PRORP", "HS", 445, "T387A"),
  feature_df("PRORP", "HS", 472, "A414V"),
  feature_df("PRORP", "DM", 207, "Y121D"),
  feature_df("PRORP", "DM", 586, "W465R")
))

classify_variant <- function(label) {
  x <- tolower(label)
  if (grepl("del|kb|bp|ex", x)) return("Deletion / structural")
  if (grepl("fs", x)) return("Frameshift")
  if (grepl("x", x)) return("Nonsense / truncating")
  if (grepl("\\[.*>.*\\]", x) || grepl("^[a-z][0-9]+[a-z]$", x)) return("Synonymous / coding")
  "Missense / residue change"
}

features_catalog$variant_type <- vapply(features_catalog$label, classify_variant, character(1))
features_catalog$span <- ifelse(
  features_catalog$aligned_start == features_catalog$aligned_end,
  as.character(features_catalog$aligned_start),
  paste0(features_catalog$aligned_start, "-", features_catalog$aligned_end)
)

rnase_summary <- list(
  RNASEH1 = list(
    title = "RNase H1",
    biological_function = "Removes RNA primers and RNA:DNA hybrids, with a major mitochondrial DNA maintenance component.",
    disease = "Chronic progressive external ophthalmoplegia with mitochondrial DNA deletions.",
    ux_hint = "Useful for inspecting whether patient variants fall in conserved regions involved in RNA:DNA substrate processing."
  ),
  RNASEH2A = list(
    title = "RNase H2 subunit A",
    biological_function = "Catalytic subunit of RNase H2, which removes ribonucleotides misincorporated into DNA.",
    disease = "Aicardi-Goutieres syndrome.",
    ux_hint = "Prioritize conserved residues near catalytic or complex-stability regions when interpreting variants."
  ),
  RNASEH2B = list(
    title = "RNase H2 subunit B",
    biological_function = "Structural subunit of RNase H2, supporting genome integrity through ribonucleotide excision repair.",
    disease = "Aicardi-Goutieres syndrome.",
    ux_hint = "Variant clustering can suggest complex-stability effects rather than direct catalytic loss."
  ),
  RNASEH2C = list(
    title = "RNase H2 subunit C",
    biological_function = "Structural subunit of RNase H2 required for efficient ribonucleotide removal from DNA.",
    disease = "Aicardi-Goutieres syndrome.",
    ux_hint = "Compare across vertebrates and yeast carefully because structural subunits can diverge in local loops."
  ),
  AGO2 = list(
    title = "Argonaute 2",
    biological_function = "Core effector of miRNA and siRNA pathways; catalyzes slicing for highly complementary targets.",
    disease = "Lessel-Kreienkamp syndrome and related neurodevelopmental phenotypes.",
    ux_hint = "Inspect whether substitutions lie in conserved PAZ/MID/PIWI-associated regions or RNA-binding surfaces."
  ),
  DICER1 = list(
    title = "Dicer 1",
    biological_function = "RNase III enzyme that processes precursor miRNAs and other double-stranded RNA substrates.",
    disease = "DICER1 tumor predisposition syndrome, GLOW syndrome, and related disorders.",
    ux_hint = "Large protein: windowed inspection is more informative than global alignment alone."
  ),
  DIS3L2 = list(
    title = "DIS3-like exonuclease 2",
    biological_function = "3' to 5' exonuclease involved in decay of uridylated RNAs.",
    disease = "Perlman syndrome and Wilms tumor predisposition.",
    ux_hint = "Deletion spans and RNB-domain variants should be interpreted differently from single-residue substitutions."
  ),
  ELAC2 = list(
    title = "ELAC2 / RNase Z",
    biological_function = "tRNA 3' trailer endonuclease with nuclear and mitochondrial RNA-processing roles.",
    disease = "Combined oxidative phosphorylation deficiency 17 and hereditary prostate cancer associations.",
    ux_hint = "Cross-species variants are particularly useful here because yeast, fly, mouse, and human models exist."
  ),
  RNASET2 = list(
    title = "RNase T2",
    biological_function = "Endolysosomal ribonuclease involved in RNA turnover and innate immune homeostasis.",
    disease = "Cystic leukoencephalopathy.",
    ux_hint = "Distinguish catalytic-site cysteine/histidine-proximal changes from broad deletion events."
  ),
  PARN = list(
    title = "PARN",
    biological_function = "Poly(A)-specific ribonuclease involved in mRNA deadenylation and telomere biology.",
    disease = "Dyskeratosis congenita, pulmonary fibrosis, and telomere biology disorders.",
    ux_hint = "Domain-aware inspection is important because catalytic, R3H, and RNA-recognition regions differ."
  ),
  PRORP = list(
    title = "PRORP / proteinaceous RNase P",
    biological_function = "Protein-only RNase P enzyme involved in tRNA 5' end maturation.",
    disease = "Candidate disease-associated RNase; variants can be explored in a comparative context.",
    ux_hint = "Use residue conservation across eukaryotic orthologs to prioritize experimental follow-up."

  )
)

clinvar_variant_url <- function(rnase, label) {
  query <- paste(rnase, gsub("\\s+", " ", label))
  paste0("https://www.ncbi.nlm.nih.gov/clinvar/?term=", utils::URLencode(query, reserved = TRUE))
}

variant_links_html <- function(rnase) {
  feats <- features_catalog[features_catalog$rnase == rnase, , drop = FALSE]
  if (nrow(feats) == 0) return("No curated variants are currently listed for this RNase.")

  links <- vapply(seq_len(nrow(feats)), function(i) {
    label <- feats$label[i]
    organism <- feats$organism_code[i]
    span <- feats$span[i]
    sprintf(
      "<a href='%s' target='_blank'><strong>%s</strong></a> <span class='small-note'>(%s, aln. %s)</span>",
      clinvar_variant_url(rnase, label), label, organism, span
    )
  }, character(1))
  paste(links, collapse = ", ")
}

rnase_info_html <- function(rnase) {
  stopifnot(rnase %in% rnase_order)
  info <- switch(
    rnase,
    RNASEH1 = paste0(
      "<strong><em>RNASEH1</em></strong> encodes ribonuclease H1, an enzyme involved in removal of RNA primers and RNA:DNA hybrids, ",
      "with an important role in mitochondrial DNA replication and maintenance. It is linked to ",
      "<a href='https://www.omim.org/entry/616479' target='_blank'><strong>progressive external ophthalmoplegia with mitochondrial DNA deletions, autosomal recessive 2</strong></a>."
    ),
    RNASEH2A = paste0(
      "<strong><em>RNASEH2A</em></strong> encodes the catalytic subunit of the RNase H2 complex, which removes ribonucleotides misincorporated into DNA. ",
      "It is associated with <a href='https://www.omim.org/entry/610333' target='_blank'><strong>Aicardi-Goutieres syndrome 4</strong></a>."
    ),
    RNASEH2B = paste0(
      "<strong><em>RNASEH2B</em></strong> encodes a non-catalytic structural subunit of RNase H2. Its variants are associated with ",
      "<a href='https://www.omim.org/entry/610181' target='_blank'><strong>Aicardi-Goutieres syndrome 2</strong></a>, often through impaired genome maintenance and chronic DNA-damage signaling."
    ),
    RNASEH2C = paste0(
      "<strong><em>RNASEH2C</em></strong> encodes a structural RNase H2 subunit required for efficient removal of ribonucleotides from DNA. ",
      "It is linked to <a href='https://www.omim.org/entry/610329' target='_blank'><strong>Aicardi-Goutieres syndrome 3</strong></a>."
    ),
    AGO2 = paste0(
      "<strong><em>AGO2</em></strong> encodes Argonaute 2, the catalytic effector of small-RNA-guided silencing. ",
      "Disease-associated variants are linked to <a href='https://www.omim.org/entry/619149' target='_blank'><strong>Lessel-Kreienkamp syndrome</strong></a> and related neurodevelopmental phenotypes."
    ),
    DICER1 = paste0(
      "<strong><em>DICER1</em></strong> encodes an RNase III enzyme required for precursor miRNA and double-stranded RNA processing. ",
      "Its variants are associated with <a href='https://www.omim.org/entry/606241' target='_blank'><strong>DICER1 syndrome</strong></a> and several tumor predisposition phenotypes."
    ),
    DIS3L2 = paste0(
      "<strong><em>DIS3L2</em></strong> encodes a cytoplasmic 3' to 5' exonuclease involved in degradation of uridylated RNAs. ",
      "Loss-of-function variants are linked to <a href='https://www.omim.org/entry/267000' target='_blank'><strong>Perlman syndrome</strong></a> and Wilms tumor predisposition."
    ),
    ELAC2 = paste0(
      "<strong><em>ELAC2</em></strong>, also known as RNase Z, encodes a tRNA 3' trailer endonuclease with nuclear and mitochondrial RNA-processing functions. ",
      "Variants have been connected to mitochondrial translation/OXPHOS defects and hereditary prostate cancer susceptibility contexts."
    ),
    RNASET2 = paste0(
      "<strong><em>RNASET2</em></strong> encodes an endolysosomal ribonuclease involved in RNA turnover and immune homeostasis. ",
      "Biallelic variants are associated with <a href='https://www.omim.org/entry/612951' target='_blank'><strong>cystic leukoencephalopathy without megalencephaly</strong></a>."
    ),
    PARN = paste0(
      "<strong><em>PARN</em></strong> encodes a poly(A)-specific ribonuclease involved in mRNA deadenylation and telomere biology. ",
      "Variants are associated with telomere biology disorders, including dyskeratosis congenita and pulmonary fibrosis."
    ),
    PRORP = paste0(
      "<strong><em>PRORP</em></strong> encodes a proteinaceous RNase P enzyme involved in tRNA 5' end maturation. ",
      "The comparative view is useful for prioritizing candidate residues by conservation across available eukaryotic orthologs."
    )
  )

  paste0(
    info,
    "<br><br><strong>Marked variants / regions in RNaseViz:</strong> ",
    variant_links_html(rnase),
    "<br><br><span class='small-note'>Click the alignment or enter a reference residue to inspect conservation at the aligned coordinate.</span>"
  )
}

rnase_config <- function(rnase) {
  stopifnot(rnase %in% rnase_order)
  feats <- features_catalog[features_catalog$rnase == rnase, , drop = FALSE]
  organisms <- unique(c("Homo sapiens", organism_lookup$organism[organism_lookup$code %in% unique(feats$organism_code)]))
  list(
    rnase = rnase,
    title = rnase_summary[[rnase]]$title,
    summary = rnase_summary[[rnase]],
    features = feats,
    features_seqvisr = lapply(seq_len(nrow(feats)), function(i) {
      row <- feats[i, ]
      if (row$aligned_start == row$aligned_end) {
        c(row$organism_code, row$aligned_start, row$label)
      } else {
        c(row$organism_code, row$aligned_start:row$aligned_end, row$label)
      }
    }),
    colors = rnase_palette[seq_len(max(1, nrow(feats)))],
    organisms = organisms,
    msa_seqvisr = file.path("data", "alignments", sprintf("reordered_for_msavisr_%sAln.fa", rnase)),
    msa_ggmsa = file.path("data", "alignments", sprintf("reordered_for_ggmsa_%sAln.fa", rnase)),
    raw_seq = file.path("data", "alignments", sprintf("%s_raw_seq.rds", rnase)),
    tree = file.path("data", "trees", sprintf("%s_tree.rds", rnase)),
    tree_image = file.path("data", "trees", rnase_asset_map[[rnase]]$tree),
    tree_image_url = paste0("rnase-trees/", rnase_asset_map[[rnase]]$tree),
    image = file.path("data", "domains", rnase_asset_map[[rnase]]$domain),
    image_url = paste0("rnase-domains/", rnase_asset_map[[rnase]]$domain),
    legacy_image = file.path("www", "images", sprintf("%s.png", rnase))
  )
}

get_organism_code <- function(organism) {
  organism_lookup$code[match(organism, organism_lookup$organism)]
}

get_organism_rds_name <- function(organism) {
  organism_lookup$rds_name[match(organism, organism_lookup$organism)]
}
