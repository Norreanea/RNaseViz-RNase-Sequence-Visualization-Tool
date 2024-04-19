# Load libraries
#install.packages("rsvg")
library(shiny)
library(seqvisr)
library(msa)
library(ggimage)
library(ggtree)
library(seqinr)
library(ape)
library(ggmsa)
library(plotly)
library(ggiraph)
library(Biostrings)
library(ggrepel)
library(rsvg)


# Load all msa and annotation
inpmsa_RNASEH1 <- "./data/reordered_for_msavisr_RNASEH1Aln.fa"
inpmsa_RNASEH1_ggmsa <- "./data/reordered_for_ggmsa_RNASEH1Aln.fa"
multfeatures_RNASEH1 <- list(c("HS", 227,"V142I"), 
                             c("HS", 242,"R157X"),
                             c("HS", 271,"A185V"))
inpmsa_RNASEH2A <- "./data/reordered_for_msavisr_RNASEH2AAln.fa"
inpmsa_RNASEH2A_ggmsa <- "./data/reordered_for_ggmsa_RNASEH2AAln.fa"
multfeatures_RNASEH2A <- list(c("HS", 93,"G37S"), 
                              c("HS", 78,"R25R (R[CGC]>R[CGT])"),
                              c("HS", 324,"R235Q"),
                              c("HS", 73,"V23V (V[GTG]>V[GTA])"),
                              c("HS", 262,"R186W"),
                              c("HS", 288,"N212I"),
                              c("SC",93,"G37S"),
                              c("MM",93,"G37S"))
inpmsa_RNASEH2B <- "./data/reordered_for_msavisr_RNASEH2BAln.fa"
inpmsa_RNASEH2B_ggmsa <- "./data/reordered_for_ggmsa_RNASEH2BAln.fa"
multfeatures_RNASEH2B <- list(c("HS", 203,"A177T"), 
                              c("HS", 211,"V185G"),
                              c("MM",203,"A177T"))
inpmsa_RNASEH2C <- "./data/reordered_for_msavisr_RNASEH2CAln.fa"
inpmsa_RNASEH2C_ggmsa <- "./data/reordered_for_ggmsa_RNASEH2CAln.fa"
multfeatures_RNASEH2C <- list(c("HS", 106,"R69W"), 
                              c("HS", 182,"K143I"))
inpmsa_AGO2 <- "./data/reordered_for_msavisr_AGO2Aln.fa"
inpmsa_AGO2_ggmsa <- "./data/reordered_for_ggmsa_AGO2Aln.fa"
multfeatures_AGO2 <- list(c("HS", 605,"L192P"), 
                          c("HS", 774,"T357M"),
                          c("HS", 781,"M364T"),
                          c("HS", 1181,"C751Y"),
                          c("HS", 1163,"G733R"),
                          c("HS", 566,"F182del"),
                          c("MM",566,"F182del"))
inpmsa_DICER1 <- "./data/reordered_for_msavisr_DICER1Aln.fa"
inpmsa_DICER1_ggmsa <- "./data/reordered_for_ggmsa_DICER1Aln.fa"
multfeatures_DICER1 <- list(c("HS", 2248,"D1713V"), 
                            c("HS", 2244,"D1709Y"),
                            c("HS", 316,"E292fs"),
                            c("HS", c(1069:1075),"2457C-G DEL"),
                            c("HS", 1095,"S839F"),
                            c("HS", 2115,"L1583R"),
                            c("HS", 599,"E503X"),
                            c("HS", 1206,"R944X"),
                            c("HS", 1054,"T798fs"),
                            c("HS", 640,"R544X "),
                            c("HS", 1835,"L1303VfsX4"),
                            c("HS", 1632,"Y1204LfsTer29"),
                            c("HS", c(1143:1193),"EX18DEL"),
                            c("DM",1205,"Q1147X"))
inpmsa_ELAC2 <- "./data/reordered_for_msavisr_ELAC2Aln.fa"
inpmsa_ELAC2_ggmsa <- "./data/reordered_for_ggmsa_ELAC2Aln.fa"
multfeatures_ELAC2 <- list(c("HS",264 ,"R211X"),
                           c("HS",678 ,"T520I"), 
                           c("HS",198,"F154L"), 
                           c("HS",532,"L423F"), 
                           c("HS",271,"S217L"), 
                           c("HS",700,"A541T"), 
                           c("HS",707,"H548Afs"), 
                           c("HS",987,"R781H"), 
                           c("HS",805,"E622V"), 
                           c("MM",700,"A541T"),
                           c("DM",198,"F154L"),
                           c("DM",678 ,"T520I"), 
                           c("SC",678 ,"T520I") 
)
inpmsa_DIS3L2 <- "./data/reordered_for_msavisr_DIS3L2Aln.fa"
inpmsa_DIS3L2_ggmsa <- "./data/reordered_for_ggmsa_DIS3L2Aln.fa"
multfeatures_DIS3L2 <- list(c("HS", c(470:580),"82.8-KB DEL EX6DEL"), 
                            c("HS", c(707:767),"22-KB DEL EX9DEL"),
                            c("HS", 891,"C489Y"),
                            c("HS", c(1167:1205),"EX19DEL"),
                            c("DM",7 ,"V7GfsX10"), 
                            c("MM",c(707:767),"22-KB DEL EX9DEL"))
inpmsa_RNASET2 <- "./data/reordered_for_msavisr_RNASET2Aln.fa"
inpmsa_RNASET2_ggmsa <- "./data/reordered_for_ggmsa_RNASET2Aln.fa"
multfeatures_RNASET2 <- list(c("HS", c(99:113),"2.5-KB DEL"), 
                             c("HS", c(36:40),"15-BP DEL"),
                             c("HS", 372,"C184R"),
                             c("HS", 377,"Q189Q Q[CAG]>Q[CAA]"),
                             c("CE", 2219,"G119E"),
                             c("CE", 136,"P55X"))
inpmsa_PARN <- "./data/reordered_for_msavisr_PARNAln.fa"
inpmsa_PARN_ggmsa <- "./data/reordered_for_ggmsa_PARNAln.fa"
multfeatures_PARN  <- list(c("HS",403,"A383V"), 
                             c("HS",c(295:320),"PARTIAL EX13DEL"),
                             c("HS",295,"G281TfsX4"),
                             c("HS",302,"N288KfsX23"),
                             c("HS",c(222:234),"PARTIAL R3H_DEL"),
                             c("HS",363,"R349W R307VfsX22"),
                             c("HS",c(321:437),"22-KB DEL"),
                             c("HS",190,"Q177X"),
                             c("HS",201,"I188IfsX7 CAF1"),
                             c("HS",445,"K421R"))
inpmsa_PRORP <- "./data/reordered_for_msavisr_PRORPAln.fa"
inpmsa_PRORP_ggmsa <- "./data/reordered_for_ggmsa_PRORPAln.fa"
multfeatures_PRORP <- list(c("HS",549,"A485V"), 
                           c("HS",470,"N412S"),
                           c("HS",498,"A434D"),
                           c("HS",509,"R445Q"),
                           c("HS",458,"S400IfsX6"),
                           c("HS",445,"T387A"),
                           c("HS",472,"A414V"),
                           c("DM",207,"Y121D"),
                           c("DM",586,"W465R"))



#Identify positions for ALL aminoacids in updated aligned seq
mapAminoAcidPositions <- function(original_seq, aligned_seq) {
  # Convert the sequences into lists of characters (amino acids)
  original_aas <- unlist(strsplit(original_seq, ""))
  aligned_aas <- unlist(strsplit(aligned_seq, ""))
  # Initialize variables for tracking positions and amino acids
  original_positions <- numeric(0)
  aligned_positions <- numeric(0)
  original_seq_amino <- character(0)
  aligned_seq_amino <- character(0)
  # Track the position in the original sequence
  original_pos_counter <- 1
  # Loop through the aligned sequence
  for (i in 1:length(aligned_aas)) {
    if (aligned_aas[i] != "-") {
      if (original_pos_counter <= length(original_aas)) {
        original_positions <- c(original_positions, original_pos_counter)
        original_seq_amino <- c(original_seq_amino, original_aas[original_pos_counter])
        aligned_seq_amino <- c(aligned_seq_amino, aligned_aas[i])
        aligned_positions <- c(aligned_positions, i)
        original_pos_counter <- original_pos_counter + 1
      }
    }
  }
  # Prepare the result as a data frame
  result <- data.frame(
    Original_seq_amino = original_seq_amino,
    Original_position = original_positions,
    Aligned_seq_amino = aligned_seq_amino,
    Aligned_position = aligned_positions
  )
  return(result)
}

# Server logic
server <- function(input, output, session) {
  rv <- reactiveValues(display = FALSE)
  observeEvent(input$visualizeBtn, {
    rnase <- input$rnase  # get the selected RNAse
    img_path <- sprintf("images/%s.png", rnase)  
    # When the button is clicked, set the reactive value to TRUE
    rv$display <- TRUE
    # Based on selected RNase, load the corresponding MSA file and features
    rnase_data <- switch(input$rnase,
                         "RNASEH1" = list(file = inpmsa_RNASEH1, features = multfeatures_RNASEH1,
                                          colors = c("MediumBlue", "DarkRed", "DarkOliveGreen")),
                         "RNASEH2A" = list(file = inpmsa_RNASEH2A, features = multfeatures_RNASEH2A,
                                           colors=c("MediumBlue", "DarkRed", "DarkOliveGreen", "DarkMagenta", "SaddleBrown", 
                                                    "FireBrick", "Navy", "DarkCyan")),
                         "RNASEH2B" = list(file = inpmsa_RNASEH2B, features = multfeatures_RNASEH2B,
                                           colors=c("MediumBlue", "DarkRed", "DarkOliveGreen", "DarkMagenta")),
                         "RNASEH2C" = list(file = inpmsa_RNASEH2C, features = multfeatures_RNASEH2C,
                                           colors=c("MediumBlue", "DarkRed")),
                         "AGO2" = list(file = inpmsa_AGO2, features = multfeatures_AGO2,
                                       colors=c("MediumBlue", "DarkRed", "DarkOliveGreen", "DarkMagenta", "SaddleBrown", 
                                                "FireBrick", "Navy", "DarkCyan", "Yellow3", "Maroon", "DarkSlateBlue")),
                         "DICER1" = list(file = inpmsa_DICER1, features = multfeatures_DICER1,
                                         colors=c("MediumBlue", "DarkRed", "DarkOliveGreen", "DarkMagenta","SaddleBrown", 
                                                  "FireBrick", "Navy", "DarkCyan", "Yellow3", "Maroon", "DarkSlateBlue", "ForestGreen","DarkGoldenrod")),
                         "ELAC2" = list(file = inpmsa_ELAC2, features = multfeatures_ELAC2,
                                        colors=c("MediumBlue", "DarkRed", "DarkOliveGreen", "DarkMagenta", "SaddleBrown", 
                                                 "FireBrick","Navy", "DarkCyan", "Yellow3", "Maroon", "DarkSlateBlue", "ForestGreen", "DarkGoldenrod") ),
                         "DIS3L2" = list(file = inpmsa_DIS3L2, features = multfeatures_DIS3L2,
                                         colors=c("MediumBlue", "DarkRed", "DarkOliveGreen", "DarkMagenta", "SaddleBrown", 
                                                  "FireBrick")),
                         "RNASET2" = list(file = inpmsa_RNASET2, features = multfeatures_RNASET2,
                                          colors=c("MediumBlue", "DarkRed", "DarkOliveGreen", "DarkMagenta", "SaddleBrown", 
                                                   "FireBrick", "Navy", "DarkCyan")),
                         "PARN" = list(file = inpmsa_PARN, features = multfeatures_PARN,
                                       colors=c("MediumBlue", "DarkRed", "DarkOliveGreen", "DarkMagenta", "SaddleBrown", 
                                                "FireBrick", "Navy", "DarkCyan", "Yellow3", "Maroon", "DarkSlateBlue")),
                         "PRORP" = list(file = inpmsa_PRORP, features = multfeatures_PRORP,
                                        colors=c("MediumBlue", "DarkRed", "DarkOliveGreen", "DarkMagenta", "SaddleBrown", 
                                                 "FireBrick", "Navy", "DarkCyan", "Yellow3", "Maroon", "DarkSlateBlue") ) )
    
    # rnase_tree <- switch(input$rnase,
    #                      "RNASEH1" = RNASEH1Tree_plot,
    #                      "RNASEH2A" = RNASEH2ATree_plot,
    #                      "RNASEH2B" = RNASEH2BTree_plot,
    #                      "RNASEH2C" = RNASEH2CTree_plot,
    #                      "AGO2" = AGO2Tree_plot,
    #                      "DICER1" = DICER1Tree_plot,
    #                      "ELAC2" = ELAC2Tree_plot,
    #                      "DIS3L2" = DIS3L2Tree_plot,
    #                      "RNASET2" = RNASET2Tree_plot,
    #                      "PARN" = PARNTree_plot,
    #                      "PRORP" = PRORPTree_plot)
    rnase_tree <- readRDS(sprintf("./data/%s_tree.rds", input$rnase))
    #rnase_info_html <- "<strong>Test</strong>: This is a <a href='https://example.com' target='_blank'>link</a>."
    
    rnase_info_html  <- switch(input$rnase,
                               "RNASEH1"=paste0("<strong>RNASEH1</strong> plays a key role in DNA replication and repair, 
                                                encoding ribonuclease H1 crucial for removing RNA primers from DNA. 
                                                It's linked to <a href='https://www.omim.org/entry/616479' target='_blank'><strong>Progressive external 
                                                ophthalmoplegia with mitochondrial DNA deletions, autosomal 
                                                recessive 2 (PEOB2)</strong></a>. Mutations include <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000412621.2/?redir=rcv' 
                                                target='_blank'><strong>V142I</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000412498.1/?redir=rcv' 
                                                target='_blank'><strong>R157X</strong></a>, and <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000412557.1/?redir=rcv' 
                                                target='_blank'><strong>A185V</strong></a>. Gene deficiency leads to mitochondrial dysfunction characterized by disorganized cristae and                    
                                                vacuolation, retained RNA primers causing stalled replication forks, and subsequent apoptosis, as observed in RNase H1-deficient mouse models. These 
                                                molecular disturbances result in decreased mtDNA and accumulation of mtDNA fragments, impaired oxidative phosphorylation (OXPHOS) system, and 
                                                disrupted ATP production, collectively contributing to the mitochondrial pathogenesis underlying PEOB2"),
                               
                               "RNASEH2A"=paste0("<strong>RNASEH2A</strong>, essential for DNA repair, encodes a catalytic subunit of the ribonuclease H2 complex. It's associated 
                                                 with <a href='https://www.omim.org/entry/610333' target='_blank'><strong>Aicardi-Goutieres syndrome 4</strong></a>. 
                                                 The gene carries mutations such as <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000004904.5/?redir=rcv' 
                                                 target='_blank'><strong>G37S</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000056302.8/?redir=rcv' 
                                                 target='_blank'><strong>R25R R[CGC]>R[CGT]</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000056303.6/?redir=rcv'
                                                 target='_blank'><strong>R235Q</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000056304.9/?redir=rcv' 
                                                 target='_blank'><strong>V23V V[GTG]>V[GTA]</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000056305.8/?redir=rcv' 
                                                 target='_blank'><strong>R186W</strong></a>, and <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001304321.4/?redir=rcv' 
                                                 target='_blank'><strong>N212I</strong></a>. Mouse models with the G37S mutation exhibit perinatal lethality and activation of the innate 
                                                 immune response through the cGAS-STING pathway, driven by gene`s failure to adequately remove misincorporated 
                                                 ribonucleotides and RNA-DNA hybrids. This leads to a DNA damage response and cellular senescence, as observed in
                                                 homozygous G37S mice. Additionally, yeast cells harboring the G42S (equivalent to G37S in humans) mutation show increased ribonucleotide incorporation into 
                                                 the nuclear genome, underscoring the enzyme's role in genome integrity."),
                               
                               "RNASEH2B"=paste0("<strong>RNASEH2B</strong>, involved in DNA replication and repair, is the non-catalytic structural subunit of the ribonuclease H2 complex. It's implicated 
                                                 in <a href='https://www.omim.org/entry/610181' target='_blank'><strong>Aicardi-Goutieres syndrome 2</strong></a>, with
                                                 mutations <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000001324.36/?redir=rcv' target='_blank'><strong>A177T</strong></a>, 
                                                 <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000001325.6/?redir=rcv' target='_blank'><strong>V185G</strong></a> among others.
                                                 Mouse models with the A174T mutation (corresponding to human A177T) show no immune phenotype, but RNASEH2B E202X mutants are characterized by impaired 
                                                 proliferation and accumulation in the G2/M phase of the cell cycle. This is attributed to chronic activation of the DNA damage response, evidenced by 
                                                 increased single-strand breaks, elevated histone H2AX phosphorylation, and activation of p53 target genes such as the cyclin-dependent kinase inhibitor 1 (p21). 
                                                 The deficiency in RNASEH2B activity leads to heightened ribonucleotide incorporation into DNA, triggering spontaneous DNA breaks and subsequent
                                                 DNA damage responses, underlying its critical role in maintaining genomic stability."),
                               
                               "RNASEH2C"=paste0("<strong>RNASEH2C</strong> is the non-catalytic structural subunit of the ribonuclease H2 complex, essential for the removal of ribonucleotides from DNA. Linked to 
                                                 <a href='https://www.omim.org/entry/610329' target='_blank'><strong>Aicardi-Goutieres syndrome 3</strong></a>, mutations 
                                                 include <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000001322.14/?redir=rcv' target='_blank'><strong>R69W</strong></a> 
                                                 and <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000001323.5/?redir=rcv' target='_blank'><strong>K143I</strong></a>. 
                                                 RNASEH2C in mice results in embryonic lethality due to severe DNA damage from ribonucleotide misincorporation, triggering an extensive DNA 
                                                 damage response and leading to significant developmental anomalies and early embryonic death. Additionally, despite the RNASEH2C ortholog 
                                                 not being designated as trusted by the Alliance for Genome Resources, our analysis utilized this yeast gene. In yeast, the K46W mutation 
                                                 (analogous to human R69W) shows increased frequencies of cytosine and guanine ribonucleotides (rC and rG), notably accumulating rCMP at 
                                                 specific genomic sites, highlighting the challenges in efficiently removing rCMP from DNA."),
                               
                               "AGO2"=paste0("<strong>AGO2</strong> is essential for miRNA and siRNA-guided mRNA cleavage. Associated with
                                             <a href='https://www.omim.org/entry/619149' target='_blank'><strong>Lessel-Kreienkamp syndrome (LESKRES)</strong></a>, mutations 
                                             include <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001289983.1/?redir=rcv' target='_blank'><strong>L192P</strong></a>, 
                                             <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001289984.2/?redir=rcv' target='_blank'><strong>T357M</strong></a>, 
                                             <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001289985.1/?redir=rcv' target='_blank'><strong>M364T</strong></a>, 
                                             <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001289986.1/?redir=rcv' target='_blank'><strong>C751Y</strong></a>, 
                                             <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001289988.1/?redir=rcv' target='_blank'><strong>G733R</strong></a>, 
                                             and <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001289987.1/?redir=rcv' target='_blank'><strong>F182del</strong></a>. 
                                             Similarity in mutation sites between human AGO2 and <em>C. elegans</em> AGO1, specifically the F180del  in AGO1 and its analogous F182del in hAGO2, allows us to infer the impact of 
                                             AGO2 mutations through studies on AGO1. This homology indicates that mutations like F182del in AGO2 could lead to enhanced mRNA binding and stalling at miRNA-mRNA hybrids, impairing 
                                             RISC loading and functionality. Such disruptions are suggestive of the molecular dysfunctions that contribute to the neurological and developmental symptoms observed in LESKRES, 
                                             including defective miRNA silencing, altered miRNA profiles, and increased sequestration of AGO2 in P-bodies, which particularly affects neural translation."),
                               "DICER1"=paste0("DICER1 is crucial for microRNA processing, linked to disorders like <a href='https://www.omim.org/entry/601200' 
                                               target='_blank'><strong>Pleuropulmonary blastoma</strong></a>, with mutations such as <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000004725.4/?redir=rcv' 
                                               target='_blank'><strong>L1583R</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000004726.5/?redir=rcv' 
                                               target='_blank'><strong>E503X</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000004727.5/?redir=rcv' 
                                               target='_blank'><strong>R944X</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000240963.4/?redir=rcv' 
                                               target='_blank'><strong>T798fs</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000004729.5/?redir=rcv' 
                                               target='_blank'><strong>R544X</strong></a>, <a href='https://www.omim.org/entry/618272' target='_blank'><strong>Glow syndrome</strong></a>, 
                                               with mutations such as <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000735852.4/?redir=rcv' target='_blank'><strong>D1713V</strong></a>, 
                                               <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000735853.4/?redir=rcv' target='_blank'><strong>D1709Y</strong></a>, 
                                               <a href='https://www.omim.org/entry/138800' target='_blank'><strong>Goiter, multinodular 1, with or without sertoli-leydig cell tumors</strong></a>, 
                                               with mutations <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000023523.4/?redir=rcv' target='_blank'><strong>E292fs</strong></a>, 
                                               <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001398192.5/?redir=rcv' target='_blank'><strong>2457C-G DEL</strong></a>, 
                                               <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000023526.4/?redir=rcv' target='_blank'><strong>S839F</strong></a>, 
                                               <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000023524.4/?redir=rcv' target='_blank'><strong>EX18DEL</strong></a> 
                                               and <a href='https://www.omim.org/entry/180295' target='_blank'><strong>Rhabdomyosarcoma, embryonal, 2</strong></a>, 
                                               with mutations including <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000056332.5/?redir=rcv' target='_blank'><strong>L1303VfsX4</strong></a>, 
                                               <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000056333.4/?redir=rcv' target='_blank'><strong>Y1204LfsTer29</strong></a>. 
                                               In several animal models, DICER1 deficiency leads to severe developmental defects and embryonic lethality, primarily through disrupted miRNA maturation. 
                                               This disruption results in the accumulation of 3p-miRNA molecules as observed in human heterozygous mutations. These miRNAs may aberrantly regulate
                                               the PI3K-AKT-mTOR signaling pathway, which is crucial for cell growth and proliferation. Alterations in this pathway are linked to the 
                                               overgrowth observed in GLOW Syndrome and the tumorigenic processes in the aforementioned conditions."),
                               "ELAC2"=paste0("ELAC2, or <em>RNase Z</em>, is involved in tRNA processing and mitochondrial function. Associated with <a href='https://www.omim.org/entry/615440' 
                                              target='_blank'><strong>Combined oxidative phosphorylation deficiency 17 (COXPD17)</strong></a> (mutations including <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000056275.3/?redir=rcv'
                                              target='_blank'><strong>R211X</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000056274.4/?redir=rcv' target='_blank'><strong>T520I</strong></a>, 
                                              <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000056276.7/?redir=rcv' target='_blank'><strong>F154L</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000056277.4/?redir=rcv' 
                                              target='_blank'><strong>L423F</strong></a>), and <a href='https://www.omim.org/entry/614731' target='_blank'><strong>Prostate cancer hereditary 2 (HPC2)</strong></a> (mutations including 
                                              <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000005358.3/?redir=rcv' target='_blank'><strong>S217L</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000005359.3/?redir=rcv' 
                                              target='_blank'><strong>A541T</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000005360.2/?redir=rcv' target='_blank'><strong>H548Afs</strong></a>, 
                                              <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000005361.4/?redir=rcv' target='_blank'><strong>R781H</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000005362.4/?redir=rcv'
                                              target='_blank'><strong>E622V</strong></a>). Mutations in ELAC2 lead to the accumulation of unprocessed mt-tRNA molecules, disrupting mitochondrial protein synthesis and respiration, 
                                              and reducing OXPHOS complex levels. This results in decreased ATP production and increased reactive oxygen species (ROS), contributing to mitochondrial disorder in COXPD17 by affecting 
                                              energy production and inducing cellular stress. In HPC2, mutations impair both mitochondrial and nuclear RNA processing, altering cellular energy states and 
                                              activating pro-inflammatory and tumorigenesis pathways. Animal models, such as the mouse prostate-specific ELAC2 knockout and A537T mutation models (human equivalent A541T), demonstrate 
                                              how ELAC2 disruptions facilitate cancer progression by modifying cellular metabolism and stress responses, particularly through the altered processing of tRNAs and miRNAs."),
                               "RNASET2"=paste0("RNASET2 is involved in RNA turnover and lysosomal degradation, associated with <a href='https://www.omim.org/entry/612951' 
                                                target='_blank'><strong>cystic leukoencephalopathy</strong></a>, with mutations such as <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000000440.2/?redir=rcv' 
                                                target='_blank'><strong>C184R</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000000441.4/?redir=rcv' target='_blank'><strong>2.5-KB DEL</strong></a>, 
                                                <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000000442/?redir=rcv' target='_blank'><strong>IVS5AS A-G -2</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000000443/?redir=rcv'
                                                target='_blank'><strong>1-BP DEL 332+1G</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000000444.2/?redir=rcv' target='_blank'><strong>15-BP DEL</strong></a>,
                                                <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000000445.5/?redir=rcv' target='_blank'><strong>Q189Q Q[CAG]&gt;Q[CAA]</strong></a>. Mutations in RNASET2 cause accumulation of uncleaved ssRNA substrates, 
                                                notably rRNA in lysosomes, leading to lysosomal dysfunction and storage disorders, as observed in zebrafish knockouts and <em>C. elegans</em> (G119E and P55X mutations) mutants. 
                                                This disrupts cellular stress responses and contributes to leukoencephalopathy by triggering innate immune pathways, increasing IFN-1 and ISGs in mouse models, resulting in neuroinflammation 
                                                and impaired phagocytosis by activated microglia. Additionally, RNASET2 deficiency in yeast stabilizes tRNAs and prevents tRNA fragment formation, impairing oxidative stress responses and 
                                                inhibiting apoptosis."),
                               "DIS3L2"=paste0("DIS3L2, important for RNA decay, is implicated in <a href='https://www.omim.org/entry/267000' target='_blank'><strong>Perlman syndrome</strong></a> primarily through its role in degrading 3’-uridylated 7SL ncRNA, a component of the signal recognition particle (SRP). 
                                               Mutations include <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000024119.2/?redir=rcv' target='_blank'><strong>82.8-KB DEL EX6DEL</strong></a>,
                                               <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000024120.2/?redir=rcv' target='_blank'><strong>22-KB DEL EX9DEL</strong></a>, 
                                               <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000024121.6/?redir=rcv' target='_blank'><strong>C489Y</strong></a>, and 
                                               <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000024122.6/?redir=rcv' target='_blank'><strong>EX19DEL</strong></a>. 
                                               Loss-of-function mutations in DIS3L2 observed in humans, lead to the accumulation of uridylated 7SL ncRNA, impairing signal recognition particle (SRP) 
                                               function and subsequently ER-targeted translation. This disruption results in calcium leakage from the ER, disturbing calcium homeostasis and contributing to 
                                               the pathology of Perlman syndrome. Additionally, Dis3l2 EX10DEL in mice (EX9DEL in human DIS3L2) and Dis3l2 V7GfsX10 in flies, the failed degradation of 3’-uridylated-poly(A) 
                                               mRNAs leads to upregulation of the oncogene Igf2, promoting tumorigenesis."),
                               "PARN"=paste0("PARN encodes a poly(A)-specific ribonuclease affecting mRNA stability and telomere maintenance, associated with <a href='https://www.omim.org/entry/616353' 
                                             target='_blank'><strong>Dyskeratosis congenita, autosomal recessive 6</strong></a> (mutations include <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000170484.4/?redir=rcv' 
                                             target='_blank'><strong>A383V</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000170485.4/?redir=rcv' target='_blank'><strong>PARTIAL EX13DEL</strong></a>, 
                                             PARTIAL R3H_DEL and <a href='https://www.ncbi.nlm.nih.gov/clinvar/variation/190291/' target='_blank'><strong>N288KfsX23</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000203540.1/?redir=rcv' 
                                             target='_blank'><strong>R349W R307VfsX22</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000203556.2/?redir=rcv' target='_blank'><strong>22-KB DEL</strong></a>) 
                                             and <a href='https://www.omim.org/entry/616371' target='_blank'><strong>Pulmonary fibrosis and/or bone marrow failure syndrome, telomere-related, 4</strong></a>(mutations include 
                                             <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000170589.3/?redir=rcv' target='_blank'><strong>IVS4AS A-G -2</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000170590.2/?redir=rcv' 
                                             target='_blank'><strong>Q177X</strong></a>, <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000170591.2/?redir=rcv' target='_blank'><strong>I188IfsX7 CAF1</strong></a>, 
                                             <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV000170592.2/?redir=rcv' target='_blank'><strong>K421R</strong></a>). In humans, mutations in the PARN gene 
                                             result in unstable TERC RNA, leading to prematurely shortened telomeres. This reduction in telomere length initiates a DNA damage response that manifests as 
                                             cell cycle arrest, apoptosis, or premature senescence—key mechanisms underlying the pathology of these conditions. In animal models, specifically mice with 
                                             PARN knock-down and zebrafish with targeted PARN disruption, there is decreased RNA stability and impaired hematopoietic RNA processing, respectively. These
                                             defects hinder effective blood cell maturation, directly linking PARN dysfunction to the clinical features observed in these diseases."),
                               "PRORP"=paste0("PRORP is essential for mitochondrial RNA processing, related to <a href='https://www.omim.org/entry/619737' 
                                              target='_blank'><strong>Combined oxidative phosphorylation deficiency 54  (COXPD54)</strong></a>, with a wide array of mutations including 
                                              <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001824882.1/?redir=rcv' target='_blank'><strong>A485V</strong></a>, 
                                              <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001825014.1/?redir=rcv' target='_blank'><strong>N412S</strong></a>, 
                                              <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001825015.1/?redir=rcv' target='_blank'><strong>A434D</strong></a>, 
                                              <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001824883.1/?redir=rcv' target='_blank'><strong>R445Q</strong></a>, 
                                              <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV001824884.1/?redir=rcv' target='_blank'><strong>S400IfsX6</strong></a>, 
                                              <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV003892100.1/?redir=rcv' target='_blank'><strong>T387A</strong></a> and 
                                              <a href='https://www.ncbi.nlm.nih.gov/clinvar/RCV002510555.2/?redir=rcv' target='_blank'><strong>A414V</strong></a>. 
                                              In both mouse (EX3DEL) and fly models (Y121D and W465R), PRORP mutations lead to the accumulation of unprocessed mt-tRNA molecules. 
                                              This results in decreased mitochondrial protein synthesis, reduced levels of OXPHOS complexes, diminished ATP production, and increased production 
                                              of reactive oxygen species (ROS), thereby inducing significant cellular stress. These molecular dysfunctions cause developmental 
                                              delays, mitochondrial dysfunction, and cardiomyopathy.")
                 
                               
    )    
    # paste("Visualizing RNase:", input$rnase, "with colors", toString(rnase_data$colors))
    # Update the UI to display RNase details
    output$rnaseInfo <- renderUI({
      HTML(paste0('<div style="margin-top:20px;">', rnase_info_html, '</div>'))
    })
    # Sequence plot
    output$seqPlot <- renderPlot({
      msavisr(mymsa = rnase_data$file, 
              myref = "HS", 
              wnon = 1.0, 
              myroi = rnase_data$features, 
              wroi = 4.0, 
              hroi = 1, 
              hmat=0.9, 
              hnon=0.9,
              cbfcols = TRUE, 
              basecolors = c("Snow1", "Snow2", "Snow3"), 
              roicolors = rnase_data$colors)+theme(text = element_text(size = 25))
    })
    output$treePlot <- renderGirafe({
      rnase_tree$data$label[is.na(rnase_tree$data$label)] <- ""
      treePlot_interactive <- rnase_tree +  ggplot2::xlim(-0.01, max(rnase_tree$data$x)+0.05)+
        ggiraph::geom_point_interactive(aes(x = ifelse(isTip==TRUE,x+0.03,0), y = y, tooltip = label), size = 15, alpha = 0) # Invisible points
      ggiraph::girafe(ggobj = treePlot_interactive)
    })
    
    # UI for RNase image
    output$rnaseImageContainer <- renderUI({
        # Display the image
        tagList(
          img(src = img_path, alt = paste("Image of", rnase), style = "width: 100%; height: auto;")#style = "max-width: 150%; max-height: 400px;"
        )
      
      
      
    })
    
  })
  output$mainPanelSegments <- renderUI({
    if (rv$display) {
      shiny.semantic::segment(
        class = "secondary",
        style = "background-color: white; padding: 20px; border: 1px solid #ccc;",
        div(
          class = "ui grid",
          div(class = "row",
              div(class = "seven wide column",
                  shiny.semantic::segment(
                    plotOutput("mutationPlot", height = "400px"),
                    class = "secondary",
                    style = "padding: 20px; border: 1px solid #ddd; background-color: white"
                  )
              ),
              div(class = "nine wide column",
                  shiny.semantic::segment(
                    girafeOutput("treePlot", height = "400px"),
                    class = "secondary",
                    style = "padding: 20px; border: 1px solid #ddd; background-color: white"
                  )
              )
          ),
          div(class = "row",
              div(class = "sixteen wide column",
                  shiny.semantic::segment(
                    plotOutput("seqPlot", height = "600px", click = "plot_click"),
                    class = "secondary",
                    style = "padding: 20px; border: 1px solid #ddd; background-color: white"
                  )
              )
          )
        )
      )
    } else {
      # Return an empty div or an instruction message before the button is clicked
      div()
    }
  })
  
  
  
  reactiveData <- reactive({
    rnaseSeq <- switch(input$rnase, 
                       RNASEH1 = inpmsa_RNASEH1_ggmsa, 
                       RNASEH2A = inpmsa_RNASEH2A_ggmsa,
                       RNASEH2B = inpmsa_RNASEH2B_ggmsa,
                       RNASEH2C = inpmsa_RNASEH2C_ggmsa,
                       AGO2 = inpmsa_AGO2_ggmsa,
                       DICER1 = inpmsa_DICER1_ggmsa,
                       RNASET2 = inpmsa_RNASET2_ggmsa,
                       PARN = inpmsa_PARN_ggmsa,
                       ELAC2 = inpmsa_ELAC2_ggmsa,
                       PRORP = inpmsa_PRORP_ggmsa,
                       DIS3L2 = inpmsa_DIS3L2_ggmsa
    )
    original_seq <- readRDS(sprintf("./data/%s_raw_seq.rds", input$rnase))
    original_seq <- as.character(original_seq$Homo_sapiens)
    # original_seq <- switch(input$rnase, 
    #                        RNASEH1 = as.character(RNASEH1$Homo_sapiens), 
    #                        RNASEH2A = as.character(RNASEH2A$Homo_sapiens),
    #                        RNASEH2B = as.character(RNASEH2B$Homo_sapiens),
    #                        RNASEH2C = as.character(RNASEH2C$Homo_sapiens),
    #                        AGO2 = as.character(AGO2$Homo_sapiens),
    #                        DICER1 = as.character(DICER1$Homo_sapiens),
    #                        RNASET2 = as.character(RNASET2$Homo_sapiens),
    #                        PARN = as.character(PARN$Homo_sapiens),
    #                        ELAC2 = as.character(ELAC2$Homo_sapiens),
    #                        PRORP = as.character(PRORP$Homo_sapiens),
    #                        DIS3L2 = as.character(DIS3L2$Homo_sapiens)
    # )
    aligned_seq <- readAAStringSet(rnaseSeq)
    aligned_seq <- as.character(aligned_seq$HS)
    list(
      rnaseSeq = rnaseSeq,
      position_mapping = mapAminoAcidPositions(original_seq, aligned_seq)
    )
  })
  # Mutation Analysis
  
  # Reactive data to store aligned positions and corresponding amino acid
  aligned_data <- eventReactive(input$highlightMutation, {
    # Ensure the mapping data is available and input is valid
    validate(
      need(reactiveData()$position_mapping, "Position mapping not available."),
      need(is.numeric(input$mutationPos), "Position must be numeric."),
      need(input$mutationPos > 0, "Position must be greater than 0.")
    )
    
    position_mapping <- reactiveData()$position_mapping
    aligned_info <- position_mapping[input$mutationPos, ]
    
    if (is.na(aligned_info$Aligned_position)) {
      return(list(aligned_position = NA, original_position = input$mutationPos, aligned_amino = NA))
    }
    
    list(
      aligned_position = aligned_info$Aligned_position,
      original_position = aligned_info$Original_position,
      aligned_amino = aligned_info$Original_seq_amino
    )
  })
  
  # Mutation Plot
  observeEvent(input$highlightMutation, {
    data <- aligned_data()  # Using the reactive aligned_data
    
    output$mutationPlot <- renderPlot({
      req(data)  # Ensure data is not NULL
      if (is.na(data$aligned_position)) {
        plot.new()
        
        return()
      }
      
      tryCatch({
        msaData <- readAAStringSet(reactiveData()$rnaseSeq, format = "fasta")
        seqLength <- width(msaData[1])
        if (seqLength >= data$aligned_position && data$aligned_position > 0) {
          ggmsa(msaData, color = "LETTER", seq_name = TRUE, by_conservation = TRUE,
                start = max(1, data$aligned_position - 2), end = min(seqLength, data$aligned_position + 2),
                font = "TimesNewRoman", position_highlight = data$aligned_position)+theme(text = element_text(size = 25))
        }
      }, error = function(e) {
        plot.new()
        text(0.5, 0.5, "An error occurred.", cex = 1.2, col = "red")
      })
    })
  })
  
  # Aligned Position Info Text
  output$alignedPosInfo <- renderText({
    req(aligned_data())  # Ensure aligned_data is not NULL
    if (is.na(aligned_data()$aligned_position)) {
      return("Selected position is outside the sequence boundaries.")
    }
    
    paste("Amino acid", aligned_data()$aligned_amino, "at reference position", aligned_data()$original_position,
          "is aligned at position", aligned_data()$aligned_position, "and is highlighted.")
  })
  
  
  
  observeEvent(input$plot_click, {
    # Ensure the reactive data exists and is not null
    req(reactiveData()$rnaseSeq)
    # Extract the necessary data from the reactive environment
    rnaseSeq <- reactiveData()$rnaseSeq
    click_pos <- round(input$plot_click$x)  # Get the x position of the click, adjust according to your plotting logic
    selectedPosition <- click_pos 
    output$mutationPlot <- renderPlot({
      tryCatch({
        msaData <- readAAStringSet(rnaseSeq, format = "fasta")
        seqLength <- width(msaData[1])
        if (!is.null(msaData) && seqLength >= selectedPosition && selectedPosition > 0) {
          ggmsa(msaData, color = "LETTER",seq_name = TRUE, by_conservation = TRUE,
                start = max(1, selectedPosition - 2), end = min(seqLength, selectedPosition + 2),
                font = "TimesNewRoman", position_highlight = selectedPosition)+theme(text = element_text(size = 25))
        }
      }, error = function(e) {
        
      })
    })
    
    
  })
  
  
  
}

