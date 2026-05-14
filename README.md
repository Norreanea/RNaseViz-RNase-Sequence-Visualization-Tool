# RNaseViz-RNase-Sequence-Visualization-Tool

RNaseViz is a web-based application built with R Shiny, providing an interactive platform for visualizing RNase sequences, their conservation across different species, and the mutations associated with Mendelian diseases.
https://norreanea.shinyapps.io/RNaseViz/

This tool was developed as part of the article [Ribonucleases in Mendelian disease: Characterization and insight from model organisms](https://www.sciencedirect.com/science/article/pii/S2352304225001023)

## Features

- Selection of RNase enzymes from a curated list of disease-relevant ribonucleases.
- Visualization of sequence conservation across orthologs using multiple sequence alignments.
- Overview cards summarizing the selected RNase, curated variants, available orthologs, scientific context, and organism abbreviations.
- Display of ortholog trees and domain organization figures for supported RNases.
- Residue-level inspection by clicking the alignment or entering a reference-sequence position manually.
- Download of whole alignments as FASTA files or annotated PDF figures with marked pathogenic variants.
- Support for multiple RNase enzymes associated with Mendelian disease.

## Usage

To use the RNaseViz application:

1. **RNase selector**: Choose the RNase family member to analyze
2. **Reference organism selector**: Select the species whose sequence will be used as the positional reference for residue mapping.
3. **Residue / Variant query box**: Enter a residue number or mutation label to jump directly to a selected site in the alignment.
4. **Alignment download panel**: Choose the export format and download either the FASTA alignment or an annotated PDF with marked pathogenic variants.
5. **Pathway overview figure**: A figure linking the selected RNase to its biological pathway and disease context.
6. **Main navigation tabs**: Switch between the two views: Overview and Alignment.
7. Summary cards**: The number of curated variants/regions and available orthologs.
8. **Scientific context panel**: Brief biological summary of the selected RNase, including function, disease relevance, and highlighted pathogenic variants.
9. **Ortholog phylogenetic tree**: Shows evolutionary relationships among orthologs of the selected RNase across species.
10. **Organism abbreviations**
11. **Curated variant table**: Interactive table listing pathogenic or otherwise marked residues/regions, searchable and filterable by organism and aligned position.

<img width="927" height="710" alt="image" src="https://github.com/user-attachments/assets/d7e7c39a-1806-4356-b8ac-b61c576104b0" />

12. **Residue inspection summary**: Displays the currently selected residue, its aligned coordinate, conservation level, and nearby curated features.
13. **Domain organization panel**: Shows protein architecture across species, with domains and marked disease-associated residues positioned along the sequence.
14. **Domain figure legend**: Explains how to read the domain map, including aligned positions, exon deletions, domain colors, and mutation markers.
15. **Multiple sequence alignment viewer**: Comparative visualization showing residue conservation, mismatches, gaps, and the positions of curated variants across orthologs.
    
<img width="851" height="912" alt="image" src="https://github.com/user-attachments/assets/4561ded4-693b-473c-b5fb-81bc283b94ad" />

Bellow is the example of DICER1 alignment report from RNaseViz. The human sequence is aligned with orthologs from other species, and disease-associated mutations are marked directly on the alignment. The colors correspond to different diseases, and each mutation label is linked to corresponding ClinVar page.

<img width="1020" height="719" alt="image" src="https://github.com/user-attachments/assets/aa217ebc-a73b-4304-9061-70dd46a4843f" />

## Contact

For questions or support, please open an issue on this [GitHub repository](https://github.com/Norreanea/RNaseViz-RNase-Sequence-Visualization-Tool/issues).

