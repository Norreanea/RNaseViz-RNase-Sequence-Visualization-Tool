# RNaseViz-RNase-Sequence-Visualization-Tool

RNaseViz is a web-based application built with R Shiny, providing an interactive platform for visualizing RNase sequences, their conservation across different species, and the mutations associated with Mendelian diseases.
https://norreanea.shinyapps.io/RNaseViz/
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

1. **Select RNase**: Choose an RNase from the "RNase family member" dropdown. The app updates automatically to the selected enzyme.
2. **Review Overview**: In the Overview tab, inspect the scientific context, curated variants / marked positions, ortholog tree, and organism abbreviations.
3. **Choose Reference Organism**: Use the "Reference organism" dropdown to set the sequence used for residue mapping and manual position queries.
4. **Inspect Residues**: In the Alignment tab, click on the global alignment or residue-window plot to inspect a selected aligned position and its conservation across orthologs.
5. **Highlight a Reference Position**: Enter a residue number or mutation-like query in the reference-position field and click "Highlight / inspect residue" to jump to that site.
6. **Explore Curated Features**: Use the curated variant table in Overview to navigate to marked disease-associated residues or regions.
7. **View Domain Organization**: For RNases with available domain figures, inspect the domain organization panel and its explanatory legend in the Alignment tab.
8. **Download Alignment Files**: Download the whole alignment as either a FASTA file or an annotated PDF figure with marked pathogenic variants for offline use.
![Figure2](https://github.com/user-attachments/assets/c43109ba-74c9-4e5b-98e7-6e06c615b9f2)

## Contact

For questions or support, please open an issue on this [GitHub repository](https://github.com/Norreanea/RNaseViz-RNase-Sequence-Visualization-Tool/issues).

