# RNaseViz refactor package

This package replaces the original  `server.R` / `ui.R` structure with a new Shiny app layout.

## What changed

- Removed runtime package installation from `server.R`.
- Split configuration, data loading, position mapping, plotting, UI, and server logic into separate scripts.
- Added cached loading for RDS/FASTA/tree files.
- Added a scientific UX layer:
  - RNase context cards.
  - searchable mutation/feature table.
  - data provenance table.
  - residue-level conservation report.
  - alignment-window plot around selected residue.
  - downloadable residue report.
