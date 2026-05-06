# RNaseViz refactor package

This package replaces the original monolithic `server.R` / `ui.R` structure with a modular Shiny app layout.

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

## How to apply

Copy the files/folders in this package into the root of your existing repository:

```text
ui.R
server.R
app.R
install_dependencies.R
R/
```

Keep your existing `data/` and `www/images/` folders unchanged.

## Install once, then deploy

```r
source("install_dependencies.R")
shiny::runApp()
```

Do not call `BiocManager::install()` from `server.R`; this was the likely reason cold starts were so slow.

## Notes

The new code expects the existing file naming convention:

```text
data/reordered_for_msavisr_<RNASE>Aln.fa
data/reordered_for_ggmsa_<RNASE>Aln.fa
data/<RNASE>_raw_seq.rds
data/<RNASE>_tree.rds
www/images/<RNASE>.png
```
