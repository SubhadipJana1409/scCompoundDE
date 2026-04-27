# scCompoundDE 0.99.0

## New features

* Initial release submitted to Bioconductor.

* `compoundDE()`: main function that decomposes pseudo-bulk DE signals
  into transcriptional and compositional components. Automatically detects
  paired vs unpaired designs. Returns a `CDEResult` S4 object.

* `filterGenesBySource()`: convenience function to extract
  transcriptional, compositional, or mixed gene lists from a `CDEResult`,
  with optional FDR filtering.

* `plotDecomposition()`: scatter plot of z-scored transcriptional (T)
  vs compositional (C) scores for all tested genes, coloured by source
  classification.

* `plotProportion()`: stacked bar chart of mean subtype proportions per
  condition with individual sample points overlaid.

* `plotTCRatio()`: histogram of TC_ratio values with classification
  threshold lines.

* `CDEResult` S4 class with slots for `deTable`, `subtypeProportions`,
  `subtypeDE`, and `params`. Accessor methods: `deTable()`,
  `subtypeProportions()`, `subtypeDE()`, `tcRatio()`.
