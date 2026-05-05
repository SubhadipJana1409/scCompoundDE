## ------------------------------------------------------------------
## Script: make_example_data.R
## Purpose: Create and save a simulated SingleCellExperiment for use
##          in inst/extdata/example_sce.rds
##
## This script documents the provenance of all data files shipped in
## inst/extdata/, as required by Bioconductor packaging guidelines.
##
## Data source:
##   Simulated count data designed to demonstrate scCompoundDE's
##   ability to decompose pseudo-bulk DE into transcriptional and
##   compositional components.
##
## Design:
##   - 6 donors (D1–D6), split evenly between ctrl and treat
##   - 1 broad cell type ("T_cell") with 2 subtypes (TypeA, TypeB)
##   - 200 genes total
##   - Genes 1–10: transcriptional signal injected (5× fold-change
##     in TypeA cells under treatment only)
##   - Unequal subtype proportions between conditions to create a
##     compositional confound (ctrl: 70% TypeA, treat: 30% TypeA)
##   - Gene11: pure compositional marker (expressed only in TypeB)
##
## Processing steps:
##   1. Generate base Poisson counts (lambda = 8)
##   2. Inject transcriptional signal in genes 1–10
##   3. Inject compositional marker in gene 11
##   4. Assign donor, cell_type, subtype, and condition metadata
##   5. Create unbalanced subtype proportions across conditions
##   6. Build SingleCellExperiment
##   7. Save as compact RDS to inst/extdata/
##
## To regenerate the data, run from the package root:
##   source("inst/scripts/make_example_data.R")
##
## Requirements:
##   BiocManager::install("SingleCellExperiment")
## ------------------------------------------------------------------

suppressPackageStartupMessages({
    library(SingleCellExperiment)
})

set.seed(2026)

## ---- Step 1: Define experimental design ----
n_genes   <- 200L
n_donors  <- 6L

## Unequal proportions: ctrl = 70% TypeA, treat = 30% TypeA
## 3 donors per condition, 40 cells per donor
n_cells_per_donor <- 40L
n_total <- n_donors * n_cells_per_donor  # 240 cells

## Ctrl donors (D1, D2, D3): 28 TypeA + 12 TypeB each
## Treat donors (D4, D5, D6): 12 TypeA + 28 TypeB each
n_a_ctrl  <- 28L;  n_b_ctrl  <- 12L
n_a_treat <- 12L;  n_b_treat <- 28L

## ---- Step 2: Generate base counts ----
message("Generating simulated count matrix...")
counts <- matrix(
    rpois(n_genes * n_total, lambda = 8L),
    nrow = n_genes, ncol = n_total
)
rownames(counts) <- paste0("Gene", seq_len(n_genes))
colnames(counts) <- paste0("Cell", seq_len(n_total))

## ---- Step 3: Assign metadata ----
## Build cell-level metadata vectors
donors   <- character(n_total)
subtypes <- character(n_total)
conds    <- character(n_total)

idx <- 1L
for (d in seq_len(3L)) {
    ## Ctrl donors: D1, D2, D3
    dn <- paste0("D", d)
    ## TypeA cells
    donors[idx:(idx + n_a_ctrl - 1L)]   <- dn
    subtypes[idx:(idx + n_a_ctrl - 1L)] <- "TypeA"
    conds[idx:(idx + n_a_ctrl - 1L)]    <- "ctrl"
    idx <- idx + n_a_ctrl
    ## TypeB cells
    donors[idx:(idx + n_b_ctrl - 1L)]   <- dn
    subtypes[idx:(idx + n_b_ctrl - 1L)] <- "TypeB"
    conds[idx:(idx + n_b_ctrl - 1L)]    <- "ctrl"
    idx <- idx + n_b_ctrl
}
for (d in seq(4L, 6L)) {
    ## Treat donors: D4, D5, D6
    dn <- paste0("D", d)
    ## TypeA cells
    donors[idx:(idx + n_a_treat - 1L)]   <- dn
    subtypes[idx:(idx + n_a_treat - 1L)] <- "TypeA"
    conds[idx:(idx + n_a_treat - 1L)]    <- "treat"
    idx <- idx + n_a_treat
    ## TypeB cells
    donors[idx:(idx + n_b_treat - 1L)]   <- dn
    subtypes[idx:(idx + n_b_treat - 1L)] <- "TypeB"
    conds[idx:(idx + n_b_treat - 1L)]    <- "treat"
    idx <- idx + n_b_treat
}

## ---- Step 4: Inject transcriptional signal ----
## Genes 1–10: 5× fold-change in TypeA treatment cells only
message("Injecting transcriptional signal in Genes 1-10...")
typeA_treat_idx <- which(subtypes == "TypeA" & conds == "treat")
counts[seq_len(10), typeA_treat_idx] <-
    counts[seq_len(10), typeA_treat_idx] * 5L

## ---- Step 5: Inject compositional marker ----
## Gene11: expressed only in TypeB cells (pure compositional signal)
message("Injecting compositional marker in Gene11...")
typeB_idx <- which(subtypes == "TypeB")
typeA_idx <- which(subtypes == "TypeA")
counts[11L, typeB_idx] <- counts[11L, typeB_idx] * 8L
counts[11L, typeA_idx] <- rpois(length(typeA_idx), lambda = 1L)

## ---- Step 6: Build SingleCellExperiment ----
message("Building SingleCellExperiment...")
sce <- SingleCellExperiment(assays = list(counts = counts))
sce$donor     <- donors
sce$cell_type <- "T_cell"
sce$subtype   <- subtypes
sce$condition <- conds

## ---- Step 7: Save ----
outdir <- file.path("inst", "extdata")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

outfile <- file.path(outdir, "example_sce.rds")
saveRDS(sce, file = outfile, compress = "xz")

fsize <- file.size(outfile)
message(
    "\n=== Data saved ===\n",
    "File:       ", outfile, "\n",
    "Size:       ", round(fsize / 1024, 1), " KB\n",
    "Dimensions: ", nrow(sce), " genes x ", ncol(sce), " cells\n",
    "Donors:     ", paste(unique(sce$donor), collapse = ", "), "\n",
    "Conditions: ", paste(unique(sce$condition), collapse = ", "), "\n",
    "Subtypes:   ", paste(unique(sce$subtype), collapse = ", "), "\n",
    "Cell type:  ", unique(sce$cell_type), "\n",
    "Source:     Simulated (Poisson counts with injected signals)\n",
    "Design:     Unbalanced subtypes + transcriptional injection\n"
)

## ---- Step 8: Verify data ----
message("\nVerifying data structure...")
message("Subtype x condition table:")
print(table(sce$subtype, sce$condition))

message("\nDonor x condition table:")
print(table(sce$donor, sce$condition))

message("\nMean expression of injected genes (treat TypeA vs ctrl TypeA):")
typeA_ctrl_idx  <- which(sce$subtype == "TypeA" & sce$condition == "ctrl")
typeA_treat_idx <- which(sce$subtype == "TypeA" & sce$condition == "treat")
for (g in paste0("Gene", seq_len(10))) {
    ctrl_mean  <- mean(counts(sce)[g, typeA_ctrl_idx])
    treat_mean <- mean(counts(sce)[g, typeA_treat_idx])
    message(sprintf("  %s: ctrl=%.1f  treat=%.1f  FC=%.1f",
                    g, ctrl_mean, treat_mean, treat_mean / ctrl_mean))
}

message("\nGene11 (compositional marker) mean by subtype:")
message(sprintf("  TypeA: %.1f", mean(counts(sce)["Gene11", typeA_idx])))
message(sprintf("  TypeB: %.1f", mean(counts(sce)["Gene11", typeB_idx])))

message("\n=== Example data generation complete ===")
