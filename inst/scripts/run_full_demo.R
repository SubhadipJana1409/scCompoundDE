## ------------------------------------------------------------------
## Script: run_full_demo.R
## Purpose: Demonstrate scCompoundDE on the Kang et al. (2018) PBMC
##          IFN-beta stimulation dataset. Shows how pseudo-bulk DE
##          signals can be decomposed into transcriptional (real
##          biology) versus compositional (proportion-driven) sources.
##
## Dataset:
##   Kang HM et al. (2018). Multiplexed droplet single-cell
##   RNA-sequencing using natural genetic variation.
##   Nature Biotechnology, 36:89-94. doi:10.1038/nbt.4042
##
##   8 lupus patients, each profiled under control and IFN-beta
##   stimulation. ~29,000 PBMCs with annotated cell types.
##
## Usage:
##   source("inst/scripts/run_full_demo.R")
##
## Requirements:
##   BiocManager::install(c("muscData", "scCompoundDE"))
## ------------------------------------------------------------------

library(scCompoundDE)
library(SingleCellExperiment)

## ---- Load Kang 2018 data ------------------------------------------------

cat("Loading Kang et al. 2018 IFN-beta PBMC dataset...\n")

if (!requireNamespace("muscData", quietly = TRUE))
    stop("Install muscData first: BiocManager::install('muscData')")

sce <- muscData::Kang18_8vs8()

## Map column names to scCompoundDE interface
sce$donor     <- sce$ind
sce$condition <- sce$stim
sce$subtype   <- sce$cell
sce$broad     <- "PBMC"
sce <- sce[, !is.na(sce$subtype)]

cat(sprintf("  %d genes, %d cells, %d donors\n",
            nrow(sce), ncol(sce), length(unique(sce$donor))))
cat(sprintf("  Subtypes: %s\n",
            paste(sort(unique(sce$subtype)), collapse = ", ")))
cat(sprintf("  Conditions: %s\n\n",
            paste(unique(sce$condition), collapse = " vs ")))

## ---- Run decomposition --------------------------------------------------

cat("Running compoundDE() on all PBMCs...\n")

t0 <- Sys.time()
res <- compoundDE(
    sce,
    broad_type   = "PBMC",
    subtype_col  = "subtype",
    broad_col    = "broad",
    donor        = "donor",
    condition    = "condition",
    min_cells    = 10L,
    min_subtypes = 2L,
    min_donors   = 2L
)
elapsed <- as.numeric(Sys.time() - t0, units = "secs")

cat(sprintf("  Done in %.1f seconds.\n\n", elapsed))
show(res)
cat("\n")

## ---- Results overview ---------------------------------------------------

dt <- as.data.frame(deTable(res))

cat(sprintf("Tested %d genes. Significant at FDR < 0.05: %d\n",
            nrow(dt),
            sum(dt$adj.P.Val < 0.05, na.rm = TRUE)))

## Breakdown by source
src <- table(dt$source[dt$adj.P.Val < 0.05])
for (s in names(src))
    cat(sprintf("  %-16s %d\n", s, src[s]))
cat("\n")

## ---- Validate against known IFN-beta genes ------------------------------
## These canonical interferon-stimulated genes should be detected as
## transcriptional — they represent real cell-intrinsic expression
## changes, not proportion shifts.

ifn_genes <- c("ISG15", "IFIT1", "IFIT3", "IFI6", "MX1",
               "OAS1", "STAT1", "IRF7", "IFI44L", "ISG20")

cat("Validation: known IFN-beta response genes\n")
cat(sprintf("%-10s %8s %10s %9s %s\n",
            "Gene", "logFC", "FDR", "TC_ratio", "Source"))
cat(strrep("-", 52), "\n")

for (g in ifn_genes) {
    if (g %in% dt$gene) {
        r <- dt[dt$gene == g, ]
        cat(sprintf("%-10s %+8.2f %10.2e %9.3f %s\n",
                    g, r$logFC, r$adj.P.Val, r$TC_ratio, r$source))
    } else {
        cat(sprintf("%-10s %s\n", g, "(not found)"))
    }
}
cat("\n")

## ---- Top compositional genes (potential artefacts) ----------------------

comp <- filterGenesBySource(res, source = "compositional", fdr_thresh = 0.05)
if (nrow(comp) > 0) {
    cat(sprintf("Top compositional genes (%d total at FDR < 0.05):\n",
                nrow(comp)))
    comp <- comp[order(comp$adj.P.Val), ]
    print(head(comp[, c("gene", "logFC", "adj.P.Val",
                         "TC_ratio", "source")], 10))
    cat("\n")
}

## ---- Per-subtype DE summary ---------------------------------------------

sde <- subtypeDE(res)
cat("Per-subtype DE models:\n")
cat(sprintf("%-25s %6s %6s\n", "Subtype", "Tested", "Sig"))
cat(strrep("-", 40), "\n")
for (nm in names(sde)) {
    n_sig <- sum(sde[[nm]]$adj.P.Val < 0.05, na.rm = TRUE)
    cat(sprintf("%-25s %6d %6d\n", nm, nrow(sde[[nm]]), n_sig))
}
cat("\n")

## ---- Subtype proportions ------------------------------------------------

cat("Mean subtype proportions per condition:\n")
for (cc in sort(unique(sce$condition))) {
    cat(sprintf("  %s:\n", cc))
    sub_cells <- sce[, sce$condition == cc]
    tbl <- table(sub_cells$subtype)
    pct <- tbl / sum(tbl) * 100
    for (i in seq_along(pct))
        cat(sprintf("    %-25s %.1f%%\n", names(pct)[i], pct[i]))
}
cat("\n")

## ---- Plots --------------------------------------------------------------

cat("Saving plots...\n")

figdir <- file.path("man", "figures")
if (!dir.exists(figdir)) dir.create(figdir, recursive = TRUE)

## Decomposition scatter
p1 <- plotDecomposition(res, fdr_thresh = 0.05, top_n = 15)
f1 <- file.path(figdir, "scCompoundDE_decomposition.png")
ggplot2::ggsave(f1, p1, width = 9, height = 7, dpi = 150)

## TC ratio histogram
p2 <- plotTCRatio(res, fdr_thresh = 0.05)
f2 <- file.path(figdir, "scCompoundDE_tc_ratio.png")
ggplot2::ggsave(f2, p2, width = 7, height = 5, dpi = 150)

## Proportion bar chart
set.seed(42)
p3 <- plotProportion(res)
f3 <- file.path(figdir, "scCompoundDE_proportions.png")
ggplot2::ggsave(f3, p3, width = 8, height = 6, dpi = 150)

cat(sprintf("  %s\n  %s\n  %s\n\n", f1, f2, f3))

## ---- Final verdict ------------------------------------------------------

n_ifn_trans <- sum(ifn_genes %in%
    dt$gene[dt$adj.P.Val < 0.05 & dt$source == "transcriptional"])

cat(strrep("=", 50), "\n")
cat(sprintf("IFN-beta genes classified as transcriptional: %d / %d\n",
            n_ifn_trans, length(ifn_genes)))

if (n_ifn_trans >= 5) {
    cat("PASS - decomposition correctly identifies IFN-beta\n")
    cat("       response as transcriptionally driven.\n")
} else {
    cat("CHECK - fewer than expected. Review parameters.\n")
}
cat(strrep("=", 50), "\n")
