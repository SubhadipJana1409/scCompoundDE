library(testthat)
library(SingleCellExperiment)
library(scCompoundDE)

# ══════════════════════════════════════════════════════════════════════════════
# Helpers
# ══════════════════════════════════════════════════════════════════════════════

#' Minimal 2-subtype, 6-donor, 2-condition SCE
#' TypeA cells are DE in treatment (first 10 genes × 4).
#' TypeB cells are unchanged → DE is purely transcriptional.
make_transcriptional_sce <- function(seed = 42, n_genes = 120,
                                     n_donors = 6,
                                     cells_per_grp = 10) {
    set.seed(seed)
    n_cells <- n_donors * 2 * cells_per_grp * 2  # 2 subtypes × 2 conds × donors
    counts  <- matrix(rpois(n_genes * n_cells, 8L),
                      nrow = n_genes, ncol = n_cells)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("Cell", seq_len(n_cells))

    # Only TypeA treatment cells are upregulated
    typeA_treat <- seq(
        n_donors * cells_per_grp + 1,
        n_donors * cells_per_grp * 2
    )
    counts[seq_len(10), typeA_treat] <-
        counts[seq_len(10), typeA_treat] * 5L

    sce <- SingleCellExperiment(assays = list(counts = counts))

    # Layout: donors cycling, alternating TypeA/TypeB, first half ctrl
    half <- n_cells / 2L
    sce$donor     <- rep(rep(paste0("D", seq_len(n_donors)),
                             each = cells_per_grp), 4L)
    sce$cell_type <- "T_cell"
    sce$subtype   <- rep(rep(c("TypeA", "TypeB"),
                             each = n_donors * cells_per_grp), 2L)
    sce$condition <- rep(c("ctrl", "treat"), each = half)
    sce
}

#' SCE where DE is driven by subtype proportion shift only
#' In ctrl: 80% TypeA (low PDCD1). In treat: 80% TypeB (high PDCD1).
#' Gene1 is a marker of TypeB (high in TypeB, low in TypeA).
make_compositional_sce <- function(seed = 99, n_genes = 120) {
    set.seed(seed)
    # We create unequal proportions of subtypes per condition
    # ctrl:  many TypeA (low marker), few TypeB (high marker)
    # treat: few TypeA,              many TypeB
    # => Gene1 appears DE but only because TypeB cells became more common

    n_typeA_ctrl  <- 60L   # 6 donors x 10 cells
    n_typeB_ctrl  <- 20L   # 6 donors x ~3 cells
    n_typeA_treat <- 20L
    n_typeB_treat <- 60L

    n_cells <- n_typeA_ctrl + n_typeB_ctrl + n_typeA_treat + n_typeB_treat
    counts  <- matrix(rpois(n_genes * n_cells, 5L),
                      nrow = n_genes, ncol = n_cells)
    rownames(counts) <- paste0("Gene", seq_len(n_genes))
    colnames(counts) <- paste0("Cell", seq_len(n_cells))

    # Gene1 highly expressed in TypeB (compositional marker)
    typeB_idx <- c(
        seq(n_typeA_ctrl + 1L, n_typeA_ctrl + n_typeB_ctrl),
        seq(n_typeA_ctrl + n_typeB_ctrl + n_typeA_treat + 1L, n_cells)
    )
    counts[1L, typeB_idx] <- counts[1L, typeB_idx] * 8L

    sce <- SingleCellExperiment(assays = list(counts = counts))

    n6 <- 6L
    donor_a_ctrl  <- rep(paste0("D", seq_len(n6)), each = 10L)
    donor_b_ctrl  <- rep(paste0("D", seq_len(n6)), length.out = n_typeB_ctrl)
    donor_a_treat <- rep(paste0("D", seq_len(n6)), length.out = n_typeA_treat)
    donor_b_treat <- rep(paste0("D", seq_len(n6)), each = 10L)

    sce$donor <- c(donor_a_ctrl, donor_b_ctrl,
                   donor_a_treat, donor_b_treat)
    sce$cell_type <- "T_cell"
    sce$subtype   <- c(rep("TypeA", n_typeA_ctrl),
                       rep("TypeB", n_typeB_ctrl),
                       rep("TypeA", n_typeA_treat),
                       rep("TypeB", n_typeB_treat))
    sce$condition <- c(rep("ctrl",  n_typeA_ctrl + n_typeB_ctrl),
                       rep("treat", n_typeA_treat + n_typeB_treat))
    sce
}

# ══════════════════════════════════════════════════════════════════════════════
# CDEResult S4 class
# ══════════════════════════════════════════════════════════════════════════════

test_that("CDEResult constructor and accessors work", {
    library(S4Vectors)
    dt <- DataFrame(
        gene      = c("G1", "G2"),
        logFC     = c(1.2, 0.1),
        AveExpr   = c(3.1, 2.0),
        t         = c(4.1, 0.5),
        P.Value   = c(0.001, 0.8),
        adj.P.Val = c(0.01, 0.9),
        B         = c(2.1, -2.0),
        T_score   = c(1.1, 0.05),
        C_score   = c(0.1, 0.09),
        T_score_z = c(1.5, 0.2),
        C_score_z = c(0.2, 0.3),
        TC_ratio  = c(0.88, 0.40),
        source    = c("transcriptional", "mixed")
    )
    pm <- matrix(c(0.6, 0.4, 0.3, 0.7), nrow = 2,
                 dimnames = list(c("D1___ctrl", "D1___treat"),
                                 c("TypeA", "TypeB")))
    obj <- CDEResult(deTable = dt, subtypeProportions = pm,
                     subtypeDE = list(), params = list(
                         broad_type = "T_cell", subtypes = c("TypeA","TypeB"),
                         condition = "condition", cond_levels = c("ctrl","treat"),
                         tc_thresh_high = 0.8, tc_thresh_low = 0.2))

    expect_s4_class(obj, "CDEResult")
    expect_equal(nrow(deTable(obj)), 2L)
    expect_equal(dim(subtypeProportions(obj)), c(2L, 2L))
    expect_length(subtypeDE(obj), 0L)
})

test_that("tcRatio accessor returns named numeric", {
    library(S4Vectors)
    dt <- DataFrame(
        gene = c("G1","G2"), logFC = c(1.2,0.1), AveExpr = c(3.1,2.0),
        t = c(4.1,0.5), P.Value = c(0.001,0.8), adj.P.Val = c(0.01,0.9),
        B = c(2.1,-2.0), T_score = c(1.1,0.05), C_score = c(0.1,0.09),
        T_score_z = c(1.5,0.2), C_score_z = c(0.2,0.3),
        TC_ratio = c(0.88,0.40), source = c("transcriptional","mixed"))
    pm <- matrix(c(0.6,0.4), nrow=1,
                 dimnames=list("D1___ctrl",c("TypeA","TypeB")))
    obj <- CDEResult(dt, pm, list(), list(broad_type="T_cell",
                     subtypes=c("TypeA","TypeB"),
                     condition="condition", cond_levels=c("ctrl","treat"),
                     tc_thresh_high=0.8, tc_thresh_low=0.2))
    r <- tcRatio(obj)
    expect_type(r, "double")
    expect_named(r)
    expect_true(all(r >= 0 & r <= 1))
})

test_that("show method prints without error", {
    library(S4Vectors)
    dt <- DataFrame(
        gene="G1", logFC=1.2, AveExpr=3.1, t=4.1,
        P.Value=0.001, adj.P.Val=0.01, B=2.1,
        T_score=1.1, C_score=0.1, T_score_z=1.5, C_score_z=0.2,
        TC_ratio=0.88, source="transcriptional")
    pm <- matrix(c(0.6,0.4), nrow=1,
                 dimnames=list("D1___ctrl",c("TypeA","TypeB")))
    obj <- CDEResult(dt, pm, list(), list(broad_type="T_cell",
                     subtypes=c("TypeA","TypeB"),
                     condition="condition", cond_levels=c("ctrl","treat"),
                     tc_thresh_high=0.8, tc_thresh_low=0.2))
    expect_output(show(obj), "CDEResult")
})

# ══════════════════════════════════════════════════════════════════════════════
# compoundDE — input validation
# ══════════════════════════════════════════════════════════════════════════════

test_that("compoundDE errors on non-SCE input", {
    expect_error(
        compoundDE(list(), broad_type = "T_cell",
                   subtype_col = "subtype", broad_col = "cell_type",
                   donor = "donor", condition = "condition"),
        regexp = "SingleCellExperiment"
    )
})

test_that("compoundDE errors on missing colData column", {
    sce <- make_transcriptional_sce()
    expect_error(
        compoundDE(sce, broad_type = "T_cell",
                   subtype_col = "MISSING", broad_col = "cell_type",
                   donor = "donor", condition = "condition",
                   min_cells = 3L, min_subtypes = 2L, min_donors = 2L),
        regexp = "not found in colData"
    )
})

test_that("compoundDE errors on missing broad_type", {
    sce <- make_transcriptional_sce()
    expect_error(
        compoundDE(sce, broad_type = "NK_cell",
                   subtype_col = "subtype", broad_col = "cell_type",
                   donor = "donor", condition = "condition",
                   min_cells = 3L, min_subtypes = 2L, min_donors = 2L),
        regexp = "not found"
    )
})

test_that("compoundDE errors when condition has more than 2 levels", {
    sce <- make_transcriptional_sce()
    sce$condition <- rep(c("a","b","c","d"), length.out = ncol(sce))
    expect_error(
        compoundDE(sce, broad_type = "T_cell",
                   subtype_col = "subtype", broad_col = "cell_type",
                   donor = "donor", condition = "condition",
                   min_cells = 1L, min_subtypes = 2L, min_donors = 2L),
        regexp = "exactly 2 levels"
    )
})

test_that("compoundDE errors on bad tc_thresh values", {
    sce <- make_transcriptional_sce()
    expect_error(
        compoundDE(sce, broad_type = "T_cell",
                   subtype_col = "subtype", broad_col = "cell_type",
                   donor = "donor", condition = "condition",
                   tc_thresh_high = 0.2, tc_thresh_low = 0.8,
                   min_cells = 3L, min_subtypes = 2L, min_donors = 2L),
        regexp = "tc_thresh_high.*tc_thresh_low"
    )
})

# ══════════════════════════════════════════════════════════════════════════════
# compoundDE — core correctness
# ══════════════════════════════════════════════════════════════════════════════

test_that("compoundDE returns a CDEResult", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    expect_s4_class(res, "CDEResult")
})

test_that("deTable has required columns", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    expected <- c("gene","logFC","P.Value","adj.P.Val",
                  "T_score","C_score","T_score_z","C_score_z",
                  "TC_ratio","source")
    expect_true(all(expected %in% names(deTable(res))))
})

test_that("TC_ratio values are in [0, 1]", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    r <- tcRatio(res)
    expect_true(all(r >= 0 & r <= 1, na.rm = TRUE))
})

test_that("source classification only contains valid values", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    src <- as.data.frame(deTable(res))[["source"]]
    expect_true(all(src %in% c("transcriptional","compositional","mixed")))
})

test_that("subtypeProportions rows sum to approximately 1", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    row_sums <- rowSums(subtypeProportions(res))
    expect_true(all(abs(row_sums - 1) < 1e-9))
})

test_that("compoundDE detects injected transcriptional signal", {
    sce <- make_transcriptional_sce(seed = 42)
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    dt  <- as.data.frame(deTable(res))
    sig <- dt[dt$adj.P.Val < 0.05 & !is.na(dt$adj.P.Val), ]
    # At least some of the first 10 injected genes should be detected
    detected <- intersect(sig$gene, paste0("Gene", seq_len(10)))
    expect_gt(length(detected), 0L)
})

test_that("injected genes lean transcriptional (TC_ratio > 0.5)", {
    sce <- make_transcriptional_sce(seed = 42)
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    dt <- as.data.frame(deTable(res))
    injected <- dt[dt$gene %in% paste0("Gene", seq_len(10)), ]
    if (nrow(injected) > 0L)
        expect_gt(median(injected$TC_ratio, na.rm = TRUE), 0.5)
})

test_that("subtypeDE list has one entry per valid subtype", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    expect_type(subtypeDE(res), "list")
    expect_gte(length(subtypeDE(res)), 2L)
})

# ══════════════════════════════════════════════════════════════════════════════
# filterGenesBySource
# ══════════════════════════════════════════════════════════════════════════════

test_that("filterGenesBySource returns a data.frame", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    out <- filterGenesBySource(res, source = "transcriptional")
    expect_s3_class(out, "data.frame")
})

test_that("filterGenesBySource returns only requested source", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    out <- filterGenesBySource(res, source = "transcriptional",
                               fdr_thresh = 1)
    expect_true(all(out$source == "transcriptional"))
})

test_that("filterGenesBySource errors on invalid source", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    expect_error(
        filterGenesBySource(res, source = "garbage"),
        regexp = "Invalid source"
    )
})

test_that("filterGenesBySource errors on non-CDEResult input", {
    expect_error(filterGenesBySource(list()), regexp = "CDEResult")
})

# ══════════════════════════════════════════════════════════════════════════════
# plotDecomposition
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotDecomposition returns a ggplot", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    p <- plotDecomposition(res)
    expect_s3_class(p, "ggplot")
})

test_that("plotDecomposition errors on non-CDEResult input", {
    expect_error(plotDecomposition(list()), regexp = "CDEResult")
})

# ══════════════════════════════════════════════════════════════════════════════
# plotProportion
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotProportion returns a ggplot", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    p <- plotProportion(res)
    expect_s3_class(p, "ggplot")
})

test_that("plotProportion errors on non-CDEResult input", {
    expect_error(plotProportion(list()), regexp = "CDEResult")
})

# ══════════════════════════════════════════════════════════════════════════════
# plotTCRatio
# ══════════════════════════════════════════════════════════════════════════════

test_that("plotTCRatio returns a ggplot", {
    sce <- make_transcriptional_sce()
    res <- compoundDE(sce, broad_type = "T_cell",
                      subtype_col = "subtype", broad_col = "cell_type",
                      donor = "donor", condition = "condition",
                      min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
    p <- plotTCRatio(res)
    expect_s3_class(p, "ggplot")
})

test_that("plotTCRatio errors on non-CDEResult input", {
    expect_error(plotTCRatio(list()), regexp = "CDEResult")
})
