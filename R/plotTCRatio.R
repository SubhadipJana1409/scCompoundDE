#' @title Histogram of TC_ratio Distribution
#'
#' @description
#' Plots the distribution of \code{TC_ratio} values across all tested
#' genes. The TC_ratio measures the fraction of each gene's DE signal
#' that is attributable to transcriptional (cell-intrinsic) versus
#' compositional (subtype proportion shift) change. Vertical dashed
#' lines show the classification thresholds.
#'
#' @param result A \code{\link{CDEResult}} object from
#'   \code{\link{compoundDE}}.
#' @param fdr_thresh A \code{numeric(1)} FDR threshold. Genes with
#'   \code{adj.P.Val >= fdr_thresh} are shown in grey.
#'   Default: \code{0.05}.
#' @param bins An \code{integer(1)} number of histogram bins.
#'   Default: \code{40L}.
#' @param colours A named \code{character} vector with colours for
#'   \code{"transcriptional"}, \code{"compositional"}, \code{"mixed"},
#'   and \code{"ns"}. Default uses a colour-blind friendly palette.
#'
#' @return A \code{ggplot2} object.
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(42)
#' n_genes <- 150
#' counts <- matrix(rpois(n_genes * 120, 8L), nrow = n_genes, ncol = 120)
#' rownames(counts) <- paste0("Gene", seq_len(n_genes))
#' colnames(counts) <- paste0("Cell", seq_len(120))
#' counts[seq_len(10), 61:90] <- counts[seq_len(10), 61:90] * 4L
#'
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$donor     <- rep(paste0("D", seq_len(6)), each = 20)
#' sce$cell_type <- "T_cell"
#' sce$subtype   <- rep(c("TypeA", "TypeB"), times = 60)
#' sce$condition <- rep(c("ctrl", "treat"), each = 60)
#'
#' result <- compoundDE(sce, broad_type = "T_cell",
#'     subtype_col = "subtype", broad_col = "cell_type",
#'     donor = "donor", condition = "condition",
#'     min_cells = 3L, min_subtypes = 2L, min_donors = 2L)
#' plotTCRatio(result)
#'
#' @seealso \code{\link{compoundDE}}, \code{\link{plotDecomposition}},
#'   \code{\link{plotProportion}}
#'
#' @importFrom ggplot2 ggplot aes geom_histogram geom_vline scale_fill_manual labs theme_bw theme element_text element_blank annotate
#' @importFrom methods is
#' @export
plotTCRatio <- function(result,
                         fdr_thresh = 0.05,
                         bins       = 40L,
                         colours    = c(
                             transcriptional = "#2196F3",
                             compositional   = "#F44336",
                             mixed           = "#FF9800",
                             ns              = "#BDBDBD"
                         )) {
    if (!is(result, "CDEResult"))
        stop("'result' must be a CDEResult object.")

    dt <- as.data.frame(deTable(result))

    dt$plot_source <- "ns"
    sig_idx <- !is.na(dt$adj.P.Val) & dt$adj.P.Val < fdr_thresh
    dt$plot_source[sig_idx] <- dt$source[sig_idx]
    dt$plot_source <- factor(dt$plot_source,
                             levels = c("transcriptional", "compositional",
                                        "mixed", "ns"))

    tc_hi <- result@params$tc_thresh_high
    tc_lo <- result@params$tc_thresh_low

    n_t <- sum(dt$plot_source == "transcriptional", na.rm = TRUE)
    n_c <- sum(dt$plot_source == "compositional",   na.rm = TRUE)
    n_m <- sum(dt$plot_source == "mixed",           na.rm = TRUE)

    ggplot(dt, aes(x    = .data[["TC_ratio"]],
                   fill = .data[["plot_source"]])) +
        geom_histogram(bins = bins, colour = "white", linewidth = 0.2) +
        geom_vline(xintercept = tc_lo, linetype = "dashed",
                   colour = "#F44336", linewidth = 0.7) +
        geom_vline(xintercept = tc_hi, linetype = "dashed",
                   colour = "#2196F3", linewidth = 0.7) +
        scale_fill_manual(
            values = colours,
            labels = c(
                transcriptional = sprintf("Transcriptional (%d)", n_t),
                compositional   = sprintf("Compositional (%d)",   n_c),
                mixed           = sprintf("Mixed (%d)",           n_m),
                ns              = sprintf("Not significant (%d)",
                                          sum(dt$plot_source == "ns"))
            )
        ) +
        annotate("text", x = tc_lo - 0.03, y = Inf, vjust = 1.5,
                 hjust = 1, size = 3, colour = "#F44336",
                 label = paste0("Compositional\n(<", tc_lo, ")")) +
        annotate("text", x = tc_hi + 0.03, y = Inf, vjust = 1.5,
                 hjust = 0, size = 3, colour = "#2196F3",
                 label = paste0("Transcriptional\n(>", tc_hi, ")")) +
        labs(
            x     = "TC_ratio  (0 = fully compositional, 1 = fully transcriptional)",
            y     = "Number of genes",
            fill  = "Gene source",
            title = paste0("TC_ratio distribution: ", result@params$broad_type),
            subtitle = paste0("FDR < ", fdr_thresh,
                              "  |  Conditions: ",
                              paste(result@params$cond_levels,
                                    collapse = " vs "))
        ) +
        theme_bw(base_size = 11) +
        theme(
            panel.grid.minor = element_blank(),
            axis.text        = element_text(size = 10),
            plot.title       = element_text(size = 12, face = "bold"),
            plot.subtitle    = element_text(size = 10, colour = "grey40")
        )
}
