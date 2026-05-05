#' @title Scatter Plot of Transcriptional vs Compositional DE Scores
#'
#' @description
#' Produces a scatter plot of the z-scored transcriptional score
#' (\code{T_score_z}) versus the compositional score (\code{C_score_z})
#' for every gene, coloured by their \code{source} classification. This
#' visualisation makes it immediately clear which genes are driven by
#' cell-intrinsic expression changes versus subtype proportion shifts.
#'
#' @param result A \code{\link{CDEResult}} object from
#'   \code{\link{compoundDE}}.
#' @param fdr_thresh A \code{numeric(1)} FDR threshold. Only genes with
#'   \code{adj.P.Val < fdr_thresh} are highlighted as significant.
#'   Default: \code{0.05}.
#' @param top_n An \code{integer(1)} number of top significant genes to
#'   label by name. Default: \code{10L}.
#' @param point_size A \code{numeric(1)} point size. Default: \code{1.2}.
#' @param point_alpha A \code{numeric(1)} point transparency. Default:
#'   \code{0.7}.
#' @param colours A named \code{character} vector with colours for
#'   \code{"transcriptional"}, \code{"compositional"}, \code{"mixed"},
#'   and \code{"ns"} (not significant). Default uses a colour-blind
#'   friendly palette.
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
#' plotDecomposition(result)
#'
#' @seealso \code{\link{compoundDE}}, \code{\link{plotTCRatio}},
#'   \code{\link{plotProportion}}
#'
#' @importFrom ggplot2 ggplot aes geom_point geom_vline geom_hline scale_color_manual labs theme_bw theme element_text element_blank geom_text annotate
#' @importFrom methods is
#' @importFrom utils head
#' @export
plotDecomposition <- function(result,
                               fdr_thresh  = 0.05,
                               top_n       = 10L,
                               point_size  = 1.2,
                               point_alpha = 0.7,
                               colours     = c(
                                   transcriptional = "#2196F3",
                                   compositional   = "#F44336",
                                   mixed           = "#FF9800",
                                   ns              = "#BDBDBD"
                               )) {
    if (!is(result, "CDEResult"))
        stop("'result' must be a CDEResult object.")

    dt <- as.data.frame(deTable(result))

    # Label only significant genes
    dt$plot_source <- "ns"
    sig_idx <- !is.na(dt$adj.P.Val) & dt$adj.P.Val < fdr_thresh
    dt$plot_source[sig_idx] <- dt$source[sig_idx]

    # Top genes to label (significant, sorted by adj.P.Val)
    sig_genes    <- dt[sig_idx, ]
    sig_genes    <- sig_genes[order(sig_genes$adj.P.Val), ]
    label_genes  <- head(sig_genes, top_n)

    p <- ggplot(dt, aes(x = .data[["T_score_z"]],
                        y = .data[["C_score_z"]],
                        colour = .data[["plot_source"]])) +
        geom_point(size = point_size, alpha = point_alpha) +
        geom_vline(xintercept = 0, linetype = "dashed",
                   colour = "grey50", linewidth = 0.4) +
        geom_hline(yintercept = 0, linetype = "dashed",
                   colour = "grey50", linewidth = 0.4) +
        scale_color_manual(
            values = colours,
            labels = c(
                transcriptional = sprintf("Transcriptional (%d)",
                    sum(dt$plot_source == "transcriptional")),
                compositional   = sprintf("Compositional (%d)",
                    sum(dt$plot_source == "compositional")),
                mixed           = sprintf("Mixed (%d)",
                    sum(dt$plot_source == "mixed")),
                ns              = sprintf("Not significant (%d)",
                    sum(dt$plot_source == "ns"))
            )
        ) +
        labs(
            x      = "Transcriptional score (z-scored T)",
            y      = "Compositional score (z-scored C)",
            colour = "Gene source",
            title  = paste0("DE decomposition: ", result@params$broad_type),
            subtitle = paste0("FDR < ", fdr_thresh,
                              " | Conditions: ",
                              paste(result@params$cond_levels, collapse=" vs "))
        ) +
        theme_bw(base_size = 11) +
        theme(
            panel.grid.minor = element_blank(),
            axis.text        = element_text(size = 10),
            plot.title       = element_text(size = 12, face = "bold"),
            plot.subtitle    = element_text(size = 10, colour = "grey40")
        )

    if (nrow(label_genes) > 0L) {
        p <- p + geom_text(
            data        = label_genes,
            aes(label   = .data[["gene"]]),
            size        = 3,
            vjust       = -0.6,
            show.legend = FALSE
        )
    }
    p
}
