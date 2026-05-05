#' @title Stacked Bar Plot of Subtype Proportions per Condition
#'
#' @description
#' Produces a stacked bar chart showing the mean proportion of each
#' subtype within the broad cell type, split by condition. Samples are
#' shown as individual points overlaid on the bars, making it easy to
#' see between-donor variability in composition. A significant
#' compositional shift across conditions is a key driver of
#' \emph{compositional} DE genes identified by \code{\link{compoundDE}}.
#'
#' @param result A \code{\link{CDEResult}} object from
#'   \code{\link{compoundDE}}.
#' @param show_points Logical. If \code{TRUE}, overlay individual sample
#'   proportions as jittered points. Default: \code{TRUE}.
#' @param colours A named or unnamed character vector of colours for
#'   the subtypes. If \code{NULL}, a default palette is used.
#'   Default: \code{NULL}.
#' @param point_size A \code{numeric(1)} point size. Default: \code{2}.
#' @param point_alpha A \code{numeric(1)} point alpha. Default:
#'   \code{0.8}.
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
#' plotProportion(result)
#'
#' @seealso \code{\link{compoundDE}}, \code{\link{plotDecomposition}},
#'   \code{\link{plotTCRatio}}
#'
#' @importFrom ggplot2 ggplot aes geom_bar geom_jitter position_jitter
#'   position_stack position_dodge scale_fill_manual scale_color_manual
#'   labs theme_bw theme element_text element_blank facet_wrap
#' @importFrom methods is
#' @importFrom stats reshape aggregate setNames
#' @export
plotProportion <- function(result,
                            show_points = TRUE,
                            colours     = NULL,
                            point_size  = 2,
                            point_alpha = 0.8) {
    if (!is(result, "CDEResult"))
        stop("'result' must be a CDEResult object.")

    prop_mat  <- subtypeProportions(result)
    subtypes  <- result@params$subtypes
    cond_lvls <- result@params$cond_levels

    # Extract condition from row names (donor___condition)
    parts    <- strsplit(rownames(prop_mat), "___", fixed = TRUE)
    cond_vec <- vapply(parts, `[`, character(1L), 2L)
    donor_vec <- vapply(parts, `[`, character(1L), 1L)

    # Keep only valid subtypes
    valid_sub <- intersect(subtypes, colnames(prop_mat))

    # Long-format data frame for plotting
    prop_long <- do.call(rbind, lapply(valid_sub, function(k) {
        data.frame(
            sample    = rownames(prop_mat),
            donor     = donor_vec,
            condition = cond_vec,
            subtype   = k,
            proportion = prop_mat[, k],
            stringsAsFactors = FALSE
        )
    }))
    prop_long$condition <- factor(prop_long$condition, levels = cond_lvls)
    prop_long$subtype   <- factor(prop_long$subtype, levels = valid_sub)

    # Mean proportion per condition x subtype
    mean_df <- aggregate(proportion ~ condition + subtype,
                         data = prop_long, FUN = mean)
    mean_df$condition <- factor(mean_df$condition, levels = cond_lvls)
    mean_df$subtype   <- factor(mean_df$subtype, levels = valid_sub)

    # Colour palette
    n_sub <- length(valid_sub)
    default_colours <- c("#2196F3","#F44336","#4CAF50","#FF9800",
                         "#9C27B0","#00BCD4","#FF5722","#8BC34A")
    if (is.null(colours)) {
        colours <- setNames(
            default_colours[seq_len(n_sub)],
            valid_sub
        )
    }

    p <- ggplot(mean_df,
                aes(x    = .data[["condition"]],
                    y    = .data[["proportion"]],
                    fill = .data[["subtype"]])) +
        geom_bar(stat = "identity", position = "stack",
                 colour = "white", linewidth = 0.3) +
        scale_fill_manual(values = colours) +
        labs(
            x     = "Condition",
            y     = "Mean proportion",
            fill  = "Subtype",
            title = paste0("Subtype composition: ", result@params$broad_type),
            subtitle = paste0("Compositional shift drives ",
                              sum(deTable(result)[["source"]] == "compositional",
                                  na.rm = TRUE),
                              " genes")
        ) +
        theme_bw(base_size = 11) +
        theme(
            panel.grid.minor = element_blank(),
            axis.text        = element_text(size = 10),
            plot.title       = element_text(size = 12, face = "bold"),
            plot.subtitle    = element_text(size = 10, colour = "grey40")
        )

    if (show_points) {
        p <- p + geom_jitter(
            data    = prop_long,
            aes(x      = .data[["condition"]],
                y      = .data[["proportion"]],
                colour = .data[["subtype"]]),
            position    = position_jitter(width = 0.08),
            size        = point_size,
            alpha       = point_alpha,
            inherit.aes = FALSE
        ) +
        scale_color_manual(values = colours, guide = "none")
    }
    p
}
