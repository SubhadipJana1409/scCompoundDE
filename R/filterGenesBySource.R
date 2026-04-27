#' @title Filter DE Genes by Source Classification
#'
#' @description
#' A convenience utility to extract subsets of the DE table from a
#' \code{\link{CDEResult}} object by \code{source} classification
#' (\code{"transcriptional"}, \code{"compositional"}, or \code{"mixed"})
#' and optionally by FDR significance.
#'
#' @param result A \code{\link{CDEResult}} object.
#' @param source A \code{character} vector of source categories to retain.
#'   Must be one or more of \code{"transcriptional"}, \code{"compositional"},
#'   \code{"mixed"}. Default: \code{"transcriptional"}.
#' @param fdr_thresh A \code{numeric(1)} FDR threshold. Only genes with
#'   \code{adj.P.Val < fdr_thresh} are returned. Set to \code{1} to
#'   return all genes of the given source. Default: \code{0.05}.
#' @param sort_by A \code{character(1)} column to sort results by.
#'   Default: \code{"adj.P.Val"}.
#'
#' @return A \code{data.frame} of filtered DE genes.
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
#'
#' # Get only transcriptionally significant genes
#' filterGenesBySource(result, source = "transcriptional")
#'
#' # Get compositional artefacts
#' filterGenesBySource(result, source = "compositional", fdr_thresh = 1)
#'
#' @seealso \code{\link{compoundDE}}, \code{\link{deTable}}
#'
#' @importFrom methods is
#' @export
filterGenesBySource <- function(result,
                                 source     = "transcriptional",
                                 fdr_thresh = 0.05,
                                 sort_by    = "adj.P.Val") {
    if (!is(result, "CDEResult"))
        stop("'result' must be a CDEResult object.")

    valid_src <- c("transcriptional", "compositional", "mixed")
    bad <- setdiff(source, valid_src)
    if (length(bad) > 0L)
        stop("Invalid source(s): ", paste(bad, collapse = ", "),
             ". Must be one or more of: ",
             paste(valid_src, collapse = ", "))

    dt <- as.data.frame(deTable(result))

    idx <- dt[["source"]] %in% source &
           !is.na(dt[["adj.P.Val"]]) &
           dt[["adj.P.Val"]] < fdr_thresh

    out <- dt[idx, ]

    if (sort_by %in% names(out)) {
        out <- out[order(out[[sort_by]]), ]
    } else {
        warning("'sort_by' column '", sort_by, "' not found. ",
                "Returning unsorted.")
    }
    out
}
