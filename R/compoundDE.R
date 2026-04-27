#' @title Compound Differential Expression: Transcriptional and Compositional
#'   Decomposition
#'
#' @description
#' \code{compoundDE} decomposes pseudo-bulk differential expression (DE)
#' signals into two orthogonal components:
#'
#' \itemize{
#'   \item \strong{Transcriptional (T):} Cell-intrinsic expression changes —
#'     the gene would be DE even if subtype proportions were held fixed.
#'   \item \strong{Compositional (C):} Signal arising from a shift in the
#'     relative abundance of subtypes — the gene appears DE only because
#'     high-expressing (or low-expressing) subtypes became more or less common.
#' }
#'
#' Standard pseudo-bulk DE tools (DESeq2, edgeR, limma-voom) confound these
#' two sources. \code{compoundDE} fits a separate limma-voom model for each
#' cell subtype within the broad population, estimates subtype proportion
#' changes across conditions, and uses a z-score-normalised decomposition to
#' assign each gene a \strong{TC_ratio} — the fraction of its DE signal
#' attributable to transcriptional change.
#'
#' @details
#' The decomposition algorithm:
#'
#' \enumerate{
#'   \item Filter sparse donor-subtype combinations (\code{min_cells}).
#'   \item Compute subtype proportion matrix \eqn{\pi_{d,k,c}} (donor ×
#'     subtype × condition).
#'   \item For each subtype \eqn{k}, aggregate pseudo-bulk and run
#'     \code{limma::voom} to obtain per-subtype log-fold-changes
#'     \eqn{\text{logFC}_{g,k}}.
#'   \item Run a broad pseudo-bulk DE (all subtypes collapsed) to obtain the
#'     observed \eqn{\text{logFC}_g}, \eqn{P}-values, and FDR.
#'   \item Compute the \strong{transcriptional component}:
#'     \deqn{T_g = \sum_k \bar\pi_k \cdot \text{logFC}_{g,k}}
#'     where \eqn{\bar\pi_k} is the mean proportion of subtype \eqn{k}
#'     across all samples.
#'   \item Compute the \strong{compositional component}:
#'     \deqn{C_g = \sum_k \Delta\pi_k \cdot \bar\mu_{g,k}}
#'     where \eqn{\Delta\pi_k = \bar\pi_{k,\text{treat}} -
#'     \bar\pi_{k,\text{ctrl}}} and \eqn{\bar\mu_{g,k}} is the mean
#'     log2 CPM of gene \eqn{g} in subtype \eqn{k}.
#'   \item Z-score normalise both \eqn{T_g} and \eqn{C_g} across genes,
#'     then compute:
#'     \deqn{\text{TC\_ratio}_g =
#'       \frac{|T_g^z|}{|T_g^z| + |C_g^z| + \varepsilon}}
#'   \item Classify each gene as \emph{transcriptional}
#'     (\eqn{\geq} \code{tc_thresh_high}), \emph{compositional}
#'     (\eqn{\leq} \code{tc_thresh_low}), or \emph{mixed}.
#' }
#'
#' @param sce A \code{\link[SingleCellExperiment]{SingleCellExperiment}}
#'   with raw counts in \code{assay(sce, "counts")}.
#' @param broad_type A \code{character(1)} value present in
#'   \code{colData(sce)[[broad_col]]} specifying the broad population to
#'   decompose (e.g. \code{"T_cell"}).
#' @param subtype_col A \code{character(1)} naming the column in
#'   \code{colData(sce)} with fine-grained subtype labels
#'   (e.g. \code{"cell_subtype"}).
#' @param broad_col A \code{character(1)} naming the column with broad cell
#'   type labels (e.g. \code{"cell_type"}).
#' @param donor A \code{character(1)} naming the donor column.
#' @param condition A \code{character(1)} naming the condition column.
#'   Must have exactly two unique values.
#' @param contrast A \code{character(1)} limma contrast string,
#'   e.g. \code{"treat - ctrl"}. If \code{NULL}, uses level 2 minus level 1.
#'   Default: \code{NULL}.
#' @param min_cells An \code{integer(1)} minimum cells per donor-condition
#'   per subtype. Samples below this are removed. Default: \code{10L}.
#' @param min_subtypes An \code{integer(1)} minimum number of subtypes
#'   required for decomposition. Default: \code{2L}.
#' @param min_cpm A \code{numeric(1)} minimum CPM for gene expression filter.
#'   Default: \code{1}.
#' @param min_donors An \code{integer(1)} minimum number of pseudo-bulk
#'   samples a gene must be expressed in. Default: \code{2L}.
#' @param tc_thresh_high A \code{numeric(1)} TC_ratio threshold above which
#'   a gene is classified as \emph{transcriptional}. Default: \code{0.8}.
#' @param tc_thresh_low A \code{numeric(1)} TC_ratio threshold below which
#'   a gene is classified as \emph{compositional}. Default: \code{0.2}.
#' @param assay_name A \code{character(1)} name of the count assay.
#'   Default: \code{"counts"}.
#' @param BPPARAM A \code{\link[BiocParallel]{BiocParallelParam}} object.
#'   Default: \code{SerialParam()}.
#'
#' @return A \code{\link{CDEResult}} object containing:
#'   \itemize{
#'     \item \code{deTable}: per-gene statistics — \code{logFC},
#'       \code{P.Value}, \code{adj.P.Val} (from the broad model) plus
#'       \code{T_score}, \code{C_score}, \code{T_score_z}, \code{C_score_z},
#'       \code{TC_ratio}, and \code{source}.
#'     \item \code{subtypeProportions}: the \eqn{\pi} matrix.
#'     \item \code{subtypeDE}: per-subtype limma results.
#'     \item \code{params}: the analysis parameters.
#'   }
#'
#' @examples
#' library(SingleCellExperiment)
#'
#' set.seed(42)
#' n_genes <- 150
#' n_cells <- 120
#'
#' counts <- matrix(rpois(n_genes * n_cells, 8L),
#'                  nrow = n_genes, ncol = n_cells)
#' rownames(counts) <- paste0("Gene", seq_len(n_genes))
#' colnames(counts) <- paste0("Cell", seq_len(n_cells))
#'
#' # Inject DE signal into first 10 genes for subtype A treatment
#' counts[seq_len(10), 61:90] <- counts[seq_len(10), 61:90] * 4L
#'
#' sce <- SingleCellExperiment(assays = list(counts = counts))
#' sce$donor     <- rep(paste0("D", seq_len(6)), each = 20)
#' sce$cell_type <- "T_cell"
#' sce$subtype   <- rep(c("TypeA", "TypeB"), times = 60)
#' sce$condition <- rep(c("ctrl", "treat"), each = 60)
#'
#' result <- compoundDE(
#'     sce,
#'     broad_type   = "T_cell",
#'     subtype_col  = "subtype",
#'     broad_col    = "cell_type",
#'     donor        = "donor",
#'     condition    = "condition",
#'     min_cells    = 3L,
#'     min_subtypes = 2L,
#'     min_donors   = 2L
#' )
#' result
#' head(as.data.frame(deTable(result)))
#'
#' @seealso \code{\link{plotDecomposition}}, \code{\link{plotProportion}},
#'   \code{\link{plotTCRatio}}, \code{\link{CDEResult}}
#'
#' @importFrom SummarizedExperiment colData assay assays
#' @importFrom S4Vectors DataFrame
#' @importFrom BiocParallel SerialParam
#' @importFrom methods is
#' @importFrom stats sd
#' @export
compoundDE <- function(sce,
                        broad_type,
                        subtype_col,
                        broad_col,
                        donor,
                        condition,
                        contrast       = NULL,
                        min_cells      = 10L,
                        min_subtypes   = 2L,
                        min_cpm        = 1,
                        min_donors     = 2L,
                        tc_thresh_high = 0.8,
                        tc_thresh_low  = 0.2,
                        assay_name     = "counts",
                        BPPARAM        = SerialParam()) {

    # ── Input validation ──────────────────────────────────────────────────────
    if (!is(sce, "SingleCellExperiment"))
        stop("'sce' must be a SingleCellExperiment object.")

    required_cols <- c(broad_col, subtype_col, donor, condition)
    missing_cols  <- required_cols[!required_cols %in% names(colData(sce))]
    if (length(missing_cols) > 0L)
        stop("Column(s) not found in colData(sce): ",
             paste(missing_cols, collapse = ", "))

    if (!assay_name %in% names(assays(sce)))
        stop("Assay '", assay_name, "' not found in sce.")

    if (!is.numeric(tc_thresh_high) ||
        !is.numeric(tc_thresh_low)  ||
        tc_thresh_high <= tc_thresh_low ||
        tc_thresh_high > 1 || tc_thresh_low < 0)
        stop("'tc_thresh_high' must be > 'tc_thresh_low', ",
             "both in [0, 1].")

    # ── Subset to broad type ──────────────────────────────────────────────────
    broad_vals <- as.character(colData(sce)[[broad_col]])
    if (!broad_type %in% broad_vals)
        stop("'broad_type' = '", broad_type, "' not found in '",
             broad_col, "' column.")

    sce_broad <- sce[, broad_vals == broad_type]

    # ── Validate condition ────────────────────────────────────────────────────
    cond_vals   <- as.character(colData(sce_broad)[[condition]])
    cond_levels <- sort(unique(cond_vals))
    if (length(cond_levels) != 2L)
        stop("'condition' must have exactly 2 levels. Found: ",
             paste(cond_levels, collapse = ", "))

    # ── Identify subtypes ─────────────────────────────────────────────────────
    subtypes <- sort(unique(as.character(
        colData(sce_broad)[[subtype_col]])))

    if (length(subtypes) < min_subtypes)
        stop("Found ", length(subtypes), " subtype(s) for '", broad_type,
             "'. Need at least ", min_subtypes, ". ",
             "Lower 'min_subtypes' or check '", subtype_col, "' column.")

    message(sprintf(
        "compoundDE: broad_type='%s' | %d subtypes: %s",
        broad_type, length(subtypes), paste(subtypes, collapse = ", ")
    ))

    # ── Subtype proportion matrix ─────────────────────────────────────────────
    message("Computing subtype proportions...")
    prop_mat  <- .computeProportions(sce_broad, donor, condition, subtype_col)
    delta_pi  <- .computeDeltaPi(prop_mat, cond_levels)

    # ── Mean log2 expression per subtype ──────────────────────────────────────
    message("Computing mean subtype expression...")
    mean_expr <- .computeMeanExpression(sce_broad, subtype_col, subtypes,
                                        assay_name)

    # ── Per-subtype DE ────────────────────────────────────────────────────────
    message("Running per-subtype differential expression...")
    subtype_de <- vector("list", length(subtypes))
    names(subtype_de) <- subtypes

    for (k in subtypes) {
        message("  Subtype: ", k)
        sce_k <- sce_broad[,
            as.character(colData(sce_broad)[[subtype_col]]) == k]

        subtype_de[[k]] <- tryCatch(
            .runSubtypeDE(sce_k,
                          donor      = donor,
                          condition  = condition,
                          cond_levels = cond_levels,
                          contrast   = contrast,
                          min_cells  = min_cells,
                          min_donors = min_donors,
                          min_cpm    = min_cpm,
                          assay_name = assay_name),
            error = function(e) {
                message("    Skipped: ", conditionMessage(e))
                NULL
            }
        )
    }

    subtype_de_valid <- Filter(Negate(is.null), subtype_de)

    if (length(subtype_de_valid) < 2L)
        stop("Fewer than 2 subtypes produced valid DE results. ",
             "Try lowering 'min_cells' or 'min_donors'.")

    message(sprintf("Valid subtypes for decomposition: %s",
                    paste(names(subtype_de_valid), collapse = ", ")))

    # ── Broad pseudo-bulk DE ──────────────────────────────────────────────────
    message("Running broad pseudo-bulk DE...")
    broad_de <- tryCatch(
        .runSubtypeDE(sce_broad,
                      donor       = donor,
                      condition   = condition,
                      cond_levels = cond_levels,
                      contrast    = contrast,
                      min_cells   = 1L,
                      min_donors  = min_donors,
                      min_cpm     = min_cpm,
                      assay_name  = assay_name),
        error = function(e) stop("Broad DE failed: ", conditionMessage(e))
    )

    message(sprintf("Broad DE: %d genes tested.", nrow(broad_de)))

    # ── Decomposition ─────────────────────────────────────────────────────────
    message("Decomposing DE signal into T and C components...")
    decomp <- .computeDecomposition(
        broad_de       = broad_de,
        subtype_de     = subtype_de_valid,
        mean_expr      = mean_expr,
        delta_pi       = delta_pi,
        prop_mat       = prop_mat,
        tc_thresh_high = tc_thresh_high,
        tc_thresh_low  = tc_thresh_low
    )

    n_t <- sum(decomp[["source"]] == "transcriptional", na.rm = TRUE)
    n_c <- sum(decomp[["source"]] == "compositional",   na.rm = TRUE)
    message(sprintf("Decomposition complete: %d transcriptional, %d compositional.",
                    n_t, n_c))

    # ── Return CDEResult ──────────────────────────────────────────────────────
    CDEResult(
        deTable            = decomp,
        subtypeProportions = prop_mat,
        subtypeDE          = subtype_de_valid,
        params = list(
            broad_type     = broad_type,
            broad_col      = broad_col,
            subtype_col    = subtype_col,
            subtypes       = names(subtype_de_valid),
            donor          = donor,
            condition      = condition,
            contrast       = contrast,
            min_cells      = min_cells,
            min_subtypes   = min_subtypes,
            min_cpm        = min_cpm,
            min_donors     = min_donors,
            tc_thresh_high = tc_thresh_high,
            tc_thresh_low  = tc_thresh_low,
            cond_levels    = cond_levels
        )
    )
}

# ── Internal assays accessor ──────────────────────────────────────────────────

#' @keywords internal
#' @noRd
assays <- function(x) SummarizedExperiment::assays(x)
