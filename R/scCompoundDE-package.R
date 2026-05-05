#' scCompoundDE: Compositional and Transcriptional Decomposition of
#' Pseudo-Bulk Differential Expression
#'
#' @description
#' scCompoundDE decomposes pseudo-bulk differential expression (DE) signals
#' into two orthogonal components: \strong{transcriptional} changes
#' (cell-intrinsic expression shifts) and \strong{compositional} changes
#' (shifts in the relative abundance of cell subtypes within a broad
#' population).
#'
#' Standard pseudo-bulk DE tools treat a pseudo-bulk sample as if it were
#' homogeneous, confounding true transcriptional change with artifactual signal
#' arising from subtype proportion shifts. For example, if T cells in disease
#' donors are predominantly exhausted (high PDCD1, TOX) while T cells in
#' healthy donors are predominantly naive (high IL7R, CCR7), a standard
#' pseudo-bulk DE analysis will report PDCD1 and IL7R as significantly DE --
#' but these genes changed because the \emph{composition} of T cells changed,
#' not because T cells themselves changed their transcriptome. scCompoundDE
#' detects and quantifies this confound for every tested gene.
#'
#' @section Key innovation:
#' For each gene, scCompoundDE computes a \strong{TC_ratio} score in
#' \code{[0, 1]}:
#' \itemize{
#'   \item \strong{TC_ratio ~= 1} -- gene is driven by cell-intrinsic
#'     transcriptional change (real biology).
#'   \item \strong{TC_ratio ~= 0} -- gene is driven by a shift in subtype
#'     proportions (compositional artefact).
#'   \item \strong{TC_ratio ~= 0.5} -- mixed signal (both components
#'     contribute).
#' }
#'
#' @section Main functions:
#' \describe{
#'   \item{\code{\link{compoundDE}}}{Run the full decomposition pipeline.
#'     Returns a \code{CDEResult} S4 object.}
#'   \item{\code{\link{filterGenesBySource}}}{Extract transcriptional,
#'     compositional, or mixed gene lists from a \code{CDEResult}.}
#'   \item{\code{\link{plotDecomposition}}}{Scatter plot of T_score vs
#'     C_score for all genes.}
#'   \item{\code{\link{plotProportion}}}{Stacked bar chart of subtype
#'     proportions per condition.}
#'   \item{\code{\link{plotTCRatio}}}{Histogram of TC_ratio distribution
#'     with classification thresholds.}
#' }
#'
#' @section How it fits into a scRNA-seq workflow:
#' \itemize{
#'   \item Use \code{scBatchQC} first to flag and remove low-quality cells.
#'   \item Use \code{scFastDE} or any pseudo-bulk tool to get DE genes.
#'   \item Use \code{scCompoundDE} to validate whether those DE genes are
#'     transcriptionally driven or compositional artefacts.
#' }
#'
#' @references
#' Crowell HL et al. (2020). muscat detects subpopulation-specific state
#' transitions from multi-sample multi-condition single-cell transcriptomics
#' data. \emph{Nature Communications}, 11, 6077.
#'
#' Then E et al. (2023). Distinguishing cell type composition and cell
#' type-specific effects in bulk tissues. \emph{bioRxiv}.
#'
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @importFrom rlang .data
#' @docType package
#' @name scCompoundDE-package
#' @aliases scCompoundDE
"_PACKAGE"
