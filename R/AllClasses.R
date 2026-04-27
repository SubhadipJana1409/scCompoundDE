#' @title CDEResult: Compound DE Result Container
#'
#' @description
#' An S4 class storing the output of \code{\link{compoundDE}}. Every
#' pseudo-bulk DE gene is assigned a \code{TC_ratio} score — the fraction
#' of its observed fold-change attributable to transcriptional (cell-intrinsic)
#' versus compositional (subtype proportion shift) signal.
#'
#' @slot deTable A \code{DataFrame} with one row per gene containing
#'   \code{logFC}, \code{AveExpr}, \code{t}, \code{P.Value},
#'   \code{adj.P.Val}, \code{B} (from the broad limma model),
#'   \code{T_score}, \code{C_score}, \code{T_score_z}, \code{C_score_z},
#'   \code{TC_ratio}, and \code{source} (transcriptional / compositional /
#'   mixed).
#' @slot subtypeProportions A \code{matrix} of subtype proportions with
#'   rows = samples (donor___condition) and columns = subtypes.
#' @slot subtypeDE A named \code{list} of per-subtype \code{DataFrame}s
#'   from the limma DE models.
#' @slot params A \code{list} of analysis parameters.
#'
#' @exportClass CDEResult
#' @importFrom methods new is
#' @importFrom S4Vectors DataFrame
setClass(
    "CDEResult",
    representation(
        deTable            = "DataFrame",
        subtypeProportions = "matrix",
        subtypeDE          = "list",
        params             = "list"
    )
)

# ── Constructor ───────────────────────────────────────────────────────────────

#' @title Constructor for CDEResult
#'
#' @description Create a new \code{CDEResult} object.
#'
#' @param deTable A \code{DataFrame} of per-gene compound DE statistics.
#' @param subtypeProportions A \code{matrix} of subtype proportions
#'   (samples x subtypes).
#' @param subtypeDE A named \code{list} of per-subtype DE \code{DataFrame}s.
#' @param params A \code{list} of analysis parameters.
#'
#' @return A \code{CDEResult} object.
#'
#' @examples
#' library(S4Vectors)
#' dt <- DataFrame(
#'     gene      = c("G1", "G2"),
#'     logFC     = c(1.2, 0.1),
#'     AveExpr   = c(3.1, 2.0),
#'     t         = c(4.1, 0.5),
#'     P.Value   = c(0.001, 0.8),
#'     adj.P.Val = c(0.01, 0.9),
#'     B         = c(2.1, -2.0),
#'     T_score   = c(1.1, 0.05),
#'     C_score   = c(0.1, 0.09),
#'     T_score_z = c(1.5, 0.2),
#'     C_score_z = c(0.2, 0.3),
#'     TC_ratio  = c(0.88, 0.40),
#'     source    = c("transcriptional", "mixed")
#' )
#' pm <- matrix(c(0.6, 0.4, 0.3, 0.7), nrow = 2,
#'              dimnames = list(c("D1___ctrl","D1___treat"),
#'                              c("TypeA","TypeB")))
#' obj <- CDEResult(deTable = dt, subtypeProportions = pm,
#'                  subtypeDE = list(), params = list(broad_type = "T_cell"))
#' obj
#'
#' @export
CDEResult <- function(deTable, subtypeProportions, subtypeDE,
                      params = list()) {
    new("CDEResult",
        deTable            = deTable,
        subtypeProportions = subtypeProportions,
        subtypeDE          = subtypeDE,
        params             = params)
}

# ── Generics ──────────────────────────────────────────────────────────────────

#' @title Accessor for the DE table in a CDEResult
#'
#' @description Returns the per-gene compound DE statistics from a
#'   \code{CDEResult} object, including the decomposition scores.
#'
#' @param x A \code{CDEResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A \code{DataFrame} with columns \code{logFC}, \code{P.Value},
#'   \code{adj.P.Val}, \code{T_score}, \code{C_score}, \code{TC_ratio},
#'   and \code{source}.
#'
#' @examples
#' library(S4Vectors)
#' dt <- DataFrame(gene = "G1", logFC = 1.2, P.Value = 0.001,
#'                 adj.P.Val = 0.01, AveExpr = 3.1, t = 4.1, B = 2.1,
#'                 T_score = 1.1, C_score = 0.1,
#'                 T_score_z = 1.5, C_score_z = 0.2,
#'                 TC_ratio = 0.88, source = "transcriptional")
#' pm <- matrix(c(0.6, 0.4), nrow = 1,
#'              dimnames = list("D1___ctrl", c("TypeA", "TypeB")))
#' obj <- CDEResult(dt, pm, list(), list(broad_type = "T_cell"))
#' deTable(obj)
#'
#' @export
setGeneric("deTable", function(x, ...) standardGeneric("deTable"))

#' @title Accessor for subtype proportions in a CDEResult
#'
#' @description Returns the subtype proportion matrix from a
#'   \code{CDEResult} object.
#'
#' @param x A \code{CDEResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A \code{matrix} with rows = samples and columns = subtypes.
#'
#' @examples
#' library(S4Vectors)
#' dt <- DataFrame(gene = "G1", logFC = 1.2, P.Value = 0.001,
#'                 adj.P.Val = 0.01, AveExpr = 3.1, t = 4.1, B = 2.1,
#'                 T_score = 1.1, C_score = 0.1,
#'                 T_score_z = 1.5, C_score_z = 0.2,
#'                 TC_ratio = 0.88, source = "transcriptional")
#' pm <- matrix(c(0.6, 0.4), nrow = 1,
#'              dimnames = list("D1___ctrl", c("TypeA", "TypeB")))
#' obj <- CDEResult(dt, pm, list(), list(broad_type = "T_cell"))
#' subtypeProportions(obj)
#'
#' @export
setGeneric("subtypeProportions",
           function(x, ...) standardGeneric("subtypeProportions"))

#' @title Accessor for per-subtype DE results in a CDEResult
#'
#' @description Returns the list of per-subtype limma DE \code{DataFrame}s.
#'
#' @param x A \code{CDEResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A named \code{list} of \code{DataFrame}s, one per subtype.
#'
#' @examples
#' library(S4Vectors)
#' dt <- DataFrame(gene = "G1", logFC = 1.2, P.Value = 0.001,
#'                 adj.P.Val = 0.01, AveExpr = 3.1, t = 4.1, B = 2.1,
#'                 T_score = 1.1, C_score = 0.1,
#'                 T_score_z = 1.5, C_score_z = 0.2,
#'                 TC_ratio = 0.88, source = "transcriptional")
#' pm <- matrix(c(0.6, 0.4), nrow = 1,
#'              dimnames = list("D1___ctrl", c("TypeA", "TypeB")))
#' obj <- CDEResult(dt, pm, list(), list(broad_type = "T_cell"))
#' subtypeDE(obj)
#'
#' @export
setGeneric("subtypeDE", function(x, ...) standardGeneric("subtypeDE"))

#' @title Accessor for TC_ratio vector in a CDEResult
#'
#' @description Returns the per-gene TC_ratio vector — the fraction of
#'   DE signal attributable to transcriptional versus compositional change.
#'   Values near 1 = transcriptional (real biology);
#'   values near 0 = compositional (artifact).
#'
#' @param x A \code{CDEResult} object.
#' @param ... Additional arguments (not used).
#'
#' @return A named \code{numeric} vector of TC_ratio values.
#'
#' @examples
#' library(S4Vectors)
#' dt <- DataFrame(gene = c("G1","G2"), logFC = c(1.2,0.1),
#'                 P.Value = c(0.001,0.8), adj.P.Val = c(0.01,0.9),
#'                 AveExpr = c(3.1,2.0), t = c(4.1,0.5), B = c(2.1,-2.0),
#'                 T_score = c(1.1,0.05), C_score = c(0.1,0.09),
#'                 T_score_z = c(1.5,0.2), C_score_z = c(0.2,0.3),
#'                 TC_ratio = c(0.88,0.40), source = c("transcriptional","mixed"))
#' pm <- matrix(c(0.6,0.4,0.3,0.7), nrow=2,
#'              dimnames=list(c("D1___ctrl","D1___treat"),c("TypeA","TypeB")))
#' obj <- CDEResult(dt, pm, list(), list(broad_type = "T_cell"))
#' tcRatio(obj)
#'
#' @export
setGeneric("tcRatio", function(x, ...) standardGeneric("tcRatio"))

# ── Methods ───────────────────────────────────────────────────────────────────

#' @title Show method for CDEResult
#'
#' @description Prints a compact summary of a \code{CDEResult} object
#'   including the decomposition breakdown.
#'
#' @param object A \code{CDEResult} object.
#'
#' @return Invisibly returns \code{object}.
#'
#' @importFrom methods show
#' @export
setMethod("show", "CDEResult", function(object) {
    dt      <- object@deTable
    n_sig   <- sum(dt[["adj.P.Val"]] < 0.05, na.rm = TRUE)
    n_trans <- sum(dt[["source"]] == "transcriptional", na.rm = TRUE)
    n_comp  <- sum(dt[["source"]] == "compositional",   na.rm = TRUE)
    n_mixed <- sum(dt[["source"]] == "mixed",           na.rm = TRUE)

    cat("CDEResult\n")
    cat("  Broad type     :", object@params$broad_type, "\n")
    cat("  Subtypes       :",
        paste(object@params$subtypes, collapse = ", "), "\n")
    cat("  Genes tested   :", nrow(dt), "\n")
    cat("  Significant    :", n_sig, "(adj.P.Val < 0.05)\n")
    cat("  ── Decomposition ──────────────────\n")
    cat("  Transcriptional:", n_trans,
        "genes (TC_ratio >=", object@params$tc_thresh_high, ")\n")
    cat("  Compositional  :", n_comp,
        "genes (TC_ratio <=", object@params$tc_thresh_low, ")\n")
    cat("  Mixed          :", n_mixed, "genes\n")
    cat("  Condition      :", object@params$condition, "\n")
    invisible(object)
})

#' @rdname deTable
#' @export
setMethod("deTable", "CDEResult", function(x, ...) x@deTable)

#' @rdname subtypeProportions
#' @export
setMethod("subtypeProportions", "CDEResult",
          function(x, ...) x@subtypeProportions)

#' @rdname subtypeDE
#' @export
setMethod("subtypeDE", "CDEResult", function(x, ...) x@subtypeDE)

#' @rdname tcRatio
#' @export
setMethod("tcRatio", "CDEResult", function(x, ...) {
    ratio <- x@deTable[["TC_ratio"]]
    names(ratio) <- as.character(x@deTable[["gene"]])
    ratio
})
