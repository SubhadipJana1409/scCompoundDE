# ── Internal helper functions for scCompoundDE ────────────────────────────────
# These are not exported. All exported API lives in compoundDE.R and
# plotDecomposition.R.

# ── Pseudo-bulk aggregation ───────────────────────────────────────────────────

#' Build pseudo-bulk matrix for one group of cells
#'
#' Aggregates raw counts per sample label using a sparse indicator matrix.
#'
#' @param sce A \code{SingleCellExperiment} subset to the target cells.
#' @param sample_ids Character vector of sample labels (length = ncol(sce)).
#' @param assay_name Name of the count assay.
#'
#' @return A list with \code{matrix} (genes x samples),
#'   \code{weights} (sqrt cell counts), and \code{ncells}.
#'
#' @keywords internal
#' @noRd
.buildPseudobulk <- function(sce, sample_ids, assay_name = "counts") {
    sample_levels <- unique(sample_ids)
    n_samples     <- length(sample_levels)
    sample_idx    <- match(sample_ids, sample_levels)

    counts_mat <- assay(sce, assay_name)

    indicator <- Matrix::sparseMatrix(
        i        = seq_along(sample_idx),
        j        = sample_idx,
        x        = 1,
        dims     = c(ncol(sce), n_samples),
        dimnames = list(colnames(sce), sample_levels)
    )

    pb_mat  <- as.matrix(counts_mat %*% indicator)
    ncells  <- tabulate(sample_idx, nbins = n_samples)
    weights <- sqrt(ncells)
    names(ncells) <- sample_levels
    names(weights) <- sample_levels

    list(matrix = pb_mat, weights = weights, ncells = ncells)
}

# ── Subtype proportion matrix ─────────────────────────────────────────────────

#' Compute subtype proportion matrix
#'
#' For each (donor, condition) sample, compute the fraction of cells
#' belonging to each subtype within the broad cell type.
#'
#' @param sce SCE subset to the broad cell type.
#' @param donor Column name for donor.
#' @param condition Column name for condition.
#' @param subtype_col Column name for fine-grained subtypes.
#'
#' @return A \code{matrix} with rows = \code{donor___condition} samples
#'   and columns = subtypes (proportions summing to 1 per row).
#'
#' @keywords internal
#' @noRd
.computeProportions <- function(sce, donor, condition, subtype_col) {
    d_vals  <- as.character(colData(sce)[[donor]])
    c_vals  <- as.character(colData(sce)[[condition]])
    s_vals  <- as.character(colData(sce)[[subtype_col]])

    sample_id    <- paste(d_vals, c_vals, sep = "___")
    sample_lvls  <- unique(sample_id)
    subtype_lvls <- sort(unique(s_vals))

    # Cell counts per sample x subtype
    cnt <- table(sample_id, s_vals)
    cnt <- cnt[sample_lvls, subtype_lvls, drop = FALSE]

    row_totals <- rowSums(cnt)
    row_totals[row_totals == 0L] <- 1L    # guard against empty samples

    prop_mat <- sweep(cnt, 1L, row_totals, FUN = "/")
    as.matrix(prop_mat)
}

# ── Delta proportion per subtype ──────────────────────────────────────────────

#' Compute mean proportion change (condition 2 minus condition 1)
#'
#' @param prop_mat Matrix from \code{.computeProportions}.
#' @param cond_levels Two-element character vector of condition levels.
#'
#' @return Named numeric vector of Δπ per subtype.
#'
#' @keywords internal
#' @noRd
.computeDeltaPi <- function(prop_mat, cond_levels) {
    parts     <- strsplit(rownames(prop_mat), "___", fixed = TRUE)
    cond_row  <- vapply(parts, `[`, character(1L), 2L)

    pi_cond1  <- prop_mat[cond_row == cond_levels[1L], , drop = FALSE]
    pi_cond2  <- prop_mat[cond_row == cond_levels[2L], , drop = FALSE]

    mean1 <- if (nrow(pi_cond1) > 0L) colMeans(pi_cond1, na.rm = TRUE) else
                  setNames(rep(0, ncol(prop_mat)), colnames(prop_mat))
    mean2 <- if (nrow(pi_cond2) > 0L) colMeans(pi_cond2, na.rm = TRUE) else
                  setNames(rep(0, ncol(prop_mat)), colnames(prop_mat))

    mean2 - mean1
}

# ── Mean log2 expression per subtype ─────────────────────────────────────────

#' Compute mean log2 CPM expression per subtype
#'
#' @param sce SCE subset to the broad cell type.
#' @param subtype_col Column name for fine-grained subtypes.
#' @param subtypes Character vector of subtype names to process.
#' @param assay_name Name of the count assay.
#'
#' @return A \code{matrix} (genes x subtypes) of mean log2 CPM.
#'
#' @keywords internal
#' @noRd
.computeMeanExpression <- function(sce, subtype_col, subtypes,
                                    assay_name = "counts") {
    s_vals <- as.character(colData(sce)[[subtype_col]])
    counts <- assay(sce, assay_name)

    expr_list <- lapply(subtypes, function(k) {
        idx <- s_vals == k
        if (sum(idx) == 0L) return(rep(NA_real_, nrow(counts)))
        sub_counts  <- counts[, idx, drop = FALSE]
        lib_sizes   <- colSums(sub_counts)
        lib_sizes[lib_sizes == 0L] <- 1L
        cpm_mat     <- sweep(sub_counts, 2L, lib_sizes / 1e6, FUN = "/")
        log2_cpm    <- log2(as.matrix(cpm_mat) + 1)
        rowMeans(log2_cpm)
    })

    expr_mat <- do.call(cbind, expr_list)
    colnames(expr_mat) <- subtypes
    rownames(expr_mat) <- rownames(counts)
    expr_mat
}

# ── Per-subtype limma-voom DE ─────────────────────────────────────────────────

#' Run limma-voom DE for one subtype
#'
#' Handles both paired (same donors in multiple conditions) and unpaired
#' designs automatically. Donors with fewer than \code{min_cells} cells
#' in any condition are removed before aggregation.
#'
#' @param sce SCE subset to one subtype.
#' @param donor Donor column name.
#' @param condition Condition column name.
#' @param cond_levels Two-level character vector (ctrl, treat).
#' @param contrast Optional contrast string in limma syntax.
#' @param min_cells Minimum cells per donor-condition sample.
#' @param min_donors Minimum samples for gene expression filter.
#' @param min_cpm Minimum CPM threshold for gene filter.
#' @param assay_name Name of count assay.
#'
#' @return A \code{DataFrame} with columns \code{gene}, \code{logFC},
#'   \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val}, \code{B}.
#'
#' @importFrom limma voom lmFit makeContrasts contrasts.fit eBayes topTable
#' @importFrom stats model.matrix
#' @keywords internal
#' @noRd
.runSubtypeDE <- function(sce, donor, condition, cond_levels,
                           contrast   = NULL,
                           min_cells  = 10L,
                           min_donors = 2L,
                           min_cpm    = 1,
                           assay_name = "counts") {

    d_vals <- as.character(colData(sce)[[donor]])
    c_vals <- as.character(colData(sce)[[condition]])

    # ── Remove sparse donor-condition combos ──────────────────────────────────
    combo        <- paste(d_vals, c_vals, sep = "___")
    combo_counts <- table(combo)
    sparse       <- names(combo_counts)[combo_counts < min_cells]
    if (length(sparse) > 0L) {
        keep   <- !combo %in% sparse
        sce    <- sce[, keep]
        d_vals <- d_vals[keep]
        c_vals <- c_vals[keep]
    }

    n_donors <- length(unique(d_vals))
    if (n_donors < min_donors)
        stop("Fewer than ", min_donors, " donors after sparse filtering.")

    # ── Detect paired design ─────────────────────────────────────────────────
    donor_conds <- tapply(c_vals, d_vals, function(x) length(unique(x)))
    is_paired   <- any(donor_conds > 1L)

    # ── Build sample labels ───────────────────────────────────────────────────
    sample_ids <- if (is_paired) {
        paste(d_vals, c_vals, sep = "___")
    } else {
        d_vals
    }

    pb <- .buildPseudobulk(sce, sample_ids, assay_name)
    pb_mat  <- pb$matrix
    weights <- pb$weights

    # ── Gene filter ───────────────────────────────────────────────────────────
    lib_sizes <- colSums(pb_mat)
    lib_sizes[lib_sizes == 0L] <- 1L
    cpm_mat <- sweep(pb_mat, 2L, lib_sizes / 1e6, FUN = "/")
    keep    <- rowSums(cpm_mat >= min_cpm) >= min(min_donors, ncol(pb_mat))
    if (sum(keep) == 0L)
        stop("No genes passed the CPM expression filter.")
    pb_filt <- pb_mat[keep, , drop = FALSE]

    sample_levels <- colnames(pb_mat)

    # ── Build design matrix ───────────────────────────────────────────────────
    if (is_paired) {
        parts   <- strsplit(sample_levels, "___", fixed = TRUE)
        cond_f  <- factor(vapply(parts, `[`, character(1L), 2L),
                          levels = cond_levels)
        donor_f <- factor(vapply(parts, `[`, character(1L), 1L))

        n_params <- length(cond_levels) + length(levels(donor_f)) - 1L
        if (ncol(pb_filt) <= n_params)
            stop("Not enough samples for paired model.")

        design <- model.matrix(~ 0 + cond_f + donor_f)
        cond_cols <- make.names(cond_levels, unique = TRUE)
        colnames(design)[seq_along(cond_cols)] <- cond_cols
        d_idx <- seq(length(cond_cols) + 1L, ncol(design))
        colnames(design)[d_idx] <- make.names(
            sub("^donor_f", "", colnames(design)[d_idx]), unique = TRUE
        )
    } else {
        # Map each donor sample to its condition
        d_to_c <- tapply(c_vals, d_vals, function(x) {
            names(sort(table(x), decreasing = TRUE))[1L]
        })
        d_to_c    <- d_to_c[sample_levels]
        cond_f    <- factor(d_to_c, levels = cond_levels)
        cond_cols <- make.names(cond_levels, unique = TRUE)
        design    <- model.matrix(~ 0 + cond_f)
        colnames(design) <- cond_cols
    }

    # ── Contrast ──────────────────────────────────────────────────────────────
    contrast_str <- if (is.null(contrast)) {
        paste(cond_cols[2L], "-", cond_cols[1L])
    } else {
        contrast
    }

    # ── limma-voom ────────────────────────────────────────────────────────────
    wt_mat <- matrix(rep(weights, each = nrow(pb_filt)),
                     nrow = nrow(pb_filt))
    v    <- voom(pb_filt, design, weights = wt_mat, plot = FALSE)
    fit  <- lmFit(v, design, weights = v$weights)
    cmat <- makeContrasts(contrasts = contrast_str, levels = design)
    fit2 <- contrasts.fit(fit, cmat)
    fit2 <- eBayes(fit2)
    tt   <- topTable(fit2, number = Inf, sort.by = "none")

    S4Vectors::DataFrame(
        gene      = rownames(tt),
        logFC     = tt[["logFC"]],
        AveExpr   = tt[["AveExpr"]],
        t         = tt[["t"]],
        P.Value   = tt[["P.Value"]],
        adj.P.Val = tt[["adj.P.Val"]],
        B         = tt[["B"]],
        row.names = rownames(tt)
    )
}

# ── Core decomposition ────────────────────────────────────────────────────────

#' Decompose broad DE into transcriptional and compositional components
#'
#' For each gene:
#'   T_g = sum_k( pi_bar_k * logFC_gk )   [proportion-weighted transcript shift]
#'   C_g = sum_k( delta_pi_k * mu_gk )    [expression-weighted proportion shift]
#' Both are z-score normalised across genes before computing TC_ratio.
#'
#' @param broad_de DataFrame from broad pseudo-bulk DE.
#' @param subtype_de Named list of per-subtype DE DataFrames.
#' @param mean_expr Matrix (genes x subtypes) of mean log2 CPM.
#' @param delta_pi Named vector of Δπ per subtype.
#' @param prop_mat Proportion matrix (samples x subtypes).
#' @param tc_thresh_high Upper TC_ratio threshold for "transcriptional".
#' @param tc_thresh_low Lower TC_ratio threshold for "compositional".
#'
#' @return A \code{DataFrame} with all DE statistics plus decomposition.
#'
#' @keywords internal
#' @noRd
.computeDecomposition <- function(broad_de, subtype_de, mean_expr,
                                   delta_pi, prop_mat,
                                   tc_thresh_high, tc_thresh_low) {

    all_genes     <- as.character(broad_de[["gene"]])
    valid_sub     <- names(subtype_de)
    n_genes       <- length(all_genes)
    n_sub         <- length(valid_sub)

    # Mean proportion of each subtype across ALL samples (pi_bar)
    pi_bar <- colMeans(prop_mat[, valid_sub, drop = FALSE], na.rm = TRUE)

    # Δπ for valid subtypes only
    delta_pi_valid <- delta_pi[valid_sub]

    # ── logFC matrix: genes x subtypes (0 if gene not in that subtype DE) ────
    logfc_mat <- matrix(0, nrow = n_genes, ncol = n_sub,
                        dimnames = list(all_genes, valid_sub))
    for (k in valid_sub) {
        de_k   <- subtype_de[[k]]
        g_k    <- as.character(de_k[["gene"]])
        common <- intersect(all_genes, g_k)
        if (length(common) > 0L)
            logfc_mat[common, k] <- de_k[match(common, g_k), "logFC"]
    }

    # ── mu matrix: genes x subtypes (mean log2 CPM) ──────────────────────────
    mu_mat <- matrix(0, nrow = n_genes, ncol = n_sub,
                     dimnames = list(all_genes, valid_sub))
    for (k in valid_sub) {
        if (k %in% colnames(mean_expr)) {
            common <- intersect(all_genes, rownames(mean_expr))
            if (length(common) > 0L)
                mu_mat[common, k] <- mean_expr[common, k]
        }
    }
    mu_mat[is.na(mu_mat)] <- 0

    # ── T_g and C_g ───────────────────────────────────────────────────────────
    # T_g = logfc_mat %*% pi_bar  (1 x n_sub) × (n_sub x 1)
    T_g <- as.vector(logfc_mat %*% pi_bar[valid_sub])
    # C_g = mu_mat   %*% delta_pi (composition-weighted expression shift)
    C_g <- as.vector(mu_mat   %*% delta_pi_valid[valid_sub])
    names(T_g) <- all_genes
    names(C_g) <- all_genes

    # ── Z-score across genes (to make T and C comparable) ────────────────────
    .zscore <- function(x) {
        s <- stats::sd(x, na.rm = TRUE)
        m <- mean(x, na.rm = TRUE)
        if (is.na(s) || s < 1e-10) return(x - m)
        (x - m) / s
    }
    T_z <- .zscore(T_g)
    C_z <- .zscore(C_g)

    # ── TC_ratio ──────────────────────────────────────────────────────────────
    abs_T    <- abs(T_z)
    abs_C    <- abs(C_z)
    TC_ratio <- abs_T / (abs_T + abs_C + 1e-10)
    names(TC_ratio) <- all_genes

    # ── Classification ────────────────────────────────────────────────────────
    source <- rep("mixed", n_genes)
    source[TC_ratio >= tc_thresh_high] <- "transcriptional"
    source[TC_ratio <= tc_thresh_low]  <- "compositional"
    names(source) <- all_genes

    S4Vectors::DataFrame(
        gene      = all_genes,
        logFC     = broad_de[["logFC"]],
        AveExpr   = broad_de[["AveExpr"]],
        t         = broad_de[["t"]],
        P.Value   = broad_de[["P.Value"]],
        adj.P.Val = broad_de[["adj.P.Val"]],
        B         = broad_de[["B"]],
        T_score   = T_g,
        C_score   = C_g,
        T_score_z = T_z,
        C_score_z = C_z,
        TC_ratio  = TC_ratio,
        source    = source,
        row.names = all_genes
    )
}
