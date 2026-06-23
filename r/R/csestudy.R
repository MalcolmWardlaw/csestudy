#' Time-series inference for a cross-sectional event study
#'
#' Estimates the event-period cross-sectional relationship between `depvar` and
#' `indepvars` and benchmarks it against the distribution of the same
#' relationship estimated on a window of pre-event ("pseudo-event") days. This
#' is the R implementation of Cohn, Johnson, Liu and Wardlaw (2026); it mirrors
#' the Stata `csestudy` command and the Python `CSEventStudy` class shipped in
#' the same repository.
#'
#' @param data A long-form panel data frame with one row per (panel, time).
#' @param depvar Character scalar naming the dependent variable column
#'   (default `"ret"`).
#' @param indepvars Character vector of independent-variable column names, or
#'   `NULL` for an intercept-only model.
#' @param event_date,pre_start,pre_end Values of `time_var` identifying the
#'   event date, the first (earliest) pre-event date, and the last (latest)
#'   pre-event date. Require `event_date > pre_end > pre_start`.
#' @param panel_var,time_var Character scalars naming the panel-identifier and
#'   time columns (defaults `"permno"` and `"date"`). `time_var` should be a
#'   sequential integer with non-trading days omitted.
#' @param sample Optional unquoted expression evaluated in `data` used to subset
#'   rows before estimation (the analogue of Stata's `[if]`), e.g.
#'   `sample = abs(prc) > 5`. Note: like the Python implementation, the filter
#'   is applied to the whole panel up front, so it also restricts the data used
#'   to build the GLS covariance matrix.
#' @param method `"ols"` (default) or `"gls"`.
#' @param npc Number of principal components for the GLS covariance matrix
#'   (default 100). Must not exceed the number of pre-event days.
#' @param solver GLS solver: `"cholesky"` (default, most precise) or
#'   `"woodbury"` (faster, slightly less precise). Requires `method = "gls"`.
#' @param verbose If `TRUE`, print progress during the (slower) GLS pre-event
#'   loop.
#'
#' @return An object of class `"csestudy"`: a list with elements `method`,
#'   `params` (named event-date coefficients), `param_names`, `p_cdf`
#'   (empirical CDF p-values), `p_parametric`, `pre_event_betas` (L x k),
#'   `n_obs_event`, `n_obs_all`, `n_pre_event_days`, and `call`. See
#'   [coef.csestudy()] and [summary.csestudy()].
#'
#' @references Cohn, J. B., Johnson, T. L., Liu, Z., and Wardlaw, M. I. (2026).
#'   Past is Prologue: Inference from the Cross Section of Returns Around an
#'   Event. *Journal of Financial Economics* 180, 104278.
#'   \doi{10.1016/j.jfineco.2026.104278}
#'
#' @examples
#' \dontrun{
#' df <- read.csv("examples/sample_data.csv")
#' tds <- sort(unique(df$date))
#' df$tdate <- match(df$date, tds) - match(17811L, tds)   # 0 at the event date
#'
#' fit <- csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
#'                 event_date = 0, pre_start = -200, pre_end = -1,
#'                 panel_var = "permno", time_var = "tdate",
#'                 sample = abs(prc) > 5, method = "gls", npc = 100)
#' summary(fit)
#' coef(fit)
#' }
#' @export
csestudy <- function(data,
                     depvar = "ret",
                     indepvars = NULL,
                     event_date,
                     pre_start,
                     pre_end,
                     panel_var = "permno",
                     time_var = "date",
                     sample = NULL,
                     method = c("ols", "gls"),
                     npc = 100L,
                     solver = c("cholesky", "woodbury"),
                     verbose = TRUE) {
    method <- match.arg(method)
    solver <- match.arg(solver)
    cl <- match.call()

    if (solver == "woodbury" && method != "gls")
        stop("solver = 'woodbury' requires method = 'gls'")

    # ---- optional row filter (non-standard evaluation, like Stata's [if]) ----
    sm <- substitute(sample)
    if (!is.null(sm)) {
        keep <- eval(sm, data, parent.frame())
        if (!is.logical(keep))
            stop("`sample` must evaluate to a logical vector")
        data <- data[!is.na(keep) & keep, , drop = FALSE]
    }

    # ---- column selection and listwise deletion ----
    needed <- c(panel_var, time_var, depvar, indepvars)
    missing_cols <- setdiff(needed, names(data))
    if (length(missing_cols))
        stop("columns not found in `data`: ", paste(missing_cols, collapse = ", "))
    data <- data[needed]
    data <- data[stats::complete.cases(data), , drop = FALSE]
    if (nrow(data) == 0L)
        stop("no observations remain after filtering and listwise deletion")

    panel_id <- data[[panel_var]]
    time_id <- data[[time_var]]
    y_all <- as.numeric(data[[depvar]])
    if (length(indepvars)) {
        X_rhs <- data.matrix(data[indepvars])
        storage.mode(X_rhs) <- "double"
    } else {
        X_rhs <- matrix(numeric(0), nrow = nrow(data), ncol = 0L)
    }

    # ---- rectangular (panel x time) lookups ----
    panels <- sort(unique(panel_id))
    times <- sort(unique(time_id))
    n_panels <- length(panels)
    n_times <- length(times)
    pcode <- match(panel_id, panels)
    tcode <- match(time_id, times)
    rect_row <- matrix(NA_integer_, n_panels, n_times)
    rect_row[cbind(pcode, tcode)] <- seq_along(panel_id)
    rect_y <- matrix(NA_real_, n_panels, n_times)
    rect_y[cbind(pcode, tcode)] <- y_all

    # ---- resolve date values to (positional) column indices ----
    event_col <- match(event_date, times)
    pre_start_col <- match(pre_start, times)
    pre_end_col <- match(pre_end, times)
    if (anyNA(c(event_col, pre_start_col, pre_end_col)))
        stop("event_date, pre_start, and pre_end must all be present in `", time_var, "`")
    if (event_col <= pre_end_col)
        stop("event_date must be after pre_end")
    if (pre_end_col <= pre_start_col)
        stop("pre_end must be after pre_start")
    n_pre <- pre_end_col - pre_start_col + 1L
    if (method == "gls" && npc > n_pre)
        stop("npc must be <= number of pre-event days (", n_pre, ")")
    pre_window_len <- event_col - pre_start_col
    if (method == "gls" && (pre_start_col - pre_window_len) < 1L)
        stop(sprintf(paste0("GLS requires %d observations before pre_start; data only ",
                            "extends %d periods before it."),
                     pre_window_len, pre_start_col - 1L))

    param_names <- c(indepvars, "_cons")
    n_params <- length(param_names)

    # ---- estimate coefficients for a single target date ----
    run_date <- function(target_col) {
        valid <- !is.na(rect_row[, target_col])

        if (method == "gls") {
            gls_pre_start <- target_col - (event_col - pre_start_col)
            gls_pre_end <- target_col - (event_col - pre_end_col)
            win <- gls_pre_start:gls_pre_end
            cols_needed <- c(win, target_col)
            present <- !is.na(rect_row[, cols_needed, drop = FALSE])
            valid <- valid & (rowSums(present) == length(cols_needed))
            # exclude (near-)all-zero return histories
            tot_abs <- rowSums(abs(rect_y[, win, drop = FALSE]), na.rm = TRUE)
            valid <- valid & (tot_abs >= 0.01)
        }

        mask <- which(valid)
        n_valid <- length(mask)
        if (n_valid == 0L)
            return(list(beta = rep(NA_real_, n_params), n = 0L))

        rows <- rect_row[mask, target_col]
        y <- y_all[rows]
        if (ncol(X_rhs) > 0L) {
            X <- cbind(X_rhs[rows, , drop = FALSE], 1)
        } else {
            X <- matrix(1, n_valid, 1L)
        }

        if (method == "gls") {
            R_pre <- t(rect_y[mask, win, drop = FALSE])    # (T_pre x N_valid)
            dec <- .pca_decompose(R_pre, npc)
            beta <- if (solver == "woodbury")
                .gls_woodbury(y, X, dec$V, dec$lam, dec$d)
            else
                .gls_cholesky(y, X, dec$V, dec$lam, dec$d)
        } else {
            beta <- .ols_fit(y, X)
        }
        list(beta = as.numeric(beta), n = n_valid)
    }

    # ---- event date ----
    ev <- run_date(event_col)
    if (ev$n == 0L)
        stop("no valid observations on the event date after filtering")

    # ---- pre-event dates (pre_end down to pre_start) ----
    all_betas <- matrix(NA_real_, n_pre + 1L, n_params,
                        dimnames = list(NULL, param_names))
    all_nobs <- integer(n_pre + 1L)
    all_betas[1L, ] <- ev$beta
    all_nobs[1L] <- ev$n

    pre_cols <- pre_end_col:pre_start_col
    for (i in seq_along(pre_cols)) {
        rd <- run_date(pre_cols[i])
        all_betas[i + 1L, ] <- rd$beta
        all_nobs[i + 1L] <- rd$n
        if (verbose && method == "gls" && (i %% 25L == 0L || i == n_pre))
            message(sprintf("  pre-event date %d / %d", i, n_pre))
    }

    sig <- .significance(all_betas, ev$beta)

    method_label <- if (method == "gls") {
        if (solver == "woodbury") "GLS (Woodbury)" else "GLS"
    } else "OLS"

    structure(
        list(method = method_label,
             params = stats::setNames(ev$beta, param_names),
             param_names = param_names,
             p_cdf = stats::setNames(sig$p_cdf, param_names),
             p_parametric = stats::setNames(sig$p_parametric, param_names),
             pre_event_betas = all_betas[-1L, , drop = FALSE],
             n_obs_event = ev$n,
             n_obs_all = all_nobs,
             n_pre_event_days = n_pre,
             call = cl),
        class = "csestudy"
    )
}
