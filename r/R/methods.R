#' @describeIn csestudy Extract the event-date coefficient vector.
#' @param object A `"csestudy"` object.
#' @param ... Ignored.
#' @export
coef.csestudy <- function(object, ...) {
    object$params
}

#' @describeIn csestudy Print a compact results table.
#' @param x A `"csestudy"` object.
#' @export
print.csestudy <- function(x, ...) {
    cat(sprintf("%s estimates with time-series-corrected errors\n", x$method))
    cat(sprintf("Number of obs              = %9d\n", x$n_obs_event))
    cat(sprintf("Number of pre-event dates  = %9d\n", x$n_pre_event_days))
    bar <- paste0(strrep("-", 13), "+", strrep("-", 38))
    cat(bar, "\n")
    cat(sprintf("%12s | %11s %9s %11s\n", "", "Coef.", "CDF p", "Param p"))
    cat(bar, "\n")
    for (j in seq_along(x$param_names))
        cat(sprintf("%12s | %11.5g %9.3f %11.3f\n",
                    x$param_names[j], x$params[j], x$p_cdf[j], x$p_parametric[j]))
    cat(bar, "\n")
    invisible(x)
}

#' @describeIn csestudy Build a summary table (a data frame of coefficients and
#'   p-values) with its own print method.
#' @export
summary.csestudy <- function(object, ...) {
    tab <- data.frame(
        coef = unname(object$params),
        p_cdf = unname(object$p_cdf),
        p_parametric = unname(object$p_parametric),
        row.names = object$param_names
    )
    structure(
        list(method = object$method,
             n_obs_event = object$n_obs_event,
             n_pre_event_days = object$n_pre_event_days,
             coefficients = tab),
        class = "summary.csestudy"
    )
}

#' @export
print.summary.csestudy <- function(x, ...) {
    cat(sprintf("%s estimates with time-series-corrected errors\n", x$method))
    cat(sprintf("Number of obs = %d;  pre-event dates = %d\n\n",
                x$n_obs_event, x$n_pre_event_days))
    print(format(x$coefficients, digits = 5))
    invisible(x)
}
