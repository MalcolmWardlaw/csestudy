# Significance statistics from the event-date and pre-event coefficient
# vectors. Mirrors `_get_significance_stats()` (csestudy.mata) and the tail of
# `fit()` (core.py).
#
# `all_betas` is an (L+1) x k matrix whose first row is the event-date estimate
# and whose remaining L rows are the pre-event estimates. `beta_event` is the
# first row.
.significance <- function(all_betas, beta_event) {
    pre <- all_betas[-1L, , drop = FALSE]
    pre_mean <- colMeans(pre, na.rm = TRUE)
    pre_std <- sqrt(apply(pre, 2L, stats::var, na.rm = TRUE))   # ddof = 1
    Lp1 <- nrow(all_betas)                                      # L + 1

    # Empirical CDF p-value: weak inequality over all L+1 dates, denominator L+1.
    dev <- abs(sweep(all_betas, 2L, pre_mean, "-"))
    dev_evt <- abs(beta_event - pre_mean)
    ge <- sweep(dev, 2L, dev_evt, FUN = ">=")
    p_cdf <- colSums(ge, na.rm = TRUE) / Lp1

    # Parametric two-tailed p-value with a (L+1)/L variance adjustment.
    t_adj <- Lp1 / (Lp1 - 1)
    z <- abs(beta_event - pre_mean) / (pre_std * sqrt(t_adj))
    p_parametric <- 2 * stats::pt(z, df = Lp1 - 2, lower.tail = FALSE)

    list(p_cdf = p_cdf, p_parametric = p_parametric,
         pre_mean = pre_mean, pre_std = pre_std)
}
