test_that("GLS (Woodbury) reproduces the Python reference and matches Cholesky", {
    df <- load_sample()
    args <- list(data = df, depvar = "ret", indepvars = "lag_LNMV",
                 event_date = 0, pre_start = -200, pre_end = -1,
                 panel_var = "permno", time_var = "tdate",
                 sample = quote(abs(prc) > 5), method = "gls", npc = 100,
                 verbose = FALSE)

    fit_w <- do.call(csestudy, c(args, list(solver = "woodbury")))
    fit_c <- do.call(csestudy, c(args, list(solver = "cholesky")))

    expect_equal(fit_w$method, "GLS (Woodbury)")
    expect_identical(fit_w$n_obs_event, 282L)

    # vs Python Woodbury reference
    expect_equal(unname(coef(fit_w)),
                 c(0.0018204128992898052, -0.01675526178588034),
                 tolerance = 1e-5)
    expect_equal(unname(fit_w$p_parametric),
                 c(0.015442603034212424, 0.11594928719763334),
                 tolerance = 1e-5)

    # Woodbury and Cholesky solve the same system -> agree closely
    expect_equal(unname(coef(fit_w)), unname(coef(fit_c)), tolerance = 1e-6)
    expect_equal(unname(fit_w$p_cdf), unname(fit_c$p_cdf), tolerance = 1e-8)
})

test_that("solver = 'woodbury' requires GLS", {
    df <- data.frame(permno = 1, tdate = 1, ret = 0)
    expect_error(
        csestudy(df, depvar = "ret", event_date = 3, pre_start = 1, pre_end = 2,
                 time_var = "tdate", method = "ols", solver = "woodbury"),
        "woodbury"
    )
})
