test_that("GLS (Cholesky) reproduces the Python reference", {
    df <- load_sample()
    fit <- csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                    event_date = 0, pre_start = -200, pre_end = -1,
                    panel_var = "permno", time_var = "tdate",
                    sample = abs(prc) > 5, method = "gls", npc = 100,
                    solver = "cholesky", verbose = FALSE)

    expect_equal(fit$method, "GLS")
    expect_identical(fit$n_obs_event, 282L)
    expect_identical(fit$n_pre_event_days, 200L)

    expect_equal(unname(coef(fit)),
                 c(0.0018204128992898106, -0.016755261785880302),
                 tolerance = 1e-6)
    expect_equal(unname(fit$p_cdf),
                 c(0.024875621890547265, 0.13930348258706468),
                 tolerance = 1e-8)
    expect_equal(unname(fit$p_parametric),
                 c(0.015442603034211867, 0.11594928719763319),
                 tolerance = 1e-6)
})
