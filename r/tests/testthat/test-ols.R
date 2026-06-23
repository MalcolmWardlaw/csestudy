# Reference values captured from the Python package (python/csestudy) on
# examples/sample_data.csv with the canonical invocation
#   depvar=ret, indepvars=lag_LNMV, sample abs(prc)>5,
#   event_date=0, pre_start=-200, pre_end=-1  (date 17811 -> sequential 0).

test_that("OLS reproduces the Python reference", {
    df <- load_sample()
    fit <- csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                    event_date = 0, pre_start = -200, pre_end = -1,
                    panel_var = "permno", time_var = "tdate",
                    sample = abs(prc) > 5, method = "ols")

    expect_s3_class(fit, "csestudy")
    expect_equal(fit$method, "OLS")
    expect_identical(fit$n_obs_event, 284L)
    expect_identical(fit$n_pre_event_days, 200L)
    expect_identical(names(coef(fit)), c("lag_LNMV", "_cons"))

    expect_equal(unname(coef(fit)),
                 c(0.0015922394478737414, -0.019274631105643605),
                 tolerance = 1e-8)
    expect_equal(unname(fit$p_cdf),
                 c(0.05970149253731343, 0.25870646766169153),
                 tolerance = 1e-8)
    expect_equal(unname(fit$p_parametric),
                 c(0.056199962498204144, 0.26060993858311726),
                 tolerance = 1e-8)
})
