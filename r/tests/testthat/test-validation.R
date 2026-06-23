# Input-validation errors. Use a small balanced synthetic panel; most of these
# errors fire before any estimation, so the data values are immaterial.

make_panel <- function(n_panel = 6, n_time = 20) {
    grid <- expand.grid(permno = seq_len(n_panel), tdate = seq_len(n_time))
    grid <- grid[order(grid$permno, grid$tdate), ]
    grid$ret <- as.numeric(seq_len(nrow(grid))) / 1000 - 0.05
    grid$lag_LNMV <- rep(seq_len(n_panel), each = n_time) / 2
    grid
}

test_that("event_date must be after pre_end", {
    df <- make_panel()
    expect_error(
        csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                 event_date = 2, pre_start = 1, pre_end = 5, time_var = "tdate"),
        "event_date must be after pre_end"
    )
})

test_that("pre_end must be after pre_start", {
    df <- make_panel()
    expect_error(
        csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                 event_date = 10, pre_start = 5, pre_end = 5, time_var = "tdate"),
        "pre_end must be after pre_start"
    )
})

test_that("npc cannot exceed the number of pre-event days", {
    df <- make_panel()
    expect_error(
        csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                 event_date = 15, pre_start = 5, pre_end = 8, time_var = "tdate",
                 method = "gls", npc = 10),
        "npc must be"
    )
})

test_that("GLS requires enough pre-window history", {
    df <- make_panel()
    expect_error(
        csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                 event_date = 15, pre_start = 5, pre_end = 8, time_var = "tdate",
                 method = "gls", npc = 3),
        "GLS requires"
    )
})

test_that("unknown dates raise an error", {
    df <- make_panel()
    expect_error(
        csestudy(df, depvar = "ret", indepvars = "lag_LNMV",
                 event_date = 999, pre_start = 1, pre_end = 5, time_var = "tdate"),
        "must all be present"
    )
})
