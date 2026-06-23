# Locate and load the shared synthetic sample panel from the repository root.
# The 6.9 MB CSV is intentionally NOT vendored into the package, so tests that
# need it skip cleanly when run from an installed copy without the repo.
sample_data_path <- function() {
    candidates <- c(
        testthat::test_path("..", "..", "..", "examples", "sample_data.csv"),
        testthat::test_path("..", "..", "examples", "sample_data.csv"),
        "examples/sample_data.csv"
    )
    for (p in candidates) if (file.exists(p)) return(normalizePath(p))
    NA_character_
}

load_sample <- function() {
    p <- sample_data_path()
    if (is.na(p)) testthat::skip("examples/sample_data.csv not available")
    df <- utils::read.csv(p)
    tds <- sort(unique(df$date))
    df$tdate <- match(df$date, tds) - match(17811L, tds)  # 0 at the event date
    df
}

# Canonical invocation arguments shared across tests.
canonical_args <- function(df) {
    list(data = df, depvar = "ret", indepvars = "lag_LNMV",
         event_date = 0, pre_start = -200, pre_end = -1,
         panel_var = "permno", time_var = "tdate")
}
