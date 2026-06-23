# Unit tests for the internal PCA decomposition on small fixed matrices.

test_that("full-rank decomposition reconstructs the sample covariance", {
    R <- matrix(c(0.01, -0.02,  0.03,  0.00,  0.015, -0.01,  0.02, -0.005,
                  0.00,  0.01, -0.02,  0.04, -0.01,   0.02, -0.03,  0.01,
                  0.02, -0.01,  0.00,  0.03,  0.01,  -0.02,  0.015, 0.00),
                nrow = 8, ncol = 3)
    dec <- csestudy:::.pca_decompose(R, npc = 3L)

    expect_equal(dim(dec$V), c(3L, 3L))
    expect_length(dec$lam, 3L)
    expect_length(dec$d, 3L)

    # With all components, V diag(lam) V' equals the (ddof = 1) sample covariance
    # of the demeaned columns, and the residual variance hits the floor.
    Omega_lowrank <- dec$V %*% diag(dec$lam) %*% t(dec$V)
    expect_equal(Omega_lowrank, stats::cov(R), tolerance = 1e-8)
    expect_true(all(dec$d <= 1e-10))

    # Factor variances are non-increasing (singular values are sorted).
    expect_false(is.unsorted(rev(dec$lam)))
})

test_that("idiosyncratic variance is floored for a rank-1 panel", {
    time_pattern <- c(1, -2, 3, -1, 2, 0.5)
    firm_scale <- c(0.01, 0.02, -0.015, 0.03)
    R <- outer(time_pattern, firm_scale)        # exactly rank 1
    dec <- csestudy:::.pca_decompose(R, npc = 1L)
    expect_true(all(dec$d == 1e-15))
})
