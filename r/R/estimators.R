# Core numerical routines. These are internal (not exported) and are written
# to reproduce, bit-for-bit where possible, the Python reference in
# python/csestudy/core.py and the Mata reference in csestudy.mata.

# PCA covariance decomposition.
#
# `R` is a (T_pre x N) pre-event return matrix (rows = pre-event days,
# columns = firms). Returns the factor loadings `V` (N x npc), the factor
# variances `lam` (npc), and the idiosyncratic variances `d` (N), such that
#   Omega = V diag(lam) V' + diag(d).
#
# Mirrors `_pca_decompose()` (core.py): demean each column (per firm), economy
# SVD, keep the first `npc` right singular vectors, ddof = 1 variances, floor
# `d` at 1e-15. Note: base `svd()` returns `v` directly (the right singular
# vectors), unlike NumPy which returns `Vt`, so no transpose is needed.
.pca_decompose <- function(R, npc) {
    A <- sweep(R, 2L, colMeans(R), "-")          # demean each firm over time
    sv <- svd(A)                                 # economy SVD; sv$v is N x min(T,N)
    V <- sv$v[, seq_len(npc), drop = FALSE]      # (N x npc) loadings
    scores <- A %*% V                            # (T_pre x npc)
    lam <- apply(scores, 2L, stats::var)         # ddof = 1
    resid <- A - tcrossprod(scores, V)           # A - scores %*% t(V)
    d <- apply(resid, 2L, stats::var)            # ddof = 1, per firm
    d <- pmax(d, 1e-15)
    list(V = V, lam = lam, d = d)
}

# OLS coefficients via the normal equations, mirroring Mata's
# `cholsolve(quadcross(X,X), quadcross(X,y))`. Column order of `X` is
# [indepvars, intercept] (intercept LAST).
.ols_fit <- function(y, X) {
    drop(solve(crossprod(X), crossprod(X, y)))
}

# GLS via Cholesky factorisation of Omega (default, most precise).
#
# Omega = V diag(lam) V' + diag(d). We factor Omega = U'U with base `chol()`
# (UPPER triangular), so the lower factor scipy uses is L = U'. Whitening
# y_w = L^{-1} y is `backsolve(U, y, transpose = TRUE)` (solves U' z = y).
.gls_cholesky <- function(y, X, V, lam, d) {
    Omega <- V %*% (diag(lam, nrow = length(lam)) %*% t(V))
    diag(Omega) <- diag(Omega) + d
    Omega <- (Omega + t(Omega)) / 2              # force symmetry (numerical)
    U <- chol(Omega)                             # upper: Omega = t(U) %*% U
    y_w <- backsolve(U, y, transpose = TRUE)     # L^{-1} y, L = t(U)
    X_w <- backsolve(U, X, transpose = TRUE)
    drop(solve(crossprod(X_w), crossprod(X_w, y_w)))
}

# GLS via the Woodbury matrix identity (faster, slightly less precise).
#
# Omega^{-1} = D^{-1} - D^{-1} V M V' D^{-1},  M = (Lambda^{-1} + V' D^{-1} V)^{-1}
# Only an (npc x npc) matrix is inverted. Mirrors `_gls_woodbury()` (core.py),
# including clipping `lam` at 1e-15 inside Lambda^{-1}.
.gls_woodbury <- function(y, X, V, lam, d) {
    d_inv <- 1 / d
    VtDinv <- sweep(t(V), 2L, d_inv, "*")        # (npc x N): V' D^{-1}
    Lambda_inv <- diag(1 / pmax(lam, 1e-15), nrow = length(lam))
    M <- solve(Lambda_inv + VtDinv %*% V)        # (npc x npc)
    omega_inv_times <- function(B) {
        # B is (N x m); d_inv recycles down rows -> row i scaled by d_inv[i]
        (d_inv * B) - (d_inv * V) %*% (M %*% (VtDinv %*% B))
    }
    Oinv_y <- omega_inv_times(y)
    Oinv_X <- omega_inv_times(X)
    drop(solve(crossprod(X, Oinv_X), crossprod(X, Oinv_y)))
}
