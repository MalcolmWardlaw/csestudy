"""
Core estimation routines for csestudy.

Mirrors the Stata csestudy.ado / csestudy.mata logic:
  - OLS or GLS cross-sectional regressions on event and pre-event dates
  - PCA-based covariance estimation for GLS (SVD, ddof=1 variances)
  - Cholesky GLS solver (default) or Woodbury identity (optional, faster)
  - Empirical CDF p-value and parametric t-statistic from pre-event distribution
"""

from __future__ import annotations

import sys
from dataclasses import dataclass, field
from typing import Literal, Sequence

import numpy as np
import pandas as pd
from scipy import linalg as sla


# ---------------------------------------------------------------------------
# Result container
# ---------------------------------------------------------------------------

@dataclass
class CSEStudyResult:
    """Results from a csestudy estimation."""

    method: str
    params: np.ndarray          # coefficient vector (event date)
    param_names: list[str]      # names matching params
    p_cdf: np.ndarray           # empirical CDF p-values
    p_parametric: np.ndarray    # parametric (t-distribution) p-values
    pre_event_betas: np.ndarray # (L x k) matrix of pre-event coefficients
    n_obs_event: int            # firms in event regression
    n_obs_all: np.ndarray       # (L+1,) obs counts per date
    n_pre_event_days: int

    def summary(self) -> str:
        """Print a formatted summary table matching Stata output style."""
        lines = []
        lines.append(f"{self.method} Estimates with Time Series Corrected Errors")
        lines.append(f"{'':>35s}Number of obs  = {self.n_obs_event:>9,d}")
        lines.append(f"{'':>23s}Number of pre-period dates = {self.n_pre_event_days:>9,d}")
        lines.append("-" * 13 + "+" + "-" * 47)
        lines.append(f"{'':>12s} | {'Coefficient':>11s}  {'CDF p-val':>9s}  {'Param p-val':>11s}")
        lines.append("-" * 13 + "+" + "-" * 47)
        for name, b, pc, pp in zip(
            self.param_names, self.params, self.p_cdf, self.p_parametric
        ):
            lines.append(f"{name:>12s} | {b:>11.5g}  {pc:>9.3f}  {pp:>11.3f}")
        lines.append("-" * 13 + "+" + "-" * 47)
        return "\n".join(lines)

    def __repr__(self) -> str:
        return self.summary()


# ---------------------------------------------------------------------------
# PCA covariance decomposition (Stata-matched)
# ---------------------------------------------------------------------------

def _pca_decompose(
    R: np.ndarray, npc: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Decompose a (T x N) return matrix into PCA loadings, factor variances,
    and idiosyncratic variances, matching Stata's csestudy conventions:
      - Demean per firm (column means).
      - Full SVD, keep first npc components.
      - ddof=1 variances throughout.

    Returns
    -------
    V : (N, npc) loadings
    lam : (npc,) factor variances
    d : (N,) idiosyncratic variances (floored at 1e-15)
    """
    A = R - R.mean(axis=0, keepdims=True)
    U, s, Vt = np.linalg.svd(A, full_matrices=False)
    V = Vt.T[:, :npc]
    scores = A @ V
    lam = np.var(scores, axis=0, ddof=1)
    resid = A - scores @ V.T
    d = np.var(resid, axis=0, ddof=1)
    d = np.clip(d, 1e-15, None)
    return V, lam, d


# ---------------------------------------------------------------------------
# GLS solvers
# ---------------------------------------------------------------------------

def _gls_cholesky(
    y: np.ndarray, X: np.ndarray, V: np.ndarray, lam: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """
    GLS via Cholesky decomposition of Omega.

    Omega = V diag(lam) V' + diag(d)

    Cholesky-factor Omega, pre-whiten y and X, then OLS.
    This is the default and most numerically precise approach.
    """
    N = len(d)
    Omega = V @ np.diag(lam) @ V.T + np.diag(d)
    # Force symmetry (numerical)
    Omega = (Omega + Omega.T) / 2.0
    L = sla.cholesky(Omega, lower=True)
    y_w = sla.solve_triangular(L, y, lower=True)
    X_w = sla.solve_triangular(L, X, lower=True)
    beta, *_ = np.linalg.lstsq(X_w, y_w, rcond=None)
    return beta


def _gls_woodbury(
    y: np.ndarray, X: np.ndarray, V: np.ndarray, lam: np.ndarray, d: np.ndarray
) -> np.ndarray:
    """
    GLS via the Woodbury matrix identity.

    Omega^{-1} = D^{-1} - D^{-1} V (Lambda^{-1} + V' D^{-1} V)^{-1} V' D^{-1}

    Only inverts a (npc x npc) matrix instead of factoring the full (N x N)
    Omega. ~50-66% faster but slightly less numerically precise.
    """
    d_inv = 1.0 / d
    VtDinv = V.T * d_inv[np.newaxis, :]  # (npc, N)
    Lambda_inv = np.diag(1.0 / np.clip(lam, 1e-15, None))
    M = sla.inv(Lambda_inv + VtDinv @ V)  # (npc, npc)

    def omega_inv_times(B: np.ndarray) -> np.ndarray:
        DinvB = d_inv[:, np.newaxis] * B
        return DinvB - (d_inv[:, np.newaxis] * V) @ (M @ (VtDinv @ B))

    Oinv_y = omega_inv_times(y[:, np.newaxis]).ravel()
    Oinv_X = omega_inv_times(X)
    XtOinvX = X.T @ Oinv_X
    XtOinvy = X.T @ Oinv_y
    beta = sla.solve(XtOinvX, XtOinvy, assume_a="pos")
    return beta


def _ols(y: np.ndarray, X: np.ndarray) -> np.ndarray:
    """OLS via QR decomposition."""
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    return beta


# ---------------------------------------------------------------------------
# Data preparation
# ---------------------------------------------------------------------------

def _build_rect_index(
    panel_id: np.ndarray, time_id: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int, int]:
    """
    Build a rectangular (n_panels x n_times) index mapping from long-form data.
    Returns row_index, unique_panels, unique_times, n_panels, n_times.
    row_index[i, j] = position in long data for panel i at time j, or -1 if missing.
    """
    unique_panels, panel_codes = np.unique(panel_id, return_inverse=True)
    unique_times = np.sort(np.unique(time_id))
    time_map = {t: j for j, t in enumerate(unique_times)}
    time_codes = np.array([time_map[t] for t in time_id])

    n_panels = len(unique_panels)
    n_times = len(unique_times)
    row_index = np.full((n_panels, n_times), -1, dtype=np.intp)
    for k in range(len(panel_id)):
        row_index[panel_codes[k], time_codes[k]] = k
    return row_index, unique_panels, unique_times, n_panels, n_times


# ---------------------------------------------------------------------------
# Main estimator class
# ---------------------------------------------------------------------------

class CSEventStudy:
    """
    Cross-sectional event study with time-series corrected inference.

    Parameters
    ----------
    df : pd.DataFrame
        Long-form panel with columns for panel id, time, dependent variable,
        and (optionally) independent variables.
    event_date : int or str
        Event date value in the time variable.
    pre_start : int or str
        First (earliest) pre-event date.
    pre_end : int or str
        Last (latest) pre-event date.
    depvar : str
        Name of the dependent variable column (default ``'ret'``).
    indepvars : list of str or None
        Names of independent variable columns. If ``None``, estimates an
        intercept-only model.
    panel_var : str
        Name of the panel identifier column (default ``'permno'``).
    time_var : str
        Name of the time variable column (default ``'date'``).
    sample : str or None
        A pandas query string applied before estimation (analogous to
        Stata's ``[if]``).
    method : ``'ols'`` or ``'gls'``
        Estimation method (default ``'ols'``).
    npc : int
        Number of principal components for GLS (default 100).
    solver : ``'cholesky'`` or ``'woodbury'``
        GLS solver. ``'cholesky'`` (default) is more numerically precise;
        ``'woodbury'`` is faster for large cross-sections.
    """

    def __init__(
        self,
        df: pd.DataFrame,
        *,
        event_date,
        pre_start,
        pre_end,
        depvar: str = "ret",
        indepvars: Sequence[str] | None = None,
        panel_var: str = "permno",
        time_var: str = "date",
        sample: str | None = None,
        method: Literal["ols", "gls"] = "ols",
        npc: int = 100,
        solver: Literal["cholesky", "woodbury"] = "cholesky",
    ):
        self.df = df
        self.event_date = event_date
        self.pre_start = pre_start
        self.pre_end = pre_end
        self.depvar = depvar
        self.indepvars = list(indepvars) if indepvars is not None else []
        self.panel_var = panel_var
        self.time_var = time_var
        self.sample = sample
        self.method = method.lower()
        self.npc = npc
        self.solver = solver.lower()

        if self.solver == "woodbury" and self.method != "gls":
            raise ValueError("solver='woodbury' requires method='gls'")

    def fit(self, verbose: bool = True) -> CSEStudyResult:
        """Run the estimation and return results."""
        # ---- data prep ----
        df = self.df.copy()
        if self.sample is not None:
            df = df.query(self.sample)

        needed_cols = [self.panel_var, self.time_var, self.depvar] + self.indepvars
        df = df[needed_cols].dropna()

        panel_id = df[self.panel_var].values
        time_id = df[self.time_var].values
        y_all = df[self.depvar].values.astype(float)
        if self.indepvars:
            X_rhs = df[self.indepvars].values.astype(float)
        else:
            X_rhs = np.empty((len(df), 0))

        # Build rectangular index
        row_index, unique_panels, unique_times, n_panels, n_times = _build_rect_index(
            panel_id, time_id
        )
        time_list = list(unique_times)

        # Resolve date values to column indices
        event_col = time_list.index(self.event_date)
        pre_start_col = time_list.index(self.pre_start)
        pre_end_col = time_list.index(self.pre_end)
        n_pre = pre_end_col - pre_start_col + 1

        if event_col <= pre_end_col:
            raise ValueError("event_date must be after pre_end")
        if pre_end_col <= pre_start_col:
            raise ValueError("pre_end must be after pre_start")
        if self.method == "gls" and self.npc > n_pre:
            raise ValueError("npc must be <= number of pre-event days")

        # For GLS, check that data extends far enough before pre_start
        pre_window_len = event_col - pre_start_col
        if self.method == "gls":
            gls_data_start = pre_start_col - pre_window_len
            if gls_data_start < 0:
                raise ValueError(
                    f"GLS requires {pre_window_len} observations before pre_start. "
                    f"Data only extends {pre_start_col} periods before pre_start."
                )

        # ---- param names ----
        param_names = self.indepvars + ["_cons"]
        n_params = len(param_names)

        # ---- helper to run one date ----
        def _run_date(target_col: int) -> tuple[np.ndarray, int]:
            """Estimate coefficients for a single date column."""
            # Identify valid panels for this date
            valid = row_index[:, target_col] >= 0  # panel has data on target date

            if self.method == "gls":
                # GLS pre-event window for this target
                offset = event_col - pre_start_col
                gls_pre_start = target_col - offset
                gls_pre_end = target_col - (event_col - pre_end_col)

                # Require non-missing y across the full GLS window + target
                cols_needed = list(range(gls_pre_start, gls_pre_end + 1)) + [target_col]
                for c in cols_needed:
                    has_data = row_index[:, c] >= 0
                    # Check y is not NaN where row_index is valid
                    y_ok = np.ones(n_panels, dtype=bool)
                    for i in range(n_panels):
                        if has_data[i]:
                            y_ok[i] = np.isfinite(y_all[row_index[i, c]])
                        else:
                            y_ok[i] = False
                    valid &= y_ok

                # Also check that returns are not all-zero across pre-event window
                for i in range(n_panels):
                    if valid[i]:
                        total = 0.0
                        for c in range(gls_pre_start, gls_pre_end + 1):
                            idx = row_index[i, c]
                            if idx >= 0:
                                total += abs(y_all[idx])
                        if total < 0.01:
                            valid[i] = False

            panel_mask = np.where(valid)[0]
            n_valid = len(panel_mask)
            if n_valid == 0:
                return np.full(n_params, np.nan), 0

            # Extract y and X for valid panels on target date
            rows = row_index[panel_mask, target_col]
            y = y_all[rows]
            if X_rhs.shape[1] > 0:
                X = np.column_stack([X_rhs[rows], np.ones(n_valid)])
            else:
                X = np.ones((n_valid, 1))

            if self.method == "gls":
                # Build pre-event return matrix (T_pre x N_valid) for PCA
                pre_cols = list(range(gls_pre_start, gls_pre_end + 1))
                R_pre = np.empty((len(pre_cols), n_valid))
                for j, pc in enumerate(pre_cols):
                    R_pre[j, :] = y_all[row_index[panel_mask, pc]]

                V, lam, d = _pca_decompose(R_pre, self.npc)

                if self.solver == "woodbury":
                    beta = _gls_woodbury(y, X, V, lam, d)
                else:
                    beta = _gls_cholesky(y, X, V, lam, d)
            else:
                beta = _ols(y, X)

            return beta, n_valid

        # ---- event date ----
        beta_event, n_event = _run_date(event_col)
        if n_event == 0:
            raise ValueError("No valid observations on event date after filtering.")

        # ---- pre-event dates ----
        all_betas = np.empty((n_pre + 1, n_params))
        all_nobs = np.empty(n_pre + 1, dtype=int)
        all_betas[0] = beta_event
        all_nobs[0] = n_event

        for i, pre_col in enumerate(range(pre_end_col, pre_start_col - 1, -1)):
            if verbose and self.method == "gls":
                pct = (i + 1) / n_pre * 100
                if pct % 10 < (i / n_pre * 100) % 10 or i == 0:
                    print(f"\rPercent complete = {pct:.0f}% ", end="", flush=True)

            b, n = _run_date(pre_col)
            all_betas[i + 1] = b
            all_nobs[i + 1] = n

        if verbose and self.method == "gls":
            print()

        # ---- significance stats ----
        pre_betas = all_betas[1:]
        pre_mean = np.nanmean(pre_betas, axis=0)
        pre_std = np.sqrt(np.nanvar(pre_betas, axis=0, ddof=1))

        # Empirical CDF p-value (Stata convention: weak inequality, denominator = L+1)
        L_plus_1 = all_betas.shape[0]
        p_cdf = np.sum(
            np.abs(all_betas - pre_mean) >= np.abs(beta_event - pre_mean), axis=0
        ) / L_plus_1

        # Parametric t-statistic and p-value
        from scipy import stats
        t_adj = L_plus_1 / (L_plus_1 - 1)
        z = np.abs(beta_event - pre_mean) / (pre_std * np.sqrt(t_adj))
        p_param = 2.0 * stats.t.sf(z, df=L_plus_1 - 2)

        method_label = "GLS" if self.method == "gls" else "OLS"
        if self.method == "gls" and self.solver == "woodbury":
            method_label += " (Woodbury)"

        return CSEStudyResult(
            method=method_label,
            params=beta_event,
            param_names=param_names,
            p_cdf=p_cdf,
            p_parametric=p_param,
            pre_event_betas=pre_betas,
            n_obs_event=n_event,
            n_obs_all=all_nobs,
            n_pre_event_days=n_pre,
        )
