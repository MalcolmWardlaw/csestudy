"""
csestudy — Time-series inference for cross-sectional event studies.

Implements the methodology of Cohn, Johnson, Liu, and Wardlaw (2026),
"Past is Prologue: Inference from the Cross Section of Returns Around
an Event," Journal of Financial Economics 180, 104278.

Usage:
    from csestudy import CSEventStudy
    result = CSEventStudy(df, event_date=0, pre_start=-200, pre_end=-1).fit()
    print(result.summary())
"""

from csestudy.core import CSEventStudy, CSEStudyResult

__all__ = ["CSEventStudy", "CSEStudyResult"]
__version__ = "0.1.0"
