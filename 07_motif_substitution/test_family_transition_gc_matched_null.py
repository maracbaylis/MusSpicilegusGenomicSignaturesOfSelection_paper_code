#!/usr/bin/env python3
"""
Shared helpers for GC/base-change-matched motif permutation analyses.

These functions are taken from the historical family-transition analysis used
in the project and are imported by the motif-analysis scripts.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def bh_fdr(pvals):
    """Benjamini-Hochberg FDR adjustment preserving missing values."""
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(pvals.shape, np.nan)

    ok = np.isfinite(pvals)
    if ok.sum() == 0:
        return out

    p = pvals[ok]
    order = np.argsort(p)
    ranked = p[order]

    q = ranked * len(ranked) / (np.arange(len(ranked)) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]

    tmp = np.empty_like(q)
    tmp[order] = np.clip(q, 0, 1)
    out[ok] = tmp

    return out


def assign_strata(df):
    """
    Assign substitution sites to GC/base-change matched strata.

    Strata are defined by:
      - substitution base class;
      - quartile of local GC content in the ancestral/non-spicilegus sequence;
      - quartile of whole-cCRE GC content.

    Duplicate quantile boundaries are dropped by pandas.qcut.
    """
    out = df.copy()

    for col in ["local_gc_other", "cre_gc"]:
        out[f"{col}_bin"] = pd.qcut(
            out[col],
            q=4,
            duplicates="drop",
        ).astype(str)

    out["stratum"] = (
        out["base_class"].astype(str)
        + "|local="
        + out["local_gc_other_bin"].astype(str)
        + "|cre="
        + out["cre_gc_bin"].astype(str)
    )

    return out.reset_index(drop=True)
