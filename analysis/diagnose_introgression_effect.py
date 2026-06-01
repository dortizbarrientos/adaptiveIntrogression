#!/usr/bin/env python3
# =============================================================================
# diagnose_introgression_effect.py   (v2 -- effect size, not just detectability)
#
# Question: does the introgression treatment ("with" vs "without") move the
# outcome by an amount that MATTERS -- not merely by an amount we can detect?
#
# WHY v1 WAS WRONG. The first version asked "are the paired values exactly
# equal?" That distinguishes a DEAD code path (byte-identical) from a LIVE one,
# but it conflates "the toggle fires" with "the toggle matters." In a paired,
# shared-seed design, introducing introgressed alleles perturbs the RNG stream,
# so the arms diverge in the 5th or 6th decimal even when the biological effect
# is nil. v1 saw that divergence and declared "real signal." It was reporting
# DETECTABILITY (which the paired design makes trivial) instead of EFFECT SIZE
# (which is the scientific question). This version fixes that.
#
# THE TWO LENSES THAT SEPARATE JITTER FROM SIGNAL
#   1. MAGNITUDE relative to the scale of the phenomenon.
#      Standardised effect = |mean(with) - mean(without)| / SD(variable across
#      all runs). "How big is the treatment effect compared with how much this
#      quantity varies across the whole experiment?" A value of 1e-5 against a
#      spread of 0.3 is ~3e-5 SD -- negligible. We also report the largest
#      per-cell effect, so a real effect confined to a few cells (as in a
#      ramp-duration sweep) is not washed out by pooling.
#   2. DIRECTIONAL CONSISTENCY.
#      Of the pairs that differ, what fraction have with > without? Real rescue
#      is DIRECTIONAL (introgression helps -> with > without most of the time,
#      fraction near 1). Pure RNG jitter is SYMMETRIC (fraction near 0.5): the
#      arms wobble up and down at random and the mean difference cancels.
#
# A variable is called MATERIAL only if BOTH hold: the standardised effect
# clears a threshold AND the direction is consistent. "Fires at all" (exact
# equality) is still reported, but only as a SECONDARY note about the code path.
#
# This is the classic significance-vs-magnitude distinction made operational:
# a paired t-test here would be hugely "significant" and scientifically empty.
#
# Read-only. Consumes the tidy CSV from aggregate_pl_sensitivity.py.
# =============================================================================

from __future__ import annotations

import argparse
import sys
import numpy as np
import pandas as pd

RESPONSES = ["persistence", "fitness_lag", "fitness_end"]
PAIR_KEYS = ["pl", "mig", "rep", "lineage"]
CELL_KEYS = ["pl", "mig", "lineage"]   # a "cell" pools the 10 reps


def analyse(df, var, min_effect, dir_band):
    """Return a dict of effect-size diagnostics for one variable."""
    w  = df[df["script"] == "with"][PAIR_KEYS + [var]].rename(columns={var: "w"})
    wo = df[df["script"] == "without"][PAIR_KEYS + [var]].rename(columns={var: "wo"})
    m = w.merge(wo, on=PAIR_KEYS, how="inner")
    n = len(m)
    if n == 0:
        return {"var": var, "n_pairs": 0}

    diff = (m["w"] - m["wo"]).to_numpy(float)
    adiff = np.abs(diff)

    # scale of the phenomenon: SD of the variable across ALL runs (both arms)
    all_vals = np.concatenate([m["w"].to_numpy(float), m["wo"].to_numpy(float)])
    sd_total = all_vals.std(ddof=1)
    grand_mean = all_vals.mean()

    # (1a) pooled standardised effect
    mean_diff = diff.mean()
    std_eff_pooled = abs(mean_diff) / sd_total if sd_total > 0 else np.nan

    # (1b) largest per-cell standardised effect (guards against pooling washout)
    m = m.assign(diff=diff)
    cell_mean = m.groupby(CELL_KEYS, observed=True)["diff"].mean()
    std_eff_max_cell = (cell_mean.abs().max() / sd_total) if sd_total > 0 else np.nan

    # (2) directional consistency among pairs that actually differ
    nz = diff[adiff > 1e-9]
    frac_with_greater = float((nz > 0).mean()) if len(nz) else np.nan
    directional = (not np.isnan(frac_with_greater)) and \
                  (abs(frac_with_greater - 0.5) > dir_band)

    # secondary: code-path-alive check
    n_exact = int((adiff <= 1e-9).sum())

    material = (std_eff_max_cell >= min_effect) and directional

    return {
        "var": var, "n_pairs": n,
        "n_exact_equal": n_exact,
        "mean_diff": mean_diff, "max_abs_diff": float(adiff.max()),
        "sd_total": sd_total, "grand_mean": grand_mean,
        "pct_of_mean": (100*abs(mean_diff)/abs(grand_mean)) if grand_mean else np.nan,
        "std_eff_pooled": std_eff_pooled,
        "std_eff_max_cell": std_eff_max_cell,
        "frac_with_greater": frac_with_greater,
        "directional": directional,
        "material": material,
    }


def main():
    ap = argparse.ArgumentParser(
        description="Effect-size diagnostic for the introgression treatment.")
    ap.add_argument("--tidy", default="pl_sensitivity_tidy.csv")
    ap.add_argument("--min-effect", type=float, default=0.1,
                    help="Min standardised effect (fraction of a SD) to call MATERIAL. Default 0.1.")
    ap.add_argument("--dir-band", type=float, default=0.15,
                    help="Directionality must fall outside 0.5 +/- this band. Default 0.15 (=> outside 0.35-0.65).")
    args = ap.parse_args()

    df = pd.read_csv(args.tidy, dtype={"pl": str, "mig": str, "rep": str})

    print("=" * 76)
    print("INTROGRESSION-EFFECT DIAGNOSTIC  (effect size, not just detectability)")
    print("=" * 76)
    print(f"  tidy frame   : {args.tidy}   ({len(df)} rows)")
    print(f"  pairing on   : {PAIR_KEYS}")
    print(f"  MATERIAL if  : standardised effect >= {args.min_effect}"
          f"  AND  direction outside 0.5 \u00b1 {args.dir_band}\n")

    res = {}
    for var in RESPONSES:
        if var not in df.columns:
            print(f"  [skip] {var} absent"); continue
        r = analyse(df, var, args.min_effect, args.dir_band)
        res[var] = r
        if r["n_pairs"] == 0:
            print(f"  {var}: no pairs"); continue
        verdict = "MATERIAL" if r["material"] else "negligible"
        fg = r["frac_with_greater"]
        fg_str = "n/a" if np.isnan(fg) else f"{100*fg:.0f}% with>without"
        print(f"  {var:<12}: {verdict:<10} | "
              f"std effect (max cell) = {r['std_eff_max_cell']:.4f} | "
              f"direction = {fg_str}")
        print(f"  {'':12}  detail: mean(with-without)={r['mean_diff']:+.3g}, "
              f"max|diff|={r['max_abs_diff']:.3g}, "
              f"={r['pct_of_mean']:.3f}% of mean, "
              f"SD(var)={r['sd_total']:.3g}")
        if r["n_exact_equal"] == r["n_pairs"]:
            print(f"  {'':12}  [code path] byte-identical in all pairs -> toggle INERT")
        elif not r["material"]:
            print(f"  {'':12}  [code path] fires (arms differ) but effect is jitter-scale")

    # ---------- overall interpretation ----------
    print("\n" + "-" * 76)
    print("INTERPRETATION")
    print("-" * 76)

    def material(v):  return res.get(v, {}).get("material", False)
    def inert(v):
        r = res.get(v); return r and r["n_pairs"] and r["n_exact_equal"] == r["n_pairs"]

    any_material = any(material(v) for v in RESPONSES)
    all_inert    = all(inert(v) for v in RESPONSES if v in res)

    if all_inert:
        print("  Every variable is byte-identical between arms.")
        print("  -> DEAD TOGGLE: the 'with' script is not enabling introgression.")
        print("     Fix the experiment before any analysis.")
    elif not any_material:
        print("  The introgression code path FIRES (arms differ in low decimals),")
        print("  but no variable shows a material, directional effect: the")
        print("  differences are stochastic jitter, not rescue.")
        print("  -> LIVE BUT INERT. Consistent with the ramp outpacing selection")
        print("     (no time for allele-frequency change). This is a DESIGN limit,")
        print("     not a code bug. The lever is ramp duration / generation time,")
        print("     not the SLiM script. Consider a ramp-duration sweep to locate")
        print("     where introgression begins to matter.")
    elif material("persistence"):
        print("  persistence shows a material, directional effect -> real")
        print("  demographic rescue. Persistence is a valid lead figure.")
    elif material("fitness_lag") or material("fitness_end"):
        print("  persistence is negligible but FITNESS is materially affected")
        print("  -> introgression works without changing census size (soft")
        print("     selection / imposed N). Lead with a fitness variable.")
    else:
        print("  Mixed pattern -- read the per-variable numbers above.")

    # exit code: 0 if any material effect, 2 if none (inert), for scripting
    sys.exit(0 if any_material else 2)


if __name__ == "__main__":
    main()
