import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from math import *
import pathlib
import pickle
import argparse
from tqdm import tqdm
import math
from scipy import optimize
from collections import OrderedDict
import pandas as pd

import pyLCIO
from pyLCIO import UTIL, EVENT
import ROOT



bib_color = "lightsteelblue"
DEFAULT_SIGNAL_CACHE = "cache/lifetimes-trackstats.pkl"
DEFAULT_BKG_CACHE    = "cache/mumu_bkg_stats_nominal_nobib.pkl"
DEFAULT_BIB_CACHE    = "/scratch/wandriscok/kate_mucoll_scripts/bib_analysis/cache/nu_bkg_stats_loose.pkl"
DEFAULT_OUT = "/scratch/miralittmann/analysis/mira_analysis_code/sig_vs_bkg_vs_bib_tracks.pdf"

DEFAULT_LIFETIME = 30
DEFAULT_SIGNAL_SAMPLES = ["1500", "2500", "3500"]
DEFAULT_WINDOW = "nominal"
DEFAULT_BKG_OPTION = "nobib"

TRACK_REQS = ["vb", "ib", "ob"]
req_to_name = {
    "vb": "≥3 VB hits",
    "ib": "≥3 VB, ≥2 IB hits",
    "ob": "≥3 VB, ≥2 IB, ≥2 OB hits"
}
bib_to_name = {
    "bib": "10% BIB",
    "nobib": "No BIB"
}

DEFAULT_FEATURES = ["pT", "hits", "velo", "mass", "beta"]

C_MM_PER_NS = 299.792458

LABELS = {
    "pT":   r"$p_T$ [GeV]",
    "hits": "Hits on track",
    "velo": "Velocity [mm/ns]",
    "mass": r"Reconstructed mass [GeV]",
    "beta": r"$\beta$ = reco v / c",
}

DEFAULT_BINS = {"pT": 30, "hits": 20, "velo": 40, "mass": 30, "beta": 60}
DEFAULT_XLIM = {
    "pT": (0, 10000),
    "velo": (0, 310),
    "mass": (0, 4500),
    "beta": (0.6, 1.01),
}
CUTS = {
    "pT": (0, 800),
    "mass": (0,500),
    "beta": (0.99, 1.05)
}

def load_pickle(path):
    path = pathlib.Path(path)
    with path.open("rb") as f:
        return pickle.load(f)

def finite(arr):
    arr = np.asarray(arr, dtype=float)
    return arr[np.isfinite(arr)]

def weights_percent(n):
    if n <= 0:
        return None
    return np.full(n, 100.0 / n, dtype=float)

def get_signal_arr(sig, lifetime, sample, req, feature):
    return finite(sig[lifetime][sample][req].get(feature, []))

def get_bkg_arr(bkg, window, req, option, feature):
    return finite(bkg[window][req][option].get(feature, []))

def get_bib_arr(bib, window, req, option, feature):
    return finite(bib[window][req][option].get(feature, []))

def get_beta(arr):
    return finite(arr) / C_MM_PER_NS

def fetch_arr(kind, *, sig=None, bkg=None, bib=None, lifetime=None, sample=None,
              window=None, req=None, option=None, feature=None):
    base = "velo" if feature == "beta" else feature

    if kind == "bib":
        if base == "velo":
            base = "velocity"
        if feature == "beta":
            base = "velocity"

    if kind == "sig":
        arr = get_signal_arr(sig, lifetime, sample, req, base)
    elif kind == "bkg":
        arr = get_bkg_arr(bkg, window, req, option, base)
    elif kind == "bib":
        arr = get_bib_arr(bib, "loose", req, "10_bib", base)

    else:
        raise ValueError("kind must be 'sig' or 'bkg'")

    if feature == "beta":
        arr = get_beta(arr)

    return arr

def hist_percent(ax, x, *, bins, xlim, color=None, label=None,
                 histtype="step", linewidth=2.0, alpha=1.0, fill=False):
    if x.size == 0:
        return
    w = weights_percent(x.size)
    if fill:
        ax.hist(
            x, bins=bins, range=xlim,
            weights=w, histtype="stepfilled",
            color=color, alpha=alpha, label=label, linewidth=0.0
        )
    ax.hist(
        x, bins=bins, range=xlim,
        weights=w, histtype=histtype,
        color=color, alpha=alpha, label=None if fill else label,
        linewidth=linewidth
    )

def integer_bins_from_data(*arrays, pad_low=0, pad_high=0):
    xs = [np.asarray(a) for a in arrays if a is not None and len(a) > 0]
    if not xs:
        return None
    x = np.concatenate(xs)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return None

    xmin = int(np.floor(x.min())) - pad_low
    xmax = int(np.ceil(x.max())) + pad_high
    return np.arange(xmin - 0.5, xmax + 1.5, 1.0)

def hist_percent_bins(ax, x, *, bins, color=None, label=None,
                      histtype="step", linewidth=2.0, alpha=1.0, fill=False):
    if x.size == 0 or bins is None:
        return
    w = weights_percent(x.size)
    if fill:
        ax.hist(
            x, bins=bins,
            weights=w, histtype="stepfilled",
            color=color, alpha=alpha, label=label, linewidth=0.0
        )
    ax.hist(
        x, bins=bins,
        weights=w, histtype=histtype,
        color=color, alpha=alpha, label=None if fill else label,
        linewidth=linewidth
    )

def to_float_array(x):
    return np.asarray(x, dtype=float)

def get_signal_raw(sig, lifetime, sample, req, feature):
    return to_float_array(sig[lifetime][sample][req].get(feature, []))

def get_bkg_raw(bkg, window, req, option, feature):
    return to_float_array(bkg[window][req][option].get(feature, []))

def get_bib_raw(bib, window, req, option, feature):
    return to_float_array(bib[window][req][option].get(feature, []))

def raw_feature(kind, *, sig=None, bkg=None, bib=None,
                lifetime=None, sample=None, window=None, req=None, option=None,
                feature=None):
    base = "velo" if feature == "beta" else feature

    if kind == "bib":
        if base in ("velo", "velocity"):
            base = "velocity"

    if kind == "sig":
        arr = get_signal_raw(sig, lifetime, sample, req, base)
    elif kind == "bkg":
        arr = get_bkg_raw(bkg, window, req, option, base)
    elif kind == "bib":
        arr = get_bib_raw(bib, "loose", req, "bib", base)
    else:
        raise ValueError("kind must be 'sig', 'bkg', or 'bib'")

    if feature == "beta":
        arr = arr / C_MM_PER_NS

    return arr

def build_cut_mask(kind, *, sig=None, bkg=None, bib=None,
                   lifetime=None, sample=None, window=None, req=None, option=None,
                   cuts=None):
    if not cuts:
        return None

    first_feat = next(iter(cuts.keys()))
    x0 = raw_feature(kind, sig=sig, bkg=bkg, bib=bib,
                     lifetime=lifetime, sample=sample, window=window, req=req, option=option,
                     feature=first_feat)
    n = len(x0)
    keep = np.ones(n, dtype=bool)

    for feat, (lo, hi) in cuts.items():
        x = raw_feature(kind, sig=sig, bkg=bkg, bib=bib,
                        lifetime=lifetime, sample=sample, window=window, req=req, option=option,
                        feature=feat)

        if len(x) != n:
            print(f"[WARN] Cannot apply track-level cuts for {kind}/{req}: feature '{feat}' has length {len(x)} != {n}")
            return None

        finite_x = np.isfinite(x)
        in_window = finite_x & (x >= lo) & (x <= hi)
        keep &= ~in_window

    return keep

def fetch_arr_with_optional_cuts(kind, *, sig=None, bkg=None, bib=None,
                                lifetime=None, sample=None, window=None, req=None, option=None,
                                feature=None, apply_cuts=False, cuts=None):
    x = raw_feature(kind, sig=sig, bkg=bkg, bib=bib,
                    lifetime=lifetime, sample=sample, window=window, req=req, option=option,
                    feature=feature)

    if apply_cuts and cuts:
        keep = build_cut_mask(kind, sig=sig, bkg=bkg, bib=bib,
                              lifetime=lifetime, sample=sample, window=window, req=req, option=option,
                              cuts=cuts)
        if keep is not None:
            x = x[keep]

    return finite(x)

def draw_cut_lines(ax, feature, cuts):
    if feature not in cuts:
        return
    lo, hi = cuts[feature]
    ax.axvline(lo, color="red", linestyle="--", linewidth=1.5, alpha=0.8)
    ax.axvline(hi, color="red", linestyle="--", linewidth=1.5, alpha=0.8)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--signal-cache", default=DEFAULT_SIGNAL_CACHE)
    parser.add_argument("--bkg-cache", default=DEFAULT_BKG_CACHE)
    parser.add_argument("--bib-cache", default=DEFAULT_BIB_CACHE)
    parser.add_argument("--out", default=DEFAULT_OUT)

    parser.add_argument("--lifetime", type=int, default=DEFAULT_LIFETIME)
    parser.add_argument("--signal-samples", nargs="+", default=DEFAULT_SIGNAL_SAMPLES)

    parser.add_argument("--window", default=DEFAULT_WINDOW)
    parser.add_argument("--bkg-option", default=DEFAULT_BKG_OPTION)

    parser.add_argument("--features", nargs="+", default=DEFAULT_FEATURES)
    parser.add_argument("--apply-cuts", action="store_true", help="Also produce a second PDF with CUTS applied to OB tracks only.")

    parser.add_argument("--only-ob", action="store_true",
                    help="Plot only the OB track requirement (single panel per feature).")
    parser.add_argument("--no-cut-lines", action="store_true",
                        help="Do not draw vertical cut lines.")

    args = parser.parse_args()
    out_cuts = None
    if args.apply_cuts:
        base, ext = os.path.splitext(args.out)
        out_cuts = f"{base}_cuts{ext}"

    sig = load_pickle(args.signal_cache)
    bkg = load_pickle(args.bkg_cache)
    bib = load_pickle(args.bib_cache)

    lifetime = args.lifetime
    samples  = args.signal_samples
    window   = args.window
    option   = args.bkg_option

    print("\n=== Track requirement efficiencies (relative to VB) ===")
    bkg_counts = {req: len(fetch_arr("bkg", bkg=bkg, window=window, req=req, option=option, feature="hits"))
                  for req in TRACK_REQS}
    bib_counts = {req: len(fetch_arr("bib", bib=bib, window=window, req=req, option=option, feature="hits"))
                  for req in TRACK_REQS}

    sig_counts = {
        s: {req: len(fetch_arr("sig", sig=sig, lifetime=lifetime, sample=s, req=req, feature="hits"))
            for req in TRACK_REQS}
        for s in samples
    }


    vb_bkg = bkg_counts["vb"]
    for req in TRACK_REQS:
        if vb_bkg > 0:
            pct = 100.0 * bkg_counts[req] / vb_bkg
            print(f"BKG {req_to_name[req]:>22}: {pct:6.2f}% ({bkg_counts[req]}/{vb_bkg})")
        else:
            print("BKG: VB count is 0; cannot compute efficiencies.")
    print(f"\n")

    vb_bib = bib_counts["vb"]
    for req in TRACK_REQS:
        if vb_bib > 0:
            pct = 100.0 * bib_counts[req] / vb_bib
            print(f"BIB {req_to_name[req]:>22}: {pct:6.2f}% ({bib_counts[req]}/{vb_bib})")
        else:
            print("BIB: VB count is 0; cannot compute efficiencies.")
    print(f"\n")


    for s in samples:
        vb_sig = sig_counts[s]["vb"]
        for req in TRACK_REQS:
            if vb_sig > 0:
                pct = 100.0 * sig_counts[s][req] / vb_sig
                print(f"SIG {s} {req_to_name[req]:>18}: {pct:6.2f}% ({sig_counts[s][req]}/{vb_sig})")
            else:
                print(f"SIG {s}: VB count is 0; cannot compute efficiencies.") 
        print(f"\n")


    print("\n=== OB efficiencies AFTER CUTS (denominator = VB uncut) ===")
    vb_bkg_den = bkg_counts["vb"]
    vb_bib_den = bib_counts["vb"]

    bkg_ob_cut = len(fetch_arr_with_optional_cuts(
        "bkg", bkg=bkg, window=window, req="ob", option=option,
        feature="hits", apply_cuts=True, cuts=CUTS
    ))
    bib_ob_cut = len(fetch_arr_with_optional_cuts(
        "bib", bib=bib, window=window, req="ob", option=option,
        feature="hits", apply_cuts=True, cuts=CUTS
    ))

    if vb_bkg_den > 0:
        print(f"BKG OB after cuts: {100.0*bkg_ob_cut/vb_bkg_den:6.2f}% ({bkg_ob_cut}/{vb_bkg_den})")
    else:
        print("BKG: VB count is 0; cannot compute cut efficiency.")

    if vb_bib_den > 0:
        print(f"BIB OB after cuts: {100.0*bib_ob_cut/vb_bib_den:6.2f}% ({bib_ob_cut}/{vb_bib_den})")
    else:
        print("BIB: VB count is 0; cannot compute cut efficiency.")

    for s in samples:
        vb_sig_den = sig_counts[s]["vb"]
        sig_ob_cut = len(fetch_arr_with_optional_cuts(
            "sig", sig=sig, lifetime=lifetime, sample=s, req="ob",
            feature="hits", apply_cuts=True, cuts=CUTS
        ))
        if vb_sig_den > 0:
            print(f"SIG {s} OB after cuts: {100.0*sig_ob_cut/vb_sig_den:6.2f}% ({sig_ob_cut}/{vb_sig_den})")
        else:
            print(f"SIG {s}: VB count is 0; cannot compute cut efficiency.")


    
    plt.style.use("seaborn-v0_8-colorblind")

    if lifetime not in sig:
        raise KeyError(f"Lifetime {lifetime} not found in signal cache. Available: {list(sig.keys())}")

    with PdfPages(args.out) as pdf:
        for feature in args.features:
            fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
            fig.suptitle(
                f"{feature}  |  signal: τ={lifetime} ns  |  bkg: {window} window, {bib_to_name[option]}  |  bib: loose window, 10%",
                fontsize=16
            )

            bins = DEFAULT_BINS.get(feature, 30)
            xlim = DEFAULT_XLIM.get(feature, None)

            for ax, req in zip(axes, TRACK_REQS):
                bkg_arr = fetch_arr("bkg", bkg=bkg, window=window, req=req, option=option, feature=feature)
                bib_arr = fetch_arr("bib", bib=bib, window=window, req=req, option=option, feature=feature)

                if feature == "hits":
                    sig_arrs = [
                        fetch_arr("sig", sig=sig, lifetime=lifetime, sample=s, req=req, feature=feature)
                        for s in samples
                    ]
                    hit_bins = integer_bins_from_data(bkg_arr, *sig_arrs, pad_low=0, pad_high=0)
                    hist_percent_bins(
                        ax, bkg_arr, bins=hit_bins,
                        color="gray", label="Mu-mu background",
                        fill=True, alpha=0.30, linewidth=0.0
                    )
                    hist_percent_bins(
                        ax, bib_arr, bins=hit_bins,
                        color=bib_color, label="BIB background",
                        fill=True, alpha=0.30, linewidth=0.0
                    )
                    hist_percent_bins(
                        ax, bkg_arr, bins=hit_bins,
                        color="gray", label=None,
                        fill=False, linewidth=1.5
                    )
                    hist_percent_bins(
                        ax, bib_arr, bins=hit_bins,
                        color=bib_color, label=None,
                        fill=False, linewidth=1.5
                    )
                    for s, sig_arr in zip(samples, sig_arrs):
                        hist_percent_bins(
                            ax, sig_arr, bins=hit_bins,
                            label=f"Signal {s}",
                            fill=False, linewidth=2.0
                        )

                    if hit_bins is not None:
                        lo = int(np.floor(hit_bins[0] + 0.5))
                        hi = int(np.ceil(hit_bins[-1] - 0.5))
                        ax.set_xticks(np.arange(lo, hi + 1, 1))

                else:
                    hist_percent(
                        ax, bkg_arr, bins=bins, xlim=xlim,
                        color="gray", label="Mu-mu background",
                        fill=True, alpha=0.30, linewidth=0.0
                    )
                    hist_percent(
                        ax, bkg_arr, bins=bins, xlim=xlim,
                        color="gray", label=None,
                        fill=False, linewidth=1.5
                    )
                    hist_percent(
                        ax, bib_arr, bins=bins, xlim=xlim,
                        color=bib_color, label="BIB background",
                        fill=True, alpha=0.30, linewidth=0.0
                    )
                    hist_percent(
                        ax, bib_arr, bins=bins, xlim=xlim,
                        color=bib_color, label=None,
                        fill=False, linewidth=1.5
                    )


                    for s in samples:
                        sig_arr = fetch_arr("sig", sig=sig, lifetime=lifetime, sample=s, req=req, feature=feature)
                        hist_percent(
                            ax, sig_arr, bins=bins, xlim=xlim,
                            label=f"Signal {s}",
                            fill=False, linewidth=2.0
                        )
                if not args.no_cut_lines: 
                    draw_cut_lines(ax, feature, CUTS)

                ax.set_title(req_to_name[req], fontsize=14)
                ax.set_xlabel(LABELS.get(feature, feature), fontsize=14)
                ax.set_ylabel("Normalized counts (%)" if req == "vb" else "", fontsize=14)
                if xlim is not None:
                    ax.set_xlim(*xlim)
                ax.grid(True, alpha=0.2)
                ax.tick_params(axis="both", which="major", labelsize=11)

                ax.text(
                    0.02, 0.98, "Muon Collider", transform=ax.transAxes,
                    ha="left", va="top", fontsize=12, fontweight="bold", style="italic"
                )
                ax.text(
                    0.02, 0.90, "MuColl_v1", transform=ax.transAxes,
                    ha="left", va="top", fontsize=10
                )

                ax.legend(fontsize=9, frameon=False, loc="upper right")

            fig.tight_layout(rect=[0, 0, 1, 0.92])
            pdf.savefig(fig)
            plt.close(fig)

    print(f"Wrote: {args.out}")

    base, ext = os.path.splitext(args.out)
    out_ob = f"{base}_ob_only{ext}"

    with PdfPages(out_ob) as pdf:
        req = "ob"
        for feature in args.features:
            fig, ax = plt.subplots(1, 1, figsize=(6.0, 4.5), sharey=True)

            bins = DEFAULT_BINS.get(feature, 30)
            xlim = DEFAULT_XLIM.get(feature, None)

            bkg_arr = fetch_arr("bkg", bkg=bkg, window=window, req=req, option=option, feature=feature)
            bib_arr = fetch_arr("bib", bib=bib, window=window, req=req, option=option, feature=feature)

            if feature == "hits":
                sig_arrs = [
                    fetch_arr("sig", sig=sig, lifetime=lifetime, sample=s, req=req, feature=feature)
                    for s in samples
                ]
                hit_bins = integer_bins_from_data(bkg_arr, *sig_arrs, pad_low=0, pad_high=0)

                hist_percent_bins(
                    ax, bkg_arr, bins=hit_bins,
                    color="gray", label="Nominal mu-mu bkg",
                    fill=True, alpha=0.30, linewidth=0.0
                )
                hist_percent_bins(
                    ax, bib_arr, bins=hit_bins,
                    color=bib_color, label=r"10% loose BIB bkg",
                    fill=True, alpha=0.30, linewidth=0.0
                )
                hist_percent_bins(
                    ax, bkg_arr, bins=hit_bins,
                    color="gray", label=None,
                    fill=False, linewidth=1.5
                )
                hist_percent_bins(
                    ax, bib_arr, bins=hit_bins,
                    color=bib_color, label=None,
                    fill=False, linewidth=1.5
                )

                for s, sig_arr in zip(samples, sig_arrs):
                    hist_percent_bins(
                        ax, sig_arr, bins=hit_bins,
                        label=f"Signal {s}, τ={lifetime} ns",
                        fill=False, linewidth=2.0
                    )

                if hit_bins is not None:
                    lo = int(np.floor(hit_bins[0] + 0.5))
                    hi = int(np.ceil(hit_bins[-1] - 0.5))
                    ax.set_xticks(np.arange(lo, hi + 1, 1))

            else:
                hist_percent(
                    ax, bkg_arr, bins=bins, xlim=xlim,
                    color="gray", label="Nominal mu-mu bkg",
                    fill=True, alpha=0.30, linewidth=0.0
                )
                hist_percent(
                    ax, bkg_arr, bins=bins, xlim=xlim,
                    color="gray", label=None,
                    fill=False, linewidth=1.5
                )
                hist_percent(
                    ax, bib_arr, bins=bins, xlim=xlim,
                    color=bib_color, label=r"10% loose BIB bkg",
                    fill=True, alpha=0.30, linewidth=0.0
                )
                hist_percent(
                    ax, bib_arr, bins=bins, xlim=xlim,
                    color=bib_color, label=None,
                    fill=False, linewidth=1.5
                )

                for s in samples:
                    sig_arr = fetch_arr("sig", sig=sig, lifetime=lifetime, sample=s, req=req, feature=feature)
                    hist_percent(
                        ax, sig_arr, bins=bins, xlim=xlim,
                        label=f"Signal {s}, τ={lifetime} ns",
                        fill=False, linewidth=2.0
                    )

            # ax.set_title(req_to_name[req], fontsize=14)
            ax.set_xlabel(LABELS.get(feature, feature), fontsize=14)
            ax.set_ylabel("Normalized counts (%)", fontsize=14)

            if xlim is not None:
                ax.set_xlim(*xlim)

            ax.grid(True, alpha=0.2)
            ax.tick_params(axis="both", which="major", labelsize=11)

            ax.text(
                0.02, 0.98, "Muon Collider", transform=ax.transAxes,
                ha="left", va="top", fontsize=12, fontweight="bold", style="italic"
            )
            ax.text(
                0.02, 0.92, "MuColl_v1", transform=ax.transAxes,
                ha="left", va="top", fontsize=10
            )
            ax.text(
                0.02, 0.86, r"$\geq$3 VB,$\geq$2 IB,$\geq$2 OB", transform=ax.transAxes,
                ha="left", va="top", fontsize=10
            )

            ax.legend(fontsize=9, frameon=False, loc="upper right")

            fig.tight_layout(rect=[0, 0, 1, 0.92])
            pdf.savefig(fig)
            plt.close(fig)

    print(f"Wrote: {out_ob}")

if __name__ == "__main__":
    main()
