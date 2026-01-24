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

DEFAULT_SIGNAL_CACHE = "cache/lifetimes-trackstats.pkl"
DEFAULT_BKG_CACHE    = "cache/mumu_bkg_stats_nominal_nobib.pkl"
DEFAULT_OUT = "/scratch/miralittmann/analysis/mira_analysis_code/sig_vs_bkg_tracks.pdf"

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

DEFAULT_BINS = {"pT": 30, "hits": 20, "velo": 40, "mass": 30, "beta": 40}
DEFAULT_XLIM = {
    "pT": (0, 10000),
    "velo": (200, 320),
    "mass": (0, 4500),
    "beta": (0.6, 1.05),
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

def get_beta(arr):
    return finite(arr) / C_MM_PER_NS

def fetch_arr(kind, *, sig=None, bkg=None, lifetime=None, sample=None,
              window=None, req=None, option=None, feature=None):
    """
    Unified array fetcher so we don't scatter feature-specific ifs everywhere.
    kind: 'sig' or 'bkg'
    """
    base = "velo" if feature == "beta" else feature

    if kind == "sig":
        arr = get_signal_arr(sig, lifetime, sample, req, base)
    elif kind == "bkg":
        arr = get_bkg_arr(bkg, window, req, option, base)
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

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--signal-cache", default=DEFAULT_SIGNAL_CACHE)
    parser.add_argument("--bkg-cache", default=DEFAULT_BKG_CACHE)
    parser.add_argument("--out", default=DEFAULT_OUT)

    parser.add_argument("--lifetime", type=int, default=DEFAULT_LIFETIME)
    parser.add_argument("--signal-samples", nargs="+", default=DEFAULT_SIGNAL_SAMPLES)

    parser.add_argument("--window", default=DEFAULT_WINDOW)
    parser.add_argument("--bkg-option", default=DEFAULT_BKG_OPTION)

    parser.add_argument("--features", nargs="+", default=DEFAULT_FEATURES)

    args = parser.parse_args()

    sig = load_pickle(args.signal_cache)
    bkg = load_pickle(args.bkg_cache)

    lifetime = args.lifetime
    samples  = args.signal_samples
    window   = args.window
    option   = args.bkg_option

    print("\n=== Track requirement efficiencies (relative to VB) ===")

    bkg_counts = {req: len(fetch_arr("bkg", bkg=bkg, window=window, req=req, option=option, feature="hits"))
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

    for s in samples:
        vb_sig = sig_counts[s]["vb"]
        for req in TRACK_REQS:
            if vb_sig > 0:
                pct = 100.0 * sig_counts[s][req] / vb_sig
                print(f"SIG {s} {req_to_name[req]:>18}: {pct:6.2f}% ({sig_counts[s][req]}/{vb_sig})")
            else:
                print(f"SIG {s}: VB count is 0; cannot compute efficiencies.") 
        print(f"\n")
    
    plt.style.use("seaborn-v0_8-colorblind")

    if lifetime not in sig:
        raise KeyError(f"Lifetime {lifetime} not found in signal cache. Available: {list(sig.keys())}")

    with PdfPages(args.out) as pdf:
        for feature in args.features:
            fig, axes = plt.subplots(1, 3, figsize=(14, 4.5), sharey=True)
            fig.suptitle(
                f"{feature}  |  signal: τ={lifetime} ns  |  bkg: {window} window, {bib_to_name[option]}",
                fontsize=16
            )

            bins = DEFAULT_BINS.get(feature, 30)
            xlim = DEFAULT_XLIM.get(feature, None)

            for ax, req in zip(axes, TRACK_REQS):
                bkg_arr = fetch_arr("bkg", bkg=bkg, window=window, req=req, option=option, feature=feature)

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
                        ax, bkg_arr, bins=hit_bins,
                        color="gray", label=None,
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

                    for s in samples:
                        sig_arr = fetch_arr("sig", sig=sig, lifetime=lifetime, sample=s, req=req, feature=feature)
                        hist_percent(
                            ax, sig_arr, bins=bins, xlim=xlim,
                            label=f"Signal {s}",
                            fill=False, linewidth=2.0
                        )

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

if __name__ == "__main__":
    main()
