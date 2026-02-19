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

sig_cache_path = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_by_event.pkl"
muon_cache_path = "/scratch/miralittmann/analysis/mira_analysis_code/cache/mumu_bkg_stats_nominal_nobib_byevent.pkl"
# coming from these cache path is event-level info with only the pT and hits cuts already applied.

m_cuts = np.arange(0, 5000 + 50, 50)
beta_cuts = np.arange(0.9, 1.01, 0.001)

signal_samples = ["1000", "1500", "2000", "2500", "3000", "3500", "4000", "4500"]
plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/scan_cuts.pdf"

XS_PB = {
    "mumu_bkg": 0.4312,
    "1000": 0.0005108,
    "1500": 0.0004715,
    "2000": 0.0004184,
    "2500": 0.0003528,
    "3000": 0.0002780,
    "3500": 0.0001980,
    "4000": 0.0001173,
    "4500": 0.0000450,
}

L_ABINV = 10.0
L_PBINV = L_ABINV * 1e6

def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)

signal = load_pickle(sig_cache_path)
muons = load_pickle(muon_cache_path)


def eff_vs_mass_threshold(masses_lead, masses_sub, denom_mask, m_cuts, mode="ge1"):
    """
    masses_lead, masses_sub: per-event arrays (NaN if no candidate)
    denom_mask: boolean array selecting denominator events
    m_cuts: array of thresholds
    mode: "ge1" or "ge2" track requirements per event
    Returns: efficiencies array same length as m_cuts
    """
    lead = np.asarray(masses_lead, dtype=float)
    sub  = np.asarray(masses_sub, dtype=float) if masses_sub is not None else None
    denom_mask = np.asarray(denom_mask, dtype=bool)

    # restrict to denominator events
    lead = lead[denom_mask]
    if sub is not None:
        sub = sub[denom_mask]

    N = lead.size
    if N == 0:
        return np.full_like(m_cuts, np.nan, dtype=float)

    effs = np.zeros_like(m_cuts, dtype=float)

    for i, mc in enumerate(m_cuts):
        if mode == "ge1":
            pass_evt = np.isfinite(lead) & (lead > mc)
        elif mode == "ge2":
            pass_evt = (
                np.isfinite(lead) & (lead > mc) &
                np.isfinite(sub)  & (sub  > mc)
            )
        else:
            raise ValueError("mode must be 'ge1' or 'ge2'")
        effs[i] = pass_evt.mean()

    return effs

def eff_vs_beta_threshold(betas_lead, betas_sub, denom_mask, beta_cuts, mode="ge1"):
    """
    betas_lead, betas_sub: per-event arrays (NaN if no candidate)
    denom_mask: boolean array selecting denominator events
    beta_cuts: array of thresholds (pass if beta < beta_cut)
    mode: "ge1" or "ge2" track requirements per event
    Returns: efficiencies array same length as beta_cuts
    """
    lead = np.asarray(betas_lead, dtype=float)
    sub  = np.asarray(betas_sub, dtype=float) if betas_sub is not None else None
    denom_mask = np.asarray(denom_mask, dtype=bool)

    # align lengths safely (avoids the mismatch you just hit)
    n = min(len(lead), len(denom_mask))
    lead = lead[:n]
    denom_mask = denom_mask[:n]
    if sub is not None:
        sub = sub[:n]

    lead = lead[denom_mask]
    if sub is not None:
        sub = sub[denom_mask]

    if lead.size == 0:
        return np.full_like(beta_cuts, np.nan, dtype=float)

    effs = np.zeros_like(beta_cuts, dtype=float)

    for i, bc in enumerate(beta_cuts):
        if mode == "ge1":
            pass_evt = np.isfinite(lead) & (lead < bc)
        elif mode == "ge2":
            pass_evt = (
                np.isfinite(lead) & (lead < bc) &
                np.isfinite(sub)  & (sub  < bc)
            )
        else:
            raise ValueError("mode must be 'ge1' or 'ge2'")
        effs[i] = pass_evt.mean()

    return effs

def style_ticks(ax, major=18, minor=16):
    ax.tick_params(axis="both", which="major", labelsize=major, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=minor, length=4, width=1.0)

def eff_grid_mass_beta(mass_lead, mass_sub, beta_lead, beta_sub,
                       denom_mask, m_cuts, beta_cuts, mode="ge1"):
    """
    Returns a 2D grid eff[j,i] for beta_cuts[j], m_cuts[i].
    Pass condition per event:
      mass > m_cut AND beta < beta_cut
    mode:
      ge1: use leading only
      ge2: require leading AND subleading both pass
    """
    mL = np.asarray(mass_lead, dtype=float)
    bL = np.asarray(beta_lead, dtype=float)
    mS = np.asarray(mass_sub, dtype=float) if mass_sub is not None else None
    bS = np.asarray(beta_sub, dtype=float) if beta_sub is not None else None
    denom_mask = np.asarray(denom_mask, dtype=bool)

    n = min(len(mL), len(bL), len(denom_mask))
    mL, bL, denom_mask = mL[:n], bL[:n], denom_mask[:n]
    if mS is not None:
        n2 = min(len(mS), n)
        mS = mS[:n2]
        bS = bS[:n2]
        mL, bL, denom_mask = mL[:n2], bL[:n2], denom_mask[:n2]

    mL = mL[denom_mask]
    bL = bL[denom_mask]
    if mS is not None:
        mS = mS[denom_mask]
        bS = bS[denom_mask]

    N = mL.size
    if N == 0:
        return np.full((len(beta_cuts), len(m_cuts)), np.nan, dtype=float)

    eff = np.zeros((len(beta_cuts), len(m_cuts)), dtype=float)

    finL = np.isfinite(mL) & np.isfinite(bL)
    if mode == "ge2":
        if mS is None or bS is None:
            return np.full((len(beta_cuts), len(m_cuts)), np.nan, dtype=float)
        finS = np.isfinite(mS) & np.isfinite(bS)

    for j, bc in enumerate(beta_cuts):
        for i, mc in enumerate(m_cuts):
            if mode == "ge1":
                pass_evt = finL & (mL > mc) & (bL < bc)
            elif mode == "ge2":
                pass_evt = (
                    finL & finS &
                    (mL > mc) & (bL < bc) &
                    (mS > mc) & (bS < bc)
                )
            else:
                raise ValueError("mode must be 'ge1' or 'ge2'")
            eff[j, i] = pass_evt.mean()

    return eff

def plot_eff_heatmap(pdf, X, Y, Z, title, xlabel, ylabel):
    fig, ax = plt.subplots(figsize=(9, 7))

    x = np.asarray(X, dtype=float)
    y = np.asarray(Y, dtype=float)

    x_edges = np.r_[x[0] - (x[1]-x[0])/2, (x[:-1] + x[1:])/2, x[-1] + (x[-1]-x[-2])/2] if len(x) > 1 else np.array([x[0]-0.5, x[0]+0.5])
    y_edges = np.r_[y[0] - (y[1]-y[0])/2, (y[:-1] + y[1:])/2, y[-1] + (y[-1]-y[-2])/2] if len(y) > 1 else np.array([y[0]-0.01, y[0]+0.01])

    im = ax.pcolormesh(x_edges, y_edges, Z, shading="auto")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label("Event efficiency", fontsize=18)

    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_title(title, fontsize=16)

    style_ticks(ax, major=18, minor=16)
    ax.grid(False)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def plot_roc_from_grids(pdf, Zb, Zs, title):
    xb = np.asarray(Zb, dtype=float).ravel()
    ys = np.asarray(Zs, dtype=float).ravel()

    m = np.isfinite(xb) & np.isfinite(ys)
    xb = xb[m]
    ys = ys[m]

    fig, ax = plt.subplots(figsize=(8.5, 7))
    ax.scatter(ys, xb, s=8)

    ax.set_xlabel(r"Expected signal events (10 ab$^{-1}$)", fontsize=20)
    ax.set_ylabel(r"Expected background events (10 ab$^{-1}$)", fontsize=20)
    ax.set_title(title, fontsize=16)
    ax.set_yscale("log")

    ax.grid(True, alpha=0.2)
    style_ticks(ax, major=18, minor=16)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def pick_best_cut_simple(Z_sig, Z_bkg, m_cuts, beta_cuts, finite_only=True):
    """
    Pick the best (m_cut, beta_cut) using:
      1) lowest background-efficiency bin (min Z_bkg)
      2) within that bin, farthest right (max signal efficiency Z_sig)

    Z_sig, Z_bkg: 2D arrays with shape (len(beta_cuts), len(m_cuts))
                  where Z[j,i] corresponds to beta_cuts[j], m_cuts[i].
    Returns dict with chosen cuts and efficiencies + indices.
    """
    Zs = np.asarray(Z_sig, dtype=float)
    Zb = np.asarray(Z_bkg, dtype=float)
    assert Zs.shape == Zb.shape, "Z_sig and Z_bkg must have same shape"
    assert Zs.shape == (len(beta_cuts), len(m_cuts)), "Grid shape mismatch"

    if finite_only:
        ok = np.isfinite(Zs) & np.isfinite(Zb)
    else:
        ok = np.ones_like(Zs, dtype=bool)

    if not np.any(ok):
        return None

    bmin = np.nanmin(np.where(ok, Zb, np.nan))

    tol = 1e-12
    best_bkg_mask = ok & (np.abs(Zb - bmin) <= tol)

    if not np.any(best_bkg_mask):
        best_bkg_mask = ok & (Zb <= bmin + 1e-9)

    Zs_in = np.where(best_bkg_mask, Zs, -np.inf)
    flat_idx = np.argmax(Zs_in)
    j, i = np.unravel_index(flat_idx, Zs.shape)

    return {
        "m_cut": float(m_cuts[i]),
        "beta_cut": float(beta_cuts[j]),
        "sig_eff": float(Zs[j, i]),
        "bkg_eff": float(Zb[j, i]),
        "i_m": int(i),
        "j_beta": int(j),
        "bkg_min": float(bmin),
    }

def best_cuts_per_masspoint(signal, muons, signal_samples, m_cuts, beta_cuts,
                           window="nominal", req="ob", option="nobib",
                           lifetime=30, mode="ge1"):
    """
    Convenience wrapper:
      - builds Z_bkg once
      - loops over signal mass points (samples)
      - returns a dict: sample -> best cut info

    Requires your existing eff_grid_mass_beta(...) function.
    """
    d_bkg = muons[window][req][option]
    lead_m_b = d_bkg.get("leading_mass", [])
    sub_m_b  = d_bkg.get("subleading_mass", [])
    lead_b_b = d_bkg.get("leading_beta", [])
    sub_b_b  = d_bkg.get("subleading_beta", [])
    denom_bkg = np.ones(len(lead_m_b), dtype=bool)

    Zb = eff_grid_mass_beta(
        lead_m_b, sub_m_b, lead_b_b, sub_b_b,
        denom_bkg, m_cuts, beta_cuts, mode=mode
    )

    out = {}
    for sample in signal_samples:
        d_sig = signal[lifetime][sample][req]
        lead_m_s = d_sig.get("leading_mass", [])
        sub_m_s  = d_sig.get("subleading_mass", [])
        lead_b_s = d_sig.get("leading_beta", [])
        sub_b_s  = d_sig.get("subleading_beta", [])

        denom_sig = np.asarray(
            d_sig.get("has1_acc_stau", np.ones(len(lead_m_s), dtype=bool)),
            dtype=bool
        )

        Zs = eff_grid_mass_beta(
            lead_m_s, sub_m_s, lead_b_s, sub_b_s,
            denom_sig, m_cuts, beta_cuts, mode=mode
        )

        best = pick_best_cut_simple(Zs, Zb, m_cuts, beta_cuts, finite_only=True)
        out[sample] = best

    return out

def expected_yield(sample_key, eff):
    return XS_PB[sample_key] * L_PBINV * eff


d_bkg = muons["nominal"]["ob"]["nobib"]

lead_bkg_mass = d_bkg.get("leading_mass", [])
sub_bkg_mass  = d_bkg.get("subleading_mass", [])
N_bkg = len(lead_bkg_mass)
print(N_bkg)
N_bkg_events = len(d_bkg.get("event_has1_all", []))
denom_bkg = np.ones(N_bkg, dtype=bool)

lead_bkg_beta = d_bkg.get("leading_beta", [])
sub_bkg_beta  = d_bkg.get("subleading_beta", [])
bkg_beta_eff_ge1 = eff_vs_beta_threshold(lead_bkg_beta, sub_bkg_beta, denom_bkg, beta_cuts, mode="ge1")
bkg_beta_eff_ge2 = eff_vs_beta_threshold(lead_bkg_beta, sub_bkg_beta, denom_bkg, beta_cuts, mode="ge2")

bkg_mass_eff_ge1 = eff_vs_mass_threshold(lead_bkg_mass, sub_bkg_mass, denom_bkg, m_cuts, mode="ge1")
bkg_mass_eff_ge2 = eff_vs_mass_threshold(lead_bkg_mass, sub_bkg_mass, denom_bkg, m_cuts, mode="ge2")


bkg_mass_yield_ge1 = XS_PB["mumu_bkg"] * L_PBINV * bkg_mass_eff_ge1
bkg_mass_yield_ge2 = XS_PB["mumu_bkg"] * L_PBINV * bkg_mass_eff_ge2

best = best_cuts_per_masspoint(
    signal=signal,
    muons=muons,
    signal_samples=signal_samples,
    m_cuts=m_cuts,
    beta_cuts=beta_cuts,
    window="nominal",
    req="ob",
    option="nobib",
    lifetime=30,
    mode="ge1",   # or "ge2"
)

for sample, info in best.items():
    print(sample, info)
    Ns = expected_yield(sample, info["sig_eff"])
    Nb = expected_yield("mumu_bkg", info["bkg_eff"])
    print(sample, "Ns=", Ns, "Nb=", Nb, "S/sqrt(B)=", Ns/np.sqrt(max(Nb,1e-12)))


with PdfPages(plot_path) as pdf:
    for mode in ["ge1", "ge2"]:
        fig, ax = plt.subplots(figsize=(9, 6))

        bkg_curve = bkg_mass_yield_ge1 if mode == "ge1" else bkg_mass_yield_ge2
        ax.plot(m_cuts, bkg_curve, linestyle="--", linewidth=3, label=f"muon bkg ({mode})")

        for sample in signal_samples:
            d_sig = signal[30][sample]["ob"]  
            lead_sig = d_sig.get("leading_mass", [])
            sub_sig  = d_sig.get("subleading_mass", [])

            denom_sig = np.asarray(d_sig.get("has1_acc_stau", np.ones(len(lead_sig), dtype=bool)), dtype=bool)

            sig_curve = eff_vs_mass_threshold(lead_sig, sub_sig, denom_sig, m_cuts, mode=mode)
            sig_mass_yield = XS_PB[sample] * L_PBINV * sig_curve
            ax.plot(m_cuts, sig_mass_yield, linewidth=2, label=f"stau {float(sample)/1000:.1f} TeV")

        ax.set_xlabel("Mass threshold cut  m > m_cut  [GeV]", fontsize=20)
        ax.set_ylabel("Event yield", fontsize=20)
        # ax.set_ylim(0, 1.02)
        ax.grid(True, alpha=0.2)
        style_ticks(ax, major=18, minor=16)
        ax.legend(fontsize=11, ncol=2, frameon=False)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    for mode in ["ge1", "ge2"]:
        fig, ax = plt.subplots(figsize=(9, 6))

        bkg_curve = bkg_beta_eff_ge1 if mode == "ge1" else bkg_beta_eff_ge2
        ax.plot(beta_cuts, bkg_curve, linestyle="--", linewidth=3, label=f"muon bkg ({mode})")

        for sample in signal_samples:
            d_sig = signal[30][sample]["ob"]  
            lead_sig = d_sig.get("leading_beta", [])
            sub_sig  = d_sig.get("subleading_beta", [])

            denom_sig = np.asarray(d_sig.get("has1_acc_stau", np.ones(len(lead_sig), dtype=bool)), dtype=bool)

            sig_curve = eff_vs_beta_threshold(lead_sig, sub_sig, denom_sig, beta_cuts, mode=mode)
            ax.plot(beta_cuts, sig_curve, linewidth=2, label=f"stau {float(sample)/1000:.1f} TeV")

        ax.set_xlabel(r"Beta threshold cut  $\beta < \beta_{\rm cut}$", fontsize=20)
        ax.set_ylabel("Event efficiency", fontsize=20)
        ax.set_ylim(0, 1.02)
        ax.grid(True, alpha=0.2)
        style_ticks(ax, major=18, minor=16)
        ax.legend(fontsize=11, ncol=2, frameon=True, loc="upper right")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

        d_bkg = muons["nominal"]["ob"]["nobib"]

        lead_m_b = d_bkg.get("leading_mass", [])
        sub_m_b  = d_bkg.get("subleading_mass", [])
        lead_b_b = d_bkg.get("leading_beta", [])
        sub_b_b  = d_bkg.get("subleading_beta", [])

        for mode in ["ge1", "ge2"]:
            Zb = eff_grid_mass_beta(
                lead_m_b, sub_m_b,
                lead_b_b, sub_b_b,
                denom_bkg, m_cuts, beta_cuts, mode=mode
            )
            Nb_grid = XS_PB["mumu_bkg"] * L_PBINV * Zb

        for sample in signal_samples:

            d_sig = signal[30][sample]["ob"]

            lead_m_s = d_sig.get("leading_mass", [])
            sub_m_s  = d_sig.get("subleading_mass", [])
            lead_b_s = d_sig.get("leading_beta", [])
            sub_b_s  = d_sig.get("subleading_beta", [])

            denom_sig = np.asarray(d_sig.get("has1_acc_stau", np.ones(len(lead_m_s), dtype=bool)), dtype=bool)

            for mode in ["ge1", "ge2"]:
                Zs = eff_grid_mass_beta(
                    lead_m_s, sub_m_s,
                    lead_b_s, sub_b_s,
                    denom_sig, m_cuts, beta_cuts, mode=mode
                )

                Ns_grid = XS_PB[sample] * L_PBINV * Zs

    lead_m_b = d_bkg.get("leading_mass", [])
    sub_m_b  = d_bkg.get("subleading_mass", [])
    lead_b_b = d_bkg.get("leading_beta", [])
    sub_b_b  = d_bkg.get("subleading_beta", [])

    for mode in ["ge1", "ge2"]:
        Zb = eff_grid_mass_beta(
            lead_m_b, sub_m_b,
            lead_b_b, sub_b_b,
            denom_bkg, m_cuts, beta_cuts, mode=mode
        )
        Nb_grid = XS_PB["mumu_bkg"] * L_PBINV * Zb

        # plot_eff_heatmap(
        #     pdf, m_cuts, beta_cuts, Zb,
        #     title=f"Muon background efficiency ({mode})",
        #     xlabel="Mass threshold cut  m > m_cut  [GeV]",
        #     ylabel=r"Beta threshold cut  $\beta < \beta_{\rm cut}$"
        # )

        for sample in signal_samples:
            d_sig = signal[30][sample]["ob"]

            lead_m_s = d_sig.get("leading_mass", [])
            sub_m_s  = d_sig.get("subleading_mass", [])
            lead_b_s = d_sig.get("leading_beta", [])
            sub_b_s  = d_sig.get("subleading_beta", [])

            denom_sig = np.asarray(
                d_sig.get("has1_acc_stau", np.ones(len(lead_m_s), dtype=bool)),
                dtype=bool
            )

            Zs = eff_grid_mass_beta(
                lead_m_s, sub_m_s,
                lead_b_s, sub_b_s,
                denom_sig, m_cuts, beta_cuts, mode=mode
            )
            Ns_grid = XS_PB[sample] * L_PBINV * Zs

            # plot_eff_heatmap(
            #     pdf, m_cuts, beta_cuts, Zs,
            #     title=f"Signal efficiency ({float(sample)/1000:.1f} TeV) ({mode})",
            #     xlabel="Mass threshold cut  m > m_cut  [GeV]",
            #     ylabel=r"Beta threshold cut  $\beta < \beta_{\rm cut}$"
            # )

            plot_roc_from_grids(
                pdf, Nb_grid, Ns_grid,
                title=f"ROC from (m_cut, beta_cut) scan: {float(sample)/1000:.1f} TeV vs muons ({mode})"
            )

print(f"Wrote scan plots to {plot_path}")


def _align_and_mask(mL, bL, mS, bS, denom_mask):
    """Align arrays to the same length and apply denom_mask."""
    mL = np.asarray(mL, dtype=float)
    bL = np.asarray(bL, dtype=float)
    denom_mask = np.asarray(denom_mask, dtype=bool)

    n = min(len(mL), len(bL), len(denom_mask))
    mL, bL, denom_mask = mL[:n], bL[:n], denom_mask[:n]

    if mS is not None and bS is not None:
        mS = np.asarray(mS, dtype=float)[:n]
        bS = np.asarray(bS, dtype=float)[:n]
    else:
        mS, bS = None, None

    mL = mL[denom_mask]
    bL = bL[denom_mask]
    if mS is not None:
        mS = mS[denom_mask]
        bS = bS[denom_mask]

    return mL, bL, mS, bS

def pick_best_cut_lowest_nonzero(
    bkg_mass_lead, bkg_beta_lead,
    bkg_mass_sub=None, bkg_beta_sub=None,
    bkg_denom_mask=None,
    sig_mass_lead=None, sig_beta_lead=None,
    sig_mass_sub=None, sig_beta_sub=None,
    sig_denom_mask=None,
    m_cuts=None, beta_cuts=None,
    mode="ge1",
    prefer_target_k=1,           
    allow_zero_background=False, 
):
    """
    Returns dict with best cuts picked at (approximately) the lowest nonzero background acceptance.
    - If allow_zero_background=False, it ignores bkg points with 0 passing events.
    - If prefer_target_k is not None (e.g. 1), it tries to pick the smallest bkg acceptance
      >= k/N (so ~1 event in sample), rather than possibly 2/N, 3/N, etc.
    """

    if bkg_denom_mask is None:
        bkg_denom_mask = np.ones(len(bkg_mass_lead), dtype=bool)
    if sig_denom_mask is None:
        sig_denom_mask = np.ones(len(sig_mass_lead), dtype=bool)

    m_cuts = np.asarray(m_cuts, dtype=float)
    beta_cuts = np.asarray(beta_cuts, dtype=float)

    b_mL, b_bL, b_mS, b_bS = _align_and_mask(
        bkg_mass_lead, bkg_beta_lead, bkg_mass_sub, bkg_beta_sub, bkg_denom_mask
    )
    s_mL, s_bL, s_mS, s_bS = _align_and_mask(
        sig_mass_lead, sig_beta_lead, sig_mass_sub, sig_beta_sub, sig_denom_mask
    )

    Nb = b_mL.size
    Ns = s_mL.size
    if Nb == 0 or Ns == 0:
        return None

    b_finL = np.isfinite(b_mL) & np.isfinite(b_bL)
    s_finL = np.isfinite(s_mL) & np.isfinite(s_bL)

    if mode == "ge2":
        if b_mS is None or b_bS is None or s_mS is None or s_bS is None:
            raise ValueError("mode='ge2' requires subleading arrays for BOTH signal and background.")
        b_finS = np.isfinite(b_mS) & np.isfinite(b_bS)
        s_finS = np.isfinite(s_mS) & np.isfinite(s_bS)

    def count_pass(mc, bc):
        if mode == "ge1":
            b_pass = b_finL & (b_mL > mc) & (b_bL < bc)
            s_pass = s_finL & (s_mL > mc) & (s_bL < bc)
        else:  # ge2
            b_pass = (b_finL & b_finS &
                      (b_mL > mc) & (b_bL < bc) &
                      (b_mS > mc) & (b_bS < bc))
            s_pass = (s_finL & s_finS &
                      (s_mL > mc) & (s_bL < bc) &
                      (s_mS > mc) & (s_bS < bc))
        return int(np.sum(b_pass)), int(np.sum(s_pass))

    best = None

    target_k = prefer_target_k if (prefer_target_k is not None) else None

    feasible_counts = set()
    for bc in beta_cuts:
        for mc in m_cuts:
            nb, _ = count_pass(mc, bc)
            if not allow_zero_background and nb == 0:
                continue
            if target_k is not None and nb < target_k:
                continue
            feasible_counts.add(nb)

    if not feasible_counts:
        if not allow_zero_background:
            return pick_best_cut_lowest_nonzero(
                bkg_mass_lead, bkg_beta_lead, bkg_mass_sub, bkg_beta_sub, bkg_denom_mask,
                sig_mass_lead, sig_beta_lead, sig_mass_sub, sig_beta_sub, sig_denom_mask,
                m_cuts, beta_cuts, mode=mode,
                prefer_target_k=prefer_target_k,
                allow_zero_background=True,
            )
        return None

    bkg_best_count = min(feasible_counts)  
    for j, bc in enumerate(beta_cuts):
        for i, mc in enumerate(m_cuts):
            nb, ns = count_pass(mc, bc)

            if not allow_zero_background and nb == 0:
                continue
            if target_k is not None and nb < target_k:
                continue
            if nb != bkg_best_count:
                continue

            sig_eff = ns / Ns
            bkg_eff = nb / Nb

            cand = {
                "m_cut": float(mc),
                "beta_cut": float(bc),
                "sig_eff": float(sig_eff),
                "bkg_eff": float(bkg_eff),
                "npass_sig": int(ns),
                "npass_bkg": int(nb),
                "N_sig": int(Ns),
                "N_bkg": int(Nb),
                "i_m": int(i),
                "j_beta": int(j),
                "mode": mode,
                "bkg_best_count": int(bkg_best_count),
            }

            if best is None:
                best = cand
            else:
                if cand["sig_eff"] > best["sig_eff"] + 1e-15:
                    best = cand
                elif abs(cand["sig_eff"] - best["sig_eff"]) <= 1e-15 and cand["m_cut"] > best["m_cut"]:
                    best = cand
                elif (abs(cand["sig_eff"] - best["sig_eff"]) <= 1e-15 and
                      cand["m_cut"] == best["m_cut"] and
                      cand["beta_cut"] < best["beta_cut"]):
                    best = cand

    return best

d_bkg = muons["nominal"]["ob"]["nobib"]
bkg_best = None

for sample in signal_samples:
    d_sig = signal[30][sample]["ob"]

    best = pick_best_cut_lowest_nonzero(
        # background
        bkg_mass_lead=d_bkg["leading_mass"],
        bkg_beta_lead=d_bkg["leading_beta"],
        bkg_mass_sub=d_bkg.get("subleading_mass", None),
        bkg_beta_sub=d_bkg.get("subleading_beta", None),
        bkg_denom_mask=np.ones(len(d_bkg["leading_mass"]), dtype=bool),

        # signal
        sig_mass_lead=d_sig["leading_mass"],
        sig_beta_lead=d_sig["leading_beta"],
        sig_mass_sub=d_sig.get("subleading_mass", None),
        sig_beta_sub=d_sig.get("subleading_beta", None),
        sig_denom_mask=np.asarray(d_sig.get("has1_acc_stau", np.ones(len(d_sig["leading_mass"]), bool)), bool),

        m_cuts=m_cuts,
        beta_cuts=beta_cuts,
        mode="ge1",

        prefer_target_k=1,            
        allow_zero_background=False,  
    )

    print(sample, best)


L_ABINV = 10.0
L_PBINV = L_ABINV * 1e6 
N_BKG = 2448

XS_PB = {
    "mumu_bkg": 0.4312,
    "1000": 0.0005108,
    "1500": 0.0004715,
    "2000": 0.0004184,
    "2500": 0.0003528,
    "3000": 0.0002780,
    "3500": 0.0001980,
    "4000": 0.0001173,
    "4500": 0.0000450,
}

B_ASSUMED = 1.0  

def expected_events(xs_pb, eff, L_pbinv=L_PBINV):
    return float(xs_pb) * float(L_pbinv) * float(eff)

def pick_cut_at_one_bkg_event(Zs, Zb, m_cuts, beta_cuts, N_bkg=N_BKG):
    """
    Choose (m_cut, beta_cut) among grid points with nbkg_pass ~= 1 event,
    maximizing signal efficiency.
    Zs, Zb: efficiencies on same grid
    Returns dict with chosen cut and efficiencies and counts.
    """
    nb = Zb * N_bkg 
    dist = np.abs(nb - 1.0)

    ok = np.isfinite(Zs) & np.isfinite(Zb)
    if not np.any(ok):
        return None

    dist = np.where(ok, dist, np.inf)

    best_dist = np.min(dist)
    tol = 1e-12
    band = ok & (dist <= best_dist + tol)

    if not np.any(band):
        band = ok & (dist <= best_dist + 1e-6)

    Zs_band = np.where(band, Zs, -np.inf)
    flat = np.argmax(Zs_band)
    j, i = np.unravel_index(flat, Zs.shape)

    return {
        "m_cut": float(m_cuts[i]),
        "beta_cut": float(beta_cuts[j]),
        "sig_eff": float(Zs[j, i]),
        "bkg_eff": float(Zb[j, i]),
        "nbkg_mc": float(nb[j, i]),
        "dist_to_1": float(dist[j, i]),
        "i_m": int(i),
        "j_beta": int(j),
    }

d_bkg = muons["nominal"]["ob"]["nobib"]
lead_m_b = d_bkg.get("leading_mass", [])
sub_m_b  = d_bkg.get("subleading_mass", None)
lead_b_b = d_bkg.get("leading_beta", [])
sub_b_b  = d_bkg.get("subleading_beta", None)
denom_bkg = np.ones(len(lead_m_b), dtype=bool)  

MODE = "ge1"  # or "ge2"

Zb = eff_grid_mass_beta(
    lead_m_b, sub_m_b, lead_b_b, sub_b_b,
    denom_bkg, m_cuts, beta_cuts, mode=MODE
)
rows = []
for sample in signal_samples:
    d_sig = signal[30][sample]["ob"]  
    lead_m_s = d_sig.get("leading_mass", [])
    sub_m_s  = d_sig.get("subleading_mass", None)
    lead_b_s = d_sig.get("leading_beta", [])
    sub_b_s  = d_sig.get("subleading_beta", None)

    denom_sig = np.asarray(
        d_sig.get("has1_acc_stau", np.ones(len(lead_m_s), dtype=bool)),
        dtype=bool
    )

    Zs = eff_grid_mass_beta(
        lead_m_s, sub_m_s, lead_b_s, sub_b_s,
        denom_sig, m_cuts, beta_cuts, mode=MODE
    )

    best = pick_cut_at_one_bkg_event(Zs, Zb, m_cuts, beta_cuts, N_bkg=N_BKG)
    if best is None:
        continue

    S_exp = expected_events(XS_PB[sample], best["sig_eff"])
    B_exp = float(B_ASSUMED)

    rows.append({
        "mass_GeV": int(sample),
        "sigma_pb": XS_PB[sample],
        "m_cut_GeV": best["m_cut"],
        "beta_cut": best["beta_cut"],
        "sig_eff": best["sig_eff"],
        "bkg_eff": best["bkg_eff"],
        "bkg_mc_pass": best["nbkg_mc"],       
        "S_10ab": S_exp,                      
        "B_assumed": B_exp,                   
        "excluded_S>=3": (S_exp >= 3.0),
        "discovery_S>=5": (S_exp >= 5.0),
    })

df = pd.DataFrame(rows).sort_values("mass_GeV")
pd.set_option("display.max_columns", None)
print(df.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
