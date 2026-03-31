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
from array import array

SAVE_EVERY_FILES = 50  

def save_cache_atomic(stats, cache_path):
    cache_path = pathlib.Path(cache_path)
    cache_path.parent.mkdir(exist_ok=True)
    tmp = str(cache_path) + ".tmp"
    with open(tmp, "wb") as f:
        pickle.dump(stats, f, protocol=pickle.HIGHEST_PROTOCOL)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp, cache_path)

# DEFINING CUTS 
ETA_MAX     = 0.8
CHI2NDF_MAX = 3.0
VRMSW_MAX   = 1.6          
OB_HITS_MIN = 2
PT_MAX      = 10_000.0   
NO_INTERCEPT = True  

def track_eta_from_tanL(tanL):
    if not np.isfinite(tanL):
        return np.nan
    return float(np.arcsinh(tanL))

def pass_ob_req(t):
    return (t["vb_hits"] >= 3) and (t["ib_hits"] >= 2) and (t["ob_hits"] >= 2)

def pass_track_level(t):
    return (
        pass_ob_req(t) and
        np.isfinite(t["eta"]) and (abs(t["eta"]) < ETA_MAX) and
        np.isfinite(t["chi2ndf"]) and (t["chi2ndf"] < CHI2NDF_MAX) and
        np.isfinite(t["vrmsw"]) and (t["vrmsw"] < VRMSW_MAX) and
        np.isfinite(t["pT"]) and (t["pT"] < PT_MAX) and
        np.isfinite(t["beta"]) and (0 < t["beta"] < 1) and
        np.isfinite(t["mass"])
    )

# sample_cuts = {
#     1000: {"pT_min" : 100, "mass_min": 950, "mass_max": 1050, "beta_max": 0.973},

# }

dirs = {"bib": "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/mumu_bkg/bib/",
        "nobib": "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/mumu_bkg/nobib/"}
# windows = ["loose", "nominal"]
# bib_options = ["bib", "nobib"]
bib_options = ["nobib"]
windows = ["nominal"]
CACHE = pathlib.Path("cache/mumu_bkg_stats_nominal_nobib_byevent.pkl")
plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/backgrounds/mumu_tracks_nominal_nobib_byevent.pdf"
print(dirs.items())
n_files = 1000000
Bfield = 3.57
speedoflight = 299792458/1000000  # mm/ns
chi2_cut = 3

SAVE_TRACK_DISPLAYS = True
MAX_TRACK_DISPLAYS = 200          
TRACK_DISPLAY_PDF = "/scratch/miralittmann/analysis/mira_analysis_code/backgrounds/mumu_track_displays.pdf"

DISPLAY_MASS_WINDOW = None #(400.0, 600.0)

track_displays = []  


parser = argparse.ArgumentParser()
parser.add_argument("--rebuild", action="store_true")
args = parser.parse_args()
rebuild = args.rebuild

def build_rel_nav(event):
    nav = {
        "VXDBarrel": UTIL.LCRelationNavigator(event.getCollection("VXDBarrelHitsRelations")),
        "VXDEndcap": UTIL.LCRelationNavigator(event.getCollection("VXDEndcapHitsRelations")),
        "ITBarrel" : UTIL.LCRelationNavigator(event.getCollection("ITBarrelHitsRelations")),
        "ITEndcap" : UTIL.LCRelationNavigator(event.getCollection("ITEndcapHitsRelations")),
        "OTBarrel" : UTIL.LCRelationNavigator(event.getCollection("OTBarrelHitsRelations")),
        "OTEndcap" : UTIL.LCRelationNavigator(event.getCollection("OTEndcapHitsRelations")),
    }
    enc = event.getCollection("ITBarrelHits").getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
    nav["_ENCODING"] = enc
    nav["_DECODER"] = UTIL.BitField64(enc)
    return nav

system_to_relname = {
    1: "VXDBarrel", 2: "VXDEndcap",
    3: "ITBarrel",  4: "ITEndcap",
    5: "OTBarrel",  6: "OTEndcap",
}

guess_velo = 299.8

track_reqs = ["vb", "ib", "ob"]

def linearfunc(p, x):
    # p[0] = velocity [mm/ns], p[1] = intercept [mm]
    return p[0] * x + p[1]

def residual(p, function_type, times, pos, spatial_unc):
    # weighted residuals
    return (function_type(p, times) - pos) / spatial_unc

def fit_rms(p, function_type, times, pos, spatial_unc):
    x = np.asarray(times, dtype=float)
    y = np.asarray(pos, dtype=float)
    s = np.asarray(spatial_unc, dtype=float)

    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(s) & (s > 0)
    x, y, s = x[m], y[m], s[m]
    if x.size == 0:
        return np.nan, np.nan

    yhat = function_type(p, x)
    r = yhat - y                
    rms_unw = float(np.sqrt(np.mean(r*r)))
    rw = r / s                  
    rms_w = float(np.sqrt(np.mean(rw*rw)))
    return rms_unw, rms_w

def time_rms_from_fit(v, t, r, time_unc, b=0.0):
    t = np.asarray(t, float)
    r = np.asarray(r, float)
    st = np.asarray(time_unc, float)

    m = np.isfinite(t) & np.isfinite(r) & np.isfinite(st) & (st > 0)
    t, r, st = t[m], r[m], st[m]
    if t.size < 3 or (not np.isfinite(v)) or abs(v) < 1e-12:
        return np.nan, np.nan

    t_pred = (r - b) / v
    dt = t - t_pred

    rms_t = float(np.sqrt(np.mean(dt * dt)))          
    rms_pull = float(np.sqrt(np.mean((dt / st) ** 2))) 
    return rms_t, rms_pull

def reco_velo(function_type, times, pos, spatial_unc):
    x = np.asarray(times, dtype=float)
    y = np.asarray(pos, dtype=float)
    s = np.asarray(spatial_unc, dtype=float)

    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(s) & (s > 0)
    x, y, s = x[m], y[m], s[m]

    if x.size < 3 or np.allclose(x, x.mean()):
        return np.nan, np.nan, np.nan, np.nan

    p0 = np.array([guess_velo, 0.0])

    fit = optimize.least_squares(
        residual, p0,
        args=(function_type, x, y, s),
        jac='2-point'
    )
    p = fit.x  
    try:
        J = fit.jac
        dof = max(1, x.size - p.size)
        chi2 = np.sum(((function_type(p, x) - y) / s) ** 2)
        sigma2 = chi2 / dof
        cov = np.linalg.inv(J.T @ J) * sigma2
        v_err = float(np.sqrt(cov[0, 0]))
        rms_unw, rms_w = fit_rms(p, function_type, x, y, s)

    except Exception:
        v_err = np.nan
        rms_unw, rms_w = np.nan, np.nan

    return float(p[0]), v_err, rms_unw, rms_w

def linearfunc_no_intercept(v, x):
    return v * x
def residual_no_intercept(v, times, pos, spatial_unc, time_unc):
    vv = float(np.atleast_1d(v)[0])
    s_eff = np.sqrt(np.asarray(spatial_unc, float)**2 + (vv * np.asarray(time_unc, float))**2)
    return (linearfunc_no_intercept(vv, times) - pos) / s_eff

def reco_velo_no_intercept(times, pos, spatial_unc, time_unc):
    x = np.asarray(times, dtype=float)
    y = np.asarray(pos, dtype=float)
    s = np.asarray(spatial_unc, dtype=float)
    st = np.asarray(time_unc, dtype=float)

    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(s) & (s > 0) & np.isfinite(st) & (st > 0)
    x, y, s, st = x[m], y[m], s[m], st[m]

    if x.size < 3 or np.allclose(x, x.mean()):
        return np.nan, np.nan, np.nan, np.nan

    v0 = np.array([guess_velo])

    def residual0(v, times, pos, spatial_unc):
        vv = float(np.atleast_1d(v)[0])
        return (vv * times - pos) / spatial_unc

    fit = optimize.least_squares(
        residual0,
        v0,
        args=(x, y, s),
        jac="2-point"
    )

    v = float(fit.x[0])

    rms_t, rms_pull = time_rms_from_fit(v, x, y, st, b=0.0)

    try:
        r = (v * x - y)
        J = fit.jac
        dof = max(1, x.size - 1)
        chi2 = np.sum((r / s) ** 2)
        sigma2 = chi2 / dof
        cov = np.linalg.inv(J.T @ J) * sigma2
        v_err = float(np.sqrt(cov[0, 0]))
    except Exception:
        v_err = np.nan

    return v, v_err, rms_t, rms_pull

stats = None
if (not rebuild) and os.path.exists(CACHE):
    with open(CACHE, "rb") as f:
        print("Loading in cached arrays...")
        stats = pickle.load(f)

if stats is None:
    stats = {
        "n_events": 0,
        "last_file": -1,

        "leading_mass":      array("f"),
        "subleading_mass":   array("f"),
        "leading_pT":        array("f"),
        "subleading_pT":     array("f"),
        "leading_beta":      array("f"),
        "subleading_beta":   array("f"),
        "leading_d0":        array("f"),
        "subleading_d0":     array("f"),
        "leading_z0":        array("f"),
        "subleading_z0":     array("f"),
        "leading_wvrms":     array("f"),
        "subleading_wvrms":  array("f"),
    }

    start_file = int(stats.get("last_file", -1)) + 1
    print(f"Resuming from file {start_file}")
         
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    
    for ifile in tqdm(range(start_file, n_files)): 
        file_name = f"mumu_bkg_reco{ifile}.slcio"
        file_path = os.path.join(dirs["nobib"], "nominal", file_name)
        if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
            print(f"couldn't open {file_path}")
            continue
        reader.open(file_path)
        for event in reader:
            tracks_by_req = {req: [] for req in track_reqs}
            all_collections = event.getCollectionNames() 
            mcp_collection = event.getCollection("MCParticle") if "MCParticle" in all_collections else None
            track_collection = event.getCollection("SiTracks") if "SiTracks" in all_collections else None 
            track_relation_collection = event.getCollection("MCParticle_SiTracks") if "MCParticle_SiTracks" in all_collections else None 
            if not (mcp_collection and track_collection and track_relation_collection):
                print("issue 1")
                continue
            test_hit_coll = event.getCollection("VXDBarrelHits")
            if test_hit_coll is None:
                continue
            encoding = test_hit_coll.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
            decoder = UTIL.BitField64(encoding)
            nav = UTIL.LCRelationNavigator(track_relation_collection)
            
            for itrack, track in enumerate(track_collection):
                chi2 = track.getChi2()
                ndf = track.getNdf()
                if (chi2/ndf) > chi2_cut:
                    continue
                track_mcps = nav.getRelatedFromObjects(track)
                track_hits = track.getTrackerHits()
                
                reco_pT = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
                if track_mcps:
                    for mcp in track_mcps:
                        if abs(mcp.getPDG()) != 13: # just looking at muon tracks 
                            continue
                        momentum = mcp.getMomentum()
                        tlv = ROOT.TLorentzVector()
                        tlv.SetPxPyPzE(momentum[0], momentum[1], momentum[2], mcp.getEnergy())
                        if abs(tlv.Eta()) > 0.8:
                            continue
                        true_pT = tlv.Perp()
                        true_beta = tlv.Beta()
                        true_velo = true_beta * speedoflight
                        true_eta = tlv.Eta()
                        true_mass = tlv.M()
                        
                        d0 = track.getD0()
                        z0 = track.getZ0()

                        vb_hits = 0
                        ib_hits = 0
                        ob_hits = 0
                        
                        track_times = []
                        track_pos = []
                        spatial_unc = []
                        time_unc = []

                        for hit in track_hits:
                            decoder.setValue(int(hit.getCellID0()))
                            system = decoder["system"].value()
                            layer = decoder["layer"].value()
                            sim_r = np.nan
                            sim_t = np.nan

                            if system in (1,2):
                                vb_hits += 0.5
                                spatial_unc.append(0.005)
                                time_unc.append(0.03)
                            elif system in (3,4):
                                ib_hits += 1
                                spatial_unc.append(0.007)
                                time_unc.append(0.06)
                            elif system in (5,6):
                                ob_hits += 1
                                spatial_unc.append(0.007)
                                time_unc.append(0.06)

                            hit_time = hit.getTime()
                            x = hit.getPosition()[0]
                            y = hit.getPosition()[1]
                            z = hit.getPosition()[2]

                            hit_pos = np.sqrt(x**2 + y**2 + z**2)
                            tof = hit_pos/speedoflight

                            resolution = 0.03
                            if system > 2:
                                resolution = 0.06

                            corrected_t = hit.getTime() - tof
                            corrected_corrected_t = 2*hit.getTime() - corrected_t

                            track_times.append(corrected_corrected_t)
                            track_pos.append(hit_pos)

                        v_fit, v_err, rms_uw, rms_w = reco_velo_no_intercept(track_times, track_pos, spatial_unc, time_unc) 
                        if np.isfinite(v_fit) and v_fit > speedoflight:
                            v_fit = speedoflight
                        beta = v_fit / speedoflight

                        total_hits = vb_hits + ib_hits + ob_hits
                        pT_res = (reco_pT - true_pT) / true_pT
                        velo_res = (v_fit - true_velo) / true_velo

                        try:
                            tanL = track.getTanLambda()
                        except Exception:
                            tanL = np.nan

                        if np.isfinite(tanL):
                            reco_pz = reco_pT * tanL
                            reco_p  = math.sqrt(reco_pT**2 + reco_pz**2)
                        else:
                            reco_p  = np.nan
                        
                        if np.isfinite(reco_p) and np.isfinite(beta) and (0 < beta <= 1):
                            m_reco = reco_p * math.sqrt(1.0/(beta*beta) - 1.0)
                        else:
                            m_reco = np.nan

                        try:
                            tanL = track.getTanLambda()
                        except Exception:
                            tanL = np.nan
                        eta = float(np.arcsinh(tanL)) if np.isfinite(tanL) else np.nan

                        track_info = {
                        "pT": float(reco_pT) if np.isfinite(reco_pT) else np.nan,
                        "beta": float(beta) if np.isfinite(beta) else np.nan,
                        "mass": float(m_reco) if np.isfinite(m_reco) else np.nan,
                        "d0": float(d0) if np.isfinite(d0) else np.nan,
                        "z0": float(z0) if np.isfinite(z0) else np.nan,
                        "vrmsw": float(rms_w) if np.isfinite(rms_w) else np.nan,
                        "chi2ndf": float(track.getChi2()/track.getNdf()) if track.getNdf() > 0 else np.nan,
                        "vb_hits": vb_hits,
                        "ib_hits": ib_hits,
                        "ob_hits": ob_hits,
                        "eta": eta,
                        } 

                        if vb_hits >= 3 and ib_hits >= 2 and ob_hits >=2:
                            tracks_by_req["ob"].append(track_info)
            
            trks = tracks_by_req["ob"]
            passing = [t for t in trks if pass_track_level(t)]
            passing.sort(key=lambda t: t["pT"], reverse=True)  

            if len(passing) >= 1:
                lead = passing[0]
                stats["leading_pT"].append(lead["pT"])
                stats["leading_beta"].append(lead["beta"])
                stats["leading_mass"].append(lead["mass"])
                stats["leading_d0"].append(lead["d0"])
                stats["leading_z0"].append(lead["z0"])
                stats["leading_wvrms"].append(lead["vrmsw"])
            else:
                stats["leading_pT"].append(np.nan)
                stats["leading_beta"].append(np.nan)
                stats["leading_mass"].append(np.nan)
                stats["leading_d0"].append(np.nan)
                stats["leading_z0"].append(np.nan)
                stats["leading_wvrms"].append(np.nan)

            if len(passing) >= 2:
                sub = passing[1]
                stats["subleading_pT"].append(sub["pT"])
                stats["subleading_beta"].append(sub["beta"])
                stats["subleading_mass"].append(sub["mass"])
                stats["subleading_d0"].append(sub["d0"])
                stats["subleading_z0"].append(sub["z0"])
                stats["subleading_wvrms"].append(sub["vrmsw"])
            else:
                stats["subleading_pT"].append(np.nan)
                stats["subleading_beta"].append(np.nan)
                stats["subleading_mass"].append(np.nan)
                stats["subleading_d0"].append(np.nan)
                stats["subleading_z0"].append(np.nan)
                stats["subleading_wvrms"].append(np.nan)

            stats["n_events"] += 1 
        
        reader.close()
        stats["last_file"] = ifile

        if (ifile % SAVE_EVERY_FILES) == 0:
            save_cache_atomic(stats, CACHE)
            print(f"[checkpoint] saved at file {ifile}, n_events={stats['n_events']}")
                    
    save_cache_atomic(stats, CACHE)
    print(f"Writing cache to {CACHE}")

labels = {"pT": r"$p_T$ [GeV]",
          "hits": "Hits on track",
          "velo": "Velocity [mm/ns]",
          "true_pT": r"True $p_T$ [GeV]",
          "true_beta": r"True $\beta$=v/c",
          "true_eta": r"True $\eta$",
          "true_mass": "True mass",
          "mass": "Reconstructed mass",
          "velo_res": r"$(v_{reco} - v_{true})/v_{true}$",
          "pT_res": r"$(p_T^{reco} - p_T^{true})/p_T^{true}$",
          "d0": r"D0 impact parameter",
          "z0": r"Z0 impact parameter",
          "vrms_mm": r"fit RMS residual [mm]",
          "vrmsw":   r"fit RMS weighted residual"}

def print_mumu_bkg_summary(stats, windows, bib_options, track_reqs):
    sep = "-" * 110

    for window in windows:
        for option in bib_options:
            print("\n" + sep)
            print(f"Muon background | window = {window} | option = {option}")
            print(sep)

            header = (
                "req   N_events | "
                "PT:   >=1pass/den  >=2pass/den  |   "
                "   MASS: >=1pass/den  >=2pass/den  |    "
                "   BETA: >=1pass/den  >=2pass/den  |    "
                "   ALL:  >=1pass/den  >=2pass/den"
            )
            print(header)
            print(sep)

            for req in track_reqs:
                d = stats[window][req][option]

                N = d["n_events"]
                if N == 0:
                    print(f"{req:>3}   0 events")
                    continue

                def fmt(n):
                    return f"{n}/{N} ({100*n/N:5.1f}%)"

                pt1   = sum(d["event_has1_pt"])
                pt2   = sum(d["event_has2_pt"])
                m1    = sum(d["event_has1_mass"])
                m2    = sum(d["event_has2_mass"])
                b1    = sum(d["event_has1_beta"])
                b2    = sum(d["event_has2_beta"])
                a1    = sum(d["event_has1_all"])
                a2    = sum(d["event_has2_all"])

                print(
                    f"{req:>3} {N:9d} | "
                    f"{fmt(pt1):>14} {fmt(pt2):>14} | "
                    f"{fmt(m1):>14} {fmt(m2):>14} | "
                    f"{fmt(b1):>14} {fmt(b2):>14} | "
                    f"{fmt(a1):>14} {fmt(a2):>14}"
                )

            print(sep)
            print("Note:")
            print(" - Denominator = all processed events (background).")
            print(" - >=1pass : event has at least one track passing the cut.")
            print(" - >=2pass : event has at least two tracks passing the cut.")
            print(sep)

# print_mumu_bkg_summary(stats, windows, bib_options, track_reqs)

# arr = np.array(stats["nominal"]["ob"]["nobib"]["leading_vrms_mm"], dtype=float)
# arr = arr[np.isfinite(arr)]
# print("N finite:", arr.size)
# print("min/median/90%/max:", np.min(arr), np.median(arr), np.quantile(arr, 0.90), np.max(arr))


def _event_norm_hist(ax, arr, N_events, bins, label):
    x = np.asarray(arr, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0 or N_events <= 0:
        return 0.0
    w = np.full_like(x, 1.0 / N_events, dtype=float)  # sums to fraction of events
    ax.hist(x, bins=bins, weights=w, histtype="step", linewidth=2, label=label)
    return x.size / N_events  # fraction of events with a finite value

def plot_mumu_leading_subleading_eventnorm(stats, windows, bib_options, track_reqs,
                                          out_pdf,
                                          bins_cfg=None,
                                          xlims_cfg=None,
                                          tick_major=18,
                                          tick_minor=16):
    if bins_cfg is None:
        bins_cfg = {
            "pT":   np.linspace(0, 7000, 51),
            "mass": np.linspace(0, 2000, 61),
            "beta": np.linspace(0.9, 1.025, 53), 
            "d0": np.linspace(-0.01,0.01, 51),
            "z0": np.linspace(-0.01,0.01, 51),
            "vrms_mm": np.linspace(0, 0.1, 51),  
            "vrmsw":   np.linspace(0, 2, 51),
        }
    if xlims_cfg is None:
        xlims_cfg = {
            "pT": (0, 7000),
            "mass": (0, 2000),
            "beta": (0.9, 1.025),
            "d0": (-0.01, 0.01),
            "z0": (-0.01, 0.01),
            "vrms_mm": (0,0.1),
            "vrmsw": (0,2)
        }

    features = [
        ("pT",   "leading_pT",   "subleading_pT",   r"$p_T$ [GeV]"),
        ("mass", "leading_mass", "subleading_mass", r"reco mass [GeV]"),
        ("beta", "leading_beta", "subleading_beta", r"$\beta$"),
        ("d0", "leading_d0", "subleading_d0", r"D0"),
        ("z0", "leading_z0", "subleading_z0", r"Z0"),
        ("vrms_mm", "leading_vrms_mm", "subleading_vrms_mm", r"Velocity fit RMS time residual"),
        ("vrmsw",   "leading_vrmsw",   "subleading_vrmsw",   r"Velocity fit RMS time weighted residual"),
    ]

    with PdfPages(out_pdf) as pdf:
        for window in windows:
            for option in bib_options:
                for req in track_reqs:
                    d = stats

                    # Prefer explicit counter; otherwise infer from array length
                    N_events = int(d.get("n_events", 0))
                    if N_events <= 0:
                        # fallback if you didn't maintain n_events
                        N_events = len(d.get("event_has1_all", []))
                    if N_events <= 0:
                        continue

                    for key, lead_key, sub_key, xlabel in features:
                        fig, ax = plt.subplots(figsize=(8, 6))

                        frac_lead = _event_norm_hist(
                            ax, d.get(lead_key, []), N_events, bins_cfg[key],
                            label="Leading"
                        )
                        frac_sub = _event_norm_hist(
                            ax, d.get(sub_key, []), N_events, bins_cfg[key],
                            label="Subleading"
                        )

                        ax.set_xlabel(xlabel, fontsize=20)
                        ax.set_ylabel("Fraction of events per bin", fontsize=20)

                        ax.tick_params(axis="both", which="major",
                                       labelsize=tick_major, length=6, width=1.5)
                        ax.tick_params(axis="both", which="minor",
                                       labelsize=tick_minor, length=4, width=1.0)

                        if key in xlims_cfg and xlims_cfg[key] is not None:
                            ax.set_xlim(*xlims_cfg[key])

                        ax.grid(True, alpha=0.2)
                        ax.legend(frameon=False, fontsize=13, loc="upper right")

                        ax.text(0.02, 0.98, "Muon Collider",
                                ha="left", va="top", transform=ax.transAxes,
                                fontsize=20, fontweight="bold", style="italic")
                        ax.text(0.02, 0.93, f"muons, {option}, {window}, req={req}",
                                ha="left", va="top", transform=ax.transAxes, fontsize=14)
                        ax.text(0.02, 0.89, f"N_events={N_events}",
                                ha="left", va="top", transform=ax.transAxes, fontsize=14)

                        ax.text(0.98, 0.02,
                                f"Frac(events w/ leading) ~ {frac_lead:.3f}\n"
                                f"Frac(events w/ subleading) ~ {frac_sub:.3f}",
                                ha="right", va="bottom", transform=ax.transAxes, fontsize=12)

                        fig.tight_layout()
                        pdf.savefig(fig)
                        plt.close(fig)

    print(f"Saved muon event-normalized leading/subleading plots to {out_pdf}")

plot_mumu_leading_subleading_eventnorm(
    stats=stats,
    windows=windows,
    bib_options=bib_options,
    track_reqs=["ob"],   
    out_pdf="/scratch/miralittmann/analysis/mira_analysis_code/backgrounds/mumu_lead_sublead.pdf",
    tick_major=20, tick_minor=18
)
