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

# DEFINING CUTS 
PT_MIN   = 800.0
MASS_MIN = 500.0
BETA_MAX = 0.99
def pass_pt(t):   return np.isfinite(t["pT"])   and (t["pT"]   > PT_MIN)
def pass_mass(t): return np.isfinite(t["mass"]) and (t["mass"] > MASS_MIN)
def pass_beta(t): return np.isfinite(t["beta"]) and (t["beta"] < BETA_MAX)
def pass_all(t):  return pass_pt(t) and pass_mass(t) and pass_beta(t)

dirs = {"bib": "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/mumu_bkg/bib/",
        "nobib": "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/mumu_bkg/nobib/"}
# windows = ["loose", "nominal"]
# bib_options = ["bib", "nobib"]
bib_options = ["nobib"]
windows = ["nominal"]
CACHE = pathlib.Path("cache/mumu_bkg_stats_nominal_nobib_byevent.pkl")
plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/backgrounds/mumu_tracks_nominal_nobib_byevent.pdf"
print(dirs.items())
n_files = 2500
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

def reco_velo(function_type, times, pos, spatial_unc):
    x = np.asarray(times, dtype=float)
    y = np.asarray(pos, dtype=float)
    s = np.asarray(spatial_unc, dtype=float)

    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(s) & (s > 0)
    x, y, s = x[m], y[m], s[m]

    if x.size < 3 or np.allclose(x, x.mean()):
        return np.nan, np.nan

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
    except Exception:
        v_err = np.nan

    return float(p[0]), v_err

def linearfunc_no_intercept(v, x):
    return v * x
def residual_no_intercept(v, times, pos, spatial_unc, time_unc):
    vv = float(np.atleast_1d(v)[0])
    s_eff = np.sqrt(np.asarray(spatial_unc, float)**2 + (vv * np.asarray(time_unc, float))**2)
    return (linearfunc_no_intercept(vv, times) - pos) / s_eff
def reco_velo_no_intercept(times, pos, spatial_unc, time_unc):
    x = np.asarray(times, dtype=float)
    y = np.asarray(pos, dtype=float)
    sr = np.asarray(spatial_unc, dtype=float)
    st = np.asarray(time_unc, dtype=float)

    m = np.isfinite(x) & np.isfinite(y) & np.isfinite(sr) & np.isfinite(st) & (sr > 0) & (st > 0)
    x, y, sr, st = x[m], y[m], sr[m], st[m]

    if x.size < 3 or np.allclose(x, x.mean()):
        return np.nan, np.nan

    v0 = np.array([guess_velo])

    fit = optimize.least_squares(
        residual_no_intercept,
        v0,
        args=(x, y, sr, st),
        jac="2-point"
    )

    v = float(fit.x[0])

    try:
        s_eff = np.sqrt(sr**2 + (v * st)**2)
        J = fit.jac
        dof = max(1, x.size - 1)
        chi2 = np.sum(((v * x - y) / s_eff) ** 2)
        sigma2 = chi2 / dof
        cov = np.linalg.inv(J.T @ J) * sigma2
        v_err = float(np.sqrt(cov[0, 0]))
    except Exception:
        v_err = np.nan

    return v, v_err



stats = None
if (not rebuild) and os.path.exists(CACHE):
    with open(CACHE, "rb") as f:
        print("Loading in cached arrays...")
        stats = pickle.load(f)

if stats is None:
    stats = {
    window: {
        req: {
        option: {
            "n_pass_pt": [], "n_pass_mass": [], "n_pass_beta": [], "n_pass_all": [],
            "event_has1_pt": [], "event_has2_pt": [],
            "event_has1_mass": [], "event_has2_mass": [],
            "event_has1_beta": [], "event_has2_beta": [],
            "event_has1_all": [], "event_has2_all": [],
            "leading_mass": [], "subleading_mass": [],
            "leading_pT": [], "subleading_pT": [],
            "leading_beta": [], "subleading_beta": [],
            "leading_hits": [], "subleading_hits": [],
            "n_events": 0,
        } for option in bib_options
        } for req in track_reqs
    } for window in windows
    }

         
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    
    for window in windows:
        total_tracks = 0
        no_eta_cut = 0
        super_lum = 0
        print(f"Analyzing {window} window...")
        for option in bib_options:
            print(f"Analyzing {option}...")
            leading_mass = []
            sub_leading_mass = []
            leading_mass_evt = []
            subleading_mass_evt = []
            stats[window]["vb"][option]["leading_mass"] = []
            stats[window]["vb"][option]["subleading_mass"] = []
            for ifile in tqdm(range(n_files)): 
                file_name = f"mumu_bkg_reco{ifile}.slcio"
                file_path = os.path.join(dirs[option], window, file_name)
                if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
                    print(f"couldn't open {file_path}")
                    continue
                reader.open(file_path)
                for event in reader:
                    tracks_by_req = {req: [] for req in track_reqs}
                    for req in track_reqs:
                        stats[window][req][option]["n_events"] += 1

                    rel_nav = build_rel_nav(event)
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
                        total_tracks += 1
                        chi2 = track.getChi2()
                        ndf = track.getNdf()
                        if (chi2/ndf) > chi2_cut:
                            continue
                        track_mcps = nav.getRelatedFromObjects(track)
                        track_hits = track.getTrackerHits()
                        
                        reco_pT = 0.3 * Bfield / fabs(track.getOmega() * 1000.)
                        if track_mcps:
                            for mcp in track_mcps:
                                momentum = mcp.getMomentum()
                                tlv = ROOT.TLorentzVector()
                                tlv.SetPxPyPzE(momentum[0], momentum[1], momentum[2], mcp.getEnergy())
                                if abs(tlv.Eta()) > 0.8:
                                    no_eta_cut += 1
                                    continue
                                true_pT = tlv.Perp()
                                true_beta = tlv.Beta()
                                true_velo = true_beta * speedoflight
                                true_eta = tlv.Eta()
                                true_mass = tlv.M()

                                vb_hits = 0
                                ib_hits = 0
                                ob_hits = 0
                                
                                track_times = []
                                track_pos = []
                                spatial_unc = []
                                time_unc = []
                                track_xyz = []

                                track_sim_r = []
                                track_sim_t = []
                                hit_is_muon = [] 


                                for hit in track_hits:
                                    decoder.setValue(int(hit.getCellID0()))
                                    system = decoder["system"].value()
                                    layer = decoder["layer"].value()
                                    sim_r = np.nan
                                    sim_t = np.nan

                                    relname = system_to_relname.get(system, None)
                                    if relname is not None:
                                        nav_rel = rel_nav[relname]

                                        simhits = nav_rel.getRelatedToObjects(hit)

                                        if simhits:
                                            sh = simhits[0] 
                                            try:
                                                mcp = sh.getMCParticle()         
                                                is_mu = (mcp is not None) and (abs(mcp.getPDG()) == 13)
                                            except Exception:
                                                is_mu = False

                                            hit_is_muon.append(is_mu)
                                            sx = sh.getPosition()[0] 
                                            sy =  sh.getPosition()[1]
                                            sz =  sh.getPosition()[2]
                                            sim_r = math.sqrt(sx*sx + sy*sy + sz*sz)
                                            sim_t = sh.getTime()

                                            track_sim_r.append(sim_r)
                                            track_sim_t.append(sim_t)

                                            is_mu = False



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
                                    track_xyz.append((x, y, z))

                                    hit_pos = np.sqrt(x**2 + y**2 + z**2)
                                    tof = hit_pos/speedoflight

                                    resolution = 0.03
                                    if system > 2:
                                        resolution = 0.06

                                    corrected_t = hit.getTime() - tof
                                    corrected_corrected_t = 2*hit.getTime() - corrected_t

                                    track_times.append(corrected_corrected_t)
                                    track_pos.append(hit_pos)

                                v_fit, v_err = reco_velo_no_intercept(track_times, track_pos, spatial_unc, time_unc) 
                                if np.isfinite(v_fit) and v_fit > speedoflight:
                                    super_lum += 1
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
                                 
                                track_info = {
                                    "pT": float(reco_pT) if np.isfinite(reco_pT) else np.nan,
                                    "beta": float(beta) if np.isfinite(beta) else np.nan,
                                    "mass": float(m_reco) if np.isfinite(m_reco) else np.nan,
                                    "hits": float(total_hits),
                                } 

                                if vb_hits >= 3 and ib_hits >= 2 and ob_hits >=2:
                                    tracks_by_req["ob"].append(track_info)
                                
                                if vb_hits >= 3 and ib_hits >= 2:
                                    tracks_by_req["ib"].append(track_info)

                                if vb_hits >= 3: 
                                    tracks_by_req["vb"].append(track_info)
                    
                    for req in track_reqs:
                        trks = tracks_by_req[req]

                        n_pt   = sum(1 for t in trks if pass_pt(t))
                        n_mass = sum(1 for t in trks if pass_mass(t))
                        n_beta = sum(1 for t in trks if pass_beta(t))
                        n_all  = sum(1 for t in trks if pass_all(t))

                        d = stats[window][req][option]

                        d["n_pass_pt"].append(n_pt)
                        d["n_pass_mass"].append(n_mass)
                        d["n_pass_beta"].append(n_beta)
                        d["n_pass_all"].append(n_all)

                        d["event_has1_pt"].append(1 if n_pt   >= 1 else 0)
                        d["event_has2_pt"].append(1 if n_pt   >= 2 else 0)
                        d["event_has1_mass"].append(1 if n_mass >= 1 else 0)
                        d["event_has2_mass"].append(1 if n_mass >= 2 else 0)
                        d["event_has1_beta"].append(1 if n_beta >= 1 else 0)
                        d["event_has2_beta"].append(1 if n_beta >= 2 else 0)
                        d["event_has1_all"].append(1 if n_all  >= 1 else 0)
                        d["event_has2_all"].append(1 if n_all  >= 2 else 0)

                        passing_all = [t for t in trks if pass_all(t)]
                        passing_all.sort(key=lambda t: t["mass"], reverse=True)

                        if len(passing_all) >= 1:
                            lead = passing_all[0]
                            d["leading_mass"].append(lead["mass"])
                            d["leading_pT"].append(lead["pT"])
                            d["leading_beta"].append(lead["beta"])
                            d["leading_hits"].append(lead["hits"])
                        else:
                            d["leading_mass"].append(float("nan"))
                            d["leading_pT"].append(float("nan"))
                            d["leading_beta"].append(float("nan"))
                            d["leading_hits"].append(float("nan"))

                        if len(passing_all) >= 2:
                            sub = passing_all[1]
                            d["subleading_mass"].append(sub["mass"])
                            d["subleading_pT"].append(sub["pT"])
                            d["subleading_beta"].append(sub["beta"])
                            d["subleading_hits"].append(sub["hits"])
                        else:
                            d["subleading_mass"].append(float("nan"))
                            d["subleading_pT"].append(float("nan"))
                            d["subleading_beta"].append(float("nan"))
                            d["subleading_hits"].append(float("nan"))



        print(f"{total_tracks} tracks")
        print(f"{no_eta_cut} tracks didn't pass eta cut")
        print(f"{super_lum} tracks with v > c")
    print(stats)
    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump(stats, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Writing cache to {CACHE}")
    print("Saved cache successfully.")

labels = {"pT": r"$p_T$ [GeV]",
          "hits": "Hits on track",
          "velo": "Velocity [mm/ns]",
          "true_pT": r"True $p_T$ [GeV]",
          "true_beta": r"True $\beta$=v/c",
          "true_eta": r"True $\eta$",
          "true_mass": "True mass",
          "mass": "Reconstructed mass",
          "velo_res": r"$(v_{reco} - v_{true})/v_{true}$",
          "pT_res": r"$(p_T^{reco} - p_T^{true})/p_T^{true}$"}

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

print_mumu_bkg_summary(stats, windows, bib_options, track_reqs)


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
        }
    if xlims_cfg is None:
        xlims_cfg = {
            "pT": (0, 7000),
            "mass": (0, 2000),
            "beta": (0.9, 1.025),
        }

    features = [
        ("pT",   "leading_pT",   "subleading_pT",   r"$p_T$ [GeV]"),
        ("mass", "leading_mass", "subleading_mass", r"reco mass [GeV]"),
        ("beta", "leading_beta", "subleading_beta", r"$\beta$"),
    ]

    with PdfPages(out_pdf) as pdf:
        for window in windows:
            for option in bib_options:
                for req in track_reqs:
                    d = stats[window][req][option]

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

                        ax.text(0.02, 0.02,
                                f"Frac(events w/ leading) ~ {frac_lead:.3f}\n"
                                f"Frac(events w/ subleading) ~ {frac_sub:.3f}",
                                ha="left", va="bottom", transform=ax.transAxes, fontsize=12)

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
