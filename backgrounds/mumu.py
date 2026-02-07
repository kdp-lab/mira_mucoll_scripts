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

dirs = {"bib": "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/mumu_bkg/bib/",
        "nobib": "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/mumu_bkg/nobib/"}
# windows = ["loose", "nominal"]
# bib_options = ["bib", "nobib"]
bib_options = ["nobib"]
windows = ["nominal"]
CACHE = pathlib.Path("cache/mumu_bkg_stats_nominal_nobib_noeta.pkl")
plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/backgrounds/mumu_tracks_nominal_nobib_noeta.pdf"
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
                "bib": {
                    "pT": [],
                    "hits": [],
                    "velo": [],
                    "mass": [],
                    "true_eta": [],
                    "true_pT": [],
                    "true_beta": [],
                    "true_mass": [],
                    "pT_res": [],
                    "velo_res": [],

                    "v_fit_err": [],
                    "full_p": [],

                    "time_span": [],
                    "radial_span": []
                },
                "nobib": {
                    "pT": [],
                    "hits": [],
                    "velo": [],
                    "mass": [],
                    "true_eta": [],
                    "true_pT": [],
                    "true_mass": [],
                    "true_beta": [],
                    "pT_res": [],
                    "velo_res": [],

                    "v_fit_err": [],
                    "full_p": [],
                    
                    "time_span": [],
                    "radial_span": []
                },
            }
            for req in track_reqs
        }
        for window in windows
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
                    event_masses = []
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
                                # if abs(tlv.Eta()) > 0.8:
                                #     no_eta_cut += 1
                                #     continue
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
                                
                                if np.isfinite(m_reco) and (vb_hits >= 3 and ib_hits >= 2 and ob_hits >= 2):
                                    event_masses.append(m_reco)

                                if SAVE_TRACK_DISPLAYS and (len(track_displays) < MAX_TRACK_DISPLAYS):
                                    ok = True
                                    if DISPLAY_MASS_WINDOW is not None:
                                        mmin, mmax = DISPLAY_MASS_WINDOW
                                        ok = np.isfinite(m_reco) and (mmin <= m_reco <= mmax)

                                    if ok and len(track_times) >= 3 and len(track_pos) >= 3 and len(track_xyz) >= 3:
                                        tt = np.asarray(track_times, dtype=float)
                                        rr = np.asarray(track_pos, dtype=float)
                                        ss = np.asarray(spatial_unc, dtype=float)

                                        mfit = np.isfinite(tt) & np.isfinite(rr) & np.isfinite(ss) & (ss > 0)
                                        tt2, rr2, ss2 = tt[mfit], rr[mfit], ss[mfit]
                                        mu_flags = np.asarray(hit_is_muon, dtype=object)
                                        mu2 = mu_flags[mfit].tolist()

                                        b_fit = np.nan
                                        if np.isfinite(v_fit) and tt2.size >= 2:
                                            w = 1.0 / (ss2 * ss2)
                                            b_fit = np.sum(w * (rr2 - v_fit * tt2)) / np.sum(w)

                                        track_displays.append({
                                            "window": window,
                                            "option": option,
                                            "req": "vb",
                                            "chi2ndf": (chi2/ndf) if (ndf != 0) else np.nan,
                                            "vb_hits": vb_hits, "ib_hits": ib_hits, "ob_hits": ob_hits,
                                            "times": tt.tolist(),
                                            "r": rr.tolist(),
                                            "xyz": track_xyz,
                                            "spatial_unc": np.asarray(spatial_unc, dtype=float).tolist(),
                                            "v_fit": float(v_fit) if np.isfinite(v_fit) else np.nan,
                                            "v_err": float(v_err) if np.isfinite(v_err) else np.nan,
                                            "b_fit": float(b_fit) if np.isfinite(b_fit) else np.nan,
                                            "beta": float(beta) if np.isfinite(beta) else np.nan,
                                            "pT": float(reco_pT) if np.isfinite(reco_pT) else np.nan,
                                            "p": float(reco_p) if np.isfinite(reco_p) else np.nan,
                                            "m_reco": float(m_reco) if np.isfinite(m_reco) else np.nan,
                                            "true_pT": float(true_pT) if np.isfinite(true_pT) else np.nan,
                                            "true_beta": float(true_beta) if np.isfinite(true_beta) else np.nan,
                                            "true_mass": float(true_mass) if np.isfinite(true_mass) else np.nan,
                                            "true_eta": float(true_eta) if np.isfinite(true_eta) else np.nan,

                                            "reco_t": tt2.tolist(),
                                            "reco_r": rr2.tolist(),
                                            "reco_unc": ss2.tolist(),

                                            "sim_t": np.asarray(track_sim_t, dtype=float)[mfit].tolist(),
                                            "sim_r": np.asarray(track_sim_r, dtype=float)[mfit].tolist(),

                                            "hit_is_muon": mu2,
                                        })

                                
                                if len(track_times) > 0:
                                    time_span = max(track_times) - min(track_times)
                                    stats[window]["vb"][option]["time_span"].append(time_span)
                                
                                if len(track_pos) > 0:
                                    radial_span = max(track_pos) - min(track_pos)
                                    stats[window]["vb"][option]["radial_span"].append(radial_span)

                                

                                if vb_hits >= 3 and ib_hits >= 2 and ob_hits >=2:
                                    stats[window]["ob"][option]["pT"].append(reco_pT)
                                    stats[window]["ob"][option]["velo"].append(v_fit)
                                    stats[window]["ob"][option]["hits"].append(total_hits)
                                    stats[window]["ob"][option]["mass"].append(m_reco)
                                
                                    stats[window]["ob"][option]["true_beta"].append(true_beta)
                                    stats[window]["ob"][option]["true_pT"].append(true_pT)
                                    stats[window]["ob"][option]["true_eta"].append(true_eta)

                                    stats[window]["ob"][option]["pT_res"].append(pT_res)
                                    stats[window]["ob"][option]["velo_res"].append(velo_res)
                                    stats[window]["ob"][option]["v_fit_err"].append(v_err)
                                    stats[window]["ob"][option]["full_p"].append(reco_p)
                                    stats[window]["ob"][option]["true_mass"].append(true_mass)
                                
                                if vb_hits >= 3 and ib_hits >= 2:
                                    stats[window]["ib"][option]["pT"].append(reco_pT)
                                    stats[window]["ib"][option]["velo"].append(v_fit)
                                    stats[window]["ib"][option]["hits"].append(total_hits)
                                    stats[window]["ib"][option]["mass"].append(m_reco)
                                
                                    stats[window]["ib"][option]["true_beta"].append(true_beta)
                                    stats[window]["ib"][option]["true_pT"].append(true_pT)
                                    stats[window]["ib"][option]["true_eta"].append(true_eta)

                                    stats[window]["ib"][option]["pT_res"].append(pT_res)
                                    stats[window]["ib"][option]["velo_res"].append(velo_res)
                                    stats[window]["ib"][option]["v_fit_err"].append(v_err)
                                    stats[window]["ib"][option]["true_mass"].append(true_mass)
                                    stats[window]["ib"][option]["full_p"].append(reco_p)

                                if vb_hits >= 3:
                                    stats[window]["vb"][option]["pT"].append(reco_pT)
                                    stats[window]["vb"][option]["velo"].append(v_fit)
                                    stats[window]["vb"][option]["hits"].append(total_hits)
                                    stats[window]["vb"][option]["mass"].append(m_reco)

                                    stats[window]["vb"][option]["true_beta"].append(true_beta)
                                    stats[window]["vb"][option]["true_pT"].append(true_pT)
                                    stats[window]["vb"][option]["true_eta"].append(true_eta)

                                    stats[window]["vb"][option]["pT_res"].append(pT_res)
                                    stats[window]["vb"][option]["velo_res"].append(velo_res)
                                    stats[window]["vb"][option]["v_fit_err"].append(v_err)
                                    stats[window]["vb"][option]["full_p"].append(reco_p)
                                    stats[window]["vb"][option]["true_mass"].append(true_mass)
                    if len(event_masses) >= 1:
                        event_masses.sort(reverse=True)
                        stats[window]["vb"][option]["leading_mass"].append(event_masses[0])
                        if len(event_masses) >= 2:
                            stats[window]["vb"][option]["subleading_mass"].append(event_masses[1])
                        else:
                            stats[window]["vb"][option]["subleading_mass"].append(float("nan"))



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

def plot_feature(feature, n_bins, x_lim=None):
    fig, ax = plt.subplots()
    feature_arr = np.asarray(stats[window]["vb"][option][feature])
    print(feature, np.min(feature_arr), np.max(feature_arr))
    if x_lim:
        mask = (feature_arr >= x_lim[0]) & (feature_arr <= x_lim[1])
        feature_arr = feature_arr[mask]
    weights = np.full_like(feature_arr, 100.0/feature_arr.size, dtype=float)
    ax.hist(feature_arr, bins=n_bins, weights=weights, histtype="step", facecolor="none")
    ax.set_xlabel(labels[feature], fontsize=20)
    ax.set_ylabel("Normalized counts", fontsize=20)
    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    if x_lim:
        ax.set_xlim(x_lim[0], x_lim[1])

    ax.text(
        0.02, 0.98,
        "Muon Collider",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=20,
        fontweight="bold",
        style="italic",
    )
    ax.text(
        0.02, 0.90,
        f"muons, {option}, {window}",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=15
    ) 
    ax.text(
        0.02, 0.83,
        "MuColl_v1",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=15
    )
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def plot_pt_vs_hits(pT_lim=(0, 10000), hits_lim=None, pT_bins=20, hits_bins=10):
    fig, ax = plt.subplots()

    pT = np.asarray(stats[window]["vb"][option]["mass"], dtype=float)
    hits = np.asarray(stats[window]["vb"][option]["velo_res"], dtype=float)

    m = np.isfinite(pT) & np.isfinite(hits)
    if pT_lim is not None:
        m &= (pT >= pT_lim[0]) & (pT <= pT_lim[1])
    if hits_lim is not None:
        m &= (hits >= hits_lim[0]) & (hits <= hits_lim[1])

    pT = pT[m]
    hits = hits[m]

    if pT.size == 0:
        plt.close(fig)
        return

    h = ax.hist2d(
        hits, pT,
        bins=[hits_bins, pT_bins],
        range=[hits_lim, pT_lim] if (hits_lim is not None and pT_lim is not None) else None,
        cmap="viridis",
    )
    cb = fig.colorbar(h[3], ax=ax)
    cb.set_label("Counts", fontsize=14)

    ax.set_xlabel("Velocity resolution", fontsize=20)
    ax.set_ylabel(r"Mass [GeV]", fontsize=20)
    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

    ax.text(0.02, 0.98, "Muon Collider", ha="left", va="top",
            transform=ax.transAxes, fontsize=20, fontweight="bold", style="italic", color="white")
    ax.text(0.02, 0.90, f"muons, {option}, {window}", ha="left", va="top",
            transform=ax.transAxes, fontsize=15, color="white")
    ax.text(0.02, 0.83, "MuColl_v1", ha="left", va="top",
            transform=ax.transAxes, fontsize=15, color="white")

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def plot_feature_with_mass_overlay(feature, n_bins,
                                   x_lim=None,
                                   m_range=(450.0, 500.0),
                                   req="vb",
                                   overlay_label=None):
    

    d = stats[window][req][option]

    feat = np.asarray(d[feature], dtype=float)
    mass = np.asarray(d["mass"], dtype=float)

    n = min(feat.size, mass.size)
    feat = feat[:n]
    mass = mass[:n]

    base = np.isfinite(feat)
    if x_lim is not None:
        base &= (feat >= x_lim[0]) & (feat <= x_lim[1])

    mmin, mmax = m_range
    sel_mass = base & np.isfinite(mass) & (mass >= mmin) & (mass <= mmax)

    feat_all = feat[base]
    feat_m   = feat[sel_mass]

    N_all = feat_all.size
    N_m   = feat_m.size

    fig, ax = plt.subplots()

    if N_all > 0:
        w_all = np.full_like(feat_all, 100.0 / N_all, dtype=float)
        ax.hist(feat_all, bins=n_bins, weights=w_all,
                histtype="step", linewidth=1.5,
                label="All tracks")

    if N_m > 0:
        w_m = np.full_like(feat_m, 100.0 / N_m, dtype=float)
        frac = 100.0 * N_m / N_all
        lab = overlay_label or f"{mmin:g} < m < {mmax:g} ({frac:.1f}%)"
        ax.hist(feat_m, bins=n_bins, weights=w_m,
                histtype="step", linewidth=2.0,
                label=lab)

    ax.set_xlabel(labels.get(feature, feature), fontsize=20)
    ax.set_ylabel("Normalized counts", fontsize=20)
    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    if x_lim is not None:
        ax.set_xlim(x_lim[0], x_lim[1])

    ax.text(0.02, 0.98, "Muon Collider", ha="left", va="top",
            transform=ax.transAxes, fontsize=20,
            fontweight="bold", style="italic")
    ax.text(0.02, 0.90, f"muons, {option}, {window}", ha="left", va="top",
            transform=ax.transAxes, fontsize=15)
    ax.text(0.02, 0.83, "MuColl_v1", ha="left", va="top",
            transform=ax.transAxes, fontsize=15)

    ax.legend(frameon=False, fontsize=12, loc="upper right")

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    print(f"[overlay] {feature}: N_all={N_all}, N_mass={N_m}, frac={N_m/N_all:.3f}")

def plot_mass_vs_one_minus_beta(req="vb",
                                m_lim=(0, 1000),
                                one_minus_beta_lim=(0, 0.02)):

    d = stats[window][req][option]

    mass = np.asarray(d["mass"], dtype=float)
    velo = np.asarray(d["velo"], dtype=float)

    n = min(mass.size, velo.size)
    mass = mass[:n]
    velo = velo[:n]

    beta = velo / speedoflight
    one_minus_beta = 1.0 - beta

    m = np.isfinite(mass) & np.isfinite(one_minus_beta)
    m &= (mass >= m_lim[0]) & (mass <= m_lim[1])
    m &= (one_minus_beta >= one_minus_beta_lim[0]) & (one_minus_beta <= one_minus_beta_lim[1])

    mass = mass[m]
    one_minus_beta = one_minus_beta[m]

    fig, ax = plt.subplots()
    ax.scatter(one_minus_beta, mass, s=4, alpha=0.3)

    ax.set_xlabel(r"$1 - \beta$", fontsize=20)
    ax.set_ylabel(r"$m_{\mathrm{reco}}$ [GeV]", fontsize=20)
    ax.set_xlim(one_minus_beta_lim)
    ax.set_ylim(m_lim)

    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

    ax.text(0.02, 0.98, "Muon Collider",
            ha="left", va="top",
            transform=ax.transAxes,
            fontsize=20, fontweight="bold", style="italic")
    ax.text(0.02, 0.90, f"{req} tracks, {option}, {window}",
            ha="left", va="top",
            transform=ax.transAxes,
            fontsize=15)
    ax.text(0.02, 0.83, "MuColl_v1",
            ha="left", va="top",
            transform=ax.transAxes,
            fontsize=15)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    print(f"[m vs 1-beta] plotted {mass.size} points")

def plot_m_over_p_vs_one_minus_beta(req="vb",
                                    m_over_p_lim=(0, 0.5),
                                    one_minus_beta_lim=(0, 0.02)):
    d = stats[window][req][option]

    mass = np.asarray(d["mass"], dtype=float)
    velo = np.asarray(d["velo"], dtype=float)
    pT   = np.asarray(d["pT"], dtype=float)

    p = pT

    n = min(mass.size, velo.size, p.size)
    mass, velo, p = mass[:n], velo[:n], p[:n]

    beta = velo / speedoflight
    one_minus_beta = 1.0 - beta

    m_over_p = mass / p

    m = np.isfinite(m_over_p) & np.isfinite(one_minus_beta)
    m &= (m_over_p >= m_over_p_lim[0]) & (m_over_p <= m_over_p_lim[1])
    m &= (one_minus_beta >= one_minus_beta_lim[0]) & (one_minus_beta <= one_minus_beta_lim[1])

    m_over_p = m_over_p[m]
    one_minus_beta = one_minus_beta[m]

    fig, ax = plt.subplots()
    ax.scatter(one_minus_beta, m_over_p, s=4, alpha=0.3)
    
    x = np.linspace(one_minus_beta_lim[0], one_minus_beta_lim[1], 400)
    y = np.sqrt(2.0 * x)

    ax.plot(x, y, color="black", linewidth=2.0,
            label=r"$\sqrt{2(1-\beta)}$")


    ax.set_xlabel(r"$1 - \beta$", fontsize=20)
    ax.set_ylabel(r"$m_{\mathrm{reco}} / p_{\mathrm{reco}}$", fontsize=20)
    ax.set_xlim(one_minus_beta_lim)
    ax.set_ylim(m_over_p_lim)

    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

    ax.text(0.02, 0.98, "Muon Collider",
            ha="left", va="top",
            transform=ax.transAxes,
            fontsize=20, fontweight="bold", style="italic")
    ax.text(0.02, 0.90, f"{req} tracks, {option}, {window}",
            ha="left", va="top",
            transform=ax.transAxes,
            fontsize=15)
    ax.text(0.02, 0.83, "MuColl_v1",
            ha="left", va="top",
            transform=ax.transAxes,
            fontsize=15)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    print(f"[m/p vs 1-beta] plotted {m_over_p.size} points")

def plot_track_display(rec, pdf):
    times = np.asarray(rec["times"], dtype=float)
    r = np.asarray(rec["r"], dtype=float)
    xyz = np.asarray(rec["xyz"], dtype=float) 
    s_unc = np.asarray(rec["spatial_unc"], dtype=float)

    m = np.isfinite(times) & np.isfinite(r)
    if xyz.ndim == 2 and xyz.shape[0] == times.size:
        m &= np.isfinite(xyz).all(axis=1)
    if s_unc.size == times.size:
        m &= np.isfinite(s_unc) & (s_unc > 0)

    times, r = times[m], r[m]
    xyz = xyz[m] if (xyz.ndim == 2 and xyz.shape[0] == m.size) else xyz
    s_unc = s_unc[m] if (s_unc.size == m.size) else None

    if times.size < 3:
        return

    order = np.argsort(times)
    times = times[order]
    r = r[order]
    if xyz.ndim == 2 and xyz.shape[0] == times.size:
        xyz = xyz[order]
    if s_unc is not None and s_unc.size == times.size:
        s_unc = s_unc[order]

    fig = plt.figure()
 
    ax = fig.add_subplot()

    if s_unc is not None:
        ax.errorbar(times, r, yerr=s_unc, fmt="o", markersize=3, linewidth=1, alpha=0.8)
    else:
        ax.scatter(times, r, s=10, alpha=0.8)

    v_fit = rec.get("v_fit", np.nan)
    b_fit = rec.get("b_fit", np.nan)
    if np.isfinite(v_fit) and np.isfinite(b_fit):
        tline = np.linspace(times.min(), times.max(), 200)
        rline = v_fit * tline + b_fit
        ax.plot(tline, rline, linewidth=2)

    ax.set_xlabel("Corrected time [ns]")
    ax.set_ylabel("r = sqrt(x^2+y^2+z^2) [mm]")
    ax.set_title("r vs time with velocity fit")

    m_reco = rec.get("m_reco", np.nan)
    pT = rec.get("pT", np.nan)
    p = rec.get("p", np.nan)
    beta = rec.get("beta", np.nan)
    true_beta = rec.get("true_beta", np.nan)
    v_err = rec.get("v_err", np.nan)

    vb = rec.get("vb_hits", np.nan)
    ib = rec.get("ib_hits", np.nan)
    ob = rec.get("ob_hits", np.nan)

    text = (
        f"m_reco = {m_reco:.2f} GeV\n"
        f"pT = {pT:.2f} GeV,  p = {p:.2f} GeV\n"
        f"beta = {beta:.6f}, {true_beta:.6f} true\n"
        f"v_fit = {v_fit:.3f} mm/ns  (err {v_err:.3f})\n"
        f"hits: vb={vb}, ib={ib}, ob={ob}\n"
        f"{rec.get('option','')} / {rec.get('window','')}"
    )
    fig.suptitle(text, fontsize=10)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def plot_track_rt_display(rec, pdf):
    reco_t = np.asarray(rec["reco_t"], dtype=float)
    reco_r = np.asarray(rec["reco_r"], dtype=float)
    reco_s = np.asarray(rec["reco_unc"], dtype=float)

    sim_t  = np.asarray(rec["sim_t"], dtype=float)
    sim_r  = np.asarray(rec["sim_r"], dtype=float)

    fig, ax = plt.subplots(figsize=(7.5, 5.5))

    m = np.isfinite(sim_t) & np.isfinite(sim_r)

    mu = np.asarray(rec.get("hit_is_muon", []), dtype=object)
    if mu.size != reco_t.size:
        mu = np.ones_like(reco_t, dtype=bool)

    mu_mask = (mu == True)
    nonmu_mask = (mu == False)

    if np.any(mu_mask):
        ax.errorbar(
            reco_t[mu_mask], reco_r[mu_mask], yerr=reco_s[mu_mask],
            fmt="o", markersize=7,
            label="Reco hits (muon)",
            alpha=0.8
        )

    if np.any(nonmu_mask):
        ax.errorbar(
            reco_t[nonmu_mask], reco_r[nonmu_mask], yerr=reco_s[nonmu_mask],
            fmt="o", markersize=7,
            label="Reco hits (non-muon)",
            alpha=0.8
        )

    if np.any(m):
        ax.scatter(
            sim_t[m], sim_r[m],
            s=50,
            color="tab:red",
            marker="x",
            label="Sim hits"
        )

    v_fit = rec.get("v_fit", np.nan)
    b_fit = rec.get("b_fit", np.nan)
    if np.isfinite(v_fit) and np.isfinite(b_fit):
        tline = np.linspace(np.min(reco_t), np.max(reco_t), 200)
        ax.plot(tline, v_fit*tline + b_fit,
                color="black", linewidth=2,
                label="Velocity fit")

    ax.set_xlabel("Time [ns]", fontsize=14)
    ax.set_ylabel(r"$r_{xyz}$ [mm]", fontsize=14)

    ax.legend(frameon=False)

    ax.set_title(
    f"m_reco = {rec['m_reco']:.1f} GeV   beta = {rec['beta']:.6f}\n"
    f"track chi2/ndf = {rec['chi2ndf']:.6f}",
    fontsize=12
    )


    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def plot_mass_logx_variable_bins(req="vb",
                                 m_range=(1e-1, 5e4),
                                 n_bins=70,
                                 normalize=False,
                                 floor=0.1,
                                 use_log_bins=True,
                                 debug=False):
    d = stats[window][req][option]
    m0 = np.asarray(d["mass"], dtype=float)

    if m0.size == 0:
        return

    finite = np.isfinite(m0)
    m = m0[finite]

    if m.size == 0:
        return

    if debug:
        print(f"[mass debug] N total={m0.size}, finite={finite.sum()}, "
              f"nan/inf={(~finite).sum()}, "
              f"==0={np.sum(finite & (m0 == 0.0))}, <0={np.sum(finite & (m0 < 0.0))}, "
              f"min_finite={np.min(m):.6g}")

    m = np.clip(m, floor, None)

    lo, hi = m_range
    m = m[(m >= lo) & (m <= hi)]
    if m.size == 0:
        return

    if use_log_bins:
        edges = np.logspace(np.log10(lo), np.log10(hi), n_bins + 1)
    else:
        edges = np.linspace(lo, hi, n_bins + 1)

    fig, ax = plt.subplots()

    if normalize:
        weights = np.full_like(m, 100.0 / m.size, dtype=float)
        ax.hist(m, bins=edges, weights=weights, histtype="step", linewidth=1.8)
        ax.set_ylabel("Normalized counts [%]", fontsize=20)
    else:
        ax.hist(m, bins=edges, histtype="step", linewidth=1.8)
        ax.set_ylabel("Counts", fontsize=20)

    ax.set_xscale("log")
    ax.set_xlabel("Reconstructed mass [GeV]", fontsize=20)

    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

    ax.text(0.02, 0.98, "Muon Collider", ha="left", va="top",
            transform=ax.transAxes, fontsize=20, fontweight="bold", style="italic")
    ax.text(0.02, 0.90, f"muons, {option}, {window}", ha="left", va="top",
            transform=ax.transAxes, fontsize=15)
    ax.text(0.02, 0.83, "MuColl_v1", ha="left", va="top",
            transform=ax.transAxes, fontsize=15)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def plot_leading_subleading_masses(leading, subleading, n_bins=60, x_lim=(0, 4500)):
    fig, ax = plt.subplots()

    lead = np.asarray(leading, dtype=float)
    sub  = np.asarray(subleading, dtype=float)

    lead = lead[np.isfinite(lead)]
    sub  = sub[np.isfinite(sub)]

    if x_lim is not None:
        lead = lead[(lead >= x_lim[0]) & (lead <= x_lim[1])]
        sub  = sub[(sub >= x_lim[0]) & (sub <= x_lim[1])]

    if lead.size > 0:
        w = np.full_like(lead, 100.0 / lead.size, dtype=float)
        ax.hist(lead, bins=n_bins, weights=w,
                histtype="step", linewidth=2.0,
                label="Leading mass")

    if sub.size > 0:
        w = np.full_like(sub, 100.0 / sub.size, dtype=float)
        ax.hist(sub, bins=n_bins, weights=w,
                histtype="step", linewidth=2.0,
                label="Subleading mass")

    ax.set_xlabel("Reconstructed mass [GeV]", fontsize=20)
    ax.set_ylabel("Normalized counts [%]", fontsize=20)
    if x_lim is not None:
        ax.set_xlim(x_lim)

    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

    ax.text(0.02, 0.98, "Muon Collider",
            ha="left", va="top",
            transform=ax.transAxes,
            fontsize=20, fontweight="bold", style="italic")
    ax.text(0.02, 0.90, f"muons, {option}, {window}",
            ha="left", va="top",
            transform=ax.transAxes,
            fontsize=15)
    ax.text(0.02, 0.83, "MuColl_v1",
            ha="left", va="top",
            transform=ax.transAxes,
            fontsize=15)

    ax.legend(frameon=False, fontsize=12)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


with PdfPages(plot_path) as pdf:
    for window in windows:
        for option in bib_options:
            pT_bins = 10
            hit_bins = 10
            velo_bins = 10
            plot_feature("velo", 20, x_lim=(290,300))
            plot_feature("true_pT", 20)
            plot_feature("true_eta", 20)
            plot_feature("true_mass", 20)
            plot_feature("mass", 20, x_lim=(0,4500))
            plot_pt_vs_hits(pT_lim=(400,600), hits_lim=(-0.5,0.5))

            plot_feature("pT_res", 20, x_lim=(-1, 1))
            plot_feature("velo_res", 20, x_lim=(-0.5, 0.2))
            plot_feature_with_mass_overlay("velo", 20, x_lim=(290,300), m_range=(400,600), req="vb")
            plot_feature_with_mass_overlay("pT", 20, x_lim=(0,10000), m_range=(400,600), req="vb")
            plot_feature_with_mass_overlay("true_eta", 20, x_lim=(-2, 2), m_range=(400,600), req="vb")
            plot_feature_with_mass_overlay("v_fit_err", 20, x_lim=(0,4), m_range=(400,600), req="ob")
            plot_feature_with_mass_overlay("full_p", 20, x_lim=(3000,7000), m_range=(400,600), req="ob")
            plot_feature_with_mass_overlay("time_span", 20, m_range=(400,600), req="vb")
            plot_feature_with_mass_overlay("radial_span", 20, m_range=(400,600), req="vb")
            plot_mass_vs_one_minus_beta(
                req="vb",
                m_lim=(0, 1000),
                one_minus_beta_lim=(0, 0.015)
            )
            plot_m_over_p_vs_one_minus_beta(
                req="vb",
                m_over_p_lim=(0, 0.4),
                one_minus_beta_lim=(0, 0.015)
            )
            plot_mass_logx_variable_bins(req="vb",m_range=(0.1, 10e3), n_bins=70, normalize=True)
            plot_leading_subleading_masses(
                stats[window]["vb"][option]["leading_mass"],
                stats[window]["vb"][option]["subleading_mass"],
                n_bins=60,
                # x_lim=(0, 1000)
            )



print(f"Saved plots to {plot_path}") 

if SAVE_TRACK_DISPLAYS and len(track_displays) > 0:
    with PdfPages(TRACK_DISPLAY_PDF) as pdf2:
        for rec in track_displays:
            plot_track_rt_display(rec, pdf2)
    print(f"Saved {len(track_displays)} track displays to {TRACK_DISPLAY_PDF}")
else:
    print("No track displays saved.")
