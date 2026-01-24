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
CACHE = pathlib.Path("cache/mumu_bkg_stats_nominal_nobib.pkl")
plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/backgrounds/mumu_tracks_nominal_nobib.pdf"
print(dirs.items())
n_files = 2500
Bfield = 3.57
speedoflight = 299792458/1000000  # mm/ns
chi2_cut = 3

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

guess_velo = 180

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
                    "mass": []
                },
                "nobib": {
                    "pT": [],
                    "hits": [],
                    "velo": [],
                    "mass": []
                }
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
            for ifile in tqdm(range(n_files)):
                file_name = f"mumu_bkg_reco{ifile}.slcio"
                file_path = os.path.join(dirs[option], window, file_name)
                if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
                    print(f"couldn't open {file_path}")
                    continue
                reader.open(file_path)
                for event in reader:
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

                                vb_hits = 0
                                ib_hits = 0
                                ob_hits = 0
                                
                                track_times = []
                                track_pos = []
                                spatial_unc = []

                                for hit in track_hits:
                                    decoder.setValue(int(hit.getCellID0()))
                                    system = decoder["system"].value()
                                    layer = decoder["layer"].value()
                                    if system in (1,2):
                                        vb_hits += 0.5
                                        spatial_unc.append(0.005)
                                    elif system in (3,4):
                                        ib_hits += 1
                                        spatial_unc.append(0.007)
                                    elif system in (5,6):
                                        ob_hits += 1
                                        spatial_unc.append(0.007)

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

                                v_fit, v_err = reco_velo(linearfunc, track_times, track_pos, spatial_unc) 
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

                                if vb_hits >= 3 and ib_hits >= 2 and ob_hits >=2:
                                    stats[window]["ob"][option]["pT"].append(reco_pT)
                                    stats[window]["ob"][option]["velo"].append(v_fit)
                                    stats[window]["ob"][option]["hits"].append(total_hits)
                                    stats[window]["ob"][option]["mass"].append(m_reco)
                                
                                if vb_hits >= 3 and ib_hits >= 2:
                                    stats[window]["ib"][option]["pT"].append(reco_pT)
                                    stats[window]["ib"][option]["velo"].append(v_fit)
                                    stats[window]["ib"][option]["hits"].append(total_hits)
                                    stats[window]["ib"][option]["mass"].append(m_reco)
                                
                                if vb_hits >= 3:
                                    stats[window]["vb"][option]["pT"].append(reco_pT)
                                    stats[window]["vb"][option]["velo"].append(v_fit)
                                    stats[window]["vb"][option]["hits"].append(total_hits)
                                    stats[window]["vb"][option]["mass"].append(m_reco)
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
          "velo": "Velocity [mm/ns]"}

def plot_feature(feature, n_bins, x_lim=None):
    fig, ax = plt.subplots()
    feature_arr = np.asarray(stats[window][option][feature])
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

def plot_pt_vs_hits(pT_lim=(0, 10000), hits_lim=None, pT_bins=40, hits_bins=10):
    fig, ax = plt.subplots()

    pT = np.asarray(stats[window][option]["pT"], dtype=float)
    hits = np.asarray(stats[window][option]["hits"], dtype=float)

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

    ax.set_xlabel("Number of hits on track", fontsize=20)
    ax.set_ylabel(r"$p_T$ [GeV]", fontsize=20)
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


with PdfPages(plot_path) as pdf:
    for window in windows:
        for option in bib_options:
            pT_bins = 10
            hit_bins = 10
            velo_bins = 10
            plot_feature("pT", pT_bins, x_lim=(0,10000))
            plot_feature("hits", hit_bins)
            plot_feature("velo", velo_bins, x_lim=(150,450))
            plot_pt_vs_hits(pT_lim=(0, 10000), hits_bins=10, pT_bins=10)
print(f"Saved plots to {plot_path}") 