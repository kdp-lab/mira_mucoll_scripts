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

dir = "/ospool/uc-shared/project/futurecolliders/wandriscok/reco/nu_background/"
windows = ["loose"]
bib_options = ["bib"]
CACHE = pathlib.Path("cache/nu_bkg_stats.pkl")
plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/backgrounds/nu_tracks.pdf"
n_files = 5
Bfield = 3.57
speedoflight = 299792458/1000000  # mm/ns

parser = argparse.ArgumentParser()
parser.add_argument("--rebuild", action="store_true")
args = parser.parse_args()
rebuild = args.rebuild

guess_velo = 180

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
        window:
        {"bib": {
            "pT": [],
            "hits": [],
            "velo": []},
        } for window in windows
        } 
    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
    
    for window in windows:
        print(f"Analyzing {window} window...")
        for option in bib_options:
            print(f"Analyzing {option}...")
            for ifile in tqdm(range(n_files)):
                file_name = f"nu_background_reco{ifile}.slcio"
                file_path = os.path.join(dir, window, option, file_name)
                if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
                    print(f"couldn't open {file_path}")
                    continue
                reader.open(file_path)
                for event in reader:
                    all_collections = event.getCollectionNames() 
                    track_collection = event.getCollection("SiTracks") if "SiTracks" in all_collections else None 
                    if not track_collection:
                        print("issue 1")
                        continue
                    test_hit_coll = event.getCollection("VXDBarrelHits")
                    if test_hit_coll is None:
                        continue
                    encoding = test_hit_coll.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                    decoder = UTIL.BitField64(encoding)
                    
                    for itrack, track in enumerate(track_collection):
                        track_hits = track.getTrackerHits()
                        
                        reco_pT = 0.3 * Bfield / fabs(track.getOmega() * 1000.)

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
                        total_hits = vb_hits + ib_hits + ob_hits
                        stats[window][option]["pT"].append(reco_pT)
                        stats[window][option]["velo"].append(v_fit)
                        stats[window][option]["hits"].append(total_hits)

    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump(stats, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Writing cache to {CACHE}")
    print("Saved cache successfully.")

print(stats["loose"]["bib"]["pT"])

def plot_feature(feature, n_bins, x_lim=None):
    fig, ax = plt.subplots()
    feature_arr = np.asarray(stats[window][option][feature])
    if x_lim:
        mask = (feature_arr >= x_lim[0]) & (feature_arr <= x_lim[1])
        feature_arr = feature_arr[mask]
    weights = np.full_like(feature_arr, 100.0/feature_arr.size, dtype=float)
    ax.hist(feature_arr, bins=n_bins, weights=weights, histtype="step", facecolor="none")
    ax.set_xlabel(feature, fontsize=20)
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
        f"BIB, {option}, {window}",
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

with PdfPages(plot_path) as pdf:
    for window in windows:
        for option in bib_options:
            pT_bins = 10
            hit_bins = 10
            velo_bins = 10
            plot_feature("pT", pT_bins, x_lim=(0,1000))
            plot_feature("hits", hit_bins)
            plot_feature("velo", velo_bins)
print(f"Saved plots to {plot_path}") 