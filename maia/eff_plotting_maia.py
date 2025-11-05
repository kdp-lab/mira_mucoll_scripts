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

import pyLCIO
from pyLCIO import UTIL, EVENT
import ROOT

hit_collections = ["VBTrackerHitsConed", "VETrackerHitsConed", "ITBarrelHitsConed", "ITEndcapHitsConed", "OTBarrelHitsConed", "OTEndcapHitsConed"]
sim_collections = ["VertexBarrelCollectionConed", "VertexEndcapCollectionConed", "InnerTrackerBarrelCollectionConed", "InnerTrackerEndcapCollectionConed", "OuterTrackerBarrelConed", "OuterTrackerEndcapCollectionConed"]
rel_collections = ["VBTrackerHitsRelationsConed", "VETrackerHitsRelationsConed", "ITBarrelHitsRelationsConed", "ITEndcapHitsRelationsConed", "OTBarrelHitsRelationsConed", "OTEndcapHitsRelationsConed"]

nominal_dir = "/ospool/uc-shared/project/futurecolliders/miralittmann/MAIA/reco/nominal/"
window_to_dir = {"nominal": nominal_dir,
                 #"medium": medium_dir,
                 #"loose": loose_dir}
                }
n_files = 2500

CACHE = pathlib.Path("cache/maia-eff.pkl")

plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/maia_efficiency.pdf"

sample_to_mass = {
    "1000_10": 1.0,
    "2500_10": 2.5,
    "3000_10": 3.0,
    "3500_10": 3.5,
    "4000_10": 4.0,
    "4500_10": 4.5
}
mass_list = [1.0, 2.5, 3.0, 3.5, 4.0, 4.5]
samples = ["1000_10", "2500_10", "3000_10", "3500_10", "4000_10", "4500_10"]
bib_options = ["bib", "nobib"]
windows = ["nominal"]

stau_ids = {1000015, -1000015, 2000015, -2000015} 
Bfield = 5
chi2_cut = 3
nhits_cut = 4
speedoflight = 299792458/1000000  # mm/ns
guess_velo = 180

def binom_se(p, n):
    # Gaussian SE on a proportion
    if n <= 0:
        return float('nan')
    return math.sqrt(p*(1-p)/n)

parser = argparse.ArgumentParser()
parser.add_argument("--rebuild", action="store_true")
args = parser.parse_args()
rebuild = args.rebuild

def acceptance(mcp):
    mcp_stau_vertex_r = np.sqrt(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2)
    travel_dist = np.sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2 + mcp.getEndpoint()[2]**2)

    mcp_stau_endpoint_r = np.sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2)
    mcp_stau_momentum = mcp.getMomentum() 
    mcp_stau_tlv = ROOT.TLorentzVector()
    mcp_stau_tlv.SetPxPyPzE(mcp_stau_momentum[0], mcp_stau_momentum[1], mcp_stau_momentum[2], mcp.getEnergy())  

#    if abs(mcp.getPDG()) not in stau_ids or travel_dist==0 or mcp_stau_vertex_r>553.0 or mcp_stau_endpoint_r < 102.0 or abs(mcp_stau_tlv.Eta())>0.8:

    if abs(mcp.getPDG()) not in stau_ids or travel_dist==0 or mcp_stau_vertex_r>553.0 or mcp_stau_endpoint_r < 102.0:
        accepted = False
    else:
        accepted = True
    return(accepted)

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


efficiencies = None
if (not rebuild) and os.path.exists(CACHE):
    with open(CACHE, "rb") as f:
        print("Loading in cached arrays...")
        efficiencies = pickle.load(f)

if efficiencies is None:

    efficiencies = {
        window: {
            option: {
                sample: {
                    "all_tracks": 0,
                    "vb_tracks": 0,
                    "ib_tracks": 0,
                    "ob_tracks": 0,
                    "accepted_staus": 0,

                    "eta": [],

                    "pT_res": {
                        "vb_hits": [],
                        "ib_hits": [],
                        "ob_hits": []
                    },

                    "velo_res": {
                        "vb_hits": [],
                        "ib_hits": [],
                        "ob_hits": [] 
                    }
                } for sample in samples
            } for option in bib_options
        } for window in windows 
    }

    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

    for window in windows:
        print(f"Analyzing {window} window...")
        for option in bib_options:
            print(f"         {option}...")
            for sample in samples: 
                print(f"                      {sample}...")
                save_info = efficiencies[window][option][sample]
                reco_dir = os.path.join(window_to_dir[window], option, sample)

                for ifile in tqdm(range(n_files)):
                    file_name = f"{sample}_reco{ifile}.slcio"
                    file_path = os.path.join(reco_dir,file_name) 
                    print(file_path)
                    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
                        print("error 0")
                        continue
                    reader.open(file_path)
                    for event in reader:
                        all_collections = event.getCollectionNames() 
                        mcp_collection = event.getCollection("MCParticle") if "MCParticle" in all_collections else None
                        track_collection = event.getCollection("SiTracks") if "SiTracks" in all_collections else None
                        track_relation_collection = event.getCollection("MCParticle_SiTracks") if "MCParticle_SiTracks" in all_collections else None

                        if not (mcp_collection and track_collection and track_relation_collection):
                            print("issue 1")
                            continue

                        test_hit_coll = event.getCollection("VBTrackerHitsConed")
                        if test_hit_coll is None:
                            continue
                        encoding = test_hit_coll.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                        decoder = UTIL.BitField64(encoding)

                        for mcp in mcp_collection:
                            accepted = acceptance(mcp)

                            if accepted == True:
                                save_info["accepted_staus"] += 1
                    
                        nav = UTIL.LCRelationNavigator(track_relation_collection)
                        
                        for itrack, track in enumerate(track_collection):
                            track_mcps = nav.getRelatedFromObjects(track)
                            track_hits = track.getTrackerHits()
                            
                            pT = 0.3 * Bfield / fabs(track.getOmega() * 1000.)

                            truth_stau = None
                            if track_mcps:
                                for mcp in track_mcps:
                                    accepted = acceptance(mcp)
                                    if abs(mcp.getPDG()) in stau_ids and accepted == True:
                                        truth_stau = mcp
                                        save_info["all_tracks"] += 1
                                    else:
                                        continue
                            
                            if truth_stau is not None:
                                mom = truth_stau.getMomentum()
                                tlv = ROOT.TLorentzVector()
                                tlv.SetPxPyPzE(mom[0], mom[1], mom[2], truth_stau.getEnergy())
                                pT_truth = tlv.Perp()    
                                beta_truth = tlv.Beta()
                                velo_truth = beta_truth * speedoflight                    
                            
                                track_chi2 = track.getChi2() / track.getNdf()
                                vb_hits = 0
                                ib_hits = 0
                                ob_hits = 0

                                if track_chi2 > chi2_cut:
                                    continue

                                track_times = []
                                track_pos = []
                                spatial_unc = []

                                for hit in track_hits:
                                    decoder.setValue(int(hit.getCellID0()))
                                    system = decoder["system"].value()
                                    layer = decoder["layer"].value()
                                    if system in (1,2):
                                        vb_hits += 1
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
                                
                                pT_res = (pT - pT_truth) / pT_truth   
                                velo_res = (v_fit - velo_truth) / velo_truth     
    
                                if vb_hits >= 3:
                                    save_info["vb_tracks"] += 1
                                    save_info["pT_res"]["vb_hits"].append(pT_res)
                                    save_info["velo_res"]["vb_hits"].append(velo_res)
                                if vb_hits >= 3 and ib_hits >= 2:
                                    save_info["ib_tracks"] += 1
                                    save_info["pT_res"]["ib_hits"].append(pT_res)
                                    save_info["velo_res"]["ib_hits"].append(velo_res)
                                if vb_hits >= 3 and ib_hits >= 2 and ob_hits >= 2:
                                    save_info["ob_tracks"] += 1
                                    save_info["pT_res"]["ob_hits"].append(pT_res)
                                    save_info["velo_res"]["ob_hits"].append(velo_res)

    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump(efficiencies, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Writing cache to {CACHE}")
    print("Saved cache successfully.")

window_to_name = {"nominal": "Nominal", 
                  "medium": "Medium", 
                  "loose": "Loose"}

def get_eff_and_err(data, which="all_tracks"):
    k = data[which]
    n = data["accepted_staus"]
    if n == 0:
        return float('nan'), float('nan')
    else:
        p = k/n
        e = binom_se(p, n)
    return 100*p, 100*e

def collect_residuals(effs, type, requirement_key, window, option):
    # requirement_key in {"vb_hits","ib_hits","ob_hits"}
    vals = []
    for sample in samples:
        vals.extend(effs[window][option][sample][type][requirement_key])
    return np.array(vals)

panels = [
    ("vb_hits", "VB requirement (≥3 VB hits)"),
    ("ib_hits", "VB+IB requirement (≥3 VB & ≥2 IB hits)"),
    ("ob_hits", "VB+IB+OB requirement (≥3 VB, ≥2 IB, ≥2 OB hits)")
]

panels0 = [
    ("vb_tracks", "≥3 VB hits"),
    ("ib_tracks", "≥3 VB, ≥2 IB hits"),
    ("ob_tracks", "≥3 VB, ≥2 IB, ≥2 OB hits")
]

print(efficiencies["nominal"]["bib"]["4500_10"]["accepted_staus"])

plt.style.use("seaborn-v0_8-colorblind")
with PdfPages(plot_path) as pdf:
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    masses = [sample_to_mass[s] for s in samples]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color_map = {"nominal": colors[0], "medium": colors[1], "loose": colors[2]}

    for ax, (req_key, title) in zip(axes, panels0):
        for window in windows: 
            trackeff_bib, trackerr_bib = [], []
            trackeff_nobib, trackerr_nobib = [], []
            
            for sample in samples: 
                d_b = efficiencies[window]["bib"][sample]
                p_b, e_b = get_eff_and_err(d_b, req_key)
                trackeff_bib.append(p_b); trackerr_bib.append(e_b)

                d_nb = efficiencies[window]["nobib"][sample]
                p_nb, e_nb = get_eff_and_err(d_nb, req_key)
                trackeff_nobib.append(p_nb); trackerr_nobib.append(e_nb)

            ax.errorbar(
                masses, trackeff_bib, yerr=trackerr_bib,
                linewidth=2, capsize=2, color=color_map[window],
                label=f"{window_to_name[window]}, 10% BIB"
            )
            ax.errorbar(
                masses, trackeff_nobib, yerr=trackerr_nobib,
                linewidth=2, capsize=2, color=color_map[window],
                label=f"{window_to_name[window]}, no BIB", fmt=":", alpha=0.7
            )

            ax.set_title(title, fontsize=20)

    for ax in axes:
        ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
        ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    
        ax.set_xlabel("Stau mass [TeV]", fontsize=20)
        axes[0].set_ylabel("Track reconstruction efficiency (%)", fontsize=22)
        ax.grid(True, alpha=0.2)
        ax.set_ylim(0,100)

    handles, labels = [], [] 
    for ax in axes:
        h,l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)

    by_label = OrderedDict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), fontsize=18, ncol=3, handlelength=1, handletextpad=0.3, columnspacing=0.5, loc="upper center",bbox_to_anchor=(0.5, 1.11))
    
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)






    bins = np.linspace(-1.0, 1.0, 81)
    fig1, axes1 = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    show_nobib_outline = True

    for ax, (req_key, title) in zip(axes1, panels):
        for w in windows:
            # data_bib = collect_residuals(efficiencies, "pT_res", req_key, w, "bib")
            # if data_bib.size > 0:
            #     weights = np.full_like(data_bib, 100.0 / data_bib.size, dtype=float)
            #     ax.hist(
            #     data_bib,
            #     bins=bins,
            #     weights=weights,
            #     histtype="step",
            #     facecolor="none",
            #     # alpha=0.35, 
            #     color=color_map[w],
            #     label=f"{window_to_name[w]} (10% BIB)",
            #     edgecolor=color_map[w],
            #     linewidth=1.8
            # )
        
            if show_nobib_outline:
                data_nobib = collect_residuals(efficiencies, "pT_res", req_key, w, "nobib")
                if data_nobib.size > 0:
                    weights = np.full_like(data_nobib, 100.0 / data_nobib.size, dtype=float)
                    ax.hist(
                        data_nobib,
                        bins=bins,
                        weights=weights,
                        histtype="step",
                        linewidth=1.8,
                        color=color_map[w],
                        label=f"{window_to_name[w]} (no BIB)",
                        linestyle="--"
                    )
        ax.axvline(0, lw=1, ls=":", color="k", alpha=0.5)
        ax.set_title(title, fontsize=15)
        ax.grid(True, alpha=0.2)

    for ax in axes1:
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r"$(p_T^{\mathrm{reco}} - p_T^{\mathrm{truth}})/p_T^{\mathrm{truth}}$", fontsize=20)
        axes1[0].set_ylabel("Tracks per bin (%) - $\Delta$=0.025", fontsize=20) 
        ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
        ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    
    handles, labels = [], []
    for ax in axes1:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)

    by_label = OrderedDict(zip(labels, handles))
    fig1.legend(
        by_label.values(), by_label.keys(),
        ncol=3, fontsize=18,
        loc="upper center", bbox_to_anchor=(0.5, 1.08),
        handlelength=1.5, handletextpad=0.5, columnspacing=1.0
    )


    fig1.tight_layout(rect=[0, 0, 1, 0.98])
    pdf.savefig(fig1, bbox_inches="tight")
    plt.close(fig1)






    bins = np.linspace(-1.0, 1.0, 81)
    fig1, axes1 = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    show_nobib_outline = True

    for ax, (req_key, title) in zip(axes1, panels):
        for w in windows:
        #     data_bib = collect_residuals(efficiencies, "velo_res", req_key, w, "bib")
        #     if len(data_bib) == 0:
        #         continue
        #     weights = np.full_like(data_bib, 100.0/data_bib.size, dtype=float)
        #     ax.hist(
        #     data_bib,
        #     bins=bins,
        #     weights=weights,
        #     histtype="step",
        #     facecolor="none",
        #     color=color_map[w],
        #     label=f"{window_to_name[w]} (10% BIB)",
        #     edgecolor=color_map[w],
        #     linewidth=1.8
        # )
        
            if show_nobib_outline:
                data_nobib = collect_residuals(efficiencies, "velo_res", req_key, w, "nobib")
                if len(data_nobib) > 0:
                    weights = np.full_like(data_nobib, 100.0/data_nobib.size, dtype=float)
                    ax.hist(
                        data_nobib,
                        bins=bins,
                        weights=weights,
                        histtype="step",
                        linewidth=1.8,
                        color=color_map[w],
                        label=f"{window_to_name[w]} (no BIB)",
                        linestyle="--"
                    )
        ax.axvline(0, lw=1, ls=":", color="k", alpha=0.5)
        ax.set_title(title, fontsize=15)
        ax.grid(True, alpha=0.2)

    for ax in axes1:
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r"$(v^{\mathrm{reco}} - v^{\mathrm{truth}})/v^{\mathrm{truth}}$", fontsize=20)
        axes1[0].set_ylabel("Tracks per bin (%) - $\Delta$ = 0.025", fontsize=20) 
        ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
        ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    
    handles, labels = [], []
    for ax in axes1:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)

    by_label = OrderedDict(zip(labels, handles))
    fig1.legend(
        by_label.values(), by_label.keys(),
        ncol=3, fontsize=18,
        loc="upper center", bbox_to_anchor=(0.5, 1.08),
        handlelength=1.5, handletextpad=0.5, columnspacing=1.0
    )


    fig1.tight_layout(rect=[0, 0, 1, 0.98])
    pdf.savefig(fig1, bbox_inches="tight")
    plt.close(fig1)





print(f"Saved plots to {plot_path}")

           
