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

hit_collections = ["VXDBarrelHits", "VXDEndcapHits", "ITBarrelHits", "ITEndcapHits", "OTBarrelHits", "OTEndcapHits"]
sim_collections = ["VertexBarrelCollection", "VertexEndcapCollection", "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection", "OuterTrackerBarrel", "OuterTrackerEndcapCollectionConed"]
rel_collections = ["VXDBarrelHitsRelations", "VXDEndcapHitsRelations", "ITBarrelHitsRelations", "ITEndcapHitsRelations", "OTBarrelHitsRelations", "OTEndcapHitsRelations"]

tight_dir = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/efficiency/seeding_10GeV/"
medium_dir = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/p26_p5_p6/seeding_10GeV/"
loose_dir =  "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/loose4/seeding_10GeV/"
window_to_dir = {"tight": tight_dir,
                 "medium": medium_dir,
                 "loose": loose_dir}
n_files = 2500

CACHE = pathlib.Path("cache/eff_from_slcio-seeding_10GeV-halfhitsvxd.pkl")

plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency-seeding_10GeV-halfhitsvxd.pdf"

sample_to_mass = {
    "1000_10": 1.0,
    "1500_10": 1.5,
    "2000_10": 2.0,
    "2500_10": 2.5,
    "3000_10": 3.0,
    "3500_10": 3.5,
    "4000_10": 4.0,
    "4500_10": 4.5
}
mass_list = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
samples = ["1000_10", "1500_10", "2000_10", "2500_10", "3000_10", "3500_10", "4000_10", "4500_10"]
bib_options = ["bib", "nobib"]
windows = ["loose", "medium", "tight"]

stau_ids = {1000015, -1000015, 2000015, -2000015} 
Bfield = 3.57
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

def iter_events(reader): 
    while True:
        try:
            evt = reader.readNextEvent() 
        except Exception as e:
            break
        if evt is None:
            break
        yield evt
def safe_get(event, name):
    try:
        return event.getCollection(name)
    except Exception:
        return None

def acceptance(mcp):
    mcp_stau_vertex_r = np.sqrt(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2)
    travel_dist = np.sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2 + mcp.getEndpoint()[2]**2)

    mcp_stau_endpoint_r = np.sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2)
    mcp_stau_momentum = mcp.getMomentum() 
    mcp_stau_tlv = ROOT.TLorentzVector()
    mcp_stau_tlv.SetPxPyPzE(mcp_stau_momentum[0], mcp_stau_momentum[1], mcp_stau_momentum[2], mcp.getEnergy())  

    if abs(mcp.getPDG()) not in stau_ids or travel_dist==0 or mcp_stau_vertex_r>553.0 or mcp_stau_endpoint_r < 1486.0 or abs(mcp_stau_tlv.Eta()) > 0.8:
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


def get_pcent_stau(track, rel_nav):
    for hit in track.getTrackerHits(): 
        for sim in rel_nav[system_to_relname[system]].getRelatedToObjects(hit):
            mcp = sim.getMCParticle()
            if mcp and abs(mcp.getPDG()) in stau_ids:
                stau_hit_count += 1
    
    track_pcent_stau = (stau_hit_count / len(track.getTrackerHits()))

    return(track_pcent_stau)
        




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

                    "pcent_stau": {
                        "vb_hits": [],
                        "ib_hits": [],
                        "ob_hits": []
                    },

                    # "pT": {"vb_hits": [],
                    #     "ib_hits": [],
                    #     "ob_hits": []}, 

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
                
                if window == "tight":
                    reco_dir = os.path.join(tight_dir, option, "tight", sample)
                else:
                    reco_dir = os.path.join(window_to_dir[window], option, sample)

                for ifile in tqdm(range(n_files)):
                    file_name = f"{sample}_reco{ifile}.slcio"
                    file_path = os.path.join(reco_dir,file_name) 
                    if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
                        print(f"couldn't open {file_path}")
                        continue
                    reader.open(file_path)
                    for event in reader:
                        rel_nav = build_rel_nav(event)
                        all_collections = event.getCollectionNames() 
                        mcp_collection = event.getCollection("MCParticle") if "MCParticle" in all_collections else None
                        track_collection = event.getCollection("SiTracks") if "SiTracks" in all_collections else None # REFITTED
                        track_relation_collection = event.getCollection("MCParticle_SiTracks") if "MCParticle_SiTracks" in all_collections else None # REFITTED

                        if not (mcp_collection and track_collection and track_relation_collection):
                            print("issue 1")
                            continue

                        test_hit_coll = event.getCollection("VXDBarrelHits")
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
                            
                            stau_hit_count = 0
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

                                    # percent stau
                                    for sim in rel_nav[system_to_relname[system]].getRelatedToObjects(hit):
                                        mcp_true = sim.getMCParticle()
                                        if mcp_true and abs(mcp_true.getPDG()) in stau_ids:
                                            stau_hit_count += 1

                                v_fit, v_err = reco_velo(linearfunc, track_times, track_pos, spatial_unc) 
                                
                                pT_res = (pT - pT_truth) / pT_truth   
                                velo_res = (v_fit - velo_truth) / velo_truth     
                                pcent_stau = stau_hit_count / len(track.getTrackerHits())
                                # print(pcent_stau, pT)
    
                                if vb_hits >= 3:
                                    save_info["vb_tracks"] += 1
                                    save_info["pT_res"]["vb_hits"].append(pT_res)
                                    save_info["velo_res"]["vb_hits"].append(velo_res)
                                    # save_info["pT"]["vb_hits"].append(pT)
                                    save_info["pcent_stau"]["vb_hits"].append(pcent_stau)
                                if vb_hits >= 3 and ib_hits >= 2:
                                    save_info["ib_tracks"] += 1
                                    save_info["pT_res"]["ib_hits"].append(pT_res)
                                    save_info["velo_res"]["ib_hits"].append(velo_res)
                                    # save_info["pT"]["ib_hits"].append(pT)
                                    save_info["pcent_stau"]["ib_hits"].append(pcent_stau)
                                if vb_hits >= 3 and ib_hits >= 2 and ob_hits >= 2:
                                    save_info["ob_tracks"] += 1
                                    save_info["pT_res"]["ob_hits"].append(pT_res)
                                    save_info["velo_res"]["ob_hits"].append(velo_res)
                                    # save_info["pT"]["ob_hits"].append(pT)
                                    save_info["pcent_stau"]["ob_hits"].append(pcent_stau)

    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump(efficiencies, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Writing cache to {CACHE}")
    print("Saved cache successfully.")

print(efficiencies["loose"]["bib"].keys())

window_to_name = {"tight": "Nominal", 
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
    ("vb_hits", "≥3 VB hits"),
    ("ib_hits", "≥3 VB & ≥2 IB hits"),
    ("ob_hits", "≥3 VB, ≥2 IB, ≥2 OB hits")
]

panels0 = [
    ("vb_tracks", "≥3 VB hits"),
    ("ib_tracks", "≥3 VB, ≥2 IB hits"),
    ("ob_tracks", "≥3 VB, ≥2 IB, ≥2 OB hits")
]

plt.style.use("seaborn-v0_8-colorblind")
with PdfPages(plot_path) as pdf:
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    masses = [sample_to_mass[s] for s in samples]
    print(masses)
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color_map = {"tight": colors[0], "medium": colors[1], "loose": colors[2]}

    for ax, (req_key, title) in zip(axes, panels0):
        for window in windows: 
            trackeff_bib, trackerr_bib = [], []
            trackeff_nobib, trackerr_nobib = [], []
            
            for sample in samples: 
                print(sample)
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

            axes[0].text(
                0.02, 0.20,
                "Muon Collider",
                ha="left", va="top",
                transform=axes[0].transAxes,
                fontsize=24,
                fontweight="bold",
                style="italic",
            )
            axes[0].text(
                0.02, 0.13,
                "Simulation 10% BIB, $\sqrt{s}$=10 TeV",
                ha="left", va="top",
                transform=axes[0].transAxes,
                fontsize=18
            ) 
            axes[0].text(
                0.02, 0.07,
                "MuColl_v1",
                ha="left", va="top",
                transform=axes[0].transAxes,
                fontsize=18
            )
            
            ax.set_title(title, fontsize=20)

    for ax in axes:
        ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
        ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    
        ax.set_xlabel("Stau mass [TeV]", fontsize=20)
        axes[0].set_ylabel("Track reconstruction efficiency (%)", fontsize=22)
        ax.grid(True, alpha=0.2)
        ax.set_ylim(0,100)

    # axes[0].text( 0.02, 0.02, r"Non-refitted tracks", transform=axes[0].transAxes, ha="left", va="bottom", fontsize=16, fontweight="bold", style="italic" ) 
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






    # bins = np.linspace(-1.0, 1.0, 81)
    # fig1, axes1 = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    # show_nobib_outline = False

    # for ax, (req_key, title) in zip(axes1, panels):
    #     for w in windows:
    #         data_bib = collect_residuals(efficiencies, "pT_res", req_key, w, "bib")
    #         if data_bib.size == 0:
    #             continue
    #         weights = np.full_like(data_bib, 100.0 / data_bib.size, dtype=float)
    #         ax.hist(
    #             data_bib, bins=bins, weights=weights,
    #             histtype="step", linewidth=1.8,
    #             facecolor="none", edgecolor=color_map[w], color=color_map[w],
    #             label=f"{window_to_name[w]} (10% BIB)"
    #         )

    #         if show_nobib_outline:
    #             data_nb = collect_residuals(efficiencies, "pT_res", req_key, w, "nobib")
    #             if data_nb.size > 0:
    #                 weights = np.full_like(data_nb, 100.0 / data_nb.size, dtype=float)
    #                 ax.hist(
    #                     data_nb, bins=bins, weights=weights,
    #                     histtype="step", linewidth=1.8,
    #                     color=color_map[w], linestyle='--'  
    #                 )

    #     ax.axvline(0, lw=1, ls=":", color="k", alpha=0.5)
    #     ax.set_title(title, fontsize=15)
    #     ax.grid(True, alpha=0.2)

    # for ax in axes1:
    #     ax.set_xlim(-1.0, 1.0)
    #     ax.set_ylim(bottom=0)
    #     ax.set_xlabel(r"$(p_T^{\mathrm{reco}} - p_T^{\mathrm{truth}})/p_T^{\mathrm{truth}}$", fontsize=20)
    #     ax.text( 0.02, 0.98, r"$m_{\rm stau}=all$", transform=ax.transAxes, ha="left", va="top", fontsize=16, fontweight="bold", style="italic" ) 

    #     ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    #     ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    # axes1[0].set_ylabel("Tracks per bin (%) - $\\Delta$=0.025", fontsize=20)

    # handles, labels = [], []
    # for ax in axes1:
    #     h, l = ax.get_legend_handles_labels()
    #     handles.extend(h); labels.extend(l)
    # by_label = OrderedDict(zip(labels, handles))
    # fig1.legend(
    #     by_label.values(), by_label.keys(),
    #     ncol=3, fontsize=18,
    #     loc="upper center", bbox_to_anchor=(0.5, 1.08),
    #     handlelength=1.5, handletextpad=0.5, columnspacing=1.0
    # )

    # fig1.tight_layout(rect=[0, 0, 1, 0.98])
    # pdf.savefig(fig1, bbox_inches="tight")
    # plt.close(fig1)




    # bins = np.linspace(-1.0, 1.0, 81)
    # fig1, axes1 = plt.subplots(1, 3, figsize=(18, 6), sharex=True, sharey=True)
    # show_nobib_outline = False

    # for ax, (req_key, title) in zip(axes1, panels):
    #     for w in windows:
    #         data_bib = collect_residuals(efficiencies, "velo_res", req_key, w, "bib")
    #         if len(data_bib) == 0:
    #             continue
    #         weights = np.full_like(data_bib, 100.0/data_bib.size, dtype=float)
    #         ax.hist(
    #         data_bib,
    #         bins=bins,
    #         weights=weights,
    #         histtype="step",
    #         facecolor="none",
    #         color=color_map[w],
    #         label=f"{window_to_name[w]} (10% BIB)",
    #         edgecolor=color_map[w],
    #         linewidth=1.8
    #     )
        
    #         if show_nobib_outline:
    #             data_nobib = collect_residuals(efficiencies, "velo_res", req_key, w, "nobib")
    #             if len(data_nobib) > 0:
    #                 weights = np.full_like(data_nobib, 100.0/data_nobib.size, dtype=float)
    #                 ax.hist(
    #                     data_nobib,
    #                     bins=bins,
    #                     weights=weights,
    #                     histtype="step",
    #                     linewidth=1.8,
    #                     color=color_map[w],
    #                     label=f"{window_to_name[w]} (no BIB)",
    #                     linestyle="--"
    #                 )
    #     ax.axvline(0, lw=1, ls=":", color="k", alpha=0.5)
    #     ax.set_title(title, fontsize=15)
    #     ax.grid(True, alpha=0.2)

    # for ax in axes1:
    #     ax.set_xlim(-0.5, 0.5)
    #     ax.set_ylim(bottom=0)
    #     ax.set_xlabel(r"$(v^{\mathrm{reco}} - v^{\mathrm{truth}})/v^{\mathrm{truth}}$", fontsize=20)
    #     axes1[0].set_ylabel("Tracks per bin (%) - $\Delta$ = 0.025", fontsize=20) 
    #     ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    #     ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    #     ax.text(
    #             0.02, 0.98, r"$m_{\rm stau}=all$",
    #             transform=ax.transAxes,
    #             ha="left", va="top",
    #             fontsize=16, fontweight="bold", style="italic"
    #     )

    
    # handles, labels = [], []
    # for ax in axes1:
    #     h, l = ax.get_legend_handles_labels()
    #     handles.extend(h)
    #     labels.extend(l)

    # by_label = OrderedDict(zip(labels, handles))
    # fig1.legend(
    #     by_label.values(), by_label.keys(),
    #     ncol=3, fontsize=18,
    #     loc="upper center", bbox_to_anchor=(0.5, 1.08),
    #     handlelength=1.5, handletextpad=0.5, columnspacing=1.0
    # )


    # fig1.tight_layout(rect=[0, 0, 1, 0.98])
    # pdf.savefig(fig1, bbox_inches="tight")
    # plt.close(fig1)



    windows_and_names = [("tight", "Nominal"), ("medium", "Medium"), ("loose", "Loose")]

    bins = np.linspace(-1, 1, 21)
    fig, axes = plt.subplots(3, 3, figsize=(18, 15), sharex=True, sharey=True)
    
    for i, (window, row_name) in enumerate(windows_and_names):
        for j,((req_key, title), ax) in enumerate(zip(panels, axes[i])):
            for sample in ["1000_10", "2500_10", "3500_10", "4000_10"]:
                pT = np.asarray(efficiencies[window]["nobib"][sample]["pT_res"][req_key])
                if pT.size==0:
                    continue
                weights = np.full_like(pT, 100.0/pT.size, dtype=float)

                ax.hist(pT, bins=bins, weights=weights, histtype="step", facecolor="none", linewidth=2.5, label=f"{sample_to_mass[sample]} TeV", linestyle="--")
                axes[0,2].text(
                    0.98, 0.98,
                    "Muon Collider",
                    ha="right", va="top",
                    transform=axes[0,2].transAxes,
                    fontsize=24,
                    fontweight="bold",
                    style="italic",
                )
                axes[0,2].text(
                    0.98, 0.90,
                    r"Simulation No BIB, $\sqrt{s}=10\ \mathrm{TeV}$",
                    ha="right", va="top",
                    transform=axes[0,2].transAxes,
                    fontsize=18
                ) 
                axes[0,2].text(
                    0.98, 0.83,
                    "MuColl_v1",
                    ha="right", va="top",
                    transform=axes[0,2].transAxes,
                    fontsize=18
                )
                # ax.set_title(title, fontsize=15)
                ax.grid(True, alpha=0.2)
                ax.grid(True, which="minor", alpha=0.1) 
                ax.minorticks_on()
                ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
                ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

                if i == 2:
                    ax.set_xlabel(r"$(p_T^{\mathrm{reco}} - p_T^{\mathrm{truth}})/p_T^{\mathrm{truth}}$",fontsize=25)
                if j == 0:
                    ax.set_ylabel("Tracks per bin (%)", fontsize=25)
            axes[i, 0].text(-0.2, 0.5, f"{row_name} window",transform=axes[i, 0].transAxes,ha="right", va="center",rotation=90, fontsize=25, fontweight="bold", style="italic")
    
    for j, (_, title) in enumerate(panels):
        axes[0, j].set_title(title, fontsize=20, fontweight="bold", pad=10)

    handles, labels = [], []
    for ax_row in axes:
        for ax in ax_row:
            h, l = ax.get_legend_handles_labels()
            handles.extend(h)
            labels.extend(l)

    by_label = OrderedDict(zip(labels, handles))
    fig.legend(
        by_label.values(), by_label.keys(),
        ncol=4, fontsize=25,
        loc="upper center", bbox_to_anchor=(0.5, 1.02),
        handlelength=1.5, handletextpad=0.5, columnspacing=1.0
    )

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)   


    bins = np.linspace(-1, 1, 41)
    fig, axes = plt.subplots(3, 3, figsize=(18, 15), sharex=True, sharey=True)
    
    for i, (window, row_name) in enumerate(windows_and_names):
        for j,((req_key, title), ax) in enumerate(zip(panels, axes[i])):
            for sample in ["1000_10", "2500_10", "3500_10", "4000_10"]:
                v = np.asarray(efficiencies[window]["nobib"][sample]["velo_res"][req_key])
                if v.size==0:
                    continue
                weights = np.full_like(v, 100.0/v.size, dtype=float)

                if req_key=="ob_hits" and window=="loose":
                    print(sample, np.mean(v)*100)
                ax.hist(v, bins=bins, weights=weights, histtype="step", facecolor="none", linewidth=2.5, label=f"{sample_to_mass[sample]} TeV", linestyle="--")
                axes[0,0].text(
                    0.02, 0.98,
                    "Muon Collider",
                    ha="left", va="top",
                    transform=axes[0,0].transAxes,
                    fontsize=24,
                    fontweight="bold",
                    style="italic",
                )
                axes[0,0].text(
                    0.02, 0.90,
                    "Simulation No BIB, $\sqrt{s}$=10 TeV",
                    ha="left", va="top",
                    transform=axes[0,0].transAxes,
                    fontsize=18
                ) 
                axes[0,0].text(
                    0.02, 0.83,
                    "MuColl_v1",
                    ha="left", va="top",
                    transform=axes[0,0].transAxes,
                    fontsize=18
                )
                
                
                # if v_nobib.size==0:
                #     continue
                # weights_nobib = np.full_like(v_nobib, 100.0/v_nobib.size, dtype=float)
                # ax.hist(v_nobib, bins=bins, weights=weights_nobib, histtype="step", facecolor="none", linewidth=2.5, label=f"{sample_to_mass[sample]} TeV, no BIB", linestyle="--")

                # ax.set_title(title, fontsize=15)
                ax.grid(True, alpha=0.2)
                ax.grid(True, which="minor", alpha=0.1) 
                ax.minorticks_on()
                ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
                ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
                ax.set_xlim(-0.5, 0.5)

                if i == 2:
                    ax.set_xlabel(r"$(v^{\mathrm{reco}} - v^{\mathrm{truth}})/v^{\mathrm{truth}}$",fontsize=25)
                if j == 0:
                    ax.set_ylabel("Tracks per bin (%)", fontsize=25)
            axes[i, 0].text(-0.17, 0.5, f"{row_name} window",transform=axes[i, 0].transAxes,ha="right", va="center",rotation=90, fontsize=25, fontweight="bold", style="italic")
    
    for j, (_, title) in enumerate(panels):
        axes[0, j].set_title(title, fontsize=20, fontweight="bold", pad=10)

    handles, labels = [], []
    for ax_row in axes:
        for ax in ax_row:
            h, l = ax.get_legend_handles_labels()
            handles.extend(h)
            labels.extend(l)

    by_label = OrderedDict(zip(labels, handles))
    fig.legend(
        by_label.values(), by_label.keys(),
        ncol=4, fontsize=25,
        loc="upper center", bbox_to_anchor=(0.5, 1.02),
        handlelength=1.5, handletextpad=0.5, columnspacing=1.0
    )

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)   




    # mass_lab = f"{sample_to_mass[sample]:g} TeV"
    # window = "loose"

    # req_panels = [
    #     ("vb_hits", "≥3 VB"),
    #     ("ib_hits", "≥3 VB & ≥2 IB"),
    #     ("ob_hits", "≥3 VB, ≥2 IB, ≥2 OB"),
    # ]

    # for sample in ("1000_10", "2500_10", "3500_10", "4000_10"):
    #     for window in ("tight", "medium", "loose"):
    #         mass_lab = f"{sample_to_mass[sample]:g} TeV"    
    #         stau_bins = np.array(np.linspace(50,101, 20), dtype=float)
    #         bin_centers  = 0.5 * (stau_bins[1:] + stau_bins[:-1])

    #         fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), sharex=True, sharey=True, constrained_layout=True)

    #         for ax, (req_key, title) in zip(axes, req_panels):
    #             pT_res = np.asarray(efficiencies[window]["bib"][sample]["pT_res"][req_key], float)
    #             pcent  = np.asarray(efficiencies[window]["bib"][sample]["pcent_stau"][req_key], float) * 100.0  # %

    #             m = np.isfinite(pT_res) & np.isfinite(pcent)
    #             pT_res = pT_res[m]
    #             pcent  = pcent[m]
    #             if pT_res.size == 0:
    #                 ax.text(0.5, 0.5, "No tracks", ha="center", va="center", transform=ax.transAxes)
    #                 continue

    #             r_lo, r_hi = np.quantile(pT_res, [0.01, 0.99])
    #             pad = 0.05 * max(1e-6, r_hi - r_lo)
    #             residual_bins = np.linspace(r_lo - pad, r_hi + pad, 25)

    #             H, xedges, yedges = np.histogram2d(pcent, pT_res, bins=[stau_bins, residual_bins])

    #             X, Y = np.meshgrid(xedges, yedges)

    #             H_plot = H.T  
    #             H_plot_masked = np.ma.masked_where(H_plot <= 0, H_plot)

    #             pcm = ax.pcolormesh(X, Y, H_plot_masked, cmap="viridis", shading="auto")
    #             ax.set_title(title, fontsize=20)
    #             ax.set_xlabel("% of true stau hits in stau tracks", fontsize=20)
    #             ax.grid(True, alpha=0.2)

    #             ax.axhline(y=0, linestyle="--", color="gray", alpha=0.5)

    #         axes[0].set_ylabel(r"$(p_T^{\mathrm{reco}} - p_T^{\mathrm{truth}})/p_T^{\mathrm{truth}}$", fontsize=20)
    #         for ax in axes:
    #             ax.set_xlim(stau_bins[0], stau_bins[-1])

    #         cbar = fig.colorbar(pcm, ax=axes.ravel().tolist(), orientation="vertical", pad=0.02)
    #         cbar.set_label("Tracks per bin", fontsize=15)

    #         fig.suptitle(f"{mass_lab} staus — {window_to_name[window]} window", fontsize=25, y=1.1)
    #         #fig.tight_layout()
    #         pdf.savefig(fig, bbox_inches="tight")
    #         plt.close(fig)



print(f"Saved plots to {plot_path}")















 
# def get_residuals_by_sample(effs, kind, requirement_key, window, option, sample):
#     #  kind in {"pT_res", "velo_res"}
#     arr = effs[window][option][sample][kind][requirement_key]
#     if isinstance(arr, list):
#         return np.asarray(arr, dtype=float)
#     return np.array([], dtype=float)

# core_range = (-0.5,0.5)          

# def gaussian_core_fit(values, core=core_range, nbins=60):
#     v = np.asarray(values, float)
#     v = v[np.isfinite(v)]
#     n_total = v.size
#     if n_total == 0:
#         return (np.nan, np.nan, np.nan, 0, 0)
#     core = v[(v >= core_range[0]) & (v <= core_range[1])]
    
#     n_core = core.size
#     frac_core = n_core / n_total if n_total > 0 else np.nan
#     if n_core < 15:
#         # not enough stats to fit robustly
#         return (np.nan, np.nan, frac_core, n_total, n_core)
    
#     hist, edges = np.histogram(core, bins=nbins, range=core_range)
#     x = 0.5 * (edges[1:] + edges[:-1])
#     if np.all(hist == 0):
#         return (np.nan, np.nan, frac_core, n_total, n_core)

#     mu0 = float(np.mean(core))
#     sig0 = float(np.std(core, ddof=1)) if n_core >= 2 else 0.1
#     if not (np.isfinite(sig0) and sig0 > 1e-4):
#         sig0 = 0.1
#     A0 = float(hist.max())

#     def _gauss(xx, A, mu, sig):
#         return A * np.exp(-0.5 * ((xx - mu) / sig) ** 2)

#     try:
#         (A, mu, sig), _ = optimize.curve_fit(
#             _gauss, x, hist, p0=[A0, mu0, sig0],
#             bounds=([0.0, core_range[0]-1.0, 1e-4],
#                     [np.inf, core_range[1]+1.0, 1.0]),
#             maxfev=10000
#         )
#         mu_fit, sigma_fit = float(mu), float(sig)
#     except Exception:
#         mu_fit = mu0
#         sigma_fit = sig0

#     return (mu_fit, sigma_fit, frac_core, n_total, n_core)


# rows = []
# req_keys = ["vb_hits", "ib_hits", "ob_hits"]  

# for window in windows:                # loose / medium / tight
#     for option in bib_options:        # "bib" / "nobib"
#         for req_key in req_keys:      # requirement category
#             for sample in samples:    # each mass sample
#                 mass_TeV = sample_to_mass[sample]

#                 # pT residuals
#                 arr_pt = get_residuals_by_sample(
#                     efficiencies, "pT_res", req_key, window, option, sample
#                 )
#                 mean_pt, std_pt, frac_pt, n_pt, n_pt_core = gaussian_core_fit(arr_pt)

#                 rows.append({
#                     "kind": "pT",
#                     "window": window,
#                     "BIB": option,
#                     "requirement": req_key,
#                     #"sample": sample,
#                     "mass_TeV": mass_TeV,
#                     "mean_core": mean_pt,
#                     "std_core": std_pt,
#                     "frac_in_core": frac_pt,
#                     #"n_total": n_pt,
#                     #"n_core": n_pt_core,
#                 })

#                 # velocity residuals
#                 arr_v = get_residuals_by_sample(
#                     efficiencies, "velo_res", req_key, window, option, sample
#                 )
#                 mean_v, std_v, frac_v, n_v, n_v_core = gaussian_core_fit(arr_v)

#                 rows.append({
#                     "kind": "velocity",
#                     "window": window,
#                     "BIB": option,
#                     "requirement": req_key,
#                     #"sample": sample,
#                     "mass_TeV": mass_TeV,
#                     "mean_core": mean_v,
#                     "std_core": std_v,
#                     "frac_in_core": frac_v,
#                     #"n_total": n_v,
#                     #"n_core": n_v_core,
#                 })

# summary_csv = "/scratch/miralittmann/analysis/mira_analysis_code/stats_for_res.csv"
# df_summary = pd.DataFrame(rows)
# df_summary = df_summary.replace([np.nan, np.inf, -np.inf], "--")
# df_summary.sort_values(["kind","mass_TeV","window","BIB","requirement"], inplace=True)
# df_summary.to_csv(summary_csv, index=False)
# print(f"Wrote stats summary (mean/std/frac_core) to {summary_csv}")
