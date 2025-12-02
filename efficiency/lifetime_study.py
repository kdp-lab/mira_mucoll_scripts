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

loose_dir =  "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/loose4/seeding_10GeV/bib"
window_to_dir = {"loose": loose_dir}
n_files = 250

CACHE = pathlib.Path("cache/lifetimes-decaypos.pkl")

plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/lifetimes-decaypos.pdf"

sample_to_mass = {
    "1000": 1.0,
    "1500": 1.5,
    "2000": 2.0,
    "2500": 2.5,
    "3000": 3.0,
    "3500": 3.5,
    "4000": 4.0,
    "4500": 4.5
}
mass_list = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
lifetimes = [3, 10, 30]

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

    mcp_stau_momentum = mcp.getMomentum() 
    mcp_stau_tlv = ROOT.TLorentzVector()
    mcp_stau_tlv.SetPxPyPzE(mcp_stau_momentum[0], mcp_stau_momentum[1], mcp_stau_momentum[2], mcp.getEnergy())  

    if abs(mcp.getPDG()) not in stau_ids or travel_dist==0 or mcp_stau_vertex_r>553.0 or abs(mcp_stau_tlv.Eta()) > 0.8: # or mcp_stau_endpoint_r < 1486.0 
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
        lifetime: {
            sample: {
                "all_tracks": 0,
                "ob_tracks": 0,
                "accepted_staus": 0,

                "pT_res": {
                    "ob_hits": []
                },

                "velo_res": {
                    "ob_hits": [] 
                },
                # "accepted": [],
                # "track": [],
                # "endpoint": []
            } for sample in sample_to_mass.keys()
        } for lifetime in lifetimes
    } 


    reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

    for lifetime in lifetimes: 
        print(f"                      {lifetime}...")
        for sample in sample_to_mass.keys():
            print(f"                        {sample}...")
            save_info = efficiencies[lifetime][sample]

            for ifile in tqdm(range(n_files)):
                file_name = f"{sample}_{lifetime}/{sample}_{lifetime}_reco{ifile}.slcio"
                file_path = os.path.join(loose_dir,file_name) 
                if not os.path.exists(file_path) or os.path.getsize(file_path) == 0:
                    print(f"couldn't find {file_path}")
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
                        print("no hit collections")
                        continue
                    encoding = test_hit_coll.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                    decoder = UTIL.BitField64(encoding)

                    for mcp in mcp_collection:
                        accepted = acceptance(mcp)
                        if accepted == True:
                            save_info["accepted"].append(1)
                        else:
                            save_info["accepted"].append(0)
                        
                        mcp_stau_endpoint_r = np.sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2)
                        save_info["endpoint"].append(mcp_stau_endpoint_r)
                
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
                                    # save_info["all_tracks"] += 1
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

                                corrected_t = hit.getTime() - tof
                                corrected_corrected_t = 2*hit.getTime() - corrected_t

                                track_times.append(corrected_corrected_t)
                                track_pos.append(hit_pos)

                                # # percent stau
                                # for sim in rel_nav[system_to_relname[system]].getRelatedToObjects(hit):
                                #     mcp_true = sim.getMCParticle()
                                #     if mcp_true and abs(mcp_true.getPDG()) in stau_ids:
                                #         stau_hit_count += 1

                            # v_fit, v_err = reco_velo(linearfunc, track_times, track_pos, spatial_unc) 
                            
                            # pT_res = (pT - pT_truth) / pT_truth   
                            # velo_res = (v_fit - velo_truth) / velo_truth     
                            
                            if vb_hits >= 3 and ib_hits >= 2 and ob_hits >= 2:
                                save_info["track"].append(1)
                                # save_info["pT_res"]["ob_hits"].append(pT_res)
                                # save_info["velo_res"]["ob_hits"].append(velo_res)
                            else:
                                save_info["track"].append(0)

    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump(efficiencies, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Writing cache to {CACHE}")
    print("Saved cache successfully.")

print(efficiencies.keys())

def get_eff_and_err(data, which="all_tracks"):
    k = data[which]
    n = data["accepted_staus"]
    # n = data["accepted"]
    if n == 0:
        return float('nan'), float('nan')
    else:
        p = k/n
        e = binom_se(p, n)
    return 100*p, 100*e

def collect_residuals(effs, type, requirement_key, window, option):
    # requirement_key in {"vb_hits","ib_hits","ob_hits"}
    vals = []
    for sample in sample_to_mass.keys():
        vals.extend(effs[window][option][sample][type][requirement_key])
    return np.array(vals)



panels = [
    ("vb_hits", "≥3 VB hits"),
    ("ib_hits", "≥3 VB & ≥2 IB hits"),
    ("ob_hits", "≥3 VB, ≥2 IB, ≥2 OB hits")
]

panels0 = [
    ("ob_tracks", "≥3 VB, ≥2 IB, ≥2 OB hits")
]

plt.style.use("seaborn-v0_8-colorblind")
with PdfPages(plot_path) as pdf:
    fig, ax = plt.subplots(figsize=(8,6))
    sorted_samples = sorted(sample_to_mass.keys(), key=lambda s: sample_to_mass[s])
    masses = [sample_to_mass[s] for s in sorted_samples]
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color_map = {3: colors[0], 10: colors[1], 30: colors[2]}
    for lifetime in lifetimes: 
        trackeff_bib, trackerr_bib = [], []
        
        for sample in sorted_samples: 
            d_b = efficiencies[lifetime][sample]
            # p_b, e_b = get_eff_and_err(d_b, "track")
            p_b, e_b = get_eff_and_err(d_b, "ob_tracks")
            trackeff_bib.append(p_b); trackerr_bib.append(e_b)

        ax.errorbar(
            masses, trackeff_bib, yerr=trackerr_bib,
            linewidth=2, capsize=2, color=color_map[lifetime],
            label=f"{lifetime}ns lifetime"
        )

    ax.text(
        0.02, 0.20,
        "Muon Collider",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=24,
        fontweight="bold",
        style="italic",
    )
    ax.text(
        0.02, 0.13,
        "Simulation 10% BIB, $\sqrt{s}$=10 TeV",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=18
    ) 
    ax.text(
        0.02, 0.07,
        "MuColl_v1",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=18
    )        

    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

    ax.set_xlabel("Stau mass [TeV]", fontsize=20)
    ax.set_ylabel("Track reconstruction efficiency (%)", fontsize=22)
    ax.grid(True, alpha=0.2)
    ax.set_ylim(0,100)

    handles, labels = [], [] 
    h,l = ax.get_legend_handles_labels()
    handles.extend(h)
    labels.extend(l)

    by_label = OrderedDict(zip(labels, handles))
    fig.legend(by_label.values(), by_label.keys(), fontsize=18, ncol=3, handlelength=1, handletextpad=0.3, columnspacing=0.5, loc="upper center",bbox_to_anchor=(0.5, 1.11))
    
    fig.tight_layout()
    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)



    # bins = np.linspace(-1, 1, 21)
    # fig, axes = plt.subplots(3, 3, figsize=(18, 15), sharex=True, sharey=True)
    
    # for i, (window, row_name) in enumerate(windows_and_names):
    #     for j,((req_key, title), ax) in enumerate(zip(panels, axes[i])):
    #         for sample in ["1000_10", "2500_10", "3500_10", "4000_10"]:
    #             pT = np.asarray(efficiencies[window]["nobib"][sample]["pT_res"][req_key])
    #             if pT.size==0:
    #                 continue
    #             weights = np.full_like(pT, 100.0/pT.size, dtype=float)

    #             ax.hist(pT, bins=bins, weights=weights, histtype="step", facecolor="none", linewidth=2.5, label=f"{sample_to_mass[sample]} TeV", linestyle="--")
    #             axes[0,2].text(
    #                 0.98, 0.98,
    #                 "Muon Collider",
    #                 ha="right", va="top",
    #                 transform=axes[0,2].transAxes,
    #                 fontsize=24,
    #                 fontweight="bold",
    #                 style="italic",
    #             )
    #             axes[0,2].text(
    #                 0.98, 0.90,
    #                 r"Simulation No BIB, $\sqrt{s}=10\ \mathrm{TeV}$",
    #                 ha="right", va="top",
    #                 transform=axes[0,2].transAxes,
    #                 fontsize=18
    #             ) 
    #             axes[0,2].text(
    #                 0.98, 0.83,
    #                 "MuColl_v1",
    #                 ha="right", va="top",
    #                 transform=axes[0,2].transAxes,
    #                 fontsize=18
    #             )
    #             # ax.set_title(title, fontsize=15)
    #             ax.grid(True, alpha=0.2)
    #             ax.grid(True, which="minor", alpha=0.1) 
    #             ax.minorticks_on()
    #             ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    #             ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

    #             if i == 2:
    #                 ax.set_xlabel(r"$(p_T^{\mathrm{reco}} - p_T^{\mathrm{truth}})/p_T^{\mathrm{truth}}$",fontsize=25)
    #             if j == 0:
    #                 ax.set_ylabel("Tracks per bin (%)", fontsize=25)
    #         axes[i, 0].text(-0.2, 0.5, f"{row_name} window",transform=axes[i, 0].transAxes,ha="right", va="center",rotation=90, fontsize=25, fontweight="bold", style="italic")
    
    # for j, (_, title) in enumerate(panels):
    #     axes[0, j].set_title(title, fontsize=20, fontweight="bold", pad=10)

    # handles, labels = [], []
    # for ax_row in axes:
    #     for ax in ax_row:
    #         h, l = ax.get_legend_handles_labels()
    #         handles.extend(h)
    #         labels.extend(l)

    # by_label = OrderedDict(zip(labels, handles))
    # fig.legend(
    #     by_label.values(), by_label.keys(),
    #     ncol=4, fontsize=25,
    #     loc="upper center", bbox_to_anchor=(0.5, 1.02),
    #     handlelength=1.5, handletextpad=0.5, columnspacing=1.0
    # )

    # fig.tight_layout(rect=[0, 0, 1, 0.96])
    # pdf.savefig(fig, bbox_inches="tight")
    # plt.close(fig)   


    # bins = np.linspace(-1, 1, 41)
    # fig, axes = plt.subplots(3, 3, figsize=(18, 15), sharex=True, sharey=True)
    
    # for i, (window, row_name) in enumerate(windows_and_names):
    #     for j,((req_key, title), ax) in enumerate(zip(panels, axes[i])):
    #         for sample in ["1000_10", "2500_10", "3500_10", "4000_10"]:
    #             v = np.asarray(efficiencies[window]["nobib"][sample]["velo_res"][req_key])
    #             if v.size==0:
    #                 continue
    #             weights = np.full_like(v, 100.0/v.size, dtype=float)
    #             ax.hist(v, bins=bins, weights=weights, histtype="step", facecolor="none", linewidth=2.5, label=f"{sample_to_mass[sample]} TeV", linestyle="--")
    #             axes[0,0].text(
    #                 0.02, 0.98,
    #                 "Muon Collider",
    #                 ha="left", va="top",
    #                 transform=axes[0,0].transAxes,
    #                 fontsize=24,
    #                 fontweight="bold",
    #                 style="italic",
    #             )
    #             axes[0,0].text(
    #                 0.02, 0.90,
    #                 "Simulation No BIB, $\sqrt{s}$=10 TeV",
    #                 ha="left", va="top",
    #                 transform=axes[0,0].transAxes,
    #                 fontsize=18
    #             ) 
    #             axes[0,0].text(
    #                 0.02, 0.83,
    #                 "MuColl_v1",
    #                 ha="left", va="top",
    #                 transform=axes[0,0].transAxes,
    #                 fontsize=18
    #             )
                

    #             ax.grid(True, alpha=0.2)
    #             ax.grid(True, which="minor", alpha=0.1) 
    #             ax.minorticks_on()
    #             ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    #             ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    #             ax.set_xlim(-0.5, 0.5)

    #             if i == 2:
    #                 ax.set_xlabel(r"$(v^{\mathrm{reco}} - v^{\mathrm{truth}})/v^{\mathrm{truth}}$",fontsize=25)
    #             if j == 0:
    #                 ax.set_ylabel("Tracks per bin (%)", fontsize=25)
    #         axes[i, 0].text(-0.17, 0.5, f"{row_name} window",transform=axes[i, 0].transAxes,ha="right", va="center",rotation=90, fontsize=25, fontweight="bold", style="italic")
    
    # for j, (_, title) in enumerate(panels):
    #     axes[0, j].set_title(title, fontsize=20, fontweight="bold", pad=10)

    # handles, labels = [], []
    # for ax_row in axes:
    #     for ax in ax_row:
    #         h, l = ax.get_legend_handles_labels()
    #         handles.extend(h)
    #         labels.extend(l)

    # by_label = OrderedDict(zip(labels, handles))
    # fig.legend(
    #     by_label.values(), by_label.keys(),
    #     ncol=4, fontsize=25,
    #     loc="upper center", bbox_to_anchor=(0.5, 1.02),
    #     handlelength=1.5, handletextpad=0.5, columnspacing=1.0
    # )

    # fig.tight_layout(rect=[0, 0, 1, 0.96])
    # pdf.savefig(fig, bbox_inches="tight")
    # plt.close(fig)   



print(f"Saved plots to {plot_path}")
