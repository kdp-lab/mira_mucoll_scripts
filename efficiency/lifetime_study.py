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
n_files = 2500

CACHE = pathlib.Path("cache/lifetimes-trackstats.pkl") #-nozero

plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/lifetimes-decaypos-full.pdf"
track_stats_plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/lifetimes-trackstats.pdf" # -nozero

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
# lifetimes = [3, 10, 30]
lifetimes = [30]

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

track_reqs = ["vb", "ib", "ob"]

efficiencies = None
if (not rebuild) and os.path.exists(CACHE):
    with open(CACHE, "rb") as f:
        print("Loading in cached arrays...")
        efficiencies = pickle.load(f)

if efficiencies is None:

    efficiencies = {
        lifetime: {
            sample: {
                track_req: {
                    "pT": [],
                    "hits": [],
                    "velo": [],
                    "mass": [],
                    "leading_mass": [],
                    "subleading_mass": []
                } for track_req in track_reqs
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
                masses = []
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

                    nav = UTIL.LCRelationNavigator(track_relation_collection)
                    rel_nav = build_rel_nav(event)
                    masses_by_req = {req: [] for req in track_reqs}

                    for mcp in mcp_collection:
                        pdg = mcp.getPDG()
                        vx = mcp.getVertex()
                        ep = mcp.getEndpoint()

                        r_decay = np.sqrt(ep[0]**2 + ep[1]**2)
                        z_decay = ep[2]

                        acc = acceptance(mcp)
                        is_stau = (abs(pdg) in stau_ids)

                        has_track = 0

                        if acc and is_stau:
                            related_tracks = nav.getRelatedToObjects(mcp)

                            for track in related_tracks:
                                track_chi2 = track.getChi2() / track.getNdf()
                                if track_chi2 > chi2_cut:
                                    continue

                                reco_pT = 0.3 * Bfield / fabs(track.getOmega() * 1000.)

                                vb_hits = 0
                                ib_hits = 0
                                ob_hits = 0

                                track_hits = track.getTrackerHits()

                                track_times = []
                                track_pos = []
                                spatial_unc = []

                                for hit in track_hits:
                                    decoder.setValue(int(hit.getCellID0()))
                                    system = decoder["system"].value()

                                    if system in (1, 2):
                                        vb_hits += 0.5
                                        spatial_unc.append(0.005)
                                    elif system in (3, 4):
                                        ib_hits += 1
                                        spatial_unc.append(0.007)
                                    elif system in (5, 6):
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
                                
                                v_fit, v_err = reco_velo(linearfunc, track_times, track_pos, spatial_unc) 
                                total_hits = vb_hits + ib_hits + ob_hits

                                try:
                                    tanL = track.getTanLambda()
                                except Exception:
                                    tanL = np.nan

                                if np.isfinite(tanL):
                                    reco_pz = reco_pT * tanL
                                    reco_p  = math.sqrt(reco_pT**2 + reco_pz**2)
                                else:
                                    reco_p  = np.nan
                                beta = v_fit / speedoflight if np.isfinite(v_fit) else np.nan
                                
                                if np.isfinite(reco_p) and np.isfinite(beta) and (0 < beta < 1):
                                    m_reco = reco_p * math.sqrt(1.0/(beta*beta) - 1.0)
                                else:
                                    m_reco = np.nan

                                if vb_hits >= 3 and ib_hits >= 2 and ob_hits >= 2:
                                    save_info["ob"]["pT"].append(reco_pT)
                                    save_info["ob"]["hits"].append(total_hits)
                                    save_info["ob"]["velo"].append(v_fit)
                                    save_info["ob"]["mass"].append(m_reco) 
                                    if np.isfinite(m_reco):
                                        masses_by_req["ob"].append(float(m_reco))
                                if vb_hits >= 3 and ib_hits >= 2:
                                    save_info["ib"]["pT"].append(reco_pT)
                                    save_info["ib"]["hits"].append(total_hits)
                                    save_info["ib"]["velo"].append(v_fit)
                                    save_info["ib"]["mass"].append(m_reco) 
                                    if np.isfinite(m_reco):
                                        masses_by_req["ib"].append(float(m_reco))
                                if vb_hits >= 3:
                                    save_info["vb"]["pT"].append(reco_pT)
                                    save_info["vb"]["hits"].append(total_hits)
                                    save_info["vb"]["velo"].append(v_fit)
                                    save_info["vb"]["mass"].append(m_reco) 
                                    if np.isfinite(m_reco):
                                        masses_by_req["vb"].append(float(m_reco))

                    for req in track_reqs:
                        arr = masses_by_req[req]
                        if len(arr) == 0:
                            continue  
                        arr.sort(reverse=True)

                        leading = arr[0]
                        subleading = arr[1] if len(arr) >= 2 else float("nan")  

                        save_info[req]["leading_mass"].append(float(leading))
                        save_info[req]["subleading_mass"].append(float(subleading))

    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump(efficiencies, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Writing cache to {CACHE}")
    print("Saved cache successfully.")

print(efficiencies.keys())

# def get_eff_and_err(data, which="has_track"):
#     accepted = np.array(data["accepted"])
#     k = np.array(data[which])

#     mask = (accepted == 1)
#     n = mask.sum()
#     if n == 0:
#         return float('nan'), float('nan')

#     # Only consider accepted staus
#     p = k[mask].mean()
#     e = binom_se(p, n)
#     return 100 * p, 100 * e

# def rz_efficiency_map(data, r_bins, z_bins):
#     r = np.array(data["pos_r"])
#     z = np.array(data["pos_z"])
#     accepted = np.array(data["accepted"])
#     has_track = np.array(data["has_track"])

#     mask = (accepted == 1)

#     r_sel = r[mask]
#     z_sel = z[mask]
#     t_sel = has_track[mask]

#     N_tot, r_edges, z_edges = np.histogram2d(
#         r_sel, z_sel,
#         bins=[r_bins, z_bins]
#     )

#     N_rec, _, _ = np.histogram2d(
#         r_sel, z_sel,
#         bins=[r_bins, z_bins],
#         weights=t_sel.astype(float)
#     )

#     with np.errstate(divide="ignore", invalid="ignore"):
#         eff = np.where(N_tot > 0, N_rec / N_tot, np.nan)

#     return eff, r_edges, z_edges

# def collect_residuals(effs, type, requirement_key, window, option):
#     # requirement_key in {"vb_hits","ib_hits","ob_hits"}
#     vals = []
#     for sample in sample_to_mass.keys():
#         vals.extend(effs[window][option][sample][type][requirement_key])
#     return np.array(vals)



# panels = [
#     ("vb_hits", "≥3 VB hits"),
#     ("ib_hits", "≥3 VB & ≥2 IB hits"),
#     ("ob_hits", "≥3 VB, ≥2 IB, ≥2 OB hits")
# ]

# panels0 = [
#     ("ob_tracks", "≥3 VB, ≥2 IB, ≥2 OB hits")
# ]

# plt.style.use("seaborn-v0_8-colorblind")
# with PdfPages(plot_path) as pdf:
#     fig, ax = plt.subplots(figsize=(8,6))
#     sorted_samples = sorted(sample_to_mass.keys(), key=lambda s: sample_to_mass[s])
#     masses = [sample_to_mass[s] for s in sorted_samples]
#     colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
#     color_map = {3: colors[0], 10: colors[1], 30: colors[2]}
#     for lifetime in lifetimes: 
#         trackeff_bib, trackerr_bib = [], []
        
#         for sample in sorted_samples: 
#             d_b = efficiencies[lifetime][sample]
#             # p_b, e_b = get_eff_and_err(d_b, "track")
#             p_b, e_b = get_eff_and_err(d_b, "has_track")
#             trackeff_bib.append(p_b); trackerr_bib.append(e_b)

#         ax.errorbar(
#             masses, trackeff_bib, yerr=trackerr_bib,
#             linewidth=2, capsize=2, color=color_map[lifetime],
#             label=f"{lifetime}ns lifetime"
#         )

#     ax.text(
#         0.02, 0.20,
#         "Muon Collider",
#         ha="left", va="top",
#         transform=ax.transAxes,
#         fontsize=24,
#         fontweight="bold",
#         style="italic",
#     )
#     ax.text(
#         0.02, 0.13,
#         "Simulation 10% BIB, $\sqrt{s}$=10 TeV",
#         ha="left", va="top",
#         transform=ax.transAxes,
#         fontsize=18
#     ) 
#     ax.text(
#         0.02, 0.07,
#         "MuColl_v1",
#         ha="left", va="top",
#         transform=ax.transAxes,
#         fontsize=18
#     )        

#     ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
#     ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

#     ax.set_xlabel("Stau mass [TeV]", fontsize=20)
#     ax.set_ylabel("Track reconstruction efficiency (%)", fontsize=22)
#     ax.grid(True, alpha=0.2)
#     ax.set_ylim(0,100)

#     handles, labels = [], [] 
#     h,l = ax.get_legend_handles_labels()
#     handles.extend(h)
#     labels.extend(l)

#     by_label = OrderedDict(zip(labels, handles))
#     fig.legend(by_label.values(), by_label.keys(), fontsize=18, ncol=3, handlelength=1, handletextpad=0.3, columnspacing=0.5, loc="upper center",bbox_to_anchor=(0.5, 1.11))
    
#     fig.tight_layout()
#     pdf.savefig(fig, bbox_inches="tight")
#     plt.close(fig)


# print(f"Saved plots to {plot_path}")


# rz_plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/lifetimes-rz-eff.pdf"
# print(np.mean(efficiencies[30]["2500"]["pos_r"]))

# r_bins = np.linspace(0, 6000, 15)     
# z_bins = np.linspace(-1500, 1500, 15) 

# vb_r = [30, 51, 74, 102]
# vb_l = 130
# ve_z = [80, 120, 200, 280]
# ve_h = [0, 15, 35, 50]

# ib_r = [127, 340, 554]
# ib_l = [963.2, 963.2, 1384.6]
# ie_z = [524, 808, 1093, 1377, 1661, 1946, 2190]
# ie_h = np.arange(100,300, (300-100)/7) 

# ob_r = [819, 1153, 1486]
# ob_l = 2528.4
# oe_z = [1310, 1617, 1883, 2190]

# plt.style.use("seaborn-v0_8-colorblind")
# with PdfPages(rz_plot_path) as pdf:
#     lifetime_to_plot = 10

#     for sample in sorted(sample_to_mass.keys(), key=lambda s: sample_to_mass[s]):
#         mass = sample_to_mass[sample]
#         data = efficiencies[lifetime_to_plot][sample]

#         eff, r_edges, z_edges = rz_efficiency_map(data, r_bins, z_bins)

#         fig, ax = plt.subplots(figsize=(8, 6))

#         R, Z = np.meshgrid(r_edges, z_edges, indexing="ij")
#         im = ax.pcolormesh(
#             Z, R, eff,  
#             shading="auto"
#         )

#         cbar = fig.colorbar(im, ax=ax)
#         cbar.set_label("Track reconstruction efficiency", fontsize=16)

#         ax.set_xlabel("z decay position [mm]", fontsize=18)
#         ax.set_ylabel("r decay position [mm]", fontsize=18)
#         ax.set_title(f"m = {mass:.1f} TeV, τ = {lifetime_to_plot} ns",
#                      fontsize=18)

#         ax.tick_params(axis="both", which="major", labelsize=14, length=6, width=1.5)
#         ax.grid(True, alpha=0.2)

#         for layer in vb_r:
#             ax.plot([-vb_l/2, vb_l/2], [layer, layer], color="black")
#             ax.plot([-vb_l/2, vb_l/2], [layer+2, layer+2], color="black")
#         for i,layer in enumerate(ib_r):
#             ax.plot([-ib_l[i]/2, ib_l[i]/2], [layer, layer], color="black", linestyle=":")
#         for layer in ob_r:
#             ax.plot([-ob_l/2, ob_l/2], [layer, layer], color="black", linestyle="--")

#         fig.tight_layout()
#         pdf.savefig(fig, bbox_inches="tight")
#         plt.close(fig)

# print(f"Saved R–z efficiency plots to {rz_plot_path}")


# decay_hist_path = "/scratch/miralittmann/analysis/mira_analysis_code/lifetimes-decaydist.pdf"
# decay_hist_path = "/scratch/miralittmann/analysis/mira_analysis_code/lifetimes-decay_r.pdf"

# plt.style.use("seaborn-v0_8-colorblind")
# with PdfPages(decay_hist_path) as pdf:
#     sorted_samples = sorted(sample_to_mass.keys(), key=lambda s: sample_to_mass[s])

#     r_bins = np.linspace(0, 6000, 21)

#     for lifetime in lifetimes:
#         for sample in sorted_samples:
#             mass = sample_to_mass[sample]
#             data = efficiencies[lifetime][sample]

#             r = np.array(data["pos_r"])
#             accepted = np.array(data["accepted"])
#             mask = (accepted == 1)

#             r_sel = r[mask]

#             if r_sel.size == 0:
#                 continue

#             fig, ax = plt.subplots(figsize=(8, 6))

#             ax.hist(
#                 r_sel,
#                 bins=r_bins,
#                 histtype="step",
#                 linewidth=2
#             )

#             ax.set_xlabel("Decay radius r [mm]", fontsize=18)
#             ax.set_ylabel("Number of staus", fontsize=18)
#             ax.set_title(
#                 f"Radial decay distribution\nm = {mass:.1f} TeV, τ = {lifetime} ns",
#                 fontsize=18
#             )

#             ax.tick_params(axis="both", which="major", labelsize=14, length=6, width=1.5)
#             ax.grid(True, alpha=0.2)

#             ax.text(
#                 0.97, 0.95,
#                 f"N = {r_sel.size}",
#                 ha="right", va="top",
#                 transform=ax.transAxes,
#                 fontsize=14
#             )

#             fig.tight_layout()
#             pdf.savefig(fig, bbox_inches="tight")
#             plt.close(fig)

# print(f"Saved radial decay histograms to {decay_hist_path}")

labels = {"pT": r"$p_T$ [GeV]",
          "hits": "Hits on track",
          "velo": "Velocity [mm/ns]",
          "leading_mass": "Leading reconstructed mass [GeV]",
          "subleading_mass": "Subleading reconstructed mass [GeV]",}

def plot_feature_all_lifetimes(sample, feature, n_bins, x_lim=None):
    fig, ax = plt.subplots()

    for lifetime in lifetimes:
        feature_arr = np.asarray(efficiencies[lifetime][sample][feature])
        if x_lim:
            mask = (feature_arr >= x_lim[0]) & (feature_arr <= x_lim[1])
            feature_arr = feature_arr[mask]

        if feature_arr.size == 0:
            continue

        weights = np.full_like(feature_arr, 100.0/feature_arr.size, dtype=float)

        ax.hist(
            feature_arr,
            bins=n_bins,
            weights=weights,
            histtype="step",
            linewidth=2,
            label=f"{lifetime} ns",
        )

    ax.set_xlabel(labels[feature], fontsize=20)
    ax.set_ylabel("Normalized counts", fontsize=20)
    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    if x_lim:
        ax.set_xlim(x_lim[0], x_lim[1])

    ax.legend(fontsize=13)

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
        f"Staus, {sample_to_mass[sample]} TeV",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=15
    )
    ax.text(
        0.02, 0.83,
        "MuColl_v1, No BIB",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=15
    )

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def plot_mass_rank_all_lifetimes(sample, track_req, n_bins=60, x_lim=None):
    fig, ax = plt.subplots()

    for lifetime in lifetimes:
        lead_arr = np.asarray(efficiencies[lifetime][sample][track_req]["leading_mass"], dtype=float)
        lead_arr = lead_arr[np.isfinite(lead_arr)]
        sub_lead_arr = np.asarray(efficiencies[lifetime][sample][track_req]["subleading_mass"], dtype=float)
        sub_lead_arr = sub_lead_arr[np.isfinite(sub_lead_arr)]

        if x_lim is not None:
            sub_lead_arr = sub_lead_arr[(sub_lead_arr >= x_lim[0]) & (sub_lead_arr <= x_lim[1])]
            lead_arr = lead_arr[(lead_arr >= x_lim[0]) & (lead_arr <= x_lim[1])]

        if sub_lead_arr.size == 0:
            continue

        if lead_arr.size == 0:
            continue

        lead_weights = np.full_like(lead_arr, 100.0 / lead_arr.size, dtype=float)
        sub_lead_weights = np.full_like(sub_lead_arr, 100.0 / sub_lead_arr.size, dtype=float)

        ax.hist(
            lead_arr,
            bins=n_bins,
            weights=lead_weights,
            histtype="step",
            linewidth=2,
            label=f"Leading mass",
        )
        ax.hist(
            sub_lead_arr,
            bins=n_bins-40,
            weights=sub_lead_weights,
            histtype="step",
            linewidth=2,
            label=f"Sub-leading mass",
        )

    ax.set_xlabel("Mass [GeV]", fontsize=20)
    ax.set_ylabel("Normalized counts", fontsize=20)
    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    if x_lim is not None:
        ax.set_xlim(x_lim[0], x_lim[1])

    ax.legend(fontsize=13)

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
        f"Staus, {sample_to_mass[sample]} TeV, req={track_req}",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=15
    )
    ax.text(
        0.02, 0.83,
        "MuColl_v1, No BIB",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=15
    )

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


with PdfPages(track_stats_plot_path) as pdf:
    for sample in sample_to_mass.keys():
        # plot_feature_all_lifetimes(sample, "pT", 10, x_lim=(0,10000))
        # plot_feature_all_lifetimes(sample, "velo", 10)
        # plot_feature_all_lifetimes(sample, "hits", 10)
        for req in track_reqs: 
            plot_mass_rank_all_lifetimes(sample, req, n_bins=60, x_lim=(0, 6000))
            plot_mass_rank_all_lifetimes(sample, req, n_bins=60, x_lim=(0, 6000))
print(f"Saved track stats to {track_stats_plot_path}")