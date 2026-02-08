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

hit_collections = ["VXDBarrelHits", "VXDEndcapHits", "ITBarrelHits", "ITEndcapHits", "OTBarrelHits", "OTEndcapHits"]
sim_collections = ["VertexBarrelCollection", "VertexEndcapCollection", "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection", "OuterTrackerBarrel", "OuterTrackerEndcapCollectionConed"]
rel_collections = ["VXDBarrelHitsRelations", "VXDEndcapHitsRelations", "ITBarrelHitsRelations", "ITEndcapHitsRelations", "OTBarrelHitsRelations", "OTEndcapHitsRelations"]

loose_dir =  "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/loose4/seeding_10GeV/nobib"
window_to_dir = {"loose": loose_dir}
n_files = 2500

CACHE = pathlib.Path("cache/sig_by_event.pkl") #-nozero

plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/lifetimes-decaypos-full.pdf"
track_stats_plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/sig_by_event.pdf" # -nozero

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
                    "leading_mass": [],
                    "subleading_mass": [],
                    "leading_pT": [],
                    "subleading_pT": [],
                    "leading_beta": [],
                    "subleading_beta": [],
                    "leading_hits": [],
                    "subleading_hits": [],

                    "n_acc_staus": [],          
                    "has1_acc_stau": [],        
                    "has2_acc_stau": [],        

                    "n_pass_pt": [],
                    "n_pass_mass": [],
                    "n_pass_beta": [],
                    "n_pass_all": [],

                    "event_has1_pt": [], "event_has2_pt": [],
                    "event_has1_mass": [], "event_has2_mass": [],
                    "event_has1_beta": [], "event_has2_beta": [],
                    "event_has1_all": [], "event_has2_all": [],
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
                    tracks_by_req = {req: [] for req in track_reqs}
                    n_acc_staus = 0

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
                            n_acc_staus += 1
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

                                track_info = {
                                    "pT": float(reco_pT) if np.isfinite(reco_pT) else np.nan,
                                    "beta": float(beta) if np.isfinite(beta) else np.nan,
                                    "mass": float(m_reco) if np.isfinite(m_reco) else np.nan,
                                }


                                if vb_hits >= 3 and ib_hits >= 2 and ob_hits >= 2:
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
                        passing_all = [t for t in trks if pass_all(t)]
                        passing_all.sort(key=lambda t: t["mass"], reverse=True)

                        save_info[req]["n_pass_pt"].append(n_pt)
                        save_info[req]["n_pass_mass"].append(n_mass)
                        save_info[req]["n_pass_beta"].append(n_beta)
                        save_info[req]["n_pass_all"].append(n_all)

                        save_info[req]["event_has1_pt"].append(1 if n_pt   >= 1 else 0)
                        save_info[req]["event_has2_pt"].append(1 if n_pt   >= 2 else 0)
                        save_info[req]["event_has1_mass"].append(1 if n_mass >= 1 else 0)
                        save_info[req]["event_has2_mass"].append(1 if n_mass >= 2 else 0)
                        save_info[req]["event_has1_beta"].append(1 if n_beta >= 1 else 0)
                        save_info[req]["event_has2_beta"].append(1 if n_beta >= 2 else 0)
                        save_info[req]["event_has1_all"].append(1 if n_all  >= 1 else 0)
                        save_info[req]["event_has2_all"].append(1 if n_all  >= 2 else 0)

                        save_info[req]["n_acc_staus"].append(n_acc_staus)
                        save_info[req]["has1_acc_stau"].append(1 if n_acc_staus >= 1 else 0)
                        save_info[req]["has2_acc_stau"].append(1 if n_acc_staus >= 2 else 0)


                        if len(passing_all) >= 1:
                            lead = passing_all[0]
                            save_info[req]["leading_mass"].append(lead["mass"])
                            save_info[req]["leading_pT"].append(lead["pT"])
                            save_info[req]["leading_beta"].append(lead["beta"])
                        else:
                            save_info[req]["leading_mass"].append(float("nan"))
                            save_info[req]["leading_pT"].append(float("nan"))
                            save_info[req]["leading_beta"].append(float("nan"))

                        if len(passing_all) >= 2:
                            sub = passing_all[1]
                            save_info[req]["subleading_mass"].append(sub["mass"])
                            save_info[req]["subleading_pT"].append(sub["pT"])
                            save_info[req]["subleading_beta"].append(sub["beta"])
                        else:
                            save_info[req]["subleading_mass"].append(float("nan"))
                            save_info[req]["subleading_pT"].append(float("nan"))
                            save_info[req]["subleading_beta"].append(float("nan"))

    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump(efficiencies, f, protocol=pickle.HIGHEST_PROTOCOL)
    print(f"Writing cache to {CACHE}")
    print("Saved cache successfully.")

print(efficiencies.keys())

labels = {"pT": r"$p_T$ [GeV]",
          "hits": "Hits on track",
          "velo": "Velocity [mm/ns]",
          "leading_mass": "Leading reconstructed mass [GeV]",
          "subleading_mass": "Subleading reconstructed mass [GeV]",}


def print_cutflow_summary(efficiencies, lifetimes, sample_to_mass, track_reqs=("ob",), use_two_stau_den_for_ge2=True):
    """
    Print event-level cutflow counts/efficiencies for each sample (and lifetime),
    using denominators based on accepted staus per event.

    Assumes for each [lifetime][sample][req] we have arrays:
      - "n_acc_staus" (0/1/2 per event)  OR ("has1_acc_stau","has2_acc_stau")
      - "event_has1_pt",   "event_has2_pt"
      - "event_has1_mass", "event_has2_mass"
      - "event_has1_beta", "event_has2_beta"
      - "event_has1_all",  "event_has2_all"

    If use_two_stau_den_for_ge2=True:
      - >=2-track efficiencies are computed with denominator = events with 2 accepted staus
      - >=1-track efficiencies are computed with denominator = events with >=1 accepted stau
    Otherwise:
      - both >=1 and >=2 efficiencies use denominator = events with >=1 accepted stau
    """
    def _safe_div(num, den):
        return float(num) / float(den) if den > 0 else float("nan")

    def _fmt(num, den):
        eff = 100.0 * _safe_div(num, den)
        if den <= 0:
            return f"{num}/{den} (nan%)"
        return f"{num}/{den} ({eff:6.2f}%)"

    # stable ordering by mass
    samples_sorted = sorted(sample_to_mass.keys(), key=lambda s: sample_to_mass[s])

    for lifetime in lifetimes:
        print("\n" + "=" * 120)
        print(f"Cutflow summary | lifetime = {lifetime} ns")
        print(f"Track-level cuts used in event_has*: pT>{PT_MIN}, mass>{MASS_MIN}, beta<{BETA_MAX}")
        print("=" * 120)

        for req in track_reqs:
            print("\n" + "-" * 120)
            print(f"Track requirement group: {req}")
            print("-" * 120)
            header = (
                "mass[TeV]  sample  N_events  N(acc>=1) N(acc=2) | "
                "PT:  >=1pass /den  >=2pass /den | "
                "MASS:>=1pass /den  >=2pass /den | "
                "BETA:>=1pass /den  >=2pass /den | "
                "ALL: >=1pass /den  >=2pass /den"
            )
            print(header)
            print("-" * 120)

            for sample in samples_sorted:
                d = efficiencies[lifetime][sample][req]

                if "n_acc_staus" in d:
                    n_acc = np.asarray(d["n_acc_staus"], dtype=int)
                    acc1 = (n_acc >= 1)
                    acc2 = (n_acc >= 2)
                else:
                    acc1 = np.asarray(d.get("has1_acc_stau", []), dtype=int) == 1
                    acc2 = np.asarray(d.get("has2_acc_stau", []), dtype=int) == 1

                N_evt = len(acc1)
                N_acc1 = int(acc1.sum())
                N_acc2 = int(acc2.sum())

                den_ge1 = N_acc1
                den_ge2 = N_acc2 if use_two_stau_den_for_ge2 else N_acc1

                def _num(key, mask, default_len):
                    arr = d.get(key, None)
                    if arr is None:
                        return 0
                    a = np.asarray(arr, dtype=int)
                    if len(a) != default_len:
                        raise ValueError(f"Length mismatch for {sample}/{req}/{key}: len={len(a)} expected {default_len}")
                    return int(a[mask].sum())

                n_pt_1   = _num("event_has1_pt",   acc1, N_evt)
                n_mass_1 = _num("event_has1_mass", acc1, N_evt)
                n_beta_1 = _num("event_has1_beta", acc1, N_evt)
                n_all_1  = _num("event_has1_all",  acc1, N_evt)

                mask2 = acc2 if use_two_stau_den_for_ge2 else acc1
                n_pt_2   = _num("event_has2_pt",   mask2, N_evt)
                n_mass_2 = _num("event_has2_mass", mask2, N_evt)
                n_beta_2 = _num("event_has2_beta", mask2, N_evt)
                n_all_2  = _num("event_has2_all",  mask2, N_evt)

                mtev = sample_to_mass[sample]

                line = (
                    f"{mtev:7.1f} {sample:>5s} {N_evt:4d} {N_acc1:4d} {N_acc2:4d} | "
                    f"pT {_fmt(n_pt_1, den_ge1):>16s} {_fmt(n_pt_2, den_ge2):>16s} | "
                    f"m {_fmt(n_mass_1, den_ge1):>16s} {_fmt(n_mass_2, den_ge2):>16s} | "
                    f"β {_fmt(n_beta_1, den_ge1):>16s} {_fmt(n_beta_2, den_ge2):>16s} | "
                    f"ALL {_fmt(n_all_1, den_ge1):>16s} {_fmt(n_all_2, den_ge2):>16s}"
                )
                print(line)

        print("\nNote:")
        print(" - '>=1pass' efficiencies use denominator = events with >=1 accepted stau.")
        if use_two_stau_den_for_ge2:
            print(" - '>=2pass' efficiencies use denominator = events with 2 accepted staus (recommended).")
        else:
            print(" - '>=2pass' efficiencies use denominator = events with >=1 accepted stau (less strict).")
        print("=" * 120 + "\n")


print_cutflow_summary(efficiencies, lifetimes, sample_to_mass, track_reqs=("ob",))

def _event_norm_hist(ax, arr, N_events, bins, label):
    """Histogram normalized to number of events (not number of entries)."""
    x = np.asarray(arr, dtype=float)
    x = x[np.isfinite(x)]
    if x.size == 0 or N_events <= 0:
        return 0.0
    w = np.full_like(x, 1.0 / N_events, dtype=float)   # sums to fraction of events
    ax.hist(x, bins=bins, weights=w, histtype="step", linewidth=2, label=label)
    return x.size / N_events

def plot_leading_subleading_eventnorm(efficiencies, lifetimes, sample_to_mass, track_reqs,
                                      out_pdf,
                                      bins_cfg=None,
                                      xlims_cfg=None):
    if bins_cfg is None:
        bins_cfg = {
            "pT":   np.linspace(0, 7000, 51),
            "mass": np.linspace(0, 6000, 61),
            "beta": np.linspace(0.6, 1.05, 53),
        }
    if xlims_cfg is None:
        xlims_cfg = {
            "pT": (0, 7000),
            "mass": (0, 6000),
            "beta": (0.6, 1.05),
        }

    features = [
        ("pT",   "leading_pT",   "subleading_pT",   r"$p_T$ [GeV]"),
        ("mass", "leading_mass", "subleading_mass", r"reco mass [GeV]"),
        ("beta", "leading_beta", "subleading_beta", r"$\beta$"),
    ]

    samples_sorted = sorted(sample_to_mass.keys(), key=lambda s: sample_to_mass[s])

    with PdfPages(out_pdf) as pdf:
        for lifetime in lifetimes:
            for sample in samples_sorted:
                mtev = sample_to_mass[sample]
                for req in track_reqs:
                    d = efficiencies[lifetime][sample][req]

                    N_events = len(d.get("leading_mass", []))
                    if N_events == 0:
                        continue

                    for key, lead_key, sub_key, xlabel in features:
                        fig, ax = plt.subplots(figsize=(8, 6))

                        frac_lead = _event_norm_hist(ax, d.get(lead_key, []), N_events, bins_cfg[key],
                                                     label="Leading")
                        frac_sub  = _event_norm_hist(ax, d.get(sub_key, []),  N_events, bins_cfg[key],
                                                     label="Subleading")

                        ax.set_xlabel(xlabel, fontsize=16)
                        ax.set_ylabel("Fraction of events per bin", fontsize=16)
                        ax.legend(frameon=True, loc="upper right", fontsize=14)

                        if key in xlims_cfg:
                            ax.set_xlim(*xlims_cfg[key])

                        ax.grid(True, alpha=0.2)

                        ax.text(0.02, 0.98, "Muon Collider", ha="left", va="top",
                                transform=ax.transAxes, fontsize=16, fontweight="bold", style="italic")
                        ax.text(0.02, 0.94, f"stau sample: {mtev:.1f} TeV, τ={lifetime} ns, req={req}",
                                ha="left", va="top", transform=ax.transAxes, fontsize=12)
                        ax.text(0.02, 0.90, f"N_events={N_events}",
                                ha="left", va="top", transform=ax.transAxes, fontsize=12)

                        ax.text(0.02, 0.02,
                                f"Frac(events w/ leading) ~ {frac_lead:.3f}\n"
                                f"Frac(events w/ subleading) ~ {frac_sub:.3f}",
                                ha="left", va="bottom", transform=ax.transAxes, fontsize=11)
                        
                        ax.tick_params(axis="both", which="major", labelsize=16)
                        ax.tick_params(axis="both", which="minor", labelsize=14)

                        fig.tight_layout()
                        pdf.savefig(fig)
                        plt.close(fig)

    print(f"Saved event-normalized leading/subleading plots to {out_pdf}")

plot_leading_subleading_eventnorm(
    efficiencies=efficiencies,
    lifetimes=lifetimes,
    sample_to_mass=sample_to_mass,
    track_reqs=["ob"],                
    out_pdf="/scratch/miralittmann/analysis/mira_analysis_code/stau_lead_sublead.pdf",
)

