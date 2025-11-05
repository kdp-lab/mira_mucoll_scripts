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

import pyLCIO
from pyLCIO import UTIL, EVENT
import ROOT

sim_collections = ["VertexBarrelCollectionConed", "VertexEndcapCollectionConed", "InnerTrackerBarrelCollectionConed", "InnerTrackerEndcapCollectionConed", "OuterTrackerBarrelConed", "OuterTrackerEndcapCollectionConed"]
CACHE = pathlib.Path("cache/truth_dists_maia.pkl")
parser = argparse.ArgumentParser()
parser.add_argument("--rebuild", action="store_true")
args = parser.parse_args()
rebuild = args.rebuild 

dir = "/ospool/uc-shared/project/futurecolliders/miralittmann/MAIA/reco/nominal/nobib"
n_files = 2500
plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/truth_dists_maia.pdf"
stau_ids = {1000015, -1000015, 2000015, -2000015} 
sample_to_mass = {
    "1000_10": 1.0,
    # "2500_10": 2.5,
    # "3000_10": 3.0,
    # "3500_10": 3.5,
    # "4000_10": 4.0,
    # "4500_10": 4.5
}
speedoflight = 299792458/1000000  # mm/ns
truths = None
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

if (not rebuild) and os.path.exists(CACHE):
    with open(CACHE, "rb") as f:
        print("Loading in cached arrays...")
        truths = pickle.load(f)

if truths is None:
    truths = {sample: {
        "pT" : [],
        "eta": [],
        "betagamma": [],
        "beta": [],
        "velocity": []
    } for sample in sample_to_mass.keys()}

    for sample in sample_to_mass.keys():
        print(sample)
        save_info = truths[sample]
        for ifile in tqdm(range(n_files)):
            file_name = f"{sample}_reco{ifile}.slcio"
            full_path = os.path.join(dir, sample, file_name)
            if not os.path.exists(full_path) or os.path.getsize(full_path)==0:
                continue

            reader.open(full_path)
            for event in reader:
                mcp_collection = event.getCollection("MCParticle")
                if not mcp_collection:
                    continue

                for mcp in mcp_collection:
                    if mcp.getPDG() not in stau_ids:
                        continue

                    mom = mcp.getMomentum()
                    tlv = ROOT.TLorentzVector()
                    tlv.SetPxPyPzE(mom[0], mom[1], mom[2], mcp.getEnergy())

                    pT = tlv.Perp()
                    beta = tlv.Beta()
                    velocity = beta * speedoflight
                    betagamma = tlv.P() / tlv.M()
                    eta = tlv.Eta()

                    save_info["pT"].append(pT)
                    save_info["eta"].append(eta)
                    save_info["betagamma"].append(betagamma)
                    save_info["velocity"].append(velocity)
                    save_info["beta"].append(beta)

            reader.close()

    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump(truths, f, protocol=pickle.HIGHEST_PROTOCOL)
        print(f"Writing cache to {CACHE}")
    print("Saved cache successfully.")

with PdfPages(plot_path) as pdf:
    plt.style.use("seaborn-v0_8-colorblind")
    
    def common_edges(arrays, nbins=50):
        data = np.concatenate([np.asarray(a, float) for a in arrays if len(a)>0])
        lo = np.nanmin(data); hi = np.nanmax(data)
        if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
            lo, hi = 0.0, 1.0
        return np.linspace(lo, hi, nbins+1)

    PERCENT = True

        # -------------------- pT --------------------
    BW = 100.0  # GeV
    fig, ax = plt.subplots(figsize=(8,6))
    all_pt = np.concatenate([np.asarray(truths[s]["pT"], float)
                            for s in sample_to_mass if len(truths[s]["pT"]) > 0])

    if all_pt.size == 0:
        edges = np.linspace(0, BW, 2)  # fallback
    else:
        lo = np.nanmin(all_pt)
        hi = np.nanmax(all_pt)
        lo_aligned = BW * np.floor(lo / BW)
        hi_aligned = BW * np.ceil(hi / BW)
        edges = np.arange(lo_aligned, hi_aligned + BW, BW)

    for sample in sample_to_mass:
        data = np.asarray(truths[sample]["pT"], float)
        if data.size == 0:
            continue
        w = np.full(data.shape, 100.0 / data.size) 
        ax.hist(data, bins=edges, weights=w, histtype="step",
                linewidth=1.5, label=sample_to_mass[sample])

    ax.set_xlabel(r"$p_T$ [GeV]", fontsize=20)
    ax.set_ylabel("Staus per 100 GeV (%)", fontsize=20)  # fixed-width label
    leg = ax.legend(frameon=True, fontsize=17, title="Stau Mass [TeV]")
    leg.get_title().set_fontsize(17)
    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    plt.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # -------------------- eta --------------------
    field, xlabel = "eta", r"$\eta$"
    fig, ax = plt.subplots(figsize=(8,6))

    edges = common_edges([truths[s][field] for s in sample_to_mass], nbins=50) 

    for sample in sample_to_mass:
        data = np.asarray(truths[sample][field], float)
        if data.size == 0:
            continue
        w = np.full(data.shape, 100.0/data.size if PERCENT else 1.0/data.size)
        ax.hist(
            data, bins=edges, weights=w, histtype="step", linewidth=1.5,
            label=sample_to_mass[sample]
        )

    ax.set_xlabel(xlabel, fontsize=20)
    ax.set_ylabel("Staus per 0.125 (%)" if PERCENT else "Staus per bin (sum=1)", fontsize=20)
    leg = ax.legend(frameon=True, fontsize=17, title="Stau Mass [TeV]")
    leg.get_title().set_fontsize(17)
    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    plt.tight_layout(); pdf.savefig(fig); plt.close(fig)

    # for field, xlabel in [
    #     ("betagamma", r"Average $\beta\gamma$"),
    #     ("velocity",  r"Average velocity [mm/ns]"),
    #     ("beta", r"Average $\beta$")
    # ]:       
    #     fig, ax = plt.subplots(figsize=(8,6))
    #     masses = [1, 2.5, 3, 3.5, 4, 4.5]
    #     all_masses = []
    #     for sample in sample_to_mass.keys():
    #         data = np.asarray(truths[sample][field])
    #         avg_data = np.mean(data)
    #         all_masses.append(avg_data)

    #     ax.plot(masses,all_masses, linewidth=4, marker="o", markersize=10, markerfacecolor="none")
    #     ax.set_xlabel("Stau mass [TeV]", fontsize=20)
    #     ax.set_ylabel(xlabel, fontsize=20)

    #     ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    #     ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0) 
        
    #     plt.tight_layout()
    #     pdf.savefig(fig)
    #     plt.close(fig)

print(f"saved plots to {plot_path}") 