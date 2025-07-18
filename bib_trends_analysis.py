import pyLCIO
from pyLCIO import UTIL
import ROOT
from pathlib import Path
import numpy as np
from math import *
import json
import os
import pickle
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.pyplot as plt

# Currently doing this for just 4TeV, 10ns staus in the medium time window but with increasing amounts of BIB
# to see how efficiency scales

bib_options = np.arange(10, 55, 5)
sim_dir = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/sim/4000_10/"
reco_bigdir = "/scratch/miralittmann/analysis/mira_analysis_code/bib_trends/4000_10/medium/"
plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/bib_trends.pdf"

def get_chunk_id(fname: str) -> int:
    base = os.path.basename(fname)
    name = base.split(".", 1)[0]
    tail = name.split("_")[-1]
    digits = ''.join(ch for ch in tail if ch.isdigit())
    return int(digits)

all_eff_data = {
    option: {
        "accepted": [],
        "total_eff": [],
        "percent_good_tracks": []
    } for option in bib_options
}
for option in bib_options: 
    reco_dir = os.path.join(reco_bigdir, f"{option}p") 
    bad_chunks = set()
    total_files_processed = 0
    sim_all_data = []
    good_reco_tracks = 0
    bad_reco_tracks = 0

    print(f"\n=================={option} percent BIB:======================")

    for reco_file in os.listdir(reco_dir):
        reco_path = os.path.join(reco_dir, reco_file)
        if not os.path.isfile(reco_path):
            continue
        with open(reco_path) as file:
            reco_data = json.load(file)
            if get_chunk_id(reco_file) in reco_data.get("bad_files", []):
                bad_chunks.add(get_chunk_id(reco_file))
                continue
        track = reco_data["match_track_info"]["n_hits"] 
        chi_sq = np.asarray(reco_data["match_track_info"]["chi_sq"])
        ndf = np.asarray(reco_data["match_track_info"]["ndf"])
        red_chi_sq = chi_sq / ndf
        good_mask = red_chi_sq < 5
        bad_mask = red_chi_sq > 5
        bad_reco_tracks += bad_mask.sum()
        good_reco_tracks += good_mask.sum() 


    print(f"found {good_reco_tracks} tracks with reduced chi^2 < 5")
    print(f"and {bad_reco_tracks} with reduced chi^2 > 5 :(")
    percent_good_track = good_reco_tracks / (good_reco_tracks + bad_reco_tracks) * 100
    all_eff_data[option]["percent_good_tracks"].append(percent_good_track)
    total_tracks = good_reco_tracks + bad_reco_tracks
    print('total tracks', total_tracks)
 
    for sim_file in os.listdir(sim_dir):
        chunk_id = get_chunk_id(sim_file)
        if chunk_id in bad_chunks:
            continue
        with open(os.path.join(sim_dir, sim_file)) as file:
            sim_data = json.load(file)
            total_files_processed += 1
            truth_staus = sim_data["mcp_stau_info"]["id"]
            hit_info = sim_data["hit_info"]
            acc_staus = sim_data["n_accepted_staus"]
            chunk_data = {"truth_staus":truth_staus, "hit_info": hit_info, "acc_staus": acc_staus}
            sim_all_data.append(chunk_data)

    total_accepted_staus = 0
    for i in range(len(sim_all_data)):
        accepted_stau_per_event = sim_all_data[i]["acc_staus"]
        total_accepted_staus += accepted_stau_per_event
    all_eff_data[option]["accepted"].append(total_accepted_staus)
    print("all accepted", total_accepted_staus)
    tot_eff = good_reco_tracks / total_accepted_staus * 100
    all_eff_data[option]["total_eff"].append(tot_eff)

    print(f"total efficiency is {tot_eff}")



with PdfPages(plot_path) as pdf:
    efficiencies = []
    percents = []
    for percent in bib_options: 
        efficiencies.append(all_eff_data[percent]["total_eff"])
        percents.append(percent)

    fig, ax = plt.subplots(figsize=(6,5))
    ax.plot(percents, efficiencies, linestyle='--', marker='.', color="tab:blue")
    ax.set_xlabel("Percent of BIB included (%)")
    ax.set_ylabel("Track reconstruction efficiency (%)")
    ax.set_title("Tracking reconstruction efficiency with increasing BIB inclusion")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    good_tracks = []
    for percent in bib_options:
        good_tracks.append(all_eff_data[percent]["percent_good_tracks"])
    fig1, ax1 = plt.subplots(figsize=(6,5))
    ax1.plot(percents, good_tracks, linestyle = '--', marker=".", color="tab:blue")
    ax1.set_xlabel("Percent of BIB included (%)")
    ax1.set_ylabel("Percent of tracks with red chi^2 < 5")
    ax1.set_title("Percent of 'good' tracks versus BIB includsion")
    fig1.tight_layout()
    pdf.savefig(fig1)
    plt.close(fig1)

    fig2, ax2 = plt.subplots(figsize=(6,5))
    ax2.plot(percents, 1/(np.asarray(good_tracks) / np.asarray(efficiencies)), linestyle = '--', marker=".", color="tab:blue")
    ax2.set_xlabel("Percent of BIB included (%)")
    ax2.set_ylabel("Percent of tracks with red chi^2 < 5 / Reconstruction efficiency")
    ax2.set_title("Percent of 'good' tracks reconstructed")
    fig2.tight_layout()
    pdf.savefig(fig2)
    plt.close(fig2)



print(f"plot saved to {plot_path}")