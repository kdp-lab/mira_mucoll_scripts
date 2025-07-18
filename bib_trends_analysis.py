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
import matplotlib as mpl

# Currently doing this for just 4TeV, 10ns staus in the medium time window but with increasing amounts of BIB
# to see how efficiency scales


# bib_options = np.arange(10, 55, 5)
bib_options = [50]
sim_dir = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/sim/4000_10/"
reco_bigdir = "/scratch/miralittmann/analysis/mira_analysis_code/bib_trends/4000_10/medium/"
plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/bib_trends.pdf"
print(len(os.listdir(sim_dir)))

def get_chunk_id(fname: str) -> int:
    base = os.path.basename(fname)
    name = base.split(".", 1)[0]
    tail = name.split("_")[-1]
    digits = ''.join(ch for ch in tail if ch.isdigit())
    return int(digits)

all_eff_data = {
    option: {
        "accepted": [],
        "eff_precut": [],
        "total_eff": [],
        "percent_good_tracks": [],
        "loss_chi2": [],
        "loss_reco": [],
        "total_matched": [],
        "eta_eff": [],
        "eta_err": [],
        "eta_centers": []
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

    eta_bins = np.linspace(-1.01, 1.01, 21)
    eta_den = np.zeros(len(eta_bins)-1, int)
    eta_num = np.zeros(len(eta_bins)-1, int) 
    accepted_staus = 0
 
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
        accepted_staus += reco_data["accepted_staus"] 

        theta_truth = np.asarray(reco_data["match_track_info"]["theta"])
        eta_truth = np.asarray(-np.log(np.tan(0.5*theta_truth)))
        eta_den += np.histogram(eta_truth, bins=eta_bins)[0]
        eta_reco = eta_truth[good_mask]
        eta_num += np.histogram(eta_reco, bins=eta_bins)[0]
    
    eta_eff = np.divide(eta_num, eta_den, out=np.zeros_like(eta_num, float), where=eta_den>0)
    eta_err = np.zeros_like(eta_eff)
    eta_m = eta_den > 0
    eta_err[eta_m] = np.sqrt(eta_eff[eta_m]* (1-eta_eff[eta_m])/eta_den[eta_m])
    eta_centers = 0.5*(eta_bins[1:] + eta_bins[:-1])

    all_eff_data[option]["eta_err"]=eta_err*100
    all_eff_data[option]["eta_eff"]=eta_eff*100
    all_eff_data[option]["eta_centers"]=eta_centers
    

    print(f"found {good_reco_tracks} tracks with reduced chi^2 < 5")
    print(f"and {bad_reco_tracks} with reduced chi^2 > 5 :(")
    percent_good_track = good_reco_tracks / (good_reco_tracks + bad_reco_tracks) * 100
    all_eff_data[option]["percent_good_tracks"]=(percent_good_track)
    total_tracks = good_reco_tracks + bad_reco_tracks
    print('total tracks', total_tracks)
    print('total accepted', accepted_staus)
     
    total_accepted_staus = accepted_staus
    all_eff_data[option]["accepted"] = total_accepted_staus
    tot_eff = good_reco_tracks / total_accepted_staus * 100
    all_eff_data[option]["total_eff"]=(tot_eff)
    all_eff_data[option]["loss_reco"]=((total_accepted_staus - (good_reco_tracks + bad_reco_tracks))/total_accepted_staus)*100
    all_eff_data[option]["eff_precut"]=((good_reco_tracks + bad_reco_tracks)/total_accepted_staus * 100)
    all_eff_data[option]["loss_chi2"]=(bad_reco_tracks / total_accepted_staus * 100)

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

    accepted_total = [all_eff_data[option]["accepted"] for option in bib_options]
    eff_precut = [all_eff_data[option]["eff_precut"] for option in bib_options]
    eff_total = [all_eff_data[option]["total_eff"] for option in bib_options]

    fig2, ax2 = plt.subplots(figsize=(6,5))
    ax2.plot(bib_options, eff_precut, label="before $\chi^2$ cut")
    ax2.plot(bib_options, eff_total, label="after $\chi^2$ cut")
    ax2.set_xlabel("Percent of BIB included (%)")
    ax2.set_ylabel("Reco efficiency (%)")
    ax2.set_title("Reco tracking efficiency before/after $\chi^2$ cut")
    ax2.legend()
    fig2.tight_layout()
    pdf.savefig(fig2)
    plt.close(fig2)

    loss_chi2_frac = [all_eff_data[p]["loss_chi2"] for p in bib_options]
    loss_reco_frac = [all_eff_data[p]["loss_reco"] for p in bib_options]
    loss_share = [lc/(lc+lr)*100 for lc, lr in zip(loss_chi2_frac, loss_reco_frac)]
        
    fig3, ax3 = plt.subplots(figsize=(6,5))
    ax3.plot(bib_options, loss_share, marker='o')
    ax3.set_xlabel("Percent of BIB (%)")
    ax3.set_ylabel("$\chi^2$-cut share of total tracking loss")
    ax3.set_title("How much of total lost reconstruction efficiency comes from $\chi^2 > 5$?")
    pdf.savefig(fig3)
    plt.close(fig3)

    fig4, ax4 = plt.subplots(figsize=(6,5))
    norm = mpl.colors.Normalize(vmin=bib_options.min(), vmax=bib_options.max())
    cmap = plt.cm.viridis 
    for option in bib_options:
        print(all_eff_data[option]["eta_err"])
        color=cmap(norm(option))
        alpha = 0.9-0.6*norm(option)
        ax4.errorbar(all_eff_data[option]["eta_centers"], all_eff_data[option]["eta_eff"], yerr=all_eff_data[option]["eta_err"], label=option, color=color, alpha=alpha)
    ax4.legend()
    ax4.set_xlabel("Truth $\eta$")
    ax4.set_ylabel("Reconstruction efficiency (%)")
    ax4.set_title("$\eta$-dependent tracking efficiency")
    pdf.savefig(fig4)
    plt.close(fig4)
    
print(f"plot saved to {plot_path}")