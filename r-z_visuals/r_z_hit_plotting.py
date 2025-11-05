import os
import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import bisect
from pathlib import Path
import matplotlib.transforms as mtrans

save_plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/4dtrackingproof.pdf"
base_dir = Path("/scratch/miralittmann/analysis/mira_analysis_code/pt_test/pt10/4000_10/")
prefix = "4000_10_reco"

num_plots = 60 
    
layer_map = [
    ("VB", 0, 30), ("VB", 1, 30), ("VB", 2, 51), ("VB", 3, 51), ("VB", 4, 74), ("VB", 5, 74), ("VB", 6, 102), ("VB", 7, 102),
    ("IB", 0, 127), ("IB", 1, 340), ("IB", 2, 554), 
    ("OB", 0, 819), ("OB", 1, 1153), ("OB", 2, 1486)
]

stau_ids = [1000015, 2000015]

layer_names = [f"{det}{lay}" for det, lay, _ in layer_map]
layer_radii = np.array([r for _,_,r in layer_map])
n_layers_tot = len(layer_map)

rad_label = {}
for det,layer,radius in layer_map:
    if radius not in rad_label:
        rad_label[radius]=f"{det} L{layer}"
match_color = "tab:green"

non_stau_match_color = "tab:pink"

def plot_file(json_path, ax):
    with open(json_path) as file:
        reco = json.load(file)

    r_all, z_all, stau_all, time_all = [], [], [], []
    for det in ("VB", "IB", "OB"):
        track_hit_info = reco["match_track_info"]["hits"]
        hitx = np.concatenate(track_hit_info[det]["x"])
        hity = np.concatenate(track_hit_info[det]["y"])
        hitz = np.concatenate(track_hit_info[det]["z"])
        time = np.concatenate(track_hit_info[det]["time"])
        hitr = np.hypot(hitx, hity)
        pdg = np.concatenate(track_hit_info[det]["pdg_id"])
        is_stau = np.isin(np.abs(pdg), stau_ids)
        print(len(hitx), len(hity))

        r_all.append(hitr)
        z_all.append(hitz)

        time_all.append(time)
        stau_all.append(is_stau)

    if not r_all:
        ax.set_axis_off()
        return
    
    r_all = np.concatenate(r_all)
    z_all = np.concatenate(z_all)    
    time_all = np.concatenate(time_all)
    stau_all = np.concatenate(stau_all)

    ax.scatter(r_all, z_all, label="matched to tracks", color=match_color, s=40)
    ax.scatter(r_all[~stau_all], z_all[~stau_all], label="non-stau matched", color=non_stau_match_color, s=40)
 
    for r in layer_radii:
        ax.axvline(r, ls='--', lw=0.8, color="grey", zorder=0)
        trans = mtrans.blended_transform_factory(ax.transData, ax.transAxes)
        ax.text(r, 1.02, rad_label[r], transform=trans, rotation=45, color="dimgray", fontsize=5)
    
    for x, y, t in zip(r_all, z_all, time_all):
        ax.annotate(f"{t:.2f}", xy=(x,y), xytext=(3,3), textcoords="offset points", fontsize=6)

    ax.set_xlabel("r [mm]")
    ax.set_ylabel("z [mm]")

with PdfPages(save_plot_path) as pdf:
    for proc_id in range(num_plots):
        paths = base_dir / f"{prefix}{proc_id}.json"            
        fig, axes = plt.subplots()
        plot_file(paths, axes)
        
        fig.suptitle(f"Process ID {proc_id}", fontsize=14, fontweight="bold") 
        pdf.savefig(fig)
        plt.close(fig)

print(f"saved plots to {save_plot_path}")