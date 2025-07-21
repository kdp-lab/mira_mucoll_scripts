import os
import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import bisect
from pathlib import Path

save_plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/rzplane.pdf"
prefix = "4000_10_reco"

num_plots = 60
    
layer_map = [
    ("VB", 0, 30), ("VB", 1, 30), ("VB", 2, 51), ("VB", 3, 51), ("VB", 4, 74), ("VB", 5, 74), ("VB", 6, 102), ("VB", 7, 102),
    ("IB", 0, 127), ("IB", 1, 340), ("IB", 2, 554), 
    ("OB", 0, 819), ("OB", 1, 1153), ("OB", 2, 1486)
]

layer_names = [f"{det}{lay}" for det, lay, _ in layer_map]
layer_radii = np.array([r for _,_,r in layer_map])
n_layers_tot = len(layer_map)
timing_windows = {"tight": (0.15, 0.3, 0.3), "medium": (0.32, 1.25, 4.0), "loose": (0.32, 3.30, 10.0)}

with PdfPages(save_plot_path) as pdf:
    for window in ("tight", "medium", "loose"):
        reco_dir = Path(f"/scratch/miralittmann/analysis/mira_analysis_code/efficiency/nobib/{window}/4000_10/hit_level")
        reco_files = []

        for i in range(num_plots):
            file_path = reco_dir / f"{prefix}{i}.json"

            if not file_path.exists():
                continue
            reco_files.append(file_path)
        
        for reco_file in reco_files:
            with open(reco_file) as file:
                reco = json.load(file)

            all_coords = {"VB": {"coords": [], "match": [], "time": []},
                        "IB": {"coords": [], "match": [], "time": []},
                        "OB": {"coords": [], "match": [], "time": []}}
        
            r_all, z_all, match_all, time_all, in_time_all = [],[],[],[],[] 
            for det in ("VB", "IB", "OB"):
                if not reco["hits_from_mcp"][det]["x"]:
                    continue
    
                hitx = (reco["hits_from_mcp"][det]["x"])
                hity = (reco["hits_from_mcp"][det]["y"])
                hitz = (reco["hits_from_mcp"][det]["z"])
                time = (reco["hits_from_mcp"][det]["time"])
                hitr = np.hypot(hitx,hity)
        
                coords = np.column_stack((hitr,hitz))
                matched = (reco["hits_from_mcp"][det]["part_of_track"])

                in_time = []

                for i in range(len(time)): 
                    if time[i] < timing_windows[window][det]:
                        in_time.append(True)
                    else:
                        in_time.append(False)
        
                r_all.append(hitr)
                z_all.append(hitz)   
                match_all.append(matched)
                time_all.append(time)
                in_time_all.append(in_time)
            
            if len(r_all) != 0:
                r_all = np.concatenate(r_all)
                z_all = np.concatenate(z_all)
                match_all = np.concatenate(match_all).astype(bool)
                time_all = np.concatenate(time_all)
                in_time_all = np.concatenate(in_time_all)
            else:
                continue
                
            match_color = "tab:green"
            unmatch_color = "tab:brown"

            fig,ax = plt.subplots(figsize=(6,5))

            rad_label = {}
            for det,layer,radius in layer_map:
                if radius not in rad_label:
                    rad_label[radius]=f"{det} L{layer}"
            radii = sorted(rad_label)
        
            ax.scatter(r_all[~match_all], z_all[~match_all], label="unmatched to tracks", color=unmatch_color, s=40)
            ax.scatter(r_all[match_all], z_all[match_all], label="matched to tracks", color=match_color, s=40)
            
            z_min, z_max = ax.get_ylim()
            for r in radii:
                ax.axvline(r, linestyle="--", linewidth=0.8, color="grey", zorder=0)
                ax.text(r, z_max + 0.05, rad_label[r], ha="center", va="bottom", fontsize=5, rotation=45, color="dimgray", clip_on=False)

            show_mask = ~match_all # only plot times for hits that weren't matched to tracks. could also do only layer x: show_mask = np.isclose(r_all, <radius of detector>, atol=1)
            for x, y, t in zip(r_all[show_mask], z_all[show_mask], time_all[show_mask]):
                ax.annotate(f"{float(t):.2f}", xy=(x,y), xytext=(4,4), textcoords="offset points", fontsize=6, color="black")
        
    
            fig.subplots_adjust(top=0.93)
            fig.suptitle(Path(reco_file).name, fontsize=12, fontweight='bold')
            ax.set_xlabel("r [mm]")
            ax.set_ylabel("z [mm]") 
            ax.legend(loc="lower right")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

print(f"plots saved to {save_plot_path}")
