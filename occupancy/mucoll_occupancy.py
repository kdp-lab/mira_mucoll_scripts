import os
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import json
import awkward as ak
from scipy import optimize
import uproot
from matplotlib.backends.backend_pdf import PdfPages

cb = plt.colormaps['tab10'].colors
occ_conv = True
to_percent = True
poisson = False

VXD_PIXEL_CM2   = (25e-4) * (25e-4)   
IT_OT_PIXEL_CM2 = (50e-4) * (1e-4)

PIXEL_AREA_CM2 = np.r_[
    np.full(8,  VXD_PIXEL_CM2),    # 0..7   VXD barrel
    np.full(16, VXD_PIXEL_CM2),    # 8..23  VXD disks
    np.full(3,  IT_OT_PIXEL_CM2),  # 24..26 IT barrel
    np.full(14, IT_OT_PIXEL_CM2),  # 27..40 IT disks
    np.full(3,  IT_OT_PIXEL_CM2),  # 41..43 OT barrel
    np.full(8,  IT_OT_PIXEL_CM2),  # 44..51 OT disks
]
assert PIXEL_AREA_CM2.shape[0] == 52

def get_occupancy(root_file):
    file = uproot.open(root_file)
    weighted_nhits = file["h_ntimehits"]
    weighted_data = weighted_nhits.to_numpy()
    edges = weighted_data[1][0:52]
    bin_contents = weighted_data[0]*10

    if occ_conv:
        mu = bin_contents * PIXEL_AREA_CM2  # mean hits per channel per event
        if poisson:
            y = 1.0 - np.exp(-mu)       # exact P(≥1) for Poisson
        else:
            y = mu                      # sparse approximation
        if to_percent:
            y *= 100.0
        return edges, y
    else:
        return edges, bin_contents

bin_edges,loose_4tev = get_occupancy("/scratch/miralittmann/analysis/mira_analysis_code/loose4tev.root")
bin_edges,nominal = get_occupancy("/scratch/miralittmann/analysis/mira_analysis_code/tight_occ_aug5.root")
bin_edges,p26_p5_6 = get_occupancy("/scratch/miralittmann/analysis/mira_analysis_code/p26_p5_p6.root")
print(np.mean(p26_p5_6[41:51]/nominal[41:51])-1)
print(p26_p5_6[41], nominal[41])
print(np.mean(loose_4tev[41:51]/nominal[41:51])-1)

# Plot the histogram
LABEL_FONTSIZE = 20      
TICK_FONTSIZE  = 16      
NOTE_FONTSIZE  = 14  

pdf_path = "/scratch/miralittmann/analysis/mira_analysis_code/bib_occupancy_9-25.pdf"
with PdfPages(pdf_path) as pdf:
    plt.style.use("seaborn-v0_8-colorblind")
    fig,ax = plt.subplots(figsize=(8,6))
    ax.hist(bin_edges, bins=bin_edges, weights=nominal, histtype="step", label="Nominal", linewidth=2, alpha=0.8)
    ax.hist(bin_edges, bins=bin_edges, weights=p26_p5_6, histtype="step", label="Medium", linewidth=2, alpha=0.8)
    ax.hist(bin_edges, bins=bin_edges, weights=loose_4tev, histtype="step", label="Loose", linewidth=2, alpha=0.8)
    ax.axvline(8, color="black", ls=":")

    ax.text(4, 5e-3, "VXD Barrel", fontsize=NOTE_FONTSIZE, rotation=90,
        ha="center", va="center", weight="bold",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7))

    ax.axvline(24, color="black", ls=":")
    ax.text(16.5, 8e-3, "VXD Disks", fontsize=NOTE_FONTSIZE, rotation=0,
        ha="center", va="center", weight="bold",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7))

    ax.axvline(27, color="black", ls=":")
    ax.text(25.5, 5e-3, "IT Barrel", fontsize=NOTE_FONTSIZE, rotation=90,
        ha="center", va="center", weight="bold",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7))

    ax.axvline(41, color="black", ls=":")
    ax.text(34, 8e-3, "IT Disks", fontsize=NOTE_FONTSIZE, rotation=0,
        ha="center", va="center", weight="bold",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7))

    ax.axvline(44, color="black", ls=":")
    ax.text(42.5, 5e-3, "OT Barrel", fontsize=NOTE_FONTSIZE, rotation=90,
        ha="center", va="center", weight="bold",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7))

    ax.text(49, 8e-3, "OT Disks", fontsize=NOTE_FONTSIZE, rotation=0,
        ha="center", va="center", weight="bold",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.7))

    #plt.title("weighted nhits")
    ax.tick_params(axis="both", which="major", labelsize=TICK_FONTSIZE, width=1.2, length=6)
    ax.tick_params(axis="both", which="minor", width=1.0, length=4)
    ax.set_xlabel("Detector Layer", fontsize=LABEL_FONTSIZE, labelpad=4)
    if occ_conv == True:
        ax.set_ylabel("Occupancy per event (%)", fontsize=LABEL_FONTSIZE, labelpad=4)
    else:
        ax.set_ylabel("Average number of hits / cm$^2$", fontsize=LABEL_FONTSIZE, labelpad=4)
    ax.legend(loc="upper right", fontsize=14, frameon=True, framealpha=1, edgecolor="black") 
    ax.set_yscale("log")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

print(f"saved plots to {pdf_path}")
