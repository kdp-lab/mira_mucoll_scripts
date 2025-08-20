import os
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt
import json
import awkward as ak
from scipy import optimize
import uproot
from matplotlib.backends.backend_pdf import PdfPages

#all samples 4 TeV, 10 ns
bib_loose = uproot.open(("/scratch/miralittmann/analysis/mira_analysis_code/loose_occ_aug5.root"))
bib_medium = uproot.open("/scratch/miralittmann/analysis/mira_analysis_code/medium_occ_aug5.root")
bib_tight = uproot.open("/scratch/miralittmann/analysis/mira_analysis_code/tight_occ_aug5.root")

weighted_nhits_tight = bib_tight["h_ntimehits"]
weighted_nhits_medium = bib_medium["h_ntimehits"]
weighted_nhits_loose = bib_loose["h_ntimehits"]

#this is number of hits divided by layer area
weighted_data_tight = weighted_nhits_tight.to_numpy()
weighted_data_medium = weighted_nhits_medium.to_numpy()
weighted_data_loose = weighted_nhits_loose.to_numpy()

bin_edges = weighted_data_tight[1][0:52]
bin_contents_tight = weighted_data_tight[0]*10
bin_contents_medium = weighted_data_medium[0]*10
bin_contents_loose = weighted_data_loose[0]*10

# Plot the histogram
LABEL_FONTSIZE = 15      
TICK_FONTSIZE  = 13      
NOTE_FONTSIZE  = 10  


pdf_path = "/scratch/miralittmann/analysis/mira_analysis_code/bib_occupancy.pdf"
with PdfPages(pdf_path) as pdf:
    fig,ax = plt.subplots()
    ax.hist(bin_edges, bins=bin_edges, weights=bin_contents_tight, histtype="step", label="tight", color="maroon", linewidth=1.5, alpha=0.8)
    ax.hist(bin_edges, bins=bin_edges, weights=bin_contents_medium, histtype="step", label="medium", color="teal", linewidth=1.5, alpha=0.8)  
    ax.hist(bin_edges, bins=bin_edges, weights=bin_contents_loose, histtype="step", label="loose", color="palevioletred", linewidth=1.5,alpha=0.8)

    ax.axvline(8, color="black", ls=":")
    ax.text(4,30,'VXD barrel',rotation=90)

    ax.axvline(24, color="black", ls=":")
    ax.text(12,50,'VXD disks',rotation=0)

    ax.axvline(27, color="black", ls=":")
    ax.text(25,30,'IT barrel',rotation=90)

    ax.axvline(41, color="black", ls=":")
    ax.text(31,50,'IT disks',rotation=0)

    ax.axvline(44, color="black", ls=":")
    ax.text(42,30,'OT barrel',rotation=90)

    ax.text(45,50,'OT disks',rotation=0)


    #plt.title("weighted nhits")
    ax.tick_params(axis='both', labelsize=TICK_FONTSIZE) 
    ax.set_xlabel("Detector Layer", fontsize=LABEL_FONTSIZE, labelpad=4)
    ax.set_ylabel("Average number of hits / cm$^2$", fontsize=LABEL_FONTSIZE, labelpad=4)
    ax.legend(loc="upper right")
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()

print(f"saved plots to {pdf_path}")
