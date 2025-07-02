import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--mass_vs_goodtracks", action= "store_true")
parser.add_argument("--mass_vs_acceptance", action="store_true")
parser.add_argument("--mass_vs_trackeff", action="store_true")

args = parser.parse_args()
mass_vs_goodtracks_arg = args.mass_vs_goodtracks
mass_vs_acceptance_arg = args.mass_vs_acceptance
mass_vs_trackeff_arg = args.mass_vs_trackeff

save_plot_path = "/scratch/miralittmann/analysis/efficiency_plots/medium_window_all_masses.pdf"

sample_to_mass = {
    "1000_10": 1.0,
    "2500_10": 2.5,
    "4000_10": 4.0,
    "4500_10": 4.5
}
mass_list = [1.0, 2.5, 4.0, 4.5]
# windows = ["tight", "medium", "loose"]
windows = ["medium"]
samples = ["1000_10", "2500_10", "4000_10", "4500_10"]
fields = ["goodtracks_bib", "goodtracks_nobib", "acceptance", "trackeff_bib", "trackeff_nobib"]

track_eff_data = {
    window: {} for window in windows
}

med_1tev = [77.28, 97.95, 38.79, 81.51, 99.74]
med_25tev = [80.65, 100.00, 60.18, 79.56, 96.31]
med_4tev = [49.03, 99.53, 69.48, 47.68, 92.05]
med_45tev = [0.00, 100.00, 73.90, 0.00, 13.67]

med_all = [med_1tev, med_25tev, med_4tev, med_45tev]

track_eff_data["medium"] = {sample: dict(zip(fields, vals)) for sample, vals in zip(samples, med_all)}

marker_map = {"tight": "v", "medium":"o", "loose":"^"}
bibcolor = "orange"
nobibcolor = "blue"

def mass_vs_goodtracks(pdf):
    fig, ax = plt.subplots()
    for window in windows:
        goodtracks_bib = []
        goodtracks_nobib = []
        for sample in samples:
            goodtracks_bib.append(track_eff_data[window][sample]["goodtracks_bib"])
            goodtracks_nobib.append(track_eff_data[window][sample]["goodtracks_nobib"]) 
        ax.plot(mass_list, goodtracks_bib, color='orange', label=f"{window} window with BIB", marker=marker_map[window], linestyle='-', markersize=6)
        ax.plot(mass_list, goodtracks_nobib, color='blue', label=f'{window} window without BIB', marker=marker_map[window], linestyle='--', markersize=6)
    ax.set_xlabel("Mass [TeV]")
    ax.set_ylabel("Percentage of tracks with Reduced $\chi^2$ < 5")
    ax.set_title("'Good' tracks")
    ax.legend(ncols=2)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def mass_vs_acceptance(pdf):
    fig, ax = plt.subplots()
    for window in windows:
        acceptance = []
        for sample in samples: 
            acceptance.append(track_eff_data[window][sample]["acceptance"])
        ax.plot(mass_list, acceptance, label=f"{window} window", marker=marker_map[window])
    ax.legend()
    ax.set_xlabel("Mass [TeV]")
    ax.set_ylabel("Percentage of accepted staus")
    ax.set_title("Acceptance")
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

def mass_vs_trackeff(pdf):
    fig, ax = plt.subplots()
    for window in windows:
        trackeff_bib = []
        trackeff_nobib = []
        for sample in samples:
            trackeff_bib.append(track_eff_data[window][sample]["trackeff_bib"])
            trackeff_nobib.append(track_eff_data[window][sample]["trackeff_nobib"]) 
        ax.plot(mass_list, trackeff_bib, color='orange', label=f"{window} window with BIB", marker=marker_map[window], linestyle='-', markersize=6)
        ax.plot(mass_list, trackeff_nobib, color='blue', label=f'{window} window without BIB', marker=marker_map[window], linestyle='--', markersize=6)
    ax.set_xlabel("Mass [TeV]")
    ax.set_ylabel("Track reconstruction efficiency")
    ax.set_title("Tracking efficiency")
    ax.legend(ncols=2)    
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


with PdfPages(save_plot_path) as pdf: 

    if mass_vs_goodtracks_arg == True:
        mass_vs_goodtracks(pdf)

    if mass_vs_acceptance_arg == True:
        mass_vs_acceptance(pdf)

    if mass_vs_trackeff_arg == True:
        mass_vs_trackeff(pdf)

    


