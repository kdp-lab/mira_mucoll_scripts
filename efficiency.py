import json
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import argparse

sim_path = "/scratch/miralittmann/analysis/efficiency_v1/sim/1000_10/"
reco_path = "/scratch/miralittmann/analysis/efficiency_v1/nobib/medium/1000_10/"

parser = argparse.ArgumentParser()
parser.add_argument("--plotting", action="store_true")
args = parser.parse_args()
plotting = args.plotting

def get_chunk_id(fname: str) -> int:
    base = os.path.basename(fname)
    name = base.split(".", 1)[0]
    tail = name.split("_")[-1]
    digits = ''.join(ch for ch in tail if ch.isdigit())
    return int(digits)

# Get bad chunks
bad_chunks = set()
for reco_file in os.listdir(reco_path):
    with open(os.path.join(reco_path, reco_file)) as file:
        chunk_reco_data = json.load(file)
        if get_chunk_id(reco_file) in chunk_reco_data.get("bad_files", []):
            bad_chunks.add(get_chunk_id(reco_file))

print(f"Bad chunks to exclude: {len(bad_chunks)}")

stau_ids = {1000015, -1000015, 2000015, -2000015}

events_data = []
total_files_processed = 0

for sim_file in os.listdir(sim_path):
    chunk_id = get_chunk_id(sim_file)
    if chunk_id in bad_chunks:
        continue
        
    with open(os.path.join(sim_path, sim_file)) as file:
        chunk_sim_data = json.load(file)
        total_files_processed += 1 
        event_data = {
            'file': sim_file,
            'truth_staus': chunk_sim_data["mcp_stau_info"]["id"],
            'hit_info': chunk_sim_data["hit_info"],
            'accepted_staus': chunk_sim_data["n_accepted_staus"]
        }
        events_data.append(event_data)

total_accepted_staus = 0
for i in range(len(events_data)):
    accepted_stau_per_event = events_data[i]["accepted_staus"]
    total_accepted_staus += accepted_stau_per_event
print(total_accepted_staus)


print(f"Files processed: {total_files_processed}")
print(f"Total events: {len(events_data)}")

total_truth_staus = sum(len(event['truth_staus']) for event in events_data)
print(f"Total truth staus: {total_truth_staus}")


print(f"\n=== (SIM) ===")
print(f"Total truth staus: {total_truth_staus}")
print(f"Total accepted staus: {total_accepted_staus}")
if total_truth_staus > 0:
    acceptance_rate = total_accepted_staus / total_truth_staus * 100
    print(f"Acceptance rate: {acceptance_rate:.2f}%")


good_reco_tracks = 0
total_reco_tracks = 0

for reco_file in os.listdir(reco_path):
    with open(os.path.join(reco_path, reco_file)) as file: 
        reco_data = json.load(file)
        if get_chunk_id(reco_file) in bad_chunks: 
            continue 
        chi_sq = reco_data["match_track_info"]["chi_sq"]
        ndf = reco_data["match_track_info"]["ndf"]
        
        for i in range(len(chi_sq)):
            red_chi_sq = chi_sq[i] / ndf[i]
            if red_chi_sq < 5:
                good_reco_tracks +=1

        total_reco_tracks += len(reco_data.get("match_stau_info", {}).get("id", []))

rate_good_tracks = (good_reco_tracks / total_reco_tracks) * 100
if total_accepted_staus > 0:
    efficiency =  good_reco_tracks / total_accepted_staus * 100
    print(f"\n=== (RESULTS) ===")
    print(f"Efficiency: {efficiency:.2f}%")
    print(f"Good tracks: {good_reco_tracks}")
    print(f"Percent of good tracks: {rate_good_tracks:.2f}%" )
else:
    print("Cannot calculate efficiency: no accepted staus")

########################################## PLOTTING ######################################
if plotting == True:
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
