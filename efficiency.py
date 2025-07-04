import json
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import pickle
import pathlib
from tqdm import tqdm

sim_dir = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/sim/"
reco_dir = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/"
save_plot_path = "/scratch/miralittmann/analysis/efficiency_plots/all_windows_all_masses.pdf"
CACHE = pathlib.Path("cache/efficiency.pkl")

sample_to_mass = {
    "1000_10": 1.0,
    "2500_10": 2.5,
    "4000_10": 4.0,
    "4500_10": 4.5
}

mass_list = [1.0, 2.5, 4.0, 4.5]
windows = ["tight", "medium", "loose"]
samples = ["1000_10", "2500_10", "4000_10", "4500_10"]
bib_options = ["bib/", "nobib/"]
fields = ["acceptance", "trackeff_bib", "trackeff_nobib"]

parser = argparse.ArgumentParser()
parser.add_argument("--plotting", action="store_true")
parser.add_argument("--rebuild", action="store_true")
args = parser.parse_args()
plotting = args.plotting
rebuild = args.rebuild

def get_chunk_id(fname: str) -> int:
    base = os.path.basename(fname)
    name = base.split(".", 1)[0]
    tail = name.split("_")[-1]
    digits = ''.join(ch for ch in tail if ch.isdigit())
    return int(digits)

stau_ids = {1000015, -1000015, 2000015, -2000015}


def build_analysis(redo=False):
    total_files_processed = 0
    all_data = {
        window: {
            option: {
                sample: [] for sample in samples
            } for option in bib_options
        } for window in windows
    }

    track_eff_data = {
        window: {
            sample: {
                field: []
                for field in fields
            } 
            for sample in samples
        }
        for window in windows
    }

    if CACHE.exists() and not redo:
        print(f"Loading previous info from {CACHE}, not redoing full analysis")
        with CACHE.open("rb") as f:
            track_eff_data, all_data = pickle.load(f)
        return track_eff_data, all_data            

    print("Rebuilding analysis arrays...")
    for window in tqdm(windows):
        for option in bib_options: 
            for sample in samples:
                sim_path = os.path.join(sim_dir, sample)
                reco_path = os.path.join(reco_dir, option, window, sample)
                
                events_data = []

                ##################################### setup + acceptance (from sim info) #############################
                for sim_file in os.listdir(sim_path):

                # Get bad chunks
                    bad_chunks = set()
                    for reco_file in os.listdir(reco_path):
                        with open(os.path.join(reco_path, reco_file)) as file:
                            chunk_reco_data = json.load(file)
                            if get_chunk_id(reco_file) in chunk_reco_data.get("bad_files", []):
                                bad_chunks.add(get_chunk_id(reco_file))
                    
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

                all_data[window][option][sample].append(events_data)

                total_accepted_staus = 0
                for i in range(len(events_data)):
                    accepted_stau_per_event = events_data[i]["accepted_staus"]
                    total_accepted_staus += accepted_stau_per_event
                 
                print(f"\n ====================================================================")
                print(f"{window} window, {sample}, {option}")  
                print(f"Total events in files: {len(events_data)}")

                total_truth_staus = sum(len(event['truth_staus']) for event in events_data)
                print(f"\n ========= sim ============")
                print(f"Total truth staus: {total_truth_staus}")
                print(f"Total accepted staus: {total_accepted_staus}")
                if total_truth_staus > 0:
                    acceptance_rate = total_accepted_staus / total_truth_staus * 100
                    print(f"Acceptance rate: {acceptance_rate:.2f}%")

                ################################### tracking efficiency ######################################### 
                good_reco_tracks = 0 

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
     
                efficiency =  good_reco_tracks / total_accepted_staus * 100
                print(f"\n ========== reco: track efficiency ==========")
                print(f"Efficiency: {efficiency:.2f}%") 
                                
                track_eff_data[window][sample]["acceptance"].append(acceptance_rate)
                if option == "bib/": 
                    track_eff_data[window][sample]["trackeff_bib"].append(efficiency)
                else:
                    track_eff_data[window][sample]["trackeff_nobib"].append(efficiency)


                ################################## hit-based efficiency ########################################
                
                # getting travel distance about the stau from sim (truth level), then calculating how many hits that should give it
                # for sim_file in os.listdir(sim_path):
                #     if get_chunk_id(sim_file) in bad_chunks:
                #         continue
                #     with open(os.path.join(sim_path, sim_file)) as file:
                #         all_sim_data = json.load(file)
                #         # for ???????
                #         distance = all_sim_data["mcp_stau_info"]["travel_distance"]
                        

                        

    
    print(f"Writing cache to {CACHE}")
    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump((track_eff_data, all_data), f, protocol=pickle.HIGHEST_PROTOCOL)
    return track_eff_data, all_data

track_eff_data, all_data = build_analysis(redo=rebuild)

########################################## PLOTTING ######################################
if plotting == True:
    print("Now making plots...")
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
        acceptance = []
        for sample in samples: 
            acceptance.append(track_eff_data["medium"][sample]["acceptance"])
        ax.plot(mass_list, acceptance)
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

    def mass_vs_eta(pdf):
        fig, ax = plt.subplots()


    with PdfPages(save_plot_path) as pdf:
        mass_vs_goodtracks(pdf)
        mass_vs_acceptance(pdf)
        mass_vs_trackeff(pdf)
    print(f"Saved plot(s) to {save_plot_path}")
