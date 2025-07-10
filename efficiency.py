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
vars = ["theta", "phi", "pt"]

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

    eff_by_variable = {
        window: {
            sample: {
                option: {
                    "theta": {
                        "eff": [],
                        "err": [],
                        "centers": []},
                    "phi": {
                        "eff": [],
                        "err": [],
                        "centers": []},
                    "pt": {
                        "eff": [],
                        "err": [],
                        "centers": []},
                } for option in bib_options
            } for sample in samples
        } for window in windows}

    if CACHE.exists() and not redo:
        print(f"Loading previous info from {CACHE}, not redoing full analysis")
        with CACHE.open("rb") as f:
            track_eff_data, all_data, eff_by_variable = pickle.load(f)
        return track_eff_data, all_data, eff_by_variable            

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
                            'accepted_staus': chunk_sim_data["n_accepted_staus"],
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

                theta_bins = np.linspace(0 - 0.1, 2.6+0.1, 21)
                theta_num = np.zeros(len(theta_bins)-1, int)
                theta_den = np.zeros(len(theta_bins)-1, int)
 
                phi_bins = np.linspace(-np.pi - 0.1, np.pi+0.1, 21)
                phi_num = np.zeros(len(phi_bins)-1, int)
                phi_den = np.zeros(len(phi_bins)-1, int)
 
                pt_bins = np.linspace(-5, 5, 21)
                pt_num = np.zeros(len(pt_bins)-1, int)
                pt_den = np.zeros(len(pt_bins)-1, int)
 
                for reco_file in os.listdir(reco_path):
                    with open(os.path.join(reco_path, reco_file)) as file: 
                        reco_data = json.load(file)
                        if get_chunk_id(reco_file) in bad_chunks: 
                            continue 

                    chi_sq = np.asarray(reco_data["match_track_info"]["chi_sq"])
                    ndf = np.asarray(reco_data["match_track_info"]["ndf"])
                    theta_truth = np.asarray(reco_data["match_track_info"]["theta"])
                    pt_truth = np.asarray(reco_data["match_track_info"]["pt"])/1000
                    phi_truth = np.asarray(reco_data['match_track_info']["phi"])

                    theta_den += np.histogram(theta_truth, bins=theta_bins)[0]
                    phi_den += np.histogram(phi_truth, bins=phi_bins)[0]
                    pt_den += np.histogram(pt_truth, bins=pt_bins)[0]
                
                    red_chi_sq = chi_sq / ndf
                    good_mask = red_chi_sq < 5
                    good_reco_tracks += good_mask.sum()

                    theta_reco = theta_truth[good_mask]
                    pt_reco = pt_truth[good_mask]
                    phi_reco = phi_truth[good_mask]

                    theta_num += np.histogram(theta_reco, bins=theta_bins)[0]
                    phi_num += np.histogram(phi_reco, bins=phi_bins)[0]
                    pt_num += np.histogram(pt_reco, bins=pt_bins)[0]

                theta_eff = np.divide(theta_num, theta_den, out=np.zeros_like(theta_num, float), where=theta_den>0)
                phi_eff = np.divide(phi_num, phi_den, out=np.zeros_like(phi_num, float), where=phi_den>0)                     
                pt_eff = np.divide(pt_num, pt_den, out=np.zeros_like(pt_num, float), where=pt_den>0)            

                theta_err = np.zeros_like(theta_eff)
                theta_m = theta_den > 0
                theta_err[theta_m] = np.sqrt(theta_eff[theta_m] * (1-theta_eff[theta_m])/theta_den[theta_m])
                theta_centers = 0.5*(theta_bins[1:] + theta_bins[:-1])
                eff_by_variable[window][sample][option]["theta"]["eff"].append(theta_eff)                
                eff_by_variable[window][sample][option]["theta"]["err"].append(theta_err)
                eff_by_variable[window][sample][option]["theta"]["centers"].append(theta_centers)

                phi_err = np.zeros_like(phi_eff)
                phi_m = phi_den > 0
                phi_err[phi_m] = np.sqrt(phi_eff[phi_m] * (1-phi_eff[phi_m])/phi_den[phi_m])
                phi_centers = 0.5*(phi_bins[1:] + phi_bins[:-1])
                eff_by_variable[window][sample][option]["phi"]["eff"].append(phi_eff)                
                eff_by_variable[window][sample][option]["phi"]["err"].append(phi_err)
                eff_by_variable[window][sample][option]["phi"]["centers"].append(phi_centers)

                pt_err = np.zeros_like(pt_eff)
                pt_m = pt_den > 0
                pt_err[pt_m] = np.sqrt(pt_eff[pt_m] * (1-pt_eff[pt_m])/pt_den[pt_m])
                pt_centers = 0.5*(pt_bins[1:] + pt_bins[:-1])
                eff_by_variable[window][sample][option]["pt"]["eff"].append(pt_eff)                
                eff_by_variable[window][sample][option]["pt"]["err"].append(pt_err)
                eff_by_variable[window][sample][option]["pt"]["centers"].append(pt_centers)


                total_efficiency =  good_reco_tracks / total_accepted_staus * 100
                print(f"\n ========== reco: track efficiency ==========")
                print(f"Efficiency: {total_efficiency:.2f}%") 
                                
                track_eff_data[window][sample]["acceptance"].append(acceptance_rate)
                if option == "bib/": 
                    track_eff_data[window][sample]["trackeff_bib"].append(total_efficiency)
                else:
                    track_eff_data[window][sample]["trackeff_nobib"].append(total_efficiency)
    
    print(f"Writing cache to {CACHE}")
    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump((track_eff_data, all_data, eff_by_variable), f, protocol=pickle.HIGHEST_PROTOCOL)
    return track_eff_data, all_data, eff_by_variable

track_eff_data, all_data, eff_by_variable = build_analysis(redo=rebuild)

print(eff_by_variable["medium"]["4000_10"]["nobib/"]["pt"])

########################################## PLOTTING ######################################
if plotting == True:
    print("Now making plots...")
    marker_map = {"tight": "v", "medium":"o", "loose":"^"}
    bibcolor = "orange"
    nobibcolor = "blue"
    window_colors = {
        'tight': 0.3,
        'medium': 0.6,
        'loose': 0.9
    } 

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
            ax.plot(mass_list, trackeff_bib, color='orange', label=f"{window}, BIB", linestyle='-', markersize=6, alpha = window_colors[window])
            ax.plot(mass_list, trackeff_nobib, color='blue', label=f'{window}, no BIB', linestyle=':', markersize=6, alpha = window_colors[window])
        ax.set_xlabel("Mass [TeV]")
        ax.set_ylabel("Track reconstruction efficiency")
        ax.set_title("Tracking efficiency")
        ax.legend(ncols=2)    
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    def eff_vs_theta(pdf):
        n_cols = len(samples)
        n_rows = len(windows)

        fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(3.5*n_cols, 2.8*n_rows)) 
        colors = {"bib/":"orange", "nobib/":"blue"}
        for r, window in enumerate(windows):
            for c, sample in enumerate(samples):
                ax = axes[r,c] if n_rows>1 else axes[c]
                for option in bib_options:
                    arr = eff_by_variable[window][sample][option]["theta"]
                    if not arr["centers"]:
                        continue
                    centers = arr["centers"][0]
                    eff = arr["eff"][0]
                    err = arr['err'][0]

                    ax.errorbar(centers, eff, yerr=err, label="BIB" if option=="bib/" else "no BIB", color=colors[option])

                    if r==0:
                        ax.set_title(f"{sample}")
                    if c==0:
                        ax.set_ylabel(f"{window}")
                    ax.grid(ls=":", alpha=0.4)

        handles, labels = axes[0,0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right")
        fig.suptitle("tracking efficiency vs truth theta")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


    def eff_vs_pt(pdf):
        n_cols = len(samples)
        n_rows = len(windows)

        fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(3.5*n_cols, 2.8*n_rows)) 
        colors = {"bib/":"orange", "nobib/":"blue"}
        for r, window in enumerate(windows):
            for c, sample in enumerate(samples):
                ax = axes[r,c] if n_rows>1 else axes[c]
                for option in bib_options:
                    arr = eff_by_variable[window][sample][option]["pt"]
                    if not arr["centers"]:
                        continue
                    centers = arr["centers"][0]
                    eff = arr["eff"][0]
                    err = arr['err'][0]

                    ax.errorbar(centers, eff, yerr=err, label="BIB" if option=="bib/" else "no BIB", color=colors[option])

                    if r==0:
                        ax.set_title(f"{sample}")
                    if c==0:
                        ax.set_ylabel(f"{window}")
                    ax.grid(ls=":", alpha=0.4)

        handles, labels = axes[0,0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right")
        fig.suptitle("tracking efficiency vs truth pt")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


    def eff_vs_phi(pdf):
        n_cols = len(samples)
        n_rows = len(windows)

        fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(3.5*n_cols, 2.8*n_rows)) 
        colors = {"bib/":"orange", "nobib/":"blue"}
        for r, window in enumerate(windows):
            for c, sample in enumerate(samples):
                ax = axes[r,c] if n_rows>1 else axes[c]
                for option in bib_options:
                    arr = eff_by_variable[window][sample][option]["phi"]
                    if not arr["centers"]:
                        continue
                    centers = arr["centers"][0]
                    eff = arr["eff"][0]
                    err = arr['err'][0]

                    ax.errorbar(centers, eff, yerr=err, label="BIB" if option=="bib/" else "no BIB", color=colors[option])

                    if r==0:
                        ax.set_title(f"{sample}")
                    if c==0:
                        ax.set_ylabel(f"{window}")
                    ax.grid(ls=":", alpha=0.4)

        handles, labels = axes[0,0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right")
        fig.suptitle("tracking efficiency vs phi")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
       

    with PdfPages(save_plot_path) as pdf:
        mass_vs_acceptance(pdf)
        mass_vs_trackeff(pdf)
        eff_vs_theta(pdf)
        eff_vs_pt(pdf)
        eff_vs_phi(pdf)
    print(f"Saved plot(s) to {save_plot_path}")
