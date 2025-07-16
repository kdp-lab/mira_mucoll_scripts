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
                "acceptance": [],
                "trackeff_bib": {
                    "eff": [],
                    "err": []
                },
                "trackeff_nobib" : {
                    "eff": [],
                    "err": []
                }
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
    
    hit_eff_data = {
        window: {
            sample: {"eff": [], "err":[]} for sample in samples
        } for window in windows
    }

    if CACHE.exists() and not redo:
        print(f"Loading previous info from {CACHE}, not redoing full analysis")
        with CACHE.open("rb") as f:
            track_eff_data, all_data, eff_by_variable, hit_eff_data = pickle.load(f)
        return track_eff_data, all_data, eff_by_variable, hit_eff_data            

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
                        reco_fpath = os.path.join(reco_path, reco_file)
                        if not os.path.isfile(reco_fpath):
                            continue
                        with open(reco_fpath) as file:
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

                mass = sample_to_mass[sample]
                pt_max = np.sqrt(10**2/4 - mass**2) 
                good_reco_tracks = 0 

                theta_bins = np.linspace(0 - 0.1, 2.6+0.1, 21)
                theta_num = np.zeros(len(theta_bins)-1, int)
                theta_den = np.zeros(len(theta_bins)-1, int)
 
                phi_bins = np.linspace(-np.pi - 0.1, np.pi+0.1, 21)
                phi_num = np.zeros(len(phi_bins)-1, int)
                phi_den = np.zeros(len(phi_bins)-1, int)
 
                pt_bins = np.linspace(0, pt_max, 21)
                pt_num = np.zeros(len(pt_bins)-1, int)
                pt_den = np.zeros(len(pt_bins)-1, int)
 
                for reco_file in os.listdir(reco_path):
                    reco_fpath =os.path.join(reco_path, reco_file)
                    if not os.path.isfile(reco_fpath):
                        continue
                    with open(reco_fpath) as file: 
                        reco_data = json.load(file)
                        if get_chunk_id(reco_file) in bad_chunks: 
                            continue 

                    chi_sq = np.asarray(reco_data["match_track_info"]["chi_sq"])
                    ndf = np.asarray(reco_data["match_track_info"]["ndf"])
                    theta_truth = np.asarray(reco_data["match_track_info"]["theta"])
                    pt_truth = np.asarray(reco_data["match_stau_info"]["pt"])/1000 
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
                total_efficiency_error = np.sqrt((total_efficiency/100) * (1-(total_efficiency/100))/total_accepted_staus)
                print(f"\n ========== reco: track efficiency ==========")
                print(f"Efficiency: {total_efficiency:.2f}%") 
                                
                track_eff_data[window][sample]["acceptance"].append(acceptance_rate)
                if option == "bib/": 
                    track_eff_data[window][sample]["trackeff_bib"]["eff"].append(total_efficiency)
                    track_eff_data[window][sample]["trackeff_bib"]["err"].append(total_efficiency_error)
                else:
                    track_eff_data[window][sample]["trackeff_nobib"]["eff"].append(total_efficiency)
                    track_eff_data[window][sample]["trackeff_nobib"]["err"].append(total_efficiency_error)


                #################################### hit efficiency ####################################
                    
                if option == "bib/":
                    continue
                            
                layer_map = [
                    ("VB", 0, 30), ("VB", 1, 30), ("VB", 2, 51), ("VB", 3, 51), ("VB", 4, 74), ("VB", 5, 74), ("VB", 6, 102), ("VB", 7, 102),
                    ("IB", 0, 127), ("IB", 1, 340), ("IB", 2, 554), 
                    ("OB", 0, 819), ("OB", 1, 1153), ("OB", 2, 1486)
                ] 
 
                layer_radii = np.array([r for _,_,r in layer_map])
                n_layers_tot = len(layer_map)

                def exp_hits_mask(dist_mm: float)->np.ndarray:
                    return layer_radii <= dist_mm

                hit_eff = {p: {"exp": np.zeros(n_layers_tot, dtype=int),
                            "obs": np.zeros(n_layers_tot, dtype=int)} for p in ("stau", "antistau")}
                
                hit_level_path = os.path.join(reco_path, "hit_level") 
                for hit_reco_file in os.listdir(hit_level_path):
                    hit_fpath = os.path.join(hit_level_path, hit_reco_file)
                    with open(hit_fpath) as hit_file:
                        reco = json.load(hit_file)

                    travel_dist = np.asarray(reco["match_stau_info"]["travel_dist"])
                    pdg = np.asarray(reco["match_stau_info"]["id"])
                    idx_stau, = np.where(pdg > 0)
                    idx_antistau, = np.where(pdg < 0)
                    truth_idx = {"stau": idx_stau, "antistau": idx_antistau}
                    
                    for species, idx in truth_idx.items():
                        if idx.size == 0:
                            continue

                        for i in idx:
                            dist = float(travel_dist[i])
                            exp_hits = exp_hits_mask(dist).astype(int)
                            hit_eff[species]["exp"] += exp_hits

                        for det in ("VB", "IB", "OB"):
                            lay = np.asarray(reco["hits_from_mcp"][det]["layer"])
                            pdg_hit = np.asarray(reco["hits_from_mcp"][det]["pdg_id"])

                            for layer in np.unique(lay[pdg_hit > 0]) if species == "stau" else np.unique(lay[pdg_hit < 0]):
                                if det=="VB":
                                    global_idx = layer
                                elif det=="IB":
                                    global_idx = 8 + layer
                                elif det=="OB":
                                    global_idx = 11 + layer
                                hit_eff[species]["obs"][global_idx] += 1

                # eff_stau = hit_eff["stau"]["obs"] / hit_eff["stau"]["exp"]
                # eff_antistau = hit_eff["antistau"]["obs"] / hit_eff["antistau"]["exp"]
                
                eff_total = (hit_eff["stau"]["obs"] + hit_eff["antistau"]["obs"]) / (hit_eff["stau"]["exp"] + hit_eff["antistau"]["exp"])
                hit_eff_data[window][sample]["eff"] = eff_total
                error = np.sqrt(eff_total * (1-eff_total)/(hit_eff["stau"]["exp"] + hit_eff["antistau"]["exp"]))
                hit_eff_data[window][sample]["err"] = error

                                

                print("\n ============ Hit level efficiency per layer: ===========")
                print(hit_eff_data[window][sample])
                  
    print(f"Writing cache to {CACHE}")
    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump((track_eff_data, all_data, eff_by_variable, hit_eff_data), f, protocol=pickle.HIGHEST_PROTOCOL)
    return track_eff_data, all_data, eff_by_variable, hit_eff_data

track_eff_data, all_data, eff_by_variable, hit_eff_data = build_analysis(redo=rebuild)

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
            err_bib = []
            err_nobib = []
            for sample in samples:
                trackeff_bib.append(track_eff_data[window][sample]["trackeff_bib"]["eff"][0])
                trackeff_nobib.append(track_eff_data[window][sample]["trackeff_nobib"]["eff"][0])
                err_bib.append(track_eff_data[window][sample]["trackeff_bib"]["err"][0]*100)
                err_nobib.append(track_eff_data[window][sample]["trackeff_nobib"]["err"][0]*100) 
            marker = marker_map[window]
            
            ax.errorbar(mass_list, trackeff_bib, yerr=(err_bib), fmt=f"{marker}-", ms=6, color='orange', alpha=window_colors[window], label=f"{window} BIB")
            ax.errorbar(mass_list, trackeff_nobib, yerr=(err_nobib), fmt=f"{marker}:", ms=6, color="blue", alpha=window_colors[window], label=f"{window} no BIB")

        ax.set_xlabel("Mass [TeV]")
        ax.set_ylabel("Track reconstruction efficiency (%)")
        ax.set_title("Total tracking efficiency vs mass")
        ax.legend()    
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

        fig, axes = plt.subplots(n_rows, n_cols, harey=True, figsize=(3.5*n_cols, 2.8*n_rows)) 
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


    def eff_by_layer(pdf): 
        layer_names = ["VB0","VB1","VB2","VB3","VB4","VB5","VB6","VB7","IB0","IB1","IB2","OB0","OB1","OB2"]
        x = np.arange(len(layer_names))
        width = 0.25
        offset = {"tight": -width, "medium": 0.0, "loose": +width}
        colors = {"tight":"tab:purple", "medium": "tab:cyan", "loose":"tab:olive"}
        for sample in samples:
            fig, ax = plt.subplots(figsize=(9,3))
            for window in windows:
                eff = hit_eff_data[window][sample]["eff"]
                err = hit_eff_data[window][sample]["err"]
                if eff is None or len(eff)==0:
                    continue
                eff = np.asarray(eff, dtype=float)
                err = np.asarray(err, dtype=float)
                err_clean = np.where(np.isfinite(err), err, 0.0)
                
                ax.bar(x + offset[window], eff, width=width, yerr=err_clean, label=window, color=colors[window])

            ax.set_xticks(x, layer_names, ha="right", fontsize=8)
            ax.set_ylim(0,1.05)
            ax.set_ylabel("Hit reconstruction efficiency")
            ax.set_title(f"Hit-level efficiency per layer for {sample}")
            for xc in (7.5, 10.5):
                ax.axvline(x=xc, color='k', ls='--', alpha=0.5)
            ax.legend(title="window", loc="upper right")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    with PdfPages(save_plot_path) as pdf:
        # mass_vs_acceptance(pdf)
        # mass_vs_trackeff(pdf)
        # eff_vs_theta(pdf)
        eff_vs_pt(pdf)
        # eff_vs_phi(pdf)
        # eff_by_layer(pdf)
        
    print(f"Saved plot(s) to {save_plot_path}")
