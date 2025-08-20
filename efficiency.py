import json
import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import pickle
import pathlib
from tqdm import tqdm

# sim_dir = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/sim/"
reco_dir = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/"
reco_dir_mira = "/scratch/miralittmann/analysis/mira_analysis_code/mira_time/v3/"
save_plot_path = "/scratch/miralittmann/analysis/efficiency_plots/new_time.pdf"
CACHE = pathlib.Path("cache/efficiency.pkl")

nhits3p5 = False
refitted = False

sample_to_mass = {
    "1000_10": 1.0,
    "2500_10": 2.5,
    "3000_10": 3.0,
    "3500_10": 3.5,
    "4000_10": 4.0,
    "4500_10": 4.5
}

mass_list = [1.0, 2.5, 3.0, 3.5, 4.0, 4.5]
windows = ["tight", "medium", "loose", "mira_time"]
# windows = ["pt10", "pt5"]
samples = ["1000_10", "2500_10", "3000_10", "3500_10", "4000_10", "4500_10"]
# bib_options = ["bib", "nobib"]
bib_options = ["nobib"]
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
                sample: {
                    "nhits": [],
                    "chi2": []
                } for sample in samples
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
                },
                "chi2tracks": [],
                "alltracks": []
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
                    "eta": {
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
    chi2bib = []
    chi2nobib = []

    if CACHE.exists() and not redo:
        print(f"Loading previous info from {CACHE}, not redoing full analysis")
        with CACHE.open("rb") as f:
            track_eff_data, all_data, eff_by_variable, hit_eff_data = pickle.load(f)
        return track_eff_data, all_data, eff_by_variable, hit_eff_data            

    print("Rebuilding analysis arrays...")
    for window in tqdm(windows):
        for option in bib_options: 
            for sample in samples:
                if refitted == True:
                    reco_path = os.path.join(reco_dir, option, window, sample, "refitted")
                elif nhits3p5 == True:
                    reco_path = os.path.join(reco_dir, option, window, sample, "3p5hits")
                elif window == "mira_time":
                    reco_path = os.path.join(reco_dir_mira, option, sample)  
                else:
                    reco_path = os.path.join(reco_dir, option, window, sample)  

                print(reco_path)
                print(f"{len(os.listdir(reco_path))} files")
                                                            
                events_data = []

                ##################################### setup + acceptance (from sim info) #############################
                bad_chunks = set()
                for reco_file in os.listdir(reco_path):
                    reco_fpath = os.path.join(reco_path, reco_file)
                    if not os.path.isfile(reco_fpath):
                        continue
                    with open(reco_fpath) as file:
                        chunk_reco_data = json.load(file)
                        if get_chunk_id(reco_file) in chunk_reco_data.get("bad_files", []):
                            bad_chunks.add(get_chunk_id(reco_file))
                
                ################################### tracking efficiency ######################################### 

                accepted_staus = 0
                mass = sample_to_mass[sample]
                pt_max = np.sqrt(10**2/4 - mass**2) 
                good_reco_tracks = 0 
                all_tracks = 0

                eta_bins = np.linspace(-0.81,0.81,51)
                eta_num = np.zeros(len(eta_bins)-1, int)
                eta_den = np.zeros(len(eta_bins)-1, int)

                theta_bins = np.linspace(0 - 0.1, 2.6+0.1, 51)
                theta_num = np.zeros(len(theta_bins)-1, int)
                theta_den = np.zeros(len(theta_bins)-1, int)

                phi_bins = np.linspace(-np.pi - 0.1, np.pi+0.1, 51)
                phi_num = np.zeros(len(phi_bins)-1, int)
                phi_den = np.zeros(len(phi_bins)-1, int)

                pt_bins = np.linspace(0, pt_max, 51)
                pt_num = np.zeros(len(pt_bins)-1, int)
                pt_den = np.zeros(len(pt_bins)-1, int)

                for reco_file in os.listdir(reco_path):
                    reco_fpath =os.path.join(reco_path, reco_file)
                    if not os.path.isfile(reco_fpath):
                        continue
                    with open(reco_fpath) as file: 
                        reco_data = json.load(file) 
                        # if get_chunk_id(reco_file) in bad_chunks: 
                        #     continue 
                    
                        accepted_staus += reco_data["accepted_staus"]
                        chi_sq = np.asarray(reco_data["match_track_info"]["chi_sq"]) 
                        ndf = np.asarray(reco_data["match_track_info"]["ndf"])
                        all_data[window][option][sample]["chi2"].append(chi_sq / ndf)
                        n_hits = np.asarray(reco_data["match_track_info"]["n_hits"])
                        all_data[window][option][sample]["nhits"].append(n_hits)
                        
                        theta_truth = np.asarray(reco_data["match_stau_info"]["theta"])
                        pt_truth = np.asarray(reco_data["match_stau_info"]["pt"])/1000 
                        phi_truth = np.asarray(reco_data['match_stau_info']["phi"])
                        eta_truth = -np.log(np.tan(theta_truth/2))

                        theta_den += np.histogram(theta_truth, bins=theta_bins)[0]
                        phi_den += np.histogram(phi_truth, bins=phi_bins)[0]
                        pt_den += np.histogram(pt_truth, bins=pt_bins)[0]
                        eta_den += np.histogram(eta_truth, bins=eta_bins)[0]
                    
                        red_chi_sq = chi_sq / ndf 
                        good_mask = red_chi_sq < 5
                        if option == "nobib":
                            nhits_mask = n_hits > 3.5
                        else:  # "bib"
                            nhits_mask = n_hits > 7

                        good_mask = (red_chi_sq < 5) & nhits_mask
                        good_reco_tracks += good_mask.sum()  

                        theta_reco = theta_truth[good_mask]
                        pt_reco = pt_truth[good_mask]
                        phi_reco = phi_truth[good_mask]
                        eta_reco = eta_truth[good_mask]

                        theta_num += np.histogram(theta_reco, bins=theta_bins)[0]
                        phi_num += np.histogram(phi_reco, bins=phi_bins)[0]
                        pt_num += np.histogram(pt_reco, bins=pt_bins)[0]
                        eta_num += np.histogram(eta_reco, bins=eta_bins)[0] 

                theta_eff = np.divide(theta_num, theta_den, out=np.zeros_like(theta_num, float), where=theta_den>0) 
                phi_eff = np.divide(phi_num, phi_den, out=np.zeros_like(phi_num, float), where=phi_den>0)                     
                pt_eff = np.divide(pt_num, pt_den, out=np.zeros_like(pt_num, float), where=pt_den>0)            
                eta_eff = np.divide(eta_num, eta_den, out=np.zeros_like(eta_num, float), where=eta_den>0)

                theta_err = np.zeros_like(theta_eff)
                theta_m = theta_den > 0
                theta_err[theta_m] = np.sqrt(theta_eff[theta_m] * (1-theta_eff[theta_m])/theta_den[theta_m])
                theta_centers = 0.5*(theta_bins[1:] + theta_bins[:-1])
                eff_by_variable[window][sample][option]["theta"]["eff"].append(theta_eff)                
                eff_by_variable[window][sample][option]["theta"]["err"].append(theta_err)
                eff_by_variable[window][sample][option]["theta"]["centers"].append(theta_centers)

                eta_err = np.zeros_like(eta_eff)
                eta_m = eta_den > 0
                eta_err[eta_m] = np.sqrt(eta_eff[eta_m] * (1-eta_eff[eta_m])/eta_den[eta_m])
                eta_centers = 0.5*(eta_bins[1:] + eta_bins[:-1])
                eff_by_variable[window][sample][option]["eta"]["eff"].append(eta_eff)                
                eff_by_variable[window][sample][option]["eta"]["err"].append(eta_err)
                eff_by_variable[window][sample][option]["eta"]["centers"].append(eta_centers)

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

                total_efficiency =  good_reco_tracks / accepted_staus * 100 
                track_eff_data[window][sample]["chi2tracks"].append(good_reco_tracks)
                track_eff_data[window][sample]["alltracks"].append(all_tracks)
                track_eff_data[window][sample]["acceptance"].append(accepted_staus) 
                total_efficiency_error = np.sqrt((total_efficiency/100) * (1-(total_efficiency/100))/accepted_staus)
                print(f"\n ========== {window} {sample}, {option}: track efficiency ==========")
                print(f"Efficiency: {total_efficiency:.2f}%") 
                print("accepted staus:", accepted_staus)
                print(good_reco_tracks, "good reco tracks")
                                
                if option == "bib": 
                    track_eff_data[window][sample]["trackeff_bib"]["eff"].append(total_efficiency)
                    track_eff_data[window][sample]["trackeff_bib"]["err"].append(total_efficiency_error)
                else:
                    track_eff_data[window][sample]["trackeff_nobib"]["eff"].append(total_efficiency)
                    track_eff_data[window][sample]["trackeff_nobib"]["err"].append(total_efficiency_error)


                #################################### hit efficiency ####################################
                    
                if option == "bib":
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

                for reco_file in os.listdir(reco_path):
                    reco_fpath =os.path.join(reco_path, reco_file)
                    if not os.path.isfile(reco_fpath):
                        continue
                    with open(reco_fpath) as file: 
                        reco = json.load(file)     
                    
                    travel_dist = np.asarray(reco["match_stau_info"]["travel_dist"])  
                    pdg = np.asarray(reco["match_stau_info"]["id"])
                    theta_truth = np.asarray(reco["match_stau_info"]["theta"])
                    r_reach     = travel_dist * np.sin(theta_truth)
                    acc_truth   = (r_reach >= 127.0)
                    
                    idx_stau, = np.where(pdg > 0)
                    idx_antistau, = np.where(pdg < 0)
                    acc_idx_stau     = np.where((pdg > 0) & acc_truth)[0]
                    acc_idx_antistau = np.where((pdg < 0) & acc_truth)[0]
                    
                    for i in acc_idx_stau:
                        rr = float(r_reach[i])
                        hit_eff["stau"]["exp"] += (layer_radii <= rr).astype(int)
                    for i in acc_idx_antistau:
                        rr = float(r_reach[i])
                        hit_eff["antistau"]["exp"] += (layer_radii <= rr).astype(int)

                    for det in ("VB", "IB", "OB"):
                        lay     = np.asarray(reco["hits_from_mcp"][det]["layer"])
                        pdg_hit = np.asarray(reco["hits_from_mcp"][det]["pdg_id"])

                        n = min(lay.size, pdg_hit.size)
                        if n == 0:
                            continue
                        lay, pdg_hit = lay[:n], pdg_hit[:n]

                        if acc_idx_stau.size > 0:
                            rr_max = float(np.max(r_reach[acc_idx_stau]))
                            m = (pdg_hit > 0)
                            if m.any():
                                for layer in np.unique(lay[m]):
                                    if det == "VB":   global_idx = layer
                                    elif det == "IB": global_idx = 8  + layer
                                    else:             global_idx = 11 + layer
                                    if layer_radii[global_idx] <= rr_max:
                                        hit_eff["stau"]["obs"][global_idx] += 1

                        if acc_idx_antistau.size > 0:
                            rr_max = float(np.max(r_reach[acc_idx_antistau]))
                            m = (pdg_hit < 0)
                            if m.any():
                                for layer in np.unique(lay[m]):
                                    if det == "VB":   global_idx = layer
                                    elif det == "IB": global_idx = 8  + layer
                                    else:             global_idx = 11 + layer
                                    if layer_radii[global_idx] <= rr_max:
                                        hit_eff["antistau"]["obs"][global_idx] += 1
                             
                eff_total = (hit_eff["stau"]["obs"] + hit_eff["antistau"]["obs"]) / (hit_eff["stau"]["exp"] + hit_eff["antistau"]["exp"])
                hit_eff_data[window][sample]["eff"] = eff_total
                error = np.sqrt(eff_total * (1-eff_total)/(hit_eff["stau"]["exp"] + hit_eff["antistau"]["exp"]))*100
                hit_eff_data[window][sample]["err"] = error

                print("\n ============ Hit level efficiency per layer: ===========")
                print(hit_eff_data[window][sample])
                
    print(f"Writing cache to {CACHE}")
    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump((track_eff_data, all_data, eff_by_variable, hit_eff_data), f, protocol=pickle.HIGHEST_PROTOCOL)
    return track_eff_data, all_data, eff_by_variable, hit_eff_data

track_eff_data, all_data, eff_by_variable, hit_eff_data = build_analysis(redo=rebuild)

#print("10:", track_eff_data["pt10"]["4000_10"]["chi2tracks"], track_eff_data["pt10"]["4000_10"]["alltracks"])
#print("5:", track_eff_data["pt5"]["4000_10"]["chi2tracks"], track_eff_data["pt5"]["4000_10"]["alltracks"])
colors = {"tight":"maroon", "medium": "teal", "loose":"palevioletred", "mira_time": "goldenrod", "pt10": "tab:blue", "pt5": "tab:green"}
########################################## PLOTTING ######################################
if plotting == True:
    print("Now making plots...")
    marker_map = {"tight": "v", "medium":"o", "loose":"^", "mira_time": "*", "pt10":"o", "pt5": "."}
    color_map = { 
        "refitted": {"bib": "lightcoral", "nobib": "tomato"},
        "raw": {"bib": "royalblue", "nobib":"skyblue"}
    }
    window_alphas = {
        'tight': 0.3,
        'medium': 0.6,
        'loose': 0.9
    }

    LABEL_FONTSIZE = 15      
    TICK_FONTSIZE  = 13      
    NOTE_FONTSIZE  = 10  

    # def chi2_hist(pdf):
    #     fig1, ax1 = plt.subplots()
    #     num_bins = 50
    #     bin_lims = np.linspace(0, 6, num_bins+1)
    #     bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
    #     bin_widths=bin_lims[1:]-bin_lims[:-1]

    #     ref_hist, _ = np.histogram(np.log10(refitted_chi2), bins=bin_lims)
    #     raw_hist, _ = np.histogram(np.log10(raw_chi2), bins=bin_lims)

    #     ref_norm = ref_hist / len(refitted_chi2)
    #     raw_norm = raw_hist / len(raw_chi2)

    #     ax1.bar(bin_centers, ref_norm, width=bin_widths, label='refitted', color='red', alpha=0.5)
    #     ax1.bar(bin_centers, raw_norm, width=bin_widths, label='raw', color='blue', alpha=0.5)
    #     ax1.axvline(x=np.log10(5), linestyle='--')
    #     ax1.set_xlabel("log_10 of track $\chi^2$ ")
    #     ax1.set_ylabel("normalized counts")
    #     ax1.set_title("normalized track $\chi^2$")


    #     plt.yscale("log")
    #     ax1.legend() 
    #     fig1.tight_layout()
    #     pdf.savefig(fig1)
    #     plt.close(fig1)
        
    # def n_hits(pdf):
    #     fig0, ax0 = plt.subplots()
    #     bin_lims = np.linspace(0, 21, 40)
    #     bin_centers = 0.5*(bin_lims[:-1]+bin_lims[1:])
    #     bin_widths = bin_lims[1:]-bin_lims[:-1]

    #     ref_hist, _ = np.histogram(refitted_nhits, bins=bin_lims)
    #     raw_hist, _ = np.histogram(raw_nhits, bins=bin_lims)
    #     ref_norm = ref_hist / len(refitted_nhits)
    #     raw_norm = raw_hist / len(raw_nhits)

    #     ax0.bar(bin_centers, ref_norm, width=bin_widths, alpha=0.5, label='refitted', color='red')
    #     ax0.bar(bin_centers, raw_norm, width=bin_widths, alpha=0.5, label='raw', color='blue')
        # # print(n_bins)
        # # ax0.hist(refitted_nhits, bins=n_bins, histtype='step', label="refitted", color='red')
        # # ax0.hist(raw_nhits, bins=n_bins, histtype='step', label='raw', color='blue')
        # ax0.legend()
        # ax0.set_title('Number of total hits per track')
        # ax0.set_xlabel("number of hits")
        # ax0.set_ylabel('track counts')
        # fig0.tight_layout()
        # pdf.savefig(fig0)
        # plt.close(fig0)

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
            
            ax.errorbar(mass_list, trackeff_bib, yerr=(err_bib), fmt=f"{marker}-", ms=9, color=colors[window], label=f"{window} BIB", linewidth=3)
            ax.errorbar(mass_list, trackeff_nobib, yerr=(err_nobib), fmt=f"{marker}:", ms=9, color=colors[window], alpha=0.5, label=f"{window} no BIB", linewidth=3)

        ax.set_xlabel("Stau mass [TeV]", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax.set_ylabel("Track reconstruction efficiency (%)", fontsize=LABEL_FONTSIZE, labelpad=4)
        ax.tick_params(axis='both', labelsize=TICK_FONTSIZE)
        # ax.set_title("Track reconstruction efficiency by stau mass")
        ax.legend()    
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    # def mass_vs_acceptance(pdf):
    #     fig, ax = plt.subplots()
    #     acceptance = []
    #     for sample in samples: 
    #         acceptance.append(track_eff_data["medium"][sample]["acceptance"])
    #     ax.plot(mass_list, acceptance)
    #     ax.legend()
    #     ax.set_xlabel("Mass [TeV]")
    #     ax.set_ylabel("Percentage of accepted staus")
    #     ax.set_title("Acceptance")
    #     fig.tight_layout()
    #     pdf.savefig(fig)
    #     plt.close(fig)


    def eff_vs_theta(pdf):
        n_cols = len(samples)
        n_rows = len(windows)

        fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(3.5*n_cols, 2.8*n_rows)) 
        colors = {"bib":"orange", "nobib":"blue"}
        for r, window in enumerate(windows):
            for c, sample in enumerate(samples):
                ax = axes[r,c]
                for option in bib_options: 
                    arr = eff_by_variable[window][sample][option]["theta"] 
                    if not arr["centers"]:
                        continue
                    centers = arr["centers"][0]
                    eff = arr["eff"][0]
                    err = arr['err'][0]

                    ax.errorbar(centers, eff, yerr=err, label="BIB" if option=="bib" else "no BIB", color=colors[option])

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

    def eff_vs_eta(pdf):
        n_cols = len(samples)
        n_rows = len(windows)
        
        fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(3.5*n_cols, 2.8*n_rows)) 
        colors = {"bib":"orange", "nobib":"blue"}
        for r, window in enumerate(windows):
            for c, sample in enumerate(samples):
                ax = axes[r,c] 
                for option in bib_options:
                    arr = eff_by_variable[window][sample][option]["eta"]
                    if not arr["centers"]:
                        continue
                    centers = arr["centers"][0]
                    eff = arr["eff"][0]
                    err = arr['err'][0]

                    ax.errorbar(centers, eff, yerr=err, label="BIB" if option=="bib" else "no BIB", color=colors[option])

                    if r==0:
                        ax.set_title(f"{sample}")
                    if c==0:
                        ax.set_ylabel(f"{window}")
                    ax.grid(ls=":", alpha=0.4)

        handles, labels = axes[0,0].get_legend_handles_labels()
        fig.legend(handles, labels, loc="upper right")
        fig.suptitle("tracking efficiency vs truth eta")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


    def eff_vs_pt(pdf):
        n_cols = len(samples)
        n_rows = len(windows)

        fig, axes = plt.subplots(n_rows, n_cols, sharey=True, figsize=(3.5*n_cols, 2.8*n_rows)) 
        colors = {"bib":"orange", "nobib":"blue"}
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

                    ax.errorbar(centers, eff, yerr=err, label="BIB" if option=="bib" else "no BIB", color=colors[option])

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


#    def eff_vs_phi(pdf):
    #     n_cols = len(samples)
    #     n_rows = len(windows)

    #     fig, axes = plt.subplots(n_rows, n_cols, sharex=True, sharey=True, figsize=(3.5*n_cols, 2.8*n_rows)) 
    #     colors = {"bib/":"orange", "nobib/":"blue"}
    #     for r, window in enumerate(windows):
    #         for c, sample in enumerate(samples):
    #             ax = axes[r,c] if n_rows>1 else axes[c]
    #             for option in bib_options:
    #                 arr = eff_by_variable[window][sample][option]["phi"]
    #                 if not arr["centers"]:
    #                     continue
    #                 centers = arr["centers"][0]
    #                 eff = arr["eff"][0]
    #                 err = arr['err'][0]

    #                 ax.errorbar(centers, eff, yerr=err, label="BIB" if option=="bib/" else "no BIB", color=colors[option])

    #                 if r==0:
    #                     ax.set_title(f"{sample}")
    #                 if c==0:
    #                     ax.set_ylabel(f"{window}")
    #                 ax.grid(ls=":", alpha=0.4)

    #     handles, labels = axes[0,0].get_legend_handles_labels()
    #     fig.legend(handles, labels, loc="upper right")
    #     fig.suptitle("tracking efficiency vs phi")
    #     fig.tight_layout()
    #     pdf.savefig(fig)
    #     plt.close(fig)


    def eff_by_layer(pdf): 
        layer_names = ["VB0","VB1","VB2","VB3","VB4","VB5","VB6","VB7","IB0","IB1","IB2","OB0","OB1","OB2"]
        x = np.arange(len(layer_names))
        width = 0.20
        offset = {"tight": -1.5*width, "medium": -width/2, "loose": width/2, "mira_time": 1.5*width, "pt10": 1.5*width, "pt5": width/2}
        for sample in samples:
            fig, ax = plt.subplots(figsize=(9,3))
            for window in windows:
                eff = np.asarray(hit_eff_data[window][sample]["eff"])*100
                err = np.asarray(hit_eff_data[window][sample]["err"])
                if eff is None or len(eff)==0:
                    continue
                eff = np.asarray(eff, dtype=float)
                err = np.asarray(err, dtype=float)
                err_clean = np.where(np.isfinite(err), err, 0.0)
                
                ax.bar(x + offset[window], eff, width=width, yerr=err_clean, label=window, color=colors[window])

            ax.set_xticks(x, layer_names, ha="right", fontsize=8, rotation=45)
            ax.set_ylim(0,105)
            #ax.set_ylabel("Hit reconstruction efficiency", fontsize=12, labelpad=4)
            ax.tick_params(axis='both', labelsize=TICK_FONTSIZE)
            # ax.set_title(f"Hit-level efficiency per layer for {sample_to_mass[sample]} TeV staus")
            for xc in (7.5, 10.5):
                ax.axvline(x=xc, color='k', ls='--', alpha=0.5)
            ax.legend(title=f"{sample},window", loc="upper right")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)

    with PdfPages(save_plot_path) as pdf:
        # mass_vs_acceptance(pdf)
        # chi2_hist(pdf)
        mass_vs_trackeff(pdf)
        # n_hits(pdf)
        # eff_vs_theta(pdf)
        # eff_vs_pt(pdf)
        # eff_vs_phi(pdf)
        eff_by_layer(pdf)
        # eff_vs_eta(pdf)
        
    print(f"Saved plot(s) to {save_plot_path}")
