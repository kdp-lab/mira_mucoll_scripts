import os
import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import bisect

sim_file = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/sim/1000_10/1000_10_sim10.json"
reco_files = ["/scratch/miralittmann/analysis/mira_analysis_code/efficiency/nobib/tight/1000_10/1000_10_reco10.json", "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/nobib/tight/1000_10/1000_10_reco0.json"]

detectors = ["VB", "IB", "OB"]
dist_to_hits_map = {
    30: 2,
    51: 4,
    74: 6,
    102: 8,
    127: 9,
    340: 10,
    554: 11,
    819: 12,
    1153: 13,
    1486: 14
}
keys = sorted(dist_to_hits_map)
values = [dist_to_hits_map[k] for k in keys]

total_all_hit_effs_per_stau = []

 
for reco_file in reco_files:
    with open(reco_file) as file:
        reco_data = json.load(file) 

        # this is per stau
        for i in range(len(reco_data["match_stau_info"]["id"])):
            
            total_n_hits_reco = reco_data["match_track_info"]["n_hits"][i]
            n_hits_vertex = reco_data["match_track_info"]['n_hits_vertex'][i] 
            n_hits_inner = reco_data["match_track_info"]['n_hits_inner'][i]
            n_hits_outer = reco_data["match_track_info"]['n_hits_outer'][i]

            travel_dist = reco_data["match_stau_info"]["travel_dist"][i]
            idx = bisect.bisect_right(keys, travel_dist) - 1
            n_hits_sim = values[idx] if idx > 0 else 0

            total_hit_eff_per_stau = total_n_hits_reco / n_hits_sim
            total_all_hit_effs_per_stau.append(total_hit_eff_per_stau)



# efficiency per variable. number of staus whos variable is in bin X and can be reconstructed / number of truth staus whos variable is in bin X 
# using eta as a first try



truth_etas = []

for reco_file in reco_files:       
    with open(reco_file) as file:
        reco_data = json.load(file)
        truth_etas.append(reco_data["match_stau_info"]["eta"])

total_good_reco_tracks = 0
truth_etas_forbins = np.concatenate(truth_etas)
eta_bins = np.linspace(np.min(truth_etas_forbins) - 0.1, np.max(truth_etas_forbins) + 0.1, 21)
eta_num = np.zeros(len(eta_bins)-1, int)
eta_den = np.zeros(len(eta_bins)-1, int)

for reco_file in reco_files:       
    with open(reco_file) as file:
        reco_data = json.load(file)
    
    all_truth_etas = np.asarray(reco_data["match_stau_info"]["eta"])
    
    eta_den += np.histogram(all_truth_etas, bins=eta_bins)[0]

                
    chi_sq = np.asarray(reco_data["match_track_info"]["chi_sq"])
    ndf = np.asarray(reco_data["match_track_info"]["ndf"]) 
    red_chi_sq = chi_sq / ndf 
    good_mask = red_chi_sq < 5

    total_good_reco_tracks += good_mask.sum()

    theta_reco = np.asarray(reco_data["match_track_info"]["theta"])[good_mask]
    eta_reco = all_truth_etas[good_mask]
    eta_num += np.histogram(eta_reco, bins=eta_bins)[0]

eff = np.divide(eta_num, eta_den, out=np.zeros_like(eta_num, float), where=eta_den>0)
err = np.zeros_like(eff)
m = eta_den > 0
err[m] = np.sqrt(eff[m] * (1-eff[m])/eta_den[m])
centers = 0.5*(eta_bins[1:] + eta_bins[:-1])
print(err)
print(total_good_reco_tracks)
    
        
print(eta_den)
print(eta_num)
print(np.concatenate(truth_etas))


    




