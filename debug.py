import os
import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import bisect

sim_file = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/sim/1000_10/1000_10_sim10.json"
reco_files = ["/scratch/miralittmann/analysis/mira_analysis_code/efficiency/nobib/loose/1000_10/hit_level/1000_10_reco10.json", "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/nobib/loose/1000_10/hit_level/1000_10_reco0.json"]

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
keys = np.array(sorted(dist_to_hits_map))
values = np.array([dist_to_hits_map[k] for k in keys])

# hit based efficiency. 

def expected_hits(dist_array):
    idx = np.searchsorted(keys, dist_array, side="right") -1
    idx = np.clip(idx, 0, None)
    return values[idx]

total_hits_eff_num = 0  # good hits from reco info
total_hits_eff_den = 0  # expected hits from sim info

total_n_hits_by_det = {
    "VB": [],
    "IB": [],
    "OB": []
}

for reco_file in reco_files:
    with open(reco_file) as file:
        reco_data = json.load(file) 

    truth_total_travel_dist = np.asarray(reco_data["match_stau_info"]["travel_dist"])
    truth_pdg = np.asarray(reco_data["match_stau_info"]["id"])
    stau_idx = np.where(truth_pdg > 0)[0]
    antistau_idx = np.where(truth_pdg< 0)[0]

    stau_total_truth_travel_dist = truth_total_travel_dist[stau_idx]
    antistau_total_truth_travel_dist = truth_total_travel_dist[antistau_idx]
    print(stau_total_truth_travel_dist, antistau_total_truth_travel_dist)

    truth_stau_nhits = expected_hits(stau_total_truth_travel_dist)
    truth_antistau_nhits = expected_hits(antistau_total_truth_travel_dist) 

    for detector in detectors: 
        pdg = np.asarray(reco_data["hits_from_mcp"][detector]["pdg_id"])
        layer_hits = np.asarray(reco_data["hits_from_mcp"][detector]["layer"]) 

        staus_idx = np.where(pdg > 0)[0]
        antistaus_idx = np.where(pdg < 0)[0]

        stau_layer_hits_bydet = layer_hits[staus_idx]
        antistau_layer_hits_bydet = layer_hits[antistaus_idx]

        total_reco_stau_hits = len(stau_layer_hits_bydet)
        total_reco_antistau_hits = len(antistau_layer_hits_bydet)  

        stau_unique_mask = np.unique(stau_layer_hits_bydet, return_index=True)[1]
        staus_oneperlayer = stau_layer_hits_bydet[np.sort(stau_unique_mask)]
        total_staus_oneperlayer = len(staus_oneperlayer)

        antistau_unique_mask = np.unique(antistau_layer_hits_bydet, return_index=True)[1]
        antistaus_oneperlayer = antistau_layer_hits_bydet[np.sort(antistau_unique_mask)]
        total_antistaus_oneperlayer = len(antistaus_oneperlayer)


        





