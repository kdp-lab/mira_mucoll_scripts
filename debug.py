import os
import json
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import bisect

sim_file = "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/sim/1000_10/1000_10_sim10.json"
reco_files = ["/scratch/miralittmann/analysis/mira_analysis_code/efficiency/nobib/loose/1000_10/hit_level/1000_10_reco200.json", "/scratch/miralittmann/analysis/mira_analysis_code/efficiency/nobib/loose/1000_10/hit_level/1000_10_reco0.json"]

# detectors = ["VB", "IB", "OB"]
# dist_to_hits_map = {
#     30: 2,
#     51: 4,
#     74: 6,
#     102: 8,
#     127: 9,
#     340: 10,
#     554: 11,
#     819: 12,
#     1153: 13,
#     1486: 14
# }

layer_map = [
    ("VB", 0, 30), ("VB", 1, 30), ("VB", 2, 51), ("VB", 3, 51), ("VB", 4, 74), ("VB", 5, 74), ("VB", 6, 102), ("VB", 7, 102),
    ("IB", 0, 127), ("IB", 1, 340), ("IB", 2, 554), 
    ("OB", 0, 819), ("OB", 1, 1153), ("OB", 2, 1486)
]

layer_names = [f"{det}{lay}" for det, lay, _ in layer_map]
layer_radii = np.array([r for _,_,r in layer_map])
n_layers_tot = len(layer_map)

def exp_hits_mask(dist_mm: float)->np.ndarray:
    return layer_radii <= dist_mm

hit_eff = {p: {"exp": np.zeros(n_layers_tot, dtype=int),
               "obs": np.zeros(n_layers_tot, dtype=int)} for p in ("stau", "antistau")}

for reco_file in reco_files:
    with open(reco_file) as file:
        reco = json.load(file)

    print(reco["match_stau_info"]["pt"])
    
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
        
eff_stau = hit_eff["stau"]["obs"] / hit_eff["stau"]["exp"]
eff_antistau = hit_eff["antistau"]["obs"] / hit_eff["antistau"]["exp"]
eff_total = (hit_eff["stau"]["obs"] + hit_eff["antistau"]["obs"]) / (hit_eff["stau"]["exp"] + hit_eff["antistau"]["exp"])
print(eff_total)



