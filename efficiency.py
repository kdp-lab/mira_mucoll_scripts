import json
import os
import numpy as np
import matplotlib.pyplot as plt

reco_path = "/scratch/miralittmann/analysis/10pbibjson/4000_10/medium_window/reco/"
sim_path = "/scratch/miralittmann/analysis/10pbibjson/4000_10/medium_window/sim/"

systems_key = ["VB", "VE", "IB", "IE", "OB", "OE"]
fields = ["x", "y", "z", "layer", "pdg"] ## TODO: why don't we care about time ? aren't we trying to make timing cuts
all_sim_data = {key: {f: [] for f in fields} for key in systems_key}

### Getting sim info ###
for sim_file in os.listdir(sim_path):
    with open(os.path.join(sim_path, sim_file)) as file:
        chunk_sim_data = json.load(file)
        for system in systems_key:
            for field in fields:
                all_sim_data[system][field].append(chunk_sim_data["hit_info"][system][field])
                # will it let me do this or it needs to be one by one because list or whatever

added_sim_data = {f: [] for f in fields} ##TODO: can't add these for layer? string value? 
for field in fields:
    added_sim_data[field].append(all_sim_data["VB"][field] + all_sim_data["IB"][field] + all_sim_data["OB"][field])

radius = (added_sim_data["x"]**2 + added_sim_data["y"]**2 + added_sim_data["z"]**2)

# ### Getting reco info ###
# for reco_file in os.listdir(reco_path):
#     with open(os.path.join(reco_path, reco_file)) as file:
#         chunk_reco_data = json.load(file)
#         for field in fields: 


       