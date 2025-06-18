import json
import os
import numpy as np
import matplotlib.pyplot as plt

reco_path = "/scratch/miralittmann/analysis/10pbibjson/4000_10/medium_window/reco/7-18/"
sim_path = "/scratch/miralittmann/analysis/10pbibjson/4000_10/medium_window/sim/7-18/"

chunk = ["hit_sim_info", "mcp_stau_info"]
detectors = ["VB", "VE", "IB", "IE", "OB", "OE"]
hit_fields_SIM = ['id', 'layer', 'pdg', 'stau', 'time', 'x', 'y', 'z']
hit_fields_RECO = ["x", "y", "z", "time", "corrected_time", "layer_hit", "side_hit", "mcp_id"]

hit_info_SIM = {detector: {hit_field: [] for hit_field in hit_fields_SIM} for detector in detectors}
mcp_stau_info_SIM = None 

match_stau_info_RECO = None
match_track_info_RECO = None

bad_chunks = set()
def get_chunk_id(fname:str) -> int:
    base = os.path.basename(fname)
    name = base.split(".", 1)[0]
    tail = name.split("_")[-1]
    digits = ''.join(ch for ch in tail if ch.isdigit())
    return int(digits)

### loading in reco info ###
for reco_file in os.listdir(reco_path):
    with open(os.path.join(reco_path, reco_file)) as file:
        chunk_reco_data = json.load(file)
        if get_chunk_id(reco_file) in chunk_reco_data.get("bad_files"):
            bad_chunks.add(get_chunk_id(reco_file))
            print(get_chunk_id(reco_file), "has a bad reco file -- skipping for analysis")
            continue

        # saving matched stau information
        if match_stau_info_RECO is None:
            match_stau_info_RECO = {field: [] for field in chunk_reco_data["match_stau_info"].keys()}
        for field, info in chunk_reco_data["match_stau_info"].items(): 
            match_stau_info_RECO[field].extend(info)
        
        # saving information of the matched tracks 
        if match_track_info_RECO is None: 
            match_track_info_RECO = {field: [] for field in chunk_reco_data["match_track_info"] if field != "hits"} 
            match_track_info_RECO["hits"] = {d:{f: [] for f in hit_fields_RECO} for d in detectors}   
        for field, info in chunk_reco_data["match_track_info"].items():
            if field == "hits":
                for d in detectors:
                    for f in hit_fields_RECO:
                        match_track_info_RECO["hits"][d][f].extend(info[d][f]) 
            else:
                match_track_info_RECO[field].extend(info)

### loading in sim info ###
for sim_file in os.listdir(sim_path):
    with open(os.path.join(sim_path, sim_file)) as file:
        if get_chunk_id(sim_file) in bad_chunks:
            continue
        chunk_sim_data = json.load(file)
        systems_key = chunk_sim_data["hit_info"].keys()
        for system in systems_key:
            for field in chunk_sim_data["hit_info"][system].keys():   
                hit_info_SIM[system][field].extend(chunk_sim_data["hit_info"][system][field])  
        if mcp_stau_info_SIM is None:
            mcp_stau_info_SIM = {field: [] for field in chunk_sim_data["mcp_stau_info"].keys()}
        for field, arr in chunk_sim_data["mcp_stau_info"].items():
            mcp_stau_info_SIM[field].extend(arr)

