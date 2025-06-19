import json
import os
import numpy as np
import matplotlib.pyplot as plt
from itertools import chain

reco_path = "/scratch/miralittmann/analysis/10pbibjson/4000_10/medium_window/reco/7-18/"
sim_path = "/scratch/miralittmann/analysis/10pbibjson/4000_10/medium_window/sim/7-18/"

chunk = ["hit_sim_info", "mcp_stau_info"]
detectors = ["VB", "VE", "IB", "IE", "OB", "OE"]
just_barrel = ["VB", "IB", "OB"]
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
            #print(get_chunk_id(reco_file), "has a bad reco file -- skipping for analysis")
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

total_reco_stau_count = len(match_stau_info_RECO["id"])
approx_pcent = (total_reco_stau_count / (2*(len(os.listdir(reco_path)) - len(bad_chunks))))*100
print("total reco stau count:", total_reco_stau_count)
print("approximate percent of staus reconstructed:", approx_pcent)
print("found ", len(bad_chunks), " files with 0 events; excluded from analysis")

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

## consolidating sim info ##
coords = ('x', 'y', 'z')
total_space_vec_SIM = {
    coord: list(chain.from_iterable(hit_info_SIM[detector][coord] for detector in just_barrel))
    for coord in coords
}

xyz = np.column_stack([total_space_vec_SIM[coord] for coord in coords])
hit_radial_distance_SIM = np.linalg.norm(xyz, axis=1)

pdg_SIM = []
for v, i, o in zip(hit_info_SIM["VB"]["pdg"], hit_info_SIM["IB"]["pdg"], hit_info_SIM["OB"]["pdg"]):
    pdg_SIM.append(v + i + o)

# layer_SIM = []
# system_SIM = []
# for v, i, o in zip(hit_info_SIM["VB"]["layer"], hit_info_SIM["IB"]["layer"], hit_info_SIM["OB"]["layer"]):
#     layer_SIM.append(v + i + o)
#     system_SIM.append(["VB"]*len(v) + ["IB"]*len(i) + ["OB"]*len(o))
# # TODO: this is exactly what tate does in her script, but i'm not sure how useful it is because doesn't every detector 
# # have a different amount of layers / numbering system for it? also not working for me

all_ids_SIM = []
print(len(hit_info_SIM["IB"]["id"]))


stau_ids = {1000015, -1000015, 2000015, -2000015}

seen_staus = set() # (event index, mcp id)
for event_idx, (pdgs, ids) in enumerate(zip(pdg_SIM, all_id_SIM)):
    for pdg, id in zip(pdgs, ids):
        if pdg in stau_ids:
            seen_staus.add((event_idx, id))
print(len(seen_staus))

# seen_staus = set() # (event index, pdg number)
# for i, pdg in enumerate(pdg_SIM):
#     for pdg in set(pdg):
#         if pdg in stau_ids:
#             seen_staus.add((i, pdg))

# cut conditions : stau, 127 <= r <340, (first?) layer == 0
