import json
import os
import numpy as np
from collections import defaultdict

sim_path = "/scratch/miralittmann/analysis/efficiency_v1/sim/1000_10/"
reco_path = "/scratch/miralittmann/analysis/efficiency_v1/nobib/medium/1000_10/"

def get_chunk_id(fname: str) -> int:
    base = os.path.basename(fname)
    name = base.split(".", 1)[0]
    tail = name.split("_")[-1]
    digits = ''.join(ch for ch in tail if ch.isdigit())
    return int(digits)

# Get bad chunks
bad_chunks = set()
for reco_file in os.listdir(reco_path):
    with open(os.path.join(reco_path, reco_file)) as file:
        chunk_reco_data = json.load(file)
        if get_chunk_id(reco_file) in chunk_reco_data.get("bad_files", []):
            bad_chunks.add(get_chunk_id(reco_file))

print(f"Bad chunks to exclude: {len(bad_chunks)}")

stau_ids = {1000015, -1000015, 2000015, -2000015}

events_data = []
total_files_processed = 0

for sim_file in os.listdir(sim_path):
    chunk_id = get_chunk_id(sim_file)
    if chunk_id in bad_chunks:
        continue
        
    with open(os.path.join(sim_path, sim_file)) as file:
        chunk_sim_data = json.load(file)
        total_files_processed += 1
        event_data = {
            'file': sim_file,
            'truth_staus': chunk_sim_data["mcp_stau_info"]["id"],
            'hit_info': chunk_sim_data["hit_info"]
        }
        events_data.append(event_data)

print(f"Files processed: {total_files_processed}")
print(f"Total events: {len(events_data)}")

total_truth_staus = sum(len(event['truth_staus']) for event in events_data)
print(f"Total truth staus: {total_truth_staus}")

total_accepted_staus = 0
total_staus_with_any_tracker_hit = 0

for event_idx, event in enumerate(events_data):
    event_stau_hits = defaultdict(list)
    for system in ["VB", "VE", "IB", "IE", "OB", "OE"]:
        hits = event['hit_info'][system]
        for i in range(len(hits['id'])):
            if hits['pdg'][i] in stau_ids:
                stau_id = hits['id'][i]
                event_stau_hits[stau_id].append({
                    'system': system,
                    'layer': hits['layer'][i],
                    'x': hits['x'][i],
                    'y': hits['y'][i],
                    'z': hits['z'][i],
                    'pdg': hits['pdg'][i]
                })
    total_staus_with_any_tracker_hit += len(event_stau_hits)

    accepted_in_this_event = set()
    for stau_id, hits in event_stau_hits.items():
        ib_hits = [hit for hit in hits if hit['system'] == 'IB']
        if ib_hits:
            min_layer = min(hit['layer'] for hit in ib_hits)
            if min_layer <= 0:  
                accepted_in_this_event.add(stau_id)
    total_accepted_staus += len(accepted_in_this_event)

print(f"\n=== (SIM) ===")
print(f"Total truth staus: {total_truth_staus}")
print(f"Total staus with any tracker hit: {total_staus_with_any_tracker_hit}")
print(f"Total accepted staus: {total_accepted_staus}")
if total_truth_staus > 0:
    acceptance_rate = total_accepted_staus / total_truth_staus * 100
    print(f"Acceptance rate: {acceptance_rate:.2f}%")

total_reco_staus = 0
print(reco_path)
for reco_file in os.listdir(reco_path):
    with open(os.path.join(reco_path, reco_file)) as file: 
        reco_data = json.load(file)
        if get_chunk_id(reco_file) in bad_chunks: 
            continue 
        total_reco_staus += len(reco_data.get("match_stau_info", {}).get("id", []))

if total_accepted_staus > 0:
    efficiency = total_reco_staus / total_accepted_staus * 100
    print(f"\n=== (RESULTS) ===")
    print(f"Efficiency: {efficiency:.2f}%")
else:
    print("Cannot calculate efficiency: no accepted staus")

