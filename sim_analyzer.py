## (attempting to) make a script that only deals with the sim files, gives you an overview of the information in them directly before digi/reco

import pyLCIO
from pyLCIO import UTIL
import ROOT
import glob
import json
from math import *
import numpy as np
import argparse

# TODO: Add option for events loop argument, in sh script as well

parser = argparse.ArgumentParser(description="Analyze SLCIO chunks.")
parser.add_argument(
    '--chunk', 
    type=int, 
    required=True, 
    help="The chunk number to analyze"
)
parser.add_argument(
    '--all-events',            
    action='store_true',       
    help='Loop over every event in the chunk instead of only the first.'
)
args = parser.parse_args()

in_path = "/ospool/uc-shared/project/futurecolliders/miralittmann/sim/efficiency/2500_10/"
out_path = "/scratch/miralittmann/analysis/efficiency_v1/sim/2500_10/"
sample = "2500_10"
chunk = args.chunk
in_file = f"{in_path}{sample}_sim{chunk}.slcio"
out_file = f"{out_path}{sample}_sim{chunk}.json"

## constants ##
stau_ids = [1000015, 2000015]

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(in_file)
def event_looper(reader, all_events):
    if all_events:
        yield from reader
    else:
        for evt in reader:
            yield evt
            break


## The collections we have: ECalBarrelCollection, ECalEndcapCollection, HCalBarrelCollection, HCalEndcapCollection, HCalRingCollection, 
## InnerTrackerBarrelCollection, InnerTrackerEndcapCollection, MCParticle, OuterTrackerBarrelCollection, OuterTrackerEndcapCollection, 
## VertexBarrelCollection, VertexEndcapCollection, YokeBarrelCollection, YokeEndcapCollection

# I'm pretty sure we only care about: MCParticle, VertexBarrelCollection, VertexEndcapCollection, InnerTrackerBarrelCollection, InnerTrackerEndcapCollection, 
# OuterTrackerBarrelCollection, OuterTrackerEndcapCollection
# These have getTypeName() "SimTrackerHit"

### empty lists for things we care about ###

systems = [
    ("VertexBarrelCollection", "VB"),
    ("VertexEndcapCollection", "VE"),
    ("InnerTrackerBarrelCollection", "IB"),
    ("InnerTrackerEndcapCollection", "IE"),
    ("OuterTrackerBarrelCollection", "OB"),
    ("OuterTrackerEndcapCollection", "OE")
]

# making empty lists for the things we want to save from each system's collection, add to like hit_info["VB"]["x"].append( ... ) etc

mcp_stau_info = {
    "evvt_idx": [],
    "id": [],
    "p_tot": [],
    "p_x": [],
    "p_y": [],
    "p_z": [],
    "energy": [],
    "m": [],
    "beta": [],
    "pt": [],
    "eta": [],
    "phi": [],
    "theta": [],
    "travel_distance": []
}
systems_key  = ["VB", "VE", "IB", "IE", "OB", "OE"]
fields = ["id", "pdg", "x", "y", "z", "time", "layer", "stau"]

hit_info = {key: {f: [] for f in fields} for key in systems_key}
n_accepted_staus = 0
n_truth_staus = 0

for evt_idx, event in enumerate(event_looper(reader, args.all_events)):   
    ### MCP information ###
    mcps = event.getCollection("MCParticle")
    for mcp in mcps: 
        # skipping staus that have stau parents, staus that decay immediately, and particles that aren't staus
        if any(abs(parent.getPDG()) in stau_ids for parent in mcp.getParents()):
            continue
        if abs(mcp.getPDG()) not in stau_ids:
            continue
        if mcp.getGeneratorStatus == 22:
            print("intermediate stau")
            continue
        
        n_truth_staus +=1

        mcp_stau_momentum = mcp.getMomentum() 
        mcp_stau_tlv = ROOT.TLorentzVector()
        mcp_stau_tlv.SetPxPyPzE(mcp_stau_momentum[0], mcp_stau_momentum[1], mcp_stau_momentum[2], mcp.getEnergy()) 
        mcp_stau_beta = mcp_stau_tlv.Beta()
        mcp_stau_endpoint_r = sqrt(mcp.getEndpoint()[0] ** 2 + mcp.getEndpoint()[1] ** 2)
        mcp_stau_endpoint_z = mcp.getEndpoint()[2] 

        mcp_stau_info["id"].append(mcp.id())        
        mcp_stau_info["p_tot"].append(mcp_stau_tlv.Mag())
        mcp_stau_info["p_x"].append(mcp_stau_momentum[0])
        mcp_stau_info["p_y"].append(mcp_stau_momentum[1])
        mcp_stau_info["p_z"].append(mcp_stau_momentum[2])
        mcp_stau_info["energy"].append(mcp.getEnergy())
        mcp_stau_info["m"].append(mcp.getMass())
        mcp_stau_info["beta"].append(mcp_stau_beta)
        mcp_stau_info["pt"].append(mcp_stau_tlv.Perp())
        mcp_stau_info["eta"].append(mcp_stau_tlv.Eta())
        mcp_stau_info["phi"].append(mcp_stau_tlv.Phi())
        mcp_stau_info["theta"].append(mcp_stau_tlv.Theta())
        mcp_stau_info["travel_distance"].append(sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2 + mcp.getEndpoint()[2]**2) - sqrt(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2 + mcp.getVertex()[2]**2))
        
        print(mcp_stau_endpoint_r)
        if mcp_stau_endpoint_r > 102.0 and abs(mcp_stau_tlv.Eta()) < 1.0:
            n_accepted_staus += 1

    ### Track information ###
    for system, key in systems:
        # looking at all of the hits in a particular system 
        hits = event.getCollection(system)
        encoding = hits.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
        decoder = pyLCIO.UTIL.BitField64(encoding)
        for hit in hits:
            # saving the info of the particle that made each hit, if it's a stau or not
            mcp = hit.getMCParticle()
            hit_info[key]["id"].append(mcp.id())
            hit_info[key]["pdg"].append(mcp.getPDG())         
            if np.abs(mcp.getPDG()) in stau_ids:
                hit_info[key]["stau"].append("yes")
            else:
                hit_info[key]["stau"].append("no")
            
            # saving the info of the actual hit
            hit_info[key]["x"].append(hit.getPosition()[0])
            hit_info[key]["y"].append(hit.getPosition()[1])
            hit_info[key]["z"].append(hit.getPosition()[2])
            hit_info[key]["time"].append(hit.getTime())
            
            decoder.setValue(hit.getCellID0())
            layer = decoder["layer"].value()
            hit_info[key]["layer"].append(layer)
reader.close()
print("Total number of truth staus per event:", n_truth_staus)
print("Accepted staus per event:", n_accepted_staus)

all_data = {"hit_info": hit_info, "mcp_stau_info": mcp_stau_info, "n_accepted_staus": n_accepted_staus}
with open(out_file, "w") as f:
    json.dump(all_data, f, indent=2, sort_keys=True)