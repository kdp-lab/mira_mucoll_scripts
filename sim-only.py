## (attempting to) make a script that only deals with the sim files, gives you an overview of the information in them directly before digi/reco

import pyLCIO
from pyLCIO import UTIL
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

in_path = "/ospool/uc-shared/project/futurecolliders/miralittmann/sim_try1/"
out_path = "/scratch/miralittmann/analysis/10pbibjson/4000_10/medium_window/sim/"
sample = "4000_10"
chunk = args.chunk
in_file = f"{in_path}{sample}_sim{chunk}.slcio"
out_file = f"{out_path}{sample}_sim{chunk}.json"

## constants ##
stau_ids = [1000015, 2000015]

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(in_file)
def event_looper(reader, all_events):
    if all_events:
        for evt in reader:
            yield evt
    else:
        evt = reader.readNextEvent()
        if evt:                        
            yield evt


## The collections we have: ECalBarrelCollection, ECalEndcapCollection, HCalBarrelCollection, HCalEndcapCollection, HCalRingCollection, 
## InnerTrackerBarrelCollection, InnerTrackerEndcapCollection, MCParticle, OuterTrackerBarrelCollection, OuterTrackerEndcapCollection, 
## VertexBarrelCollection, VertexEndcapCollection, YokeBarrelCollection, YokeEndcapCollection

# I'm pretty sure we only care about: MCParticle, VertexBarrelCollection, VertexEndcapCollection, InnerTrackerBarrelCollection, InnerTrackerEndcapCollection, 
# OuterTrackerBarrelCollection, OuterTrackerEndcapCollection
# These have getTypeName() "SimTrackerHit"

### empty lists for things we care about ###

hit_info = {
    "pdg_id": [],
    "x": [],
    "y": [],
    "z": [],
    "layer": []
}

systems = [
    ("VertexBarrelCollection", "VB"),
    ("VertexEndcapCollection", "VE"),
    ("InnerTrackerBarrelCollection", "IB"),
    ("InnerTrackerEndcapCollection", "IE"),
    ("OuterTrackerBarrelCollection", "OB"),
    ("OuterTrackerEndcapCollection", "OE")
]

# making empty lists for the things we want to save from each system's collection, add to like hit_info["VB"]["x"].append( ... ) etc
systems_key  = ["VB", "VE", "IB", "IE", "OB", "OE"]
fields = ["pdg", "x", "y", "z", "time", "layer", "stau"]

hit_info = {key: {f: [] for f in fields} for key in systems_key}

for event in event_looper(reader, args.all_events): 
    # total stau count
    n_truth_staus = 0
    mcps = event.getCollection("MCParticle")
    for mcp in mcps:
        # skipping staus that have stau parents (intermediate staus)
        if any(abs(parent.getPDG()) in stau_ids for parent in mcp.getParents()):
            continue
        if abs(mcp.getPDG()) in stau_ids:
            n_truth_staus +=1

    for system, key in systems:
        # looking at all of the hits in a particular system 
        hits = event.getCollection(system)
        encoding = hits.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
        decoder = pyLCIO.UTIL.BitField64(encoding)
        for hit in hits:
            # saving the info of the particle that made each hit, if it's a stau or not
            mcp = hit.getMCParticle()
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

with open(out_file, "w") as f:
    json.dump(hit_info, f, indent=2)