import pyLCIO
from pyLCIO import UTIL
import ROOT
from pathlib import Path
import glob
import json
from math import *
import numpy as np
import argparse

# Trying to rewrite Tate's analysis code for slcio files

# argument parser stuff (i think this is just for the sh script, finding chunks etc)
ROOT.gROOT.SetBatch()


# Set up the argument parser
parser = argparse.ArgumentParser(description="Analyze SLCIO chunks.")
parser.add_argument(
    '--chunk', 
    type=int, 
    required=True, 
    help="The chunk number to analyze"
)

parser.add_argument(
    '--eventsPerChunk', 
    type=int, 
    required=False, 
    default=10,
    help="The chunk number to analyze"
)

# Parse the arguments
args = parser.parse_args()

### constants ###
Bfield = 3.57
c = 299792458/1000000  # mm/ns 
stau_ids = [1000015, 2000015]

### input files ##
in_path = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco_try1/"
save_path = "/scratch/miralittmann/analysis/10pbibjson/4000_10/medium_window/sim/"
file_name = "4000_10" 
chunk = args.chunk
in_file = f"{in_path}{file_name}_reco{chunk}.slcio"
out_file = f"{save_path}{file_name}_reco{chunk}.json"
print("looking for chunk", chunk)
file_name = glob.glob(in_file)
print("found %i files:"%len(file_name), in_file)

#### MAKE BIG EMPTY LISTS #### fill later like stau["val"].append(val)
# make lists of info we want to eventually save for MCParticles (truth level) 
mcp_stau_info = {
    "beta": [],
    "m": [], 
    "p_tot": [],
    "p_x": [],
    "p_y": [],
    "p_z": [],
    "energy": [],
    "pt": [],
    "eta": [],
    "phi": [],
    "theta": [],
    "travel_distance": []
}

# make lists of info we want to eventually save for MATCHED reconstructed stau tracks -- i.e. 
# truth level information that corresponds to a reconstructed track (truth level but only if there is also reco)
match_stau_info = {
    "pt": [],
    "eta": [],
    "phi": [], 
    "theta": []
}

# make lists of info about stau tracks (exclusively reco type information)
match_track_info = {
    "pt": [],
    "theta": [],
    "phi": [],
    "ndf": [],
    "chi_sq": [],
    "d0": [],
    "z0": [],
    "n_hits": [],
    "pt_res": [],
    "n_hits_vertex": [],
    "n_hits_inner": [],
    "n_hits_outer": [],
    "x": [],
    "y": [],
    "z": [],
    "time": [],
    "corrected_time": [],
    "detector_hit": [],
    "layer_hit": [],
    "side_hit": []
}

# make lists of info for hits (reco)
# these are from sim information, accounts of how they hit the detectors -- truth level. and the reco ones are predictions by the sim files?? 
# NOTE check this because I'm confused
layers = ["VB", "VE", "IB", "IE", "OB", "OE"]
sim_fields = ["x", "y", "z", "time", "pdg", "mcpid", "layer"]
reco_fields = ["x", "y", "z", "time", "layer"]

sim_hit_info = {layer: {f: [] for f in sim_fields} for layer in layers}
reco_hit_info = {layer: {f: [] for f in reco_fields}  for layer in layers}
# then later fill like: sim_hit_info["VB"]["x"].append(sim_x_value)

# make counter variables: number of: mcp staus, inner hits, i, matched tracks, dupes (?), fake tracks, truth staus
n_reco_stau = 0
n_inner_hits = 0
n_event = 0
n_match_tracks = 0
n_fake_tracks = 0 # TODO: add fake tracks tracking
n_truth_stau = 0

# tate loops over events, but at least the way i've been running samples is only with one event per file. so not doing that right now
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(file_name)
event = reader.readNextEvent()

### Taking the collections we want from the event ###
# core physics level collections
mcp_collection, track_collection, rel_collection = (event.getCollection(n) for n in ("MCParticle", "SiTracks_Refitted", "MCParticle_SiTracks_Refitted"))
relation = pyLCIO.UTIL.LCRelationNavigator(rel_collection)

# hit collections
layer_names_for_hits = ("ITBarrelHits", "ITEndcapHits", "OTBarrelHits", "OTEndcapHits", "VXDBarrelHits", "VXDEndcapHits")
hit_collections = {name: event.getCollection(name) for name in layer_names_for_hits}

# relations
hit_relations = {name: UTIL.LCRelationNavigator(event.getCollection(f"{name}Relations")) for name in layer_names_for_hits}
 
# use like: n_itb_hits  = hit_collections["ITBarrelHits"].getNumberOfElements()

mapping = [
    ("VertexBarrelCollection", hit_relations["VXDBarrelHits"], sim_hit_info["VB"], reco_hit_info["VB"]),
    ("VertexEndcapCollection", hit_relations["VXDEndcapHits"], sim_hit_info["VE"], reco_hit_info["VE"]),
    ("InnerTrackerBarrelCollection", hit_relations["ITBarrelHits"], sim_hit_info["IB"], reco_hit_info["IB"]),
    ("InnerTrackerEndcapCollection", hit_relations["ITEndcapHits"], sim_hit_info["IE"], reco_hit_info["IE"]),
    ("OuterTrackerBarrelCollection", hit_relations["OTBarrelHits"], sim_hit_info["OB"], reco_hit_info["OB"]), 
    ("OuterTrackerEndcapCollection", hit_relations["OTEndcapHits"], sim_hit_info["OE"], reco_hit_info["OE"])
]

# Big loop 1: here we are saving the hit information for each detector layer 1 by one from sim to the lists we made above -- sim and what sim predicts reco is?
# probably a way to combine this loop with the second loop below? or this is what can be done in a separate file
for coll_name, relation, sim_info, reco_info in mapping:
    try:
        collection = event.getCollection(coll_name) # we are already doing this above to make hit_collections list, do we need to do it again? or not do it there?
        for hit in collection: 
            sim_info["x"].append(hit.getPosition()[0]) 
            sim_info["y"].append(hit.getPosition()[1]) 
            sim_info["z"].append(hit.getPosition()[2]) 
            sim_info["time"].append(hit.getTime()) 

            # Get the info about the MCParticle -- what kind of particle it is, and add to the sim info list
            mcp = hit.getMCParticle()  
            sim_info["pdg"].append(mcp.getPDG()) if mcp else None
            sim_info["mcpid"].append(mcp.id()) if mcp else None

            # Get the layer info -- not exactly sure how this works
            encoding = collection.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
            decoder = pyLCIO.UTIL.BitField64(encoding)
            cellID = int(hit.getCellID0())
            sim_info["layer"].append(decoder["system"].value())

            # then reco info. maybe this is what the problem was with tate's script, should i move it out of here?
            reco_hit = relation.getRelatedFromObjects(hit)
            if len(reco_hit) > 0:
                reco_hit = reco_hit[0]
                rel_position = reco_hit.getPosition()
                reco_info["x"].append(rel_position[0])
                reco_info["y"].append(rel_position[1])
                reco_info["z"].append(rel_position[2])
                reco_info["time"].append(reco_hit.getTime())
                reco_info["layer"].append(decoder["system"].value())

    except Exception as e:
        print("Couldn't access {coll_name}: {e}")

# Loop 2: looping through particles in the relation lists, checking if they are staus, saving hit info from them (reco info)
for mcp in mcp_collection:
    # selecting for only mcps that are staus and not intermediate staus (created by another stau)
    if abs(mcp.getPDG()) not in stau_ids:
        continue
    if any(abs(parent.getPDG()) in stau_ids for parent in mcp.getParents()):
        continue
    n_truth_stau += 1 # counting number of true staus even if we can't reconstruct them ()
    stau_tracks = relation.getRelatedToObjects(mcp) # getting the stau's tracks (so this is reco info)

    for track in stau_tracks:
        # here again we are making lists and counter variables again specifically for the tracks, i think i don't need to do this?
        n_layers_crossed = 0.0
        seen_layers = set() # making an empty set that we can add the layers to as we go through them so that we aren't double counting hits
        n_pix_hits = 0
        n_inner_hits = 0
        n_outer_hits = 0 
        # TODO: are we overwriting these lists? we have them earlier. also made them at the top in the big st_... lists, so can just add to those?
        
        # setting up tracker mapping: goes like pixel barrel, pixel endcap, inner barrel, inner endcap, etc
        layer_map = {1: "pix", 2: "pix", 3: "inner", 4: "inner", 5: "outer", 6: "outer"}
        encoding = hit_collections["ITBarrelHits"].getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)  
        decoder = pyLCIO.UTIL.BitField64(encoding)
        # we're indexing only the inner tracker collection to get the shape for the encoder, then use the same one for all layers 

        # in tate's script, the following (for mcp_stau) is inside the "for hit" loop. but i don't think that's necessary because the mcp lists are just truth-level and 
        # don't depend on the reco'd hit information? but TODO: check this.
        mcp_stau_momentum = mcp.getMomentum()
        mcp_stau_mass = mcp.getMass()
        mcp_stau_tlv = ROOT.TLorentzVector()
        mcp_stau_tlv.setPxPyPzE(mcp_stau_momentum[0], mcp_stau_momentum[1], mcp_stau_momentum[2], mcp.getEnergy()) 
        mcp_stau_beta = mcp_stau_tlv.Beta()
        
        mcp_stau_info["p_tot"].append(mcp_stau_tlv.Mag())
        mcp_stau_info["p_x"].append(mcp_stau_momentum[0])
        mcp_stau_info["p_y"].append(mcp_stau_momentum[1])
        mcp_stau_info["p_z"].append(mcp_stau_momentum[2])
        mcp_stau_info["energy"].append(mcp.getEnergy())
        mcp_stau_info["m"].append(mcp_stau_mass)
        mcp_stau_info["beta"].append(mcp_stau_beta)
        mcp_stau_info["pt"].append(mcp_stau_tlv.Perp())
        mcp_stau_info["eta"].append(mcp_stau_tlv.Eta())
        mcp_stau_info["phi"].append(mcp_stau_tlv.Phi())
        mcp_stau_info["theta"].append(mcp_stau_tlv.Theta())
        mcp_stau_info["travel_distance"].append(sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2 + mcp.getEndpoint()[2]**2) - sqrt(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2 + mcp.getVertex()[2]**2))

        for hit in track.getTrackerHits():
            decoder.setValue(int(hit.getCellID0()))
            system = decoder["system"].value()
            layer = decoder["layer"].value()
            side = decoder["side"].value() # i don't know what this means/if we're actually using it
            
            if layer in seen_layers: # TODO: I think this will work the same as before??? Check
                continue
            
            n_layers_crossed += 0.5 if system in (1, 2) else 1.0
            seen_layers.add(layer)
            if layer_map["system"] == "pix": 
                n_pix_hits += 1
            elif layer_map["system"] == "inner": 
                n_inner_hits +=1
            else: 
                n_outer_hits += 1
            # TODO: do we want to save info about if it's barrel/endcap? could be helpful?

            track_x_pos = hit.getPosition()[0]
            track_y_pos = hit.getPosition()[1]
            track_z_pos = hit.getPosition()[2]
            
            # going to ignore tracks that don't have 3.5+ hits since they aren't reliable for reco
            n_total_hits = (n_pix_hits)/2.0 + n_inner_hits + n_outer_hits
            if n_total_hits < 3.5:
                continue 
            n_reco_stau += 1
            #### so, only doing the following if there are 3.5 hits or more, meaning the track and mcp can be matched. now looking more at reco info for that 

            distance = sqrt(track_x_pos**2 + track_y_pos**2 + track_z_pos**2)
            flight_time = distance/c

            if system > 2:
                resolution = 0.06
            else: 
                resolution = 0.03
                
            track_corrected_time = hit.getTime()*(1.+ROOT.TRandom3(0).Gaus(0., resolution)) - flight_time # this is complicated messed up timing i am not going to mess with rn
            # CHECK: set TRandom(3) to 0 because it was previously set to ievt, but we again only have one event per file right now. change? because i think this means new time every run
            # now for reco'd track information for these matched things
                               
            track_tlv = ROOT.TLorentzVector()
            theta = np.pi/2 - np.arctan(track.getTanLambda())
            phi = track.getPhi()
            eta = -np.log(np.tan(theta/2))
            p_t = 0.3 * Bfield /fabs(track.getOmega() * 1000.)
            p_t_res = (mcp_stau_tlv.Perp() - p_t) / mcp_stau_tlv.Perp()
            track_tlv.SetPtEtaPhiE(p_t, eta, phi, 0)
            n_track_hits = track.getTrackerHits().size() # how is this different than the other n hits, what is size
            
            # add to match_track_info lists. TODO: Tate has these in the form of lists (like [p_t] instead of just pt), do i need to do this
            match_track_info["pt"].append(p_t)
            match_track_info["theta"].append(theta)
            match_track_info["phi"].append(phi)
            match_track_info["ndf"].append(track.getNdf()) # what is this
            match_track_info["chi_sq"].append(track.getChi2())
            match_track_info["d0"].append(track.getD0())
            match_track_info["z0"].append(track.getZ0())
            match_track_info["pt_res"].append(p_t_res) 

            match_track_info["n_hits"].append(n_track_hits)
            match_track_info["n_hits_vertex"].append(n_pix_hits)
            match_track_info["n_hits_inner"].append(n_inner_hits)
            match_track_info["n_hits_outer"].append(n_outer_hits)

            match_track_info["x"].append(track_x_pos)
            match_track_info["y"].append(track_y_pos)
            match_track_info["z"].append(track_z_pos)
            match_track_info["time"].append(hit.getTime())
            match_track_info["corrected_time"].append(track_corrected_time)
            match_track_info["detector_hit"].append(system)
            match_track_info["layer_hit"].append(layer)
            match_track_info["side_hit"].append(side) 

            # add to match_stau_info lists -- this is the same sim info as before but now that we have confirmed the track for it we can add it to these lists as well
            # TODO: I think? but then why don't we save all the other information too, or mark it somehow with like "matched stau"
            # i guess because the reco only gives us this information so thisis the part we can compare with sim? and keep it separate from the others?
            match_stau_info["pt"].append(mcp_stau_info["pt"])
            match_stau_info["eta"].append(mcp_stau_info["eta"])
            match_stau_info["phi"].append(mcp_stau_info["phi"])
            match_stau_info["theta"].append(mcp_stau_info["theta"])

reader.close()

# phew.
# NOTE: n_reco_stau is number of staus that have been reconstructed?? or number of reconstructed staus that correspond to a truth one? 

print("\nSummary Statistics:")
print("Found:")
print(f"Total truth-level staus: {n_truth_stau}")
print(f"Reconstructable staus: {n_reco_stau}")
print(f"Non-reconstructable staus: {n_truth_stau - n_reco_stau}")

# TODO: save everything to a json file