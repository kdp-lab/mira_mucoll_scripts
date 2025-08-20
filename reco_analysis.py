import pyLCIO
from pyLCIO import UTIL
import ROOT
from pathlib import Path
import glob
import json
from math import *
import numpy as np
import argparse
import pickle
import pathlib
from collections import Counter

# Set up the argument parser
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
parser.add_argument(
    "--sample",
    type=str,
    required=True,
    help="sample to analyze"
)

# def event_looper(reader, all_events):
#     while True:
#         try: 
#             evt = reader.readNextEvent()
#         except pyLCIO.LcioException as err: 
#             continue
#         if evt is None:
#             break
#         yield evt
#         if not all_events:
#             break
#     # else:
#     #     for evt in reader:
#     #         yield evt
#     #         break

def event_looper(reader, all_events):
    for evt in reader:
        yield evt
        if not all_events:
            break
# Parse the arguments
args = parser.parse_args()
sample = args.sample
chunk = args.chunk
n_hits_cut = 0
#n_hits_cut = 0

### constants ###
Bfield = 3.57
c = 299792458/1000000  # mm/ns 
stau_ids = [1000015, 2000015]
rand = ROOT.TRandom3(0)    

### input files ###

in_file = f"{sample}_reco{chunk}.slcio"
out_file = f"{sample}_reco{chunk}.json"
print("looking for chunk", chunk)

#### MAKE BIG EMPTY LISTS #### 

# make lists of info we want to eventually save for MATCHED reconstructed stau tracks -- i.e. 
# truth level information that corresponds to a reconstructed track (truth level but only if there is also reco)
# TODO: why are these the only parameters we're saving for the matched staus? 
match_stau_info = {
    "pt": [],
    "eta": [],
    "phi": [], 
    "theta": [],
    "id": [],
    "travel_dist": []
}

systems = ["VB", "VE", "IB", "IE", "OB", "OE"]
hit_fields = ["x", "y", "z", "time", "corrected_time", "layer_hit", "side_hit", "mcp_id", "pdg_id"]
hit_collection_names = ["ITBarrelHits", "ITEndcapHits", "OTBarrelHits", "OTEndcapHits", "VXDBarrelHits", "VXDEndcapHits"]

# make lists of info about stau tracks (exclusively reco type information)
match_track_info = {
    "hits": {system: {f: [] for f in hit_fields} for system in systems}, # fill like match_track_info["hits"][system]["<field_val>""].append(<val>)
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
}

hits_from_mcp = {
    system: {
        "x": [],
        "y": [],
        "z": [],
        "time": [],
        "pdg_id": [],
        "layer": [],
        "mcp_id": [],
        "is_stau_hit": [],
        "vertex_x": [],
        "vertex_y": [],
        "vertex_z": []
    } for system in systems
}

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
print(in_file)
reader.open(in_file)

seen_staus = set()
n_truth_staus = 0
n_reco_staus = 0
accepted_staus = 0
# Looping through events but only if we need to (I've been running with one event per chunk, if this isn't the case add --all-events to the sh script line 7)
for event in event_looper(reader, args.all_events): 
    decoders = {}
    for col_name in ("VBHits", "VEHits","ITBarrelHits", "ITEndcapHits","OTBarrelHits", "OTEndcapHits"):
        try:
            col = event.getCollection(col_name)
            enc = col.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
            decoders[col_name] = pyLCIO.UTIL.BitField64(enc)
        except pyLCIO.EVENT.DataNotAvailableException:
            continue            
    mcp_collection = event.getCollection("MCParticle")
    rel_collection = event.getCollection("MCParticle_SiTracks") # OR REFITTED
    relation = pyLCIO.UTIL.LCRelationNavigator(rel_collection)


    trk_coll = event.getCollection("SiTracks")

    reco_hit_colls = {
        "VXDBarrel" : event.getCollection("VXDBarrelHits"),
        "VXDEndcap" : event.getCollection("VXDEndcapHits"),
        "ITBarrel"  : event.getCollection("ITBarrelHits"),
        "ITEndcap"  : event.getCollection("ITEndcapHits"),
        "OTBarrel"  : event.getCollection("OTBarrelHits"),
        "OTEndcap"  : event.getCollection("OTEndcapHits"),
    }

    sim_hit_colls  = {
        "VXDBarrel" : event.getCollection("VertexBarrelCollection"),
        "VXDEndcap" : event.getCollection("VertexEndcapCollection"),
        "ITBarrel"  : event.getCollection("InnerTrackerBarrelCollection"),
        "ITEndcap"  : event.getCollection("InnerTrackerEndcapCollection"),
        "OTBarrel"  : event.getCollection("OuterTrackerBarrelCollection"),
        "OTEndcap"  : event.getCollection("OuterTrackerEndcapCollection"),
    }
 
    rel_nav = {
        "VXDBarrel" : pyLCIO.UTIL.LCRelationNavigator(
                        event.getCollection("VXDBarrelHitsRelations")),
        "VXDEndcap" : pyLCIO.UTIL.LCRelationNavigator(
                        event.getCollection("VXDEndcapHitsRelations")),
        "ITBarrel"  : pyLCIO.UTIL.LCRelationNavigator(
                        event.getCollection("ITBarrelHitsRelations")),
        "ITEndcap"  : pyLCIO.UTIL.LCRelationNavigator(
                        event.getCollection("ITEndcapHitsRelations")),
        "OTBarrel"  : pyLCIO.UTIL.LCRelationNavigator(
                        event.getCollection("OTBarrelHitsRelations")),
        "OTEndcap"  : pyLCIO.UTIL.LCRelationNavigator(
                        event.getCollection("OTEndcapHitsRelations")),
    }
    system_to_det = {1:'VXDBarrel', 2:'VXDEndcap',
                    3:'ITBarrel',  4:'ITEndcap',
                    5:'OTBarrel',  6:'OTEndcap'}
    fallback_decoder = next(iter(decoders.values()))
    
    stau_track_set = set() 
    for track in trk_coll:
        for parent in relation.getRelatedFromObjects(track):
            if abs(parent.getPDG()) in stau_ids:
                stau_track_set.add(track)
                break

    hits_in_stau_tracks = set()
    stau_reco_hits = []
    for track in stau_track_set:
        stau_reco_hits.extend(track.getTrackerHits())
        for hit in track.getTrackerHits():
            hits_in_stau_tracks.add(hit)

    for reco_hit in stau_reco_hits:
        cellID = int(reco_hit.getCellID0())
        decoder = decoders["ITBarrelHits"]
        decoder.setValue(cellID)
        system = decoder["system"].value()
        layer = decoder["layer"].value()
        det_key = {1:'VB',2:'VE',3:'IB',4:'IE',5:'OB',6:'OE'}[system]
        #if det_key in ("VE", "IE", "OE"):
        #    continue
        
        det_name = {1:'VXDBarrel',2:'VXDEndcap', 3:'ITBarrel', 4:'ITEndcap',5:'OTBarrel', 6:'OTEndcap'}[system]
        for sim in rel_nav[det_name].getRelatedToObjects(reco_hit):
            mcp = sim.getMCParticle()
            if not mcp:
                continue
            is_stau_hit = abs(mcp.getPDG()) in stau_ids
            pdg = mcp.getPDG()   
            #if pdg not in stau_ids:
            #    continue

            print("GENERATOR STATUS:",mcp.getGeneratorStatus())
            hits_from_mcp[det_key]["x"].append(reco_hit.getPosition()[0])
            hits_from_mcp[det_key]["y"].append(reco_hit.getPosition()[1])
            hits_from_mcp[det_key]["z"].append(reco_hit.getPosition()[2])
            hits_from_mcp[det_key]["time"].append(reco_hit.getTime())
            hits_from_mcp[det_key]["pdg_id"].append(mcp.getPDG())
            hits_from_mcp[det_key]["layer"].append(layer)
            hits_from_mcp[det_key]["mcp_id"].append(mcp.id())
            hits_from_mcp[det_key]["vertex_x"].append(mcp.getVertex()[0])
            hits_from_mcp[det_key]["vertex_y"].append(mcp.getVertex()[1])
            hits_from_mcp[det_key]["vertex_z"].append(mcp.getVertex()[2])
            hits_from_mcp[det_key]["is_stau_hit"].append(is_stau_hit)

    for mcp in mcp_collection:
 
        mcp_stau_vertex_r = sqrt(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2)
        travel_dist = sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2 + mcp.getEndpoint()[2]**2) - sqrt(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2 + mcp.getVertex()[2]**2)

        if abs(mcp.getPDG()) not in stau_ids:
            continue 
        if mcp.id() in seen_staus:
            continue 
        if travel_dist == 0:
            continue
        if mcp_stau_vertex_r > 553.0:
            continue

        n_truth_staus +=1 # already did this in sim-only analysis, but maybe can serve as sanity check
        reco_success = False

        stau_tracks = relation.getRelatedToObjects(mcp) # getting the stau's tracks (so this is reco info now) TODO: would it be better if we got this from alltracks??
        # or no because we want to have these tracks be reco stau hits

        # truth level information we want to save about the staus if they can be reconstructed (matched staus) 
        # NOTE: this is also done in the sim-only file but I'm redoing it here because I don't know how to get the relations saved 
        mcp_stau_momentum = mcp.getMomentum() 
        mcp_stau_tlv = ROOT.TLorentzVector()
        mcp_stau_tlv.SetPxPyPzE(mcp_stau_momentum[0], mcp_stau_momentum[1], mcp_stau_momentum[2], mcp.getEnergy()) 

        mcp_pt = mcp_stau_tlv.Perp()
        mcp_eta = mcp_stau_tlv.Eta()
        mcp_phi = mcp_stau_tlv.Phi()
        mcp_theta = mcp_stau_tlv.Theta()
        mcp_stau_endpoint_r = sqrt(mcp.getEndpoint()[0] ** 2 + mcp.getEndpoint()[1] ** 2)
        
        print(len(stau_tracks))
        if mcp_stau_endpoint_r < 102.0 or abs(mcp_stau_tlv.Eta()) > 0.8:
            continue
        if len(stau_tracks)>1:
            print("HEY more than one track per stau here")
            continue

        accepted_staus += 1
        for track in stau_tracks:
            track_hits = {s:{f:[] for f in hit_fields} for s in systems}
            seen_layers = set() # making an empty set that we can add the layers to as we go through them so that we aren't double counting hits
            n_pix_hits = 0
            n_inner_hits = 0
            n_outer_hits = 0 
            
            # setting up tracker mapping: goes like pixel barrel, pixel endcap, inner barrel, inner endcap, etc
            system_map = {1: "VB", 2: "VE", 3: "IB", 4: "IE", 5: "OB", 6: "OE"}

            for hit in track.getTrackerHits():
                encoding = event.getCollection("ITBarrelHits").getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                ## TODO: assuming that every collection has the same decoder parameters?
                decoder = pyLCIO.UTIL.BitField64(encoding)
                cellID = int(hit.getCellID0())
                decoder.setValue(cellID) 
                
                system = decoder["system"].value()
                layer = decoder["layer"].value()
                side = decoder["side"].value() # NOTE: are we actually using this info?
                
                if (system,layer) in seen_layers: 
                    continue
                seen_layers.add((system,layer))

                if system_map[system] in ("VB","VE"): 
                    n_pix_hits += 1
                elif system_map[system] in ("IB", "IE"): 
                    n_inner_hits +=1
                else: 
                    n_outer_hits += 1 

                hit_x_pos = hit.getPosition()[0] 
                hit_y_pos = hit.getPosition()[1] 
                hit_z_pos = hit.getPosition()[2] 

                distance = sqrt(hit_x_pos**2 + hit_y_pos**2 + hit_z_pos**2)
                flight_time = distance/c
                
                if system > 2:
                    resolution = 0.06
                else: 
                    resolution = 0.03
                    
                track_corrected_time = hit.getTime() # based on 1st principles calculations this is what we want
                # now for reco'd track information for these matched things
                det_key = system_map[system] 

                track_hits[det_key]["x"].append(hit_x_pos)
                track_hits[det_key]["y"].append(hit_y_pos)
                track_hits[det_key]["z"].append(hit_z_pos)
                track_hits[det_key]["time"].append(hit.getTime())
                track_hits[det_key]["corrected_time"].append(track_corrected_time)
                track_hits[det_key]["layer_hit"].append(layer)
                track_hits[det_key]["side_hit"].append(side)
                track_hits[det_key]["mcp_id"].append(mcp.id()) 
                track_hits[det_key]["pdg_id"].append(mcp.getPDG())
      
            # going to ignore tracks that don't have 3.5+ hits since they aren't reliable for reco              
            n_total_hits = (n_pix_hits)/2.0 + n_inner_hits + n_outer_hits 
            # NOTE: dividing pix hits by two because current geometry uses doublet layers, this will not be true soon
            seen_staus.add(mcp.id())
            #if n_total_hits < n_hits_cut:
            #    print("Not enough hits to reconstruct stau.")
            #    continue  
            reco_success = True 
            
            for s in systems:
                for f in hit_fields:
                    match_track_info["hits"][s][f].append(track_hits[s][f])
            # reference like: match_track_info["hits"][<system>][<field>][<track index>]

            track_tlv = ROOT.TLorentzVector()
            theta = np.pi/2 - np.arctan(track.getTanLambda())
            phi = track.getPhi()
            eta = -np.log(np.tan(theta/2))
            p_t = 0.3 * Bfield /fabs(track.getOmega() * 1000.)
            p_t_res = (mcp_stau_tlv.Perp() - p_t) / mcp_stau_tlv.Perp()
            track_tlv.SetPtEtaPhiE(p_t, eta, phi, 0)
            n_track_hits = track.getTrackerHits().size() # how is this different than the other n hits, what is size
                 
            match_track_info["pt"].append(p_t)
            match_track_info["theta"].append(theta)
            match_track_info["phi"].append(phi)
            match_track_info["ndf"].append(track.getNdf()) # NOTE: what is this
            match_track_info["chi_sq"].append(track.getChi2())
            match_track_info["d0"].append(track.getD0())
            match_track_info["z0"].append(track.getZ0())
            match_track_info["pt_res"].append(p_t_res) 

            match_track_info["n_hits"].append(n_track_hits)
            match_track_info["n_hits_vertex"].append(n_pix_hits)
            match_track_info["n_hits_inner"].append(n_inner_hits)
            match_track_info["n_hits_outer"].append(n_outer_hits)
 
            # add to match_stau_info lists -- this is the same sim info as before but now that we have confirmed the track for it we can add it to these lists as well
            # TODO: I think? but then why don't we save all the other information too, or mark it somehow with like "matched stau"
            # i guess because the reco only gives us this information so thisis the part we can compare with sim? and keep it separate from the others?
        if reco_success:
            n_reco_staus += 1
            match_stau_info["pt"].append(mcp_pt)
            match_stau_info["eta"].append(mcp_eta)
            match_stau_info["phi"].append(mcp_phi)
            match_stau_info["theta"].append(mcp_theta)
            match_stau_info["id"].append(mcp.getPDG())
            match_stau_info["travel_dist"].append(travel_dist)

reader.close()

bad_files = []
print("\nSummary Statistics:")
print("Found:")
print(f"Total truth-level staus: {n_truth_staus}")
print(f"Reconstructable staus: {n_reco_staus}")
print(f"Non-reconstructable staus: {n_truth_staus - n_reco_staus}")
if n_truth_staus == 0:
    print(chunk, "has no events?")
    bad_files.append(chunk)


all_data = {"match_stau_info": match_stau_info, "match_track_info": match_track_info, "bad_files": bad_files, "hits_from_mcp": hits_from_mcp, "accepted_staus": accepted_staus}
with open(out_file, "w") as f:
    json.dump(all_data, f, indent=2, sort_keys=True)