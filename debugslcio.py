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

reco_file = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/efficiency/nobib/medium/4000_10/4000_10_reco100.slcio"

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reader.open(reco_file)

stau_ids = [1000015, 2000015]
seen_staus = set()
decoders = {}

for ievt, event in enumerate(reader):
    trk_coll = event.getCollection("SiTracks_Refitted")

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

    hits_in_tracks = set()
    for track in trk_coll:
        for hit in track.getTrackerHits():
            hits_in_tracks.add(hit)

    system_map = {1: "VB", 2: "VE", 3: "IB", 4: "IE", 5: "OB", 6: "OE"}

    for det_name, simhits in sim_hit_colls.items(): 
        encoding = simhits.getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
        decoder = pyLCIO.UTIL.BitField64(encoding)
        rel = rel_nav[det_name]

        for simhit in simhits:
            mcp = simhit.getMCParticle()
            pdgid = mcp.getPDG()
            mcpid = mcp.id()
                       
            decoder.setValue(int(simhit.getCellID0()))
            system = system_map[decoder["system"].value()] 
            layer = decoder["layer"].value()
            side = decoder["side"].value()

            reco_hits = rel.getRelatedFromObjects(simhit) 
            for reco_hit in reco_hits:
                in_track = (reco_hit in hits_in_tracks)
 