import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from math import *
import pathlib
import pickle
from dataclasses import dataclass

import pyLCIO
from pyLCIO import UTIL, EVENT
import ROOT

## trying to determine when the latest hits are for each sub-detector in order to set "loose" timing windows for full hit acceptance

loose4_dir = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/loose4/nobib/4000_10"
loose4p5_dir = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/loose4p5/4500_10"

reco_hit_collections = {"VB": "VXDBarrelHits",
                   "IB": "ITBarrelHits",
                   "OB": "OTBarrelHits"}
sim_hit_collections = {"VB": "VertexBarrelCollection",
                   "IB": "InnerTrackerBarrelCollection",
                   "OB": "OuterTrackerBarrelCollection"}

stau_ids = {1000015, -1000015, 2000015, -2000015}

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

def get_latest_hit(subdet, sample, dir, option):
    hits = []
    for ifile in range(100):
        file_name = f"{sample}_reco{ifile}.slcio"
        file_path = os.path.join(dir,file_name) 
        reader.open(file_path)
       
        for event in reader:
            if option == "reco":
                hit_collection = reco_hit_collections[subdet] 
            elif option == "sim":
                hit_collection = sim_hit_collections[subdet]
            coll = event.getCollection(hit_collection) 

            encoding = coll.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
            decoder = UTIL.BitField64(encoding)
            
            for hit in coll:
                decoder.setValue(int(hit.getCellID0()))
                if option == "sim":
                    mcp = hit.getMCParticle() 
                    if not mcp or mcp.getPDG() not in stau_ids:
                        continue
                layer = decoder["layer"].value()
                
                if subdet == "VB":
                    lastlayer = 7
                else:
                    lastlayer = 2  
                if layer == lastlayer:
                    hits.append(hit.getTime()) 
    print(np.mean(hits))
    print(np.max(hits))

# get_latest_hit("VB", "4500_10", loose4p5_dir, "sim")        
# get_latest_hit("IB", "4500_10", loose4p5_dir, "sim")        
# get_latest_hit("OB", "4500_10", loose4p5_dir, "sim")       
# print(" ") 

# get_latest_hit("VB", "4000_10", loose4_dir, "sim")        
# get_latest_hit("IB", "4000_10", loose4_dir, "sim")        
get_latest_hit("OB", "4000_10", loose4_dir, "sim")   