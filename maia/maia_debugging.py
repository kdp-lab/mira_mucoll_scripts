import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import pyLCIO
from pyLCIO import UTIL, EVENT
import ROOT

unconed_hit_collections = ["IBTrackerHits", "IETrackerHits", "VBTrackerHits", "VETrackerHits", "OBTrackerHits", "OETrackerHits"]
coned_hit_collections = ["IBTrackerHitsConed", "IETrackerHitsConed", "VBTrackerHitsConed", "VETrackerHitsConed", "OBTrackerHitsConed", "OETrackerHitsConed"]

coned_sim_collections = ["VertexBarrelCollectionConed", "VertexEndcapCollectionConed", "InnerTrackerBarrelCollectionConed", "InnerTrackerEndcapCollectionConed", "OuterTrackerBarrelCollectionConed", "OuterTrackerEndcapCollectionConed"]

mcp_collection = "MCParticle"
stau_ids = {1000015, -1000015, 2000015, -2000015}

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

def analyze_slcio(path, coned=True):
    reader.open(path)
    print(path)
    for event in reader:
        names = event.getCollectionNames()
        mcp_collection = event.getCollection("MCParticle")
        track_collection = event.getCollection("SiTracks") if "SiTracks" in names else None
        track_relation_collection = event.getCollection("MCParticle_SiTracks") if "MCParticle_SiTracks" in names else None
        
        print(len(mcp_collection), " mcps")             

        if track_collection:
            print(len(track_collection), " tracks")
        else:
            print("no tracks") 
        if track_relation_collection:
            print(len(track_relation_collection), " track relations")
        else:
            print("no track relations")
        if coned==True:
            sim_collections = ["VertexBarrelCollectionConed", "VertexEndcapCollectionConed", "InnerTrackerBarrelCollectionConed", "InnerTrackerEndcapCollectionConed", "OuterTrackerBarrelCollectionConed", "OuterTrackerEndcapCollectionConed"]
        else:
            sim_collections = ["VertexBarrelCollection", "VertexEndcapCollection", "InnerTrackerBarrelCollection", "InnerTrackerEndcapCollection", "OuterTrackerBarrelCollection", "OuterTrackerEndcapCollection"]

        for coll_name in sim_collections:
            coll = event.getCollection(coll_name)
            print(coll_name)
            encoding = coll.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
            decoder = UTIL.BitField64(encoding)
            other_hits = 0
            for hit in coll:
                mcp = hit.getMCParticle() 
                if not mcp or mcp.getPDG() not in stau_ids: 
                    other_hits +=1
                    continue
                decoder.setValue(int(hit.getCellID0()))
                system = decoder["system"].value()
                layer = decoder["layer"].value()

                print(mcp.getPDG(), " hit layer ", layer)
            print(f"+ {other_hits} non-stau hits")

        nav = UTIL.LCRelationNavigator(track_relation_collection)
        for i,track in enumerate(track_collection):
            track_mcps = nav.getRelatedFromObjects(track)
            hits = track.getTrackerHits()

            if track_mcps or coned:
                print(f"\ntrack {i}:")
                print(f"Track reduced chi^2: {(track.getChi2() / track.getNdf()):.2f}")
                for mcp in track_mcps:
                    print("PDG: ", mcp.getPDG())
                    print("Generator status: ", mcp.getGeneratorStatus())

                    print(len(hits), " hits in track ")
                for hit in hits:
                    decoder.setValue(int(hit.getCellID0()))
                    system = decoder["system"].value()
                    layer = decoder["layer"].value()
                    x = hit.getPosition()[0]
                    y = hit.getPosition()[1]
                    z = hit.getPosition()[2]
                    print(f"hit system {system}, layer {layer} at {hit.getTime():.2f} ns (pos {x:.2f}, {y:.2f}, {z:.2f})")

                if coned == False:
                    rel_names = ("VBTrackerHitsRelations","IBTrackerHitsRelations","OBTrackerHitsRelations",
                    "VETrackerHitsRelations","IETrackerHitsRelations","OETrackerHitsRelations")
                else:
                    rel_names = ("VBTrackerHitsRelationsConed","IBTrackerHitsRelationsConed","OBTrackerHitsRelationsConed",
                    "VETrackerHitsRelationsConed","IETrackerHitsRelationsConed","OETrackerHitsRelationsConed")

                navs = {r: UTIL.LCRelationNavigator(event.getCollection(r)) for r in rel_names if r in names}
                pdgs = []
                for h in track.getTrackerHits():
                    pdg = next((sh.getMCParticle().getPDG()
                                for nav in navs.values()
                                for sh in (nav.getRelatedToObjects(h) or nav.getRelatedFromObjects(h) or [])
                                if sh.getMCParticle()), None)
                    pdgs.append(pdg)
                print("hit PDGs:", pdgs)



            
    reader.close()

# print(f" coned, no bib")
# analyze_slcio("/ospool/uc-shared/project/futurecolliders/miralittmann/maia/09-04_reco_coned_nobib.slcio", coned=True)
print(f"\n coned, bib")
analyze_slcio("/ospool/uc-shared/project/futurecolliders/miralittmann/maia/09-04_reco_coned_bib.slcio", coned=True)
# print(f"\n unconed, no bib")
# analyze_slcio("/ospool/uc-shared/project/futurecolliders/miralittmann/maia/09-04_reco_unconed_nobib.slcio", coned=False)
# print(f"\n unconed, bib")
# analyze_slcio("/ospool/uc-shared/project/futurecolliders/miralittmann/maia/09-04_reco_unconed_bib.slcio", coned=False)
print(f"\n ====================== fix attempt 10% BIB ================== ") 
analyze_slcio("/ospool/uc-shared/project/futurecolliders/miralittmann/maia/09-05_reco_10pbib.slcio")
