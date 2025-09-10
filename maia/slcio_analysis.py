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

hit_collections = ["IBTrackerHitsConed", "IETrackerHitsConed", "VBTrackerHitsConed", "VETrackerHitsConed", "OBTrackerHitsConed", "OETrackerHitsConed"]
sim_collections = ["VertexBarrelCollectionConed", "VertexEndcapCollectionConed", "InnerTrackerBarrelCollectionConed", "InnerTrackerEndcapCollectionConed", "OuterTrackerBarrelCollectionConed", "OuterTrackerEndcapCollectionConed"]
rel_collections = ("VBTrackerHitsRelationsConed","IBTrackerHitsRelationsConed","OBTrackerHitsRelationsConed", "VETrackerHitsRelationsConed","IETrackerHitsRelationsConed","OETrackerHitsRelationsConed")
sim_to_name = {
    "VertexBarrelCollectionConed": "VB",
    "VertexEndcapCollectionConed": "VE",
    "InnerTrackerBarrelCollectionConed": "IB",
    "InnerTrackerEndcapCollectionConed": "IE",
    "OuterTrackerBarrelCollectionConed": "OB", 
    "OuterTrackerEndcapCollectionConed": "OE",
}

stau_ids = {1000015, -1000015, 2000015, -2000015} 
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
reco_slcio_path = "/ospool/uc-shared/project/futurecolliders/miralittmann/maia/reco/10pbib"
n_files = 50

Bfield = 5.0
chi2_cut = 3
nhits_cut = 4

CACHE = pathlib.Path("cache/slcio_analysis.pkl")

@dataclass
class AnalyzeResults:
    total_staus: int
    total_accepted_staus: int
    total_tracks: int
    total_good_tracks: int
    total_contaminated_tracks: int
    total_track_pT: list
    total_pT_resolution: list
    yes_ob_pT_res: list
    yes_ob_pTs: list
    no_ob_pT_res: list
    no_ob_pTs: list

# returns:  total_staus, total_accepted_staus, total_tracks, total_good_tracks, total_contaminated_tracks, total_track_pT, total_pT_resolution
def analyze_single_slcio(path, event_print=False, track_print=False, redo=False):

    if CACHE.exists() and not redo:
        print(f"Loading previous arrays from {CACHE}, not re-calculating entire analysis")
        with CACHE.open("rb") as f:
            obj = pickle.load(f)
        if isinstance(obj, AnalyzeResults):
            return obj
        else:
            return AnalyzeResults(*obj)

    print(f"Looking at files in {path}")
    total_staus = 0
    total_tracks = 0
    total_good_tracks = 0 
    total_accepted_staus = 0
    total_contaminated_tracks = 0
    total_tracks_morethanone_bib = 0
    total_track_pT = []
    total_pT_resolution = []
    yes_ob_pTs = []
    no_ob_pTs = []
    yes_ob_pT_res = []
    no_ob_pT_res = []

    for ifile in range(n_files):
        file_name = f"4000_10_reco{ifile}.slcio"
        file_path = os.path.join(path,file_name)
        reader.open(file_path) 
        event_tracks = 0
        event_good_tracks = 0
        event_bib_hits = 0
        event_hits_by_detector = {
            "VB": {
                "stau_hits": 0,
                "bib_hits": 0,
            },
            "VE": {
                "stau_hits": 0,
                "bib_hits": 0,
            },
            "IB": {
                "stau_hits": 0,
                "bib_hits": 0,
            },
            "IE": {
                "stau_hits": 0,
                "bib_hits": 0,
            },
            "OB": {
                "stau_hits": 0,
                "bib_hits": 0,
            },
            "OE": {
                "stau_hits": 0,
                "bib_hits": 0,
            },
        }
        event_track_info = {}


        for event in reader:
            all_collections = event.getCollectionNames()
            mcp_collection = event.getCollection("MCParticle")
            track_collection = event.getCollection("SiTracks") if "SiTracks" in all_collections else None
            track_relation_collection = event.getCollection("MCParticle_SiTracks") if "MCParticle_SiTracks" in all_collections else None

            if track_collection:
                event_tracks += len(track_collection)

            # mcp by mcp acceptance
            seen_staus = ()
            for mcp in mcp_collection:
                mcp_stau_vertex_r = np.sqrt(mcp.getVertex()[0]**2 + mcp.getVertex()[1]**2)
                travel_dist = np.sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2 + mcp.getEndpoint()[2]**2)
                
                if abs(mcp.getPDG()) not in stau_ids or mcp.id() in seen_staus or travel_dist==0 or mcp_stau_vertex_r>553.0:
                    continue
                
                total_staus += 1
                mcp_stau_endpoint_r = np.sqrt(mcp.getEndpoint()[0]**2 + mcp.getEndpoint()[1]**2)
                mcp_stau_momentum = mcp.getMomentum() 
                mcp_stau_tlv = ROOT.TLorentzVector()
                mcp_stau_tlv.SetPxPyPzE(mcp_stau_momentum[0], mcp_stau_momentum[1], mcp_stau_momentum[2], mcp.getEnergy()) 
                
                if mcp_stau_endpoint_r < 102.0 or abs(mcp_stau_tlv.Eta()) > 0.8:
                    continue
                total_accepted_staus += 1

            # detector by detector hit counts            
            for collection in sim_collections:
                coll = event.getCollection(collection)
                det_name = sim_to_name[collection]
                encoding = coll.getParameters().getStringVal(EVENT.LCIO.CellIDEncoding)
                decoder = UTIL.BitField64(encoding)
                for hit in coll:
                    mcp = hit.getMCParticle()
                    if not mcp or mcp.getPDG() not in stau_ids:
                        event_bib_hits += 1
                        event_hits_by_detector[det_name[:2]]["bib_hits"] += 1 
                        continue
                    event_hits_by_detector[det_name[:2]]["stau_hits"] += 1 

            # track by track hit information
            nav = UTIL.LCRelationNavigator(track_relation_collection)
            for itrack, track in enumerate(track_collection):
                track_mcps = nav.getRelatedFromObjects(track)
                track_hits = track.getTrackerHits()
                total_tracks += 1

                pT = 0.3 * Bfield / fabs(track.getOmega() * 1000.) 

                truth_stau = None
                if track_mcps:
                    for mcp in track_mcps:
                        if abs(mcp.getPDG()) in stau_ids:
                            truth_stau = mcp
                            break

                pT_truth = None
                if truth_stau is not None:
                    mom = truth_stau.getMomentum()
                    tlv = ROOT.TLorentzVector()
                    tlv.SetPxPyPzE(mom[0], mom[1], mom[2], truth_stau.getEnergy())
                    pT_truth = tlv.Perp()

                if track_mcps:
                    track_chi2 = track.getChi2() / track.getNdf()
                    track_hit_count = len(track_hits)
                    vb_hits = 0
                    ib_hits = 0
                    ob_hits = 0

                    if track_chi2 < chi2_cut and track_hit_count > nhits_cut:
                        good=True
                        total_good_tracks += 1
                        event_good_tracks += 1
                    else:
                        good=False
                    
                    for mcp in track_mcps:
                        track_pdg = mcp.getPDG()

                    
                    hit_in_ob = False
                    hit_positions = []
                    for hit in track_hits:
                        decoder.setValue(int(hit.getCellID0()))
                        system = decoder["system"].value()
                        if system==1:
                            vb_hits += 1
                        elif system==3:
                            ib_hits += 1
                        elif system==5:
                            ob_hits += 1
                            hit_in_ob = True

                        layer = decoder["layer"].value()                            
                        x = hit.getPosition()[0]
                        y = hit.getPosition()[1]
                        z = hit.getPosition()[2]

                        hit_positions.append((x,y,z, system,layer))
 
                    navs = {r: UTIL.LCRelationNavigator(event.getCollection(r)) for r in rel_collections}
                    pdgs = []
                    for h in track.getTrackerHits():
                        pdg = next((sh.getMCParticle().getPDG()
                                    for nav in navs.values()
                                    for sh in (nav.getRelatedToObjects(h) or nav.getRelatedFromObjects(h) or [])
                                    if sh.getMCParticle()), None)
                        pdgs.append(pdg)

                    num_bib_hits = pdgs.count(None)
                    if num_bib_hits > 0:
                        total_contaminated_tracks += 1
                    if num_bib_hits > 1:
                        total_tracks_morethanone_bib += 1

                    event_track_info[itrack] = {
                        "pdg": track_pdg,
                        "pos": hit_positions,
                        "chi2": track_chi2,
                        "hits_by_det": (vb_hits, ib_hits, ob_hits),
                        "good": good,
                        "pT": pT,
                        "hit_in_ob": hit_in_ob
                    }

                    total_track_pT.append(pT)
                    total_pT_resolution.append((pT - pT_truth)/pT_truth)
                    if hit_in_ob:
                        yes_ob_pT_res.append((pT - pT_truth)/pT_truth)
                        yes_ob_pTs.append(pT)
                    else:
                        no_ob_pT_res.append((pT - pT_truth)/pT_truth)
                        no_ob_pTs.append(pT)
                    
                    if track_print == True:
                        print(f"=============== Track {itrack} ==============")
                        print(f"PDG: {track_pdg}")
                        for ihit,(x, y, z, system, layer) in enumerate(hit_positions):
                            print(f"{pdgs[ihit]} Hit: System {system}, Layer {layer} at ({x:.2f}, {y:.2f}, {z:.2f})")
                        print(f"Reduced Chi^2: {track_chi2:.2f}")
                        print(f"{vb_hits} VB hits, {ib_hits} IB hits, {ob_hits} OB hits")
                        print(f"Good track? {good}")                
                        print(f"Track pT value: {event_track_info[itrack]['pT']:.2f} GeV")
                        if pT_truth is not None:
                            print(f"pT resolution: {((pT - pT_truth)/pT_truth)*100:.2f}%")

            if event_print==True:
                print(f"============================= Event {ifile} Summary =============================")
                print(f"{len(track_collection)} tracks")
                print(f"{len(mcp_collection)} MCParticles")
                print(f"{len(track_relation_collection)} track relations")
                print(f"{event_bib_hits} total BIB hits")
    
    print(f"\n ====================================== Complete summary from {n_files} events ====================================")
    print(f"{total_staus} staus")
    print(f"{total_accepted_staus} accepted staus")
    print(f"{total_tracks} tracks")
    print(f"{total_good_tracks} good stau tracks")
    print(f"{(total_good_tracks / total_accepted_staus)*100:.2f}% efficiency")
    print(f"{total_contaminated_tracks} stau tracks contaminated by BIB hits where {total_tracks_morethanone_bib} have more than one BIB hit")

    res = AnalyzeResults(
        total_staus,
        total_accepted_staus,
        total_tracks,
        total_good_tracks,
        total_contaminated_tracks,
        total_track_pT,
        total_pT_resolution,
        yes_ob_pT_res,
        yes_ob_pTs,
        no_ob_pT_res,
        no_ob_pTs,
    )
    print(f"Writing cache --> {CACHE}")
    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump(res, f, protocol=pickle.HIGHEST_PROTOCOL)
    return res

plot_path = "/scratch/miralittmann/analysis/mira_analysis_code/maia/pT_hists.pdf"
# making histograms of total_track_pT and total_pT_resolution
def make_pt_plots(res):
    with PdfPages(plot_path) as pdf:
        fig, ax = plt.subplots(figsize=(7,6))
        n, bins, _ = ax.hist(res.total_track_pT, bins=50, histtype="stepfilled", alpha=0.25)
        centers = 0.5*(bins[1:] + bins[:-1])
        ax.errorbar(centers, n, yerr=np.sqrt(n), fmt=".", ms=3)
        ax.set_xlabel("Track $p_T$ [GeV]", fontsize=15)
        ax.set_ylabel("Counts", fontsize=15)
        ax.tick_params(axis="both", labelsize=14)
        pdf.savefig(fig)
        plt.close(fig)

        reso = np.array(res.total_pT_resolution) * 100 
        fig, ax = plt.subplots(figsize=(7,6))
        n, bins, _ = ax.hist(reso, bins=50, histtype="stepfilled", alpha=0.25)
        centers=0.5*(bins[1:] + bins[:-1])
        ax.errorbar(centers, n, yerr=np.sqrt(n), fmt='.', ms=3)
        ax.set_xlabel(r"$(p_T^{\mathrm{reco}}-p_T^{\mathrm{truth}})/p_T^{\mathrm{truth}}$ [%]", fontsize=15)
        ax.set_ylabel("Counts", fontsize=15)
        ax.tick_params(axis='both', labelsize=13)
        pdf.savefig(fig)
        plt.close(fig)

        bins_pt = np.linspace(0,3500,51)
        n_yes, _ = np.histogram(res.yes_ob_pTs, bins=bins_pt)
        n_no, _ = np.histogram(res.no_ob_pTs, bins=bins_pt)
        centers = 0.5 * (bins_pt[1:] + bins_pt[:-1])
        fig, ax = plt.subplots(figsize=(7,6))
        ax.step(bins_pt[:-1], n_yes, where="post", label="Has OB hit", color="C0")
        ax.step(bins_pt[:-1], n_no,  where="post", label="No OB hit",  color="C1")
        ax.errorbar(centers, n_yes, yerr=np.sqrt(n_yes), fmt='.', ms=3, color="C0")
        ax.errorbar(centers, n_no,  yerr=np.sqrt(n_no),  fmt='.', ms=3, color="C1")
        ax.set_xlabel(r"Track $p_T$ [GeV]", fontsize=15)
        ax.set_ylabel("Counts", fontsize=15)
        ax.tick_params(axis="both", labelsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend()
        pdf.savefig(fig); plt.close(fig)

        bins_res = np.linspace(-100, 10, 51)
        yes_ob_pT_res = np.array(res.yes_ob_pT_res) * 100
        no_ob_pT_res = np.array(res.no_ob_pT_res) * 100
        n_yes, _ = np.histogram(yes_ob_pT_res, bins=bins_res)
        n_no, _ = np.histogram(no_ob_pT_res, bins=bins_res)
        centers = 0.5 * (bins_res[1:] + bins_res[:-1])
        fig, ax = plt.subplots(figsize=(7,6))
        ax.step(bins_res[:-1], n_yes, where="post", label="Has OB hit", color="C0")
        ax.step(bins_res[:-1], n_no,  where="post", label="No OB hit",  color="C1")
        ax.errorbar(centers, n_yes, yerr=np.sqrt(n_yes), fmt='.', ms=3, color="C0")
        ax.errorbar(centers, n_no,  yerr=np.sqrt(n_no),  fmt='.', ms=3, color="C1")
        ax.set_xlabel(r"$(p_T^{\rm reco}-p_T^{\rm truth})/p_T^{\rm truth}$ [%]", fontsize=15)
        ax.set_ylabel("Counts", fontsize=15)
        ax.tick_params(axis="both", labelsize=14)
        ax.grid(True, alpha=0.3)
        ax.legend()
        pdf.savefig(fig); plt.close(fig)

    print(f"Plots saved to {plot_path}")    


res = analyze_single_slcio(reco_slcio_path, event_print=True, track_print=True, redo=True)
make_pt_plots(res)
print(f"For the {len(res.yes_ob_pTs)} tracks that have at least one hit in OB, mean resolution is {np.mean(res.yes_ob_pT_res)*100:.2f}% and mean pT is {np.mean(res.yes_ob_pTs):.2f} GeV")
print(f"For the {len(res.no_ob_pTs)} tracks that DON'T have at least one hit in OB, mean resolution is {np.mean(res.no_ob_pT_res)*100:.2f}% and mean pT is {np.mean(res.no_ob_pTs):.2f} GeV")