import numpy as np
import pandas as pd
import os 
import pyLCIO
from pyLCIO import UTIL
from math import *
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from tqdm import tqdm
import pickle 
import pathlib
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--rebuild", action="store_true", help="rescan the slcio files")
parser.add_argument("--get_avgs", action="store_true", help="get the average hit time values from all the files, print out a chart")
parser.add_argument("--bib_reco", action="store_true", help="analyze all signals from bib reco files") # have to add this if you want to plot bib stuff 

parser.add_argument("--histograms", action="store_true", help="plot time distributions for each layer, detector signal")
parser.add_argument("--signal_vs_theta", action="store_true", help="plot times vs theta for each layer, detector")
parser.add_argument("--signal_vs_z", action="store_true", help="plot times vs z coord for each layer, detector")

parser.add_argument("--bib_vs_theta", action="store_true", help="plot bib times versus theta")
parser.add_argument("--bib_vs_z", action="store_true", help="plot bib times versus z")
parser.add_argument("--bibstograms", action="store_true", help="plot all hits time distributions with signal windows marked")

args = parser.parse_args()

bib = args.bib_reco
get_avgs = args.get_avgs
redo = args.rebuild
histograms = args.histograms
signal_vs_theta = args.signal_vs_theta
signal_vs_z = args.signal_vs_z
bibstograms = args.bibstograms
bib_vs_theta = args.bib_vs_theta
bib_vs_z = args.bib_vs_z

############# ESTIMATES ####################
stau_masses = np.array([1, 2.5, 4, 4.5]) # in TeV. just using what we have rn in the tbl files, could do more
com_energy = 10 # also in TeV, CoM energy in muon collider
stau_energy = com_energy / 2 # assuming it's split evenly (is this fair?)
c = 299792458/1000000  # mm/ns 

# E^2 = m^2c^4 + p^2c^2

# locations of subdetector layers (in mm) (just barrel)
subdet_locs = {
    "VB_0": 30,
    "VB_1": 32,
    "VB_2": 51,
    "VB_3": 53,
    "VB_4": 74,
    "VB_5": 76,
    "VB_6": 102,
    "VB_7": 103,

    "IB_0": 127,
    "IB_1": 340,
    "IB_2": 554,

    "OB_0": 819,
    "OB_1": 1153,
    "OB_2": 1486
}

def get_tof_reg():
    p = np.sqrt(stau_energy**2 - stau_masses**2)
    beta = p / stau_energy
    v = beta * c

    df = pd.DataFrame(index=stau_masses)
    df.index.name = "mass [TeV]"

    for subdet, distance in subdet_locs.items():
        stau_tof_raw = distance / v
        
        # for corrected time: computing time it would take particle traveling at speed of light to get to each subdet
        ref_particle_tof = distance / c

        corr_tof = stau_tof_raw - ref_particle_tof
        df[subdet] = corr_tof
 
    return df

stau_tof_df = get_tof_reg()
print("\n===============================================")
print("Mira's calculations")
print(stau_tof_df)
print(" ")

################ Want to plot hit time distributions for different subdetectors, layers (vs z or theta) ################################
reco_dir = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/first_principles/"
sample_names = ["1000_10", "2500_10", "4000_10", "4500_10"] 
subdetectors = ["VB", "IB", "OB"]
# rand = ROOT.TRandom3(0) 
system_map = {1: "VB", 2: "VE", 3: "IB", 4: "IE", 5: "OB", 6: "OE"}
stau_ids = [1000015, 2000015]
num_files = 9

CACHE = pathlib.Path("cache/reco_bib.pkl")

def build_analysis(sample_names, num_files, redo=False):

    if CACHE.exists() and not redo:
        print(f"Loading previous arrays from {CACHE}, not redoing full analysis")
        with CACHE.open("rb") as f:
            reco_info, all_hits = pickle.load(f)
        return reco_info, all_hits
    
    
    # collecting info. Sample --> Detector --> Layer --> (Hit times, Hit z, Hit theta)
    reco_info = {sample: {
        "VB": {
            0: [],
            1: [],
            2: [],
            3: [],
            4: [],
            5: [],
            6: [],
            7: [] 
        },
        "IB": {
            0: [],
            1: [],
            2: [],
    },
        "OB": { 
            0: [],
            1: [],
            2: [],
    }
    } for sample in sample_names}

    all_hits = {sample: {
        "VB": {
            0: [],
            1: [],
            2: [],
            3: [],
            4: [],
            5: [],
            6: [],
            7: []
        },
        "IB": {
            0: [],
            1: [],
            2: []
        },
        "OB": {
            0: [],
            1: [],
            2: []
        }
    } for sample in sample_names}

    collection_mapping = {"VB": "VXDBarrelHits", 
                        "VE": "VXDEndcapHits", 
                        "IB": "ITBarrelHits", 
                        "IE": "ITEndcapHits",
                        "OB": "OTBarrelHits",
                        "OE": "OTEndcapHits"}

    print("Analyzing reco files...")
    # for each of the masses we're looking at 
    for sample in sample_names:
        reco_path = os.path.join(reco_dir, sample)
        file_prefix = f"{sample}_reco"
        # for each of the 500 files
        for i in tqdm(range(num_files)): 
            file_name = f"{file_prefix}{i}.slcio"
            file_path = os.path.join(reco_path, "no_bib/", file_name)

            if not os.path.isfile(file_path):
                print(f"File not found: {file_path}")
                continue
            
            reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
            reader.open(file_path)

            event = reader.readNextEvent()
            if event is None:
                print(f"No event in: {file_path}")
                reader.close()
                continue
        
            encoding = event.getCollection("ITBarrelHits").getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
            decoder = pyLCIO.UTIL.BitField64(encoding)
    
            # print(f"Analyzing {file_path}")
            seen_staus = set()
            mcp_collection = event.getCollection("MCParticle")
            rel_collection = event.getCollection("MCParticle_SiTracks_Refitted")
            relation = pyLCIO.UTIL.LCRelationNavigator(rel_collection)

            # for one of the staus in the event. filtering for staus that are non-intermediate and unique
            for mcp in mcp_collection:
                if abs(mcp.getPDG()) not in stau_ids:
                    continue 
                if any(abs(parent.getPDG()) in stau_ids for parent in mcp.getParents()):
                    continue
                if mcp.id() in seen_staus:
                    continue

                stau_tracks = relation.getRelatedToObjects(mcp)
                for track in stau_tracks:
                    seen_layers = set()
                    for hit in track.getTrackerHits():
                        cellID = int(hit.getCellID0())
                        decoder.setValue(cellID)

                        detector = decoder["system"].value()
                        detector_key = system_map[detector]
                        layer = decoder["layer"].value()
                        if (detector,layer) in seen_layers: 
                            continue
                        if detector in (2, 4, 6):
                            continue
                        seen_layers.add((detector,layer))

                        hit_x_pos = hit.getPosition()[0] 
                        hit_y_pos = hit.getPosition()[1] 
                        hit_z_pos = hit.getPosition()[2] 

                        distance = sqrt(hit_x_pos**2 + hit_y_pos**2 + hit_z_pos**2) 

                        theta = np.arccos(hit_z_pos / distance)

                        reco_info[sample][detector_key][layer].append((hit.getTime(), hit_z_pos, theta))
                seen_staus.add(mcp.id())                                    
            reader.close()


        ###################################### BIB hits ##########################################
        if bib == True: 
            print("Analyzing all hits from BIB files...")
            for i in tqdm(range(num_files)): 
                file_name = f"{file_prefix}{i}.slcio"
        
                bib_file_path = os.path.join(reco_path, "bib/", file_name)
                if not os.path.isfile(bib_file_path):
                    print(f"File not found: {bib_file_path}")
                    continue
    
                reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
                reader.open(file_path)

                event = reader.readNextEvent()
                if event is None:
                    print(f"No event in: {file_path}")
                    reader.close()
                    continue
        
                encoding = event.getCollection("ITBarrelHits").getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                decoder = pyLCIO.UTIL.BitField64(encoding)
    
                reader.open(bib_file_path)
                event = reader.readNextEvent()
                if event is None:
                    print(f"No event in {bib_file_path}")
                    continue
                
                for detector, collection_name in collection_mapping.items():
                    if detector in ("VE", "IE", "OE"):
                        continue

                    hit_collection = event.getCollection(collection_name)
                    for i in range(hit_collection.getNumberOfElements()):
                        hit = hit_collection.getElementAt(i)
                        hit_time = hit.getTime()
                        cellID = int(hit.getCellID0())
                        decoder.setValue(cellID)

                        layer = decoder["layer"].value()

                        x = hit.getPosition()[0]
                        y = hit.getPosition()[1]
                        z = hit.getPosition()[2]
                        distance = np.sqrt(x**2 + y**2 + z**2)
                        theta = np.arccos(z / distance)

                        all_hits[sample][detector][layer].append((hit.getTime(), z, theta))
                reader.close()
        
    print(f"Writing cache -> {CACHE}")
    CACHE.parent.mkdir(exist_ok=True)
    with CACHE.open("wb") as f:
        pickle.dump((reco_info, all_hits), f, protocol=pickle.HIGHEST_PROTOCOL)
    return reco_info, all_hits

reco_info, all_hits = build_analysis(sample_names, num_files, redo=args.rebuild)


################################# converting averages to dataframe to check before making plots ##########################

if get_avgs == True:
    rows = []
    for sample, det_dict in reco_info.items():
        for subdet, layer_dict in det_dict.items():
            for layer, hits in layer_dict.items():
                times = [hit[0] for hit in hits]
                
                avg_time = sum(times) / len(times) if times else math.nan

                rows.append(
                    {
                        "mass": sample,
                        "subdet": subdet,
                        "layer": layer,
                        "avg_time": avg_time
                    })

    df = pd.DataFrame(rows)
    pivot = df.pivot_table(index=["mass", "subdet"], columns="layer", values="avg_time")
    df["mass_TeV"] = (df["mass"].str.split("_", expand=True)[0] .astype(int).div(1000.0))
    df["det_layer"] = df["subdet"] + "_" + df["layer"].astype(str)
    reco_df = (df.pivot_table(index="mass_TeV", columns="det_layer",values="avg_time").sort_index(axis=1, key=lambda idx: [("VB","IB","OB").index(col.split("_")[0]) * 10 + int(col.split("_")[1]) for col in idx]))
    reco_df.index.name = "mass [TeV]"

    pd.set_option("display.max_columns", None)
    print("========================================")
    print("Average Reconstructed Times")
    print(" ")
    print(reco_df.to_string(float_format="{:0.7f}".format))
    print(" ")

    colors = {1.0: 'r', 2.5: 'g', 4.0: 'b', 4.5: 'c'}
    layers_order = reco_df.columns
    x_pos = range(len(layers_order))
    avg_pdf_path = "/scratch/miralittmann/analysis/mira_analysis_code/first_principles/plots/avgs_vs_layer.pdf"
    avg_pdf = PdfPages(avg_pdf_path)
    with PdfPages(avg_pdf_path) as pdf:
        fig, ax = plt.subplots(figsize=(9,4.5))
        for mass in reco_df.index:
            y = reco_df.loc[mass].values
            ax.plot(x_pos, y, marker='o',linewidth=1.6, label=f"{mass:.1f} TeV", color=colors.get(mass, 'k'))
            ax.set_xticks(x_pos)
        ax.set_xticklabels(layers_order, rotation=45, ha='right')
        ax.set_xlabel("Tracker layer")
        ax.set_ylabel("⟨t_corr⟩  [ns]")
        ax.set_title("Average corrected time vs. detector layer")
        ax.grid(alpha=0.3, linestyle='--', linewidth=0.5)
        ax.legend(title="Stau mass")

        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
            
#################################################### Plotting ##########################################################
plot_dir = "/scratch/miralittmann/analysis/mira_analysis_code/first_principles/plots/"
colors = {1.0: 'r', 2.5: 'g', 4.0: 'b', 4.5: 'c'}

pdf_path = os.path.join(plot_dir, "bib_vs_z.pdf") # CHANGE NAME HERE
pdf = PdfPages(pdf_path)

with PdfPages(pdf_path) as pdf:
    for sample in sample_names:
        mass = int(sample.split("_")[0]) / 1000.0 
        for subdet in subdetectors:
            if subdet == "VB":
                num_layers = 8
            else:
                num_layers = 3

            for layer in range(num_layers):
                triples = reco_info[sample][subdet][layer]
                if not triples:
                    continue

                signal_times, signal_z, signal_theta = zip(*triples)
                signal_times = list(signal_times)
                signal_z = list(signal_z)
                signal_theta = list(signal_theta)
                
                fig, ax = plt.subplots()
                if histograms == True:
                    ax.hist(signal_times, bins=100, histtype="stepfilled")
                    ax.axvline(x=stau_tof_df.loc[mass, f"{subdet}_{layer}"], color='r', linestyle = '--')
                    ax.set_title(f"{mass} TeV: {subdet} layer {layer}")
                    ax.set_ylabel("hits")
                    ax.set_xlabel("corrected time [ns]") 

                if signal_vs_z == True:
                    ax.scatter(signal_z, signal_times)
                    ax.set_title(f"{mass} TeV: {subdet} layer {layer}")
                    ax.set_ylabel("corrected time [ns]")
                    ax.set_xlabel("z [mm]")
                
                if signal_vs_theta == True:
                    ax.scatter(signal_theta, signal_times)
                    ax.set_title(f"{mass} TeV: {subdet} layer {layer}")
                    ax.set_ylabel("corrected time [ns]")
                    ax.set_xlabel("Theta [rad]")


                if bib == True:
                    bib_triples = all_hits[sample][subdet][layer]
                    if not bib_triples:
                        continue

                    bib_times, bib_z, bib_theta = zip(*bib_triples)
                    bib_times = list(bib_times)
                    bib_z = list(bib_z)
                    bib_theta = list(bib_theta)

 
                    if bibstograms == True:
                        weights = np.full(len(bib_times), 1.0 /num_files)

                        n, bins, _ = ax.hist(bib_times, bins=100, weights=weights)
                        max_bin = n.max()
                        signal_min = min(signal_times)
                        signal_max = max(signal_times)
                        
                        ax.axvline(x=signal_min, color='r', linestyle='--', label='signal time region')
                        ax.axvline(x=signal_max, color='r', linestyle='--')
                        ax.set_title(f"{mass} TeV: {subdet} layer {layer} [ALL HITS, BIB]")
                        ax.legend()
                        ax.set_ylabel("average hits per event")
                        ax.set_xlabel("corrected time [ns]")
                        ax.text(0.97,0.03, f"max hits {max_bin:.0f}", transform=ax.transAxes, ha="right", va="bottom", bbox=dict(boxstyle="round,pad=0.25",facecolor="white", alpha=0.6))

                        if subdet == "VB":
                            ax.set_xlim(0,0.32)
                        elif subdet == "IB":
                            ax.set_xlim(0,3)
                            ax.set_ylim(top=1000)
                        else: 
                            ax.set_xlim(0,8)
                            ax.set_ylim(top=1000)

                    if bib_vs_theta == True:
                        if sample != "1000_10": # because this should be the same for every mass ?
                            continue
                        ax.scatter(bib_theta, bib_times)
                        ax.set_title(f"{mass} TeV: {subdet} layer {layer} [ALL HITS, BIB]")
                        ax.set_ylabel("corrected time [ns]")
                        ax.set_xlabel("Theta [rad]")

                    if bib_vs_z == True:
                        if sample != "1000_10": # because this should be the same for every mass ?
                            continue
                        ax.scatter(bib_z, bib_times)
                        ax.set_title(f"{mass} TeV: {subdet} layer {layer} [ALL HITS, BIB]")
                        ax.set_ylabel("corrected time [ns]")
                        ax.set_xlabel("z coordinate [mm]")

                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

print(f"Wrote plot(s) to {pdf_path}")