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

print("Analyzing reco files...")
# for each of the masses we're looking at 
for sample in sample_names:
    reco_path = os.path.join(reco_dir, sample)
    file_prefix = f"{sample}_reco"
    # for each of the 500 files
    for i in tqdm(range(9)): # TODO: cahnge back to 500 
        file_name = f"{file_prefix}{i}.slcio"
        file_path = os.path.join(reco_path, "bib/", file_name)

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
                    encoding = event.getCollection("ITBarrelHits").getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
                    decoder = pyLCIO.UTIL.BitField64(encoding)
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
                    flight_time = distance/c

                    theta = np.arccos(hit_z_pos / distance)

                    reco_info[sample][detector_key][layer].append((hit.getTime(), hit_z_pos, theta))
            seen_staus.add(mcp.id())
                    
        reader.close()

################################# converting averages to dataframe to check before making plots ##########################
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
        




#################### Plotting ############################
plot_dir = "/scratch/miralittmann/analysis/mira_analysis_code/first_principles/plots/"

histograms = False
vs_theta = False
vs_z = False 
colors = {1.0: 'r', 2.5: 'g', 4.0: 'b', 4.5: 'c'}

pdf_path = os.path.join(plot_dir, "vs_theta.pdf") # CHANGE NAME HERE
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

                times, z, theta = zip(*triples)
                times = list(times)
                z = list(z)
                theta = list(theta)
                
                fig, ax = plt.subplots()
                if histograms == True:
                    ax.hist(times, bins=100, histtype="stepfilled")
                    ax.axvline(x=stau_tof_df.loc[mass, f"{subdet}_{layer}"], color='r', linestyle = '--')
                    ax.set_title(f"{mass} TeV: {subdet} layer {layer}")
                    ax.set_ylabel("hits")
                    ax.set_xlabel("corrected time [ns]") 

                if vs_z == True:
                    ax.scatter(z, times)
                    ax.set_title(f"{mass} TeV: {subdet} layer {layer}")
                    ax.set_ylabel("corrected time [ns]")
                    ax.set_xlabel("z [mm]")
                
                if vs_theta == True:
                    ax.scatter(theta, times)
                    ax.set_title(f"{mass} TeV: {subdet} layer {layer}")
                    ax.set_ylabel("corrected time [ns]")
                    ax.set_xlabel("Theta [rad]")

                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)
