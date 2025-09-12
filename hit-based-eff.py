import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pyLCIO
from pyLCIO import UTIL
import ROOT

from tqdm import tqdm

samples = ["1000_10", "2500_10", "3000_10", "3500_10", "4000_10", "4500_10"]
windows = ["loose", "medium", "tight"]

reg_window_path = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/efficiency/nobib/"
# mira_path = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/mira_time/v2/nobib/"
pdf_path = "/scratch/miralittmann/analysis/mira_analysis_code/hit-based-eff_track_agnostic.pdf"

n_files = 50
stau_ids = {1000015, -1000015, 2000015, -2000015}
layer_map = [
    ("VB", 0, 30), ("VB", 1, 30), ("VB", 2, 51), ("VB", 3, 51),
    ("VB", 4, 74), ("VB", 5, 74), ("VB", 6, 102), ("VB", 7, 102),
    ("IB", 0, 127), ("IB", 1, 340), ("IB", 2, 554),
    ("OB", 0, 819), ("OB", 1, 1153), ("OB", 2, 1486),
]
system_to_det = {
    1: "VB", 2: "VE",
    3: "IB",  4: "IE",
    5: "OB",  6: "OE",
}

L = len(layer_map)
eff_stats = {}

# for each sample, window, file: 
# for each stau in the event: get how long it's traveled from SIM (r distance in mm) --> assume it's passed thru every layer 
# --> get boolean array of if we expect it to have passed through each layer 

# then, for same stau (use mcp id), get the RECO hits it's made in each subdetector, layer, and make another boolean array for this

def build_rel_nav(event):
    nav = {
        "VXDBarrel": UTIL.LCRelationNavigator(event.getCollection("VXDBarrelHitsRelations")),
        "VXDEndcap": UTIL.LCRelationNavigator(event.getCollection("VXDEndcapHitsRelations")),
        "ITBarrel" : UTIL.LCRelationNavigator(event.getCollection("ITBarrelHitsRelations")),
        "ITEndcap" : UTIL.LCRelationNavigator(event.getCollection("ITEndcapHitsRelations")),
        "OTBarrel" : UTIL.LCRelationNavigator(event.getCollection("OTBarrelHitsRelations")),
        "OTEndcap" : UTIL.LCRelationNavigator(event.getCollection("OTEndcapHitsRelations")),
    }
    nav["_ENCODING"] = event.getCollection("ITBarrelHits").getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
    return nav

layer_radii = np.array([r for _,_,r in layer_map])
def exp_hits_mask(dist_mm: float)->np.ndarray:
                    return layer_radii <= dist_mm
    
def decode_system_layer(hit, encoding_str):
    dec = UTIL.BitField64(encoding_str)
    dec.setValue(int(hit.getCellID0()))
    system = dec["system"].value()
    layer  = dec["layer"].value()
    return system, layer

layer_index = {(det, layer): i for i, (det, layer, _) in enumerate(layer_map)}

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
for i in tqdm(range(n_files)):
    for sample in samples: 
        for window in windows:
            key_sw = (sample, window)
            if key_sw not in eff_stats:
                eff_stats[key_sw] = {
                     "exp": np.zeros(L, dtype=np.int64),
                     "obs": np.zeros(L, dtype=np.int64),
                }
            
            # if window == "mira_time":
            #     fname = os.path.join(mira_path, sample, f"{sample}_reco{i}.slcio")
            # else:
            fname = os.path.join(reg_window_path,window,sample,f"{sample}_reco{i}.slcio")
            if not os.path.exists(fname):
                continue

            reader.open(fname)

            for event in reader:
                rel_nav = build_rel_nav(event)

                # rel_collection   = event.getCollection("MCParticle_SiTracks")
                # track_collection = event.getCollection("SiTracks")
                # sim_particle_collection = event.getCollection("MCParticle")

                # # SIM: expected hits
                # mcp_reach = {}
                # for mcp in sim_particle_collection:
                #     if abs(mcp.getPDG()) not in stau_ids:
                #         continue
                #     r_end = np.hypot(mcp.getEndpoint()[0], mcp.getEndpoint()[1])
                #     if r_end >= 127:
                #         m = exp_hits_mask(r_end)
                #         mcp_reach[mcp.id()] = m
                #         eff_stats[key_sw]["exp"] += m.astype(np.int64)

                # RECO: observed hits
                # seen = set()
                # for idx in range(rel_collection.getNumberOfElements()):
                #     rel = rel_collection.getElementAt(idx)
                #     mcp_from_relation = rel.getFrom()
                #     track = rel.getTo()
                                    
                #     if not mcp_from_relation or abs(mcp_from_relation.getPDG()) not in stau_ids:
                #         continue
                #     if mcp_from_relation.id() in seen:
                #         continue
                #     seen.add(mcp_from_relation.id())
                    
                #     exp_mask = mcp_reach.get(mcp_from_relation.id())
                #     if exp_mask is None:
                #         continue
                                                            
                enc = rel_nav["_ENCODING"]

                sim_specs = [
                    ("VXDBarrel", "VXDBarrelCollection", "VXDBarrel"),
                    ("ITBarrel",  "ITBarrelCollection",  "ITBarrel"),
                    ("OTBarrel",  "OTBarrelCollection",  "OTBarrel"),
                ]

                for det_key, sim_coll_name, nav_key in sim_specs:
                    if sim_coll_name not in event.getCollectionNames():
                        continue
                    sim_coll = event.getCollection(sim_coll_name)
                    nav = rel_nav[nav_key]

                    for simhit in sim_coll:
                        system, layer = decode_system_layer(simhit, enc)
                        det_name = system_to_det.get(system)
                        if det_name in ("VE", "IE", "OE"):
                            continue

                        idx = layer_index.get((det_name[:2], layer))
                        if idx is None:
                            continue
                        
                        eff_stats[key_sw]["exp"][idx] += 1

                        related_reco_hits = nav.getRelatedFromObjects(simhit) or []
                        if related_reco_hits:
                            eff_stats[key_sw]["obs"][idx] += 1
                    

                    # obs_hits = np.zeros(len(layer_map),dtype=bool)
                    # for hit in track.getTrackerHits():
                    #     system, layer = decode_system_layer(hit, enc) 
                    #     det_name = system_to_det.get(system)
                    #     if det_name in ("VE", "IE", "OE"):
                    #         continue
                    #     key = (det_name[:2], layer) 
                    #     obs_hits[layer_index[key]] = True 

                    # eff_stats[key_sw]["obs"] += (obs_hits & exp_mask).astype(np.int64) 
                    
            reader.close()

print(eff_stats)

# plotting
hit_eff_data = {}
for (sample, window), counts in eff_stats.items():
    exp = counts["exp"]
    obs = counts["obs"]
    eff = np.divide(obs, exp, out=np.full(L, np.nan, dtype=float), where=exp>0)
    err = np.sqrt(eff*(1-eff)/exp)
    err = np.where(exp>0, err, np.nan)

    if window not in hit_eff_data:
        hit_eff_data[window] = {}
    hit_eff_data[window][sample] = {"eff": eff, "err": err}
    
layer_names = ["VB0","VB1","VB2","VB3","VB4","VB5","VB6","VB7","IB0","IB1","IB2","OB0","OB1","OB2"]
x = np.arange(len(layer_names))
width = 0.30
offset = {"tight": -width, "medium": 0, "loose": width}
colors = {"tight":"maroon", "medium": "teal", "loose":"palevioletred", "mira_time": "goldenrod", "pt10": "tab:blue", "pt5": "tab:green"}

sample_to_mass = {"1000_10": "1 TeV Staus", 
                  "2500_10": "2.5 TeV Staus",
                  "3000_10": "3 TeV Staus",
                  "3500_10": "3.5 TeV Staus", 
                  "4000_10": "4 TeV Staus",
                  "4500_10": "4.5 TeV Staus"}
window_labels = {"loose": "Loose", 
                 "medium": "Medium",
                 "tight": "Nominal", 
                 "mira_time": "Newest"}
with PdfPages(pdf_path) as pdf:
    for sample in samples:
        fig, ax = plt.subplots(figsize=(9,4))
        for window in windows:
            eff = np.asarray(hit_eff_data[window][sample]["eff"], dtype=float) * 100.0
            err = np.asarray(hit_eff_data[window][sample]["err"], dtype=float) * 100.0
            
            eff_clean = np.where(np.isfinite(eff), eff, 0.0)
            err_clean = np.where(np.isfinite(err), err, 0.0)
            
            ax.bar(x + offset.get(window, 0.0), eff_clean, width=width,
                   yerr=err_clean, label=window_labels[window],
                   color=colors.get(window, None), capsize=2, linewidth=0)
        
        ax.set_xticks(x, layer_names, ha="right", fontsize=14, rotation=45)
        ax.set_ylim(0, 105)
        ax.set_ylabel("Hit reconstruction efficiency (%)", fontsize=14)
        ax.set_xlabel("Detector Subsystem, Layer", fontsize=14)
        ax.tick_params(axis='y', labelsize=14)

        for xc in (7.5, 10.5):
            ax.axvline(x=xc, color='k', ls='--', alpha=0.5)
        ax.legend(title=f"{sample_to_mass[sample]}", loc="upper right")
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

print(f"saved plots to {pdf_path}")
