import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from collections import OrderedDict

import pyLCIO
from pyLCIO import UTIL
import ROOT

from tqdm import tqdm
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--rebuild", action="store_true")
args = parser.parse_args()

samples = ["1000_10", "2500_10", "3000_10", "3500_10", "4000_10", "4500_10"]
windows = ["tight", "medium", "loose"]

tight_path = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/efficiency/seeding_10GeV/refit/bib/tight/"
loose_path = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/loose4/seeding_10GeV/refit/bib/"
medium_path = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/p26_p5_p6/seeding_10GeV/refit/bib/"
window_paths = [tight_path, loose_path, medium_path]
window_to_path = {"tight":tight_path,
                  "medium":medium_path,
                  "loose":loose_path}

pdf_path = "/scratch/miralittmann/analysis/mira_analysis_code/n-hits-study-bib-refit.pdf"

n_files = 500
stau_ids = {1000015, -1000015, 2000015, -2000015}

system_to_det = {
    1: "VB", 2: "VE",
    3: "IB",  4: "IE",
    5: "OB",  6: "OE",
}
system_to_relname = {
    1: "VXDBarrel", 2: "VXDEndcap",
    3: "ITBarrel",  4: "ITEndcap",
    5: "OTBarrel",  6: "OTEndcap",
}

rebuild = args.rebuild # REMEMBER TO CHANGE THIS IF YOU WANT TO RESAVE THE INFO
cache_path = "/scratch/miralittmann/analysis/mira_analysis_code/track_n_hits_cache-refit.pkl"

def build_rel_nav(event):
    nav = {
        "VXDBarrel": UTIL.LCRelationNavigator(event.getCollection("VXDBarrelHitsRelations")),
        "VXDEndcap": UTIL.LCRelationNavigator(event.getCollection("VXDEndcapHitsRelations")),
        "ITBarrel" : UTIL.LCRelationNavigator(event.getCollection("ITBarrelHitsRelations")),
        "ITEndcap" : UTIL.LCRelationNavigator(event.getCollection("ITEndcapHitsRelations")),
        "OTBarrel" : UTIL.LCRelationNavigator(event.getCollection("OTBarrelHitsRelations")),
        "OTEndcap" : UTIL.LCRelationNavigator(event.getCollection("OTEndcapHitsRelations")),
    }
    enc = event.getCollection("ITBarrelHits").getParameters().getStringVal(pyLCIO.EVENT.LCIO.CellIDEncoding)
    nav["_ENCODING"] = enc
    nav["_DECODER"] = UTIL.BitField64(enc)
    return nav

def decode_system_layer(hit, decoder):
    decoder.setValue(int(hit.getCellID0())) 
    system = decoder["system"].value()
    layer  = decoder["layer"].value()
    return system, layer

def is_stau_hit(hit, rel_nav, system):
    det_rel = rel_nav.get(system_to_relname.get(system))
    if det_rel is None:
        return False
    related = det_rel.getRelatedToObjects(hit)  
    if not related:
        return False  
    for sim in related:
        mcp = sim.getMCParticle()
        if mcp and abs(mcp.getPDG()) in stau_ids:
            return True
    return False

reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()
track_n_hits = None
if (not rebuild) and os.path.exists(cache_path):
    with open(cache_path, "rb") as f:
        print("loading in cached arrays...")
        track_n_hits = pickle.load(f)

if track_n_hits is None:
    print("rebuilding analysis arrays...")
    track_n_hits = {
        window: {
            sample: {
                "n_vb": [],
                "n_ib": [],
                "n_ob": [],
                "n_tot": [],
                "pcent_stau": [],
                "chi2": []
            }
            for sample in samples
        }
        for window in windows
    }
    for i in tqdm(range(n_files)):
        for sample in samples: 
            for window in windows:                 
                fname = os.path.join(window_to_path[window],sample,f"{sample}_reco{i}.slcio")
                if not os.path.exists(fname):
                    continue
                
                
                reader.open(fname)
                for event in reader:
                    rel_nav = build_rel_nav(event)

                    rel_collection = event.getCollection("MCParticle_SiTracks_Refitted") if "MCParticle_SiTracks_Refitted" in event.getCollectionNames() else None # _Refitted
                    if not rel_collection:
                        continue

                    for idx in range(rel_collection.getNumberOfElements()):
                        rel = rel_collection.getElementAt(idx)
                        mcp = rel.getFrom()
                        track = rel.getTo()
                    
                        # only looking at stau tracks
                        if not mcp or abs(mcp.getPDG()) not in stau_ids:
                                continue
                        enc = rel_nav["_ENCODING"]

                        n_vb_hits = 0
                        n_ib_hits = 0
                        n_ob_hits = 0
                        
                        stau_hit_count = 0
                        hit_count = 0

                        for hit in track.getTrackerHits(): 
                            system, layer = decode_system_layer(hit, rel_nav["_DECODER"])
                            det_name = system_to_det.get(system)
                            if det_name == "VB":
                                n_vb_hits += 1
                            elif det_name == "IB":
                                n_ib_hits += 1
                            elif det_name == "OB":
                                n_ob_hits += 1

                            for sim in rel_nav[system_to_relname[system]].getRelatedToObjects(hit):
                                mcp = sim.getMCParticle()
                                if mcp and abs(mcp.getPDG()) in stau_ids:
                                    stau_hit_count += 1
                            
                            hit_count += 1                                                             
                            
                        total_track_hits = n_vb_hits + n_ib_hits + n_ob_hits
                        track_pcent_stau = (stau_hit_count / len(track.getTrackerHits())) if len(track.getTrackerHits()) else np.nan
                        red_chi2 = track.getChi2()/track.getNdf() 

                        track_n_hits[window][sample]["n_tot"].append(total_track_hits)
                        track_n_hits[window][sample]["n_vb"].append(n_vb_hits)
                        track_n_hits[window][sample]["n_ib"].append(n_ib_hits)
                        track_n_hits[window][sample]["n_ob"].append(n_ob_hits)
                        track_n_hits[window][sample]["pcent_stau"].append(track_pcent_stau)
                        track_n_hits[window][sample]["chi2"].append(red_chi2)

                reader.close()
    os.makedirs(os.path.dirname(cache_path), exist_ok=True)
    with open(cache_path, "wb") as f:
        print("saving to cache...")
        pickle.dump(track_n_hits, f, protocol=pickle.HIGHEST_PROTOCOL)
        print("cache saved.")

markers = {"VB": "v", "IB": "o", "OB": "^"}
sample_to_mass_labels = {"1000_10": "1 TeV Staus", 
                  "2500_10": "2.5 TeV Staus",
                  "3000_10": "3 TeV Staus",
                  "3500_10": "3.5 TeV Staus", 
                  "4000_10": "4 TeV Staus",
                  "4500_10": "4.5 TeV Staus"}
window_labels = {"loose": "Loose", 
                 "medium": "Medium",
                 "tight": "Nominal", 
                 "mira_time": "Newest"}
sample_to_mass = {"1000_10": 1,
                  "2500_10": 2.5, 
                  "3000_10": 3,
                  "3500_10": 3.5,
                  "4000_10": 4,
                  "4500_10": 4.5}


with PdfPages(pdf_path) as pdf:
    plt.style.use("seaborn-v0_8-colorblind")
    dets  = ["VB", "IB", "OB", "Total"]
    x = np.arange(len(dets))
    width= 0.24 

    fig, axes = plt.subplots(1, 3, figsize=(18,6), sharex=True, sharey=True)
    for window in windows:
        masses = []
        vb_pcent_stau = []
        ib_pcent_stau = []
        ob_pcent_stau = []
        for sample in samples:
            masses.append(sample_to_mass[sample])
            vb = np.asarray(track_n_hits[window][sample]["n_vb"], dtype=float)
            ib = np.asarray(track_n_hits[window][sample]["n_ib"], dtype=float)
            ob = np.asarray(track_n_hits[window][sample]["n_ob"], dtype=float)
            chi2 = np.asarray(track_n_hits[window][sample]["chi2"], dtype=float)

            pcent = np.asarray(track_n_hits[window][sample]["pcent_stau"], dtype=float)*100

            vb_mask = (vb>=1.5) #& (chi2<5)
            ib_mask = (vb>=1.5) & (ib>=2) #& (chi2<5) 
            ob_mask = (vb>=1.5) & (ib>=2) & (ob>=2) #& (chi2<5)

            vb_values = pcent[vb_mask]
            ib_values = pcent[ib_mask]
            ob_values = pcent[ob_mask]
            vb_pcent_stau.append(np.mean(vb_values) if vb_values.size else np.nan) 
            ib_pcent_stau.append(np.mean(ib_values) if ib_values.size else np.nan)
            ob_pcent_stau.append(np.mean(ob_values) if ob_values.size else np.nan)

        order = np.argsort(masses)
        masses = np.array(masses)[order]
        vb_pcent_stau = np.array(vb_pcent_stau)[order]
        ib_pcent_stau = np.array(ib_pcent_stau)[order]
        ob_pcent_stau = np.array(ob_pcent_stau)[order]
        axes[0].plot(masses, vb_pcent_stau, marker="o", linewidth=3, label=window_labels[window], ms=7)
        axes[1].plot(masses, ib_pcent_stau, marker="o", linewidth=3, label=window_labels[window], ms=7)
        axes[2].plot(masses, ob_pcent_stau, marker="o", linewidth=3, label=window_labels[window], ms=7)
        
        axes[0].set_title("$\geq$ 3 VB hits", fontsize=20)
        axes[1].set_title("$\geq$ 3 VB, $\geq$ 2 IB hits", fontsize=20)
        axes[2].set_title("$\geq$ 3 VB, $\geq$ 2 IB, $\geq$ 2 OB hits", fontsize=20)
    for ax in axes: 
        ax.set_xlabel("Stau mass [TeV]", fontsize=20)
        ax.grid(True, alpha=0.2)
        ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
        ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

    axes[0].set_ylabel("Average % of stau hits in stau tracks", fontsize=20)
    
    handles, labels = [], []
    for ax in axes:
        h, l = ax.get_legend_handles_labels()
        handles.extend(h)
        labels.extend(l)

    by_label = OrderedDict(zip(labels,handles))
    fig.legend(by_label.values(), by_label.keys(), ncol=3, fontsize=20, 
                loc="upper center", bbox_to_anchor=(0.5, 0.98), handlelength=1.5,
                handletextpad=0.5, columnspacing=1.0) 

    fig.tight_layout(rect=[0,0,1,0.86]) 
    pdf.savefig(fig)
    plt.close(fig)

    fig3, ax3 = plt.subplots(figsize=(10,7.5))
    x_new = np.arange(len(samples))
    width = 0.22

    base_colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    window_colors = {w: base_colors[i] for i, w in enumerate(windows)}  # loose, medium, tight
    
    det_alphas = {"VB": 0.35, "IB": 0.65, "OB": 0.95}
    det_order  = [("n_vb", "VB"), ("n_ib", "IB"), ("n_ob", "OB")]

    avg_by = {w: {s: {} for s in samples} for w in windows}
    max_tot = 0.0
    for w in windows:
        for s in samples:
            vb = np.mean(track_n_hits[w][s]["n_vb"]) if track_n_hits[w][s]["n_vb"] else 0.0
            ib = np.mean(track_n_hits[w][s]["n_ib"]) if track_n_hits[w][s]["n_ib"] else 0.0
            ob = np.mean(track_n_hits[w][s]["n_ob"]) if track_n_hits[w][s]["n_ob"] else 0.0
            tot = vb + ib + ob
            avg_by[w][s] = {"VB": vb, "IB": ib, "OB": ob, "TOT": tot}
            max_tot = max(max_tot, tot)

    # Plot stacks: same base color per window, different alpha per detector
    for iw, w in enumerate(windows):
        col = window_colors[w]
        xw = x_new + (iw - (len(windows) - 1) / 2) * width
        bottoms = np.zeros_like(x_new, dtype=float)
        for (det_key, det_tag) in det_order:
            heights = np.array([avg_by[w][s][det_tag] for s in samples], dtype=float)
            ax3.bar(
                xw, heights, width=width, bottom=bottoms,
                color=col, alpha=det_alphas[det_tag],
                edgecolor=col, linewidth=1.0
            )
            bottoms += heights

    # Legends: windows use the base colors; detector types use gray with matching alphas
    win_handles = [
        Line2D([0], [0], marker='s', linestyle='', color=window_colors[w],
            markersize=10, label=window_labels[w])
        for w in windows
    ]
    det_handles = [
        Patch(facecolor='gray', alpha=det_alphas["VB"], label="Vertex Barrel"),
        Patch(facecolor='gray', alpha=det_alphas["IB"], label="Inner Barrel"),
        Patch(facecolor='gray', alpha=det_alphas["OB"], label="Outer Barrel"),
    ]

    leg1 = ax3.legend(handles=win_handles, loc="upper center",
                    bbox_to_anchor=(0.5, 1.15), ncol=len(win_handles), frameon=True, fontsize=14,  handlelength=1, handletextpad=0.3, columnspacing=0.5)
    ax3.add_artist(leg1)
    leg2 = ax3.legend(handles=det_handles, loc="upper center",
                    bbox_to_anchor=(0.5, 1.08), ncol=len(det_handles), frameon=True, fontsize=14, handlelength=1, handletextpad=0.3, columnspacing=0.5)
    ax3.add_artist(leg2)

    ax3.set_xticks(x_new)
    ax3.set_xticklabels([sample_to_mass_labels[s] for s in samples], rotation=18, fontsize=16)
    ax3.set_ylabel("Average hits per track", fontsize=25)
    ax3.set_ylim(0, max_tot * 1.15)
    ax3.margins(x=0.02)

    pdf.savefig(fig3)
    plt.close(fig3)




    fig4, ax4 = plt.subplots()
    fig5, ax5 = plt.subplots()
    all_mean_lt = []
    all_mean_ge = []
    all_nmean_lt = []
    all_nmean_ge = []
    for window in windows:
        masses = [1, 2.5, 3, 3.5, 4, 4.5]
        mean_lt50 = []
        mean_ge50 = []
        nmean_lt = []
        nmean_ge = []

        for sample in samples: 
            pcent = np.array(track_n_hits[window][sample]["pcent_stau"])*100
            nhits = np.array(track_n_hits[window][sample]["n_ib"]) + np.array(track_n_hits[window][sample]["n_ob"])
            chi2 = np.array(track_n_hits[window][sample]["chi2"])

            chi2_lt = chi2[pcent<50]
            chi2_ge = chi2[pcent>=50]

            chi2_n_lt = chi2[nhits<2]
            chi2_n_ge = chi2[nhits>=2]

            mean_lt50.append(float(np.mean(chi2_lt)) if chi2_lt.size else 0.0)
            mean_ge50.append(float(np.mean(chi2_ge)) if chi2_ge.size else 0.0)

            nmean_lt.append(float(np.mean(chi2_n_lt)) if chi2_n_lt.size else 0.0)
            nmean_ge.append(float(np.mean(chi2_n_ge)) if chi2_n_ge.size else 0.0)

        ax4.plot(masses, mean_lt50, label=f"{window_labels[window]}, <50% stau", linestyle="--")
        ax4.plot(masses, mean_ge50, label=f"{window_labels[window]}, $\geq$50% stau", linestyle="-")

        ax5.plot(masses, nmean_lt, label=f"{window_labels[window]}, <2 outer hits", linestyle="--")
        ax5.plot(masses, nmean_ge, label=f"{window_labels[window]}, $\geq$2 outer hits", linestyle="-")
        
        print(f"\n{window}")
        print("avg red. $\chi^2$ for tracks that are mostly staus:", np.mean(mean_ge50))
        print("avg red. $\chi^2$ for tracks that are mostly BIB", np.mean(mean_lt50))
        all_mean_lt.append(np.mean(mean_lt50))
        all_mean_ge.append(np.mean(mean_ge50))
        all_nmean_ge.append(np.mean(chi2_n_ge))
        all_nmean_lt.append(np.mean(chi2_n_lt))

    print(f"\n regardless of mass or window")
    print("avg red. $\chi^2$ for tracks that are mostly staus:", np.mean(all_mean_ge))
    print("avg red. $\chi^2$ for tracks that are mostly BIB", np.mean(all_mean_lt))
    
    print(f"\n regardless of mass or window")
    print("avg red. $\chi^2$ for tracks with less than 2 inner/outer hits:", np.mean(all_nmean_lt))
    print("avg red. $\chi^2$ for tracks with more than 2 inner/outer hits:", np.mean(all_nmean_ge))

    ax4.set_xlabel("Stau mass (TeV)", fontsize=14)
    ax4.set_ylabel("Average track reduced $\chi^2$", fontsize=14)
    ax4.legend(ncol=3, fontsize=7)
    pdf.savefig(fig4)
    plt.close(fig4)
    pdf.savefig(fig5)
    plt.close(fig5)



    all_chi2 = np.array([
        val
        for window in windows
        for sample in samples
        for val in track_n_hits[window][sample]["chi2"]
    ], dtype=float)

    edges = np.r_[ 
        np.logspace(-2, 1, 61),          
        np.logspace(1, 4, 31)[1:],       
        np.logspace(4, 6.5, 16)[1:]      
    ]
 
    counts, edges = np.histogram(all_chi2, bins=edges)

    fig, ax = plt.subplots()
    ax.stairs(counts, edges, label=r'$\chi^2$')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'No Refit $\chi^2$', fontsize=20)
    ax.set_ylabel('Count', fontsize=20)
    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

    print(np.percentile(all_chi2, [0,25,50,75,95,99,99.9]))





print(f"saved plots to {pdf_path}")           