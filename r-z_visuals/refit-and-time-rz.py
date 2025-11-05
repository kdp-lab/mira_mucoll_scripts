import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

import pyLCIO
from pyLCIO import UTIL
import ROOT

in_dir   = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/bib_trends/4000_10/si-tracks_test/7-25/"
n_files  = 50
outfile  = "/scratch/miralittmann/analysis/mira_analysis_code/stau_prerefit_rz.pdf"

STAU_COLOR = "tab:green"
BIB_COLOR  = "tab:red"

system_to_det = {
    1: "VXDBarrel", 2: "VXDEndcap",
    3: "ITBarrel",  4: "ITEndcap",
    5: "OTBarrel",  6: "OTEndcap",
}
stau_ids = {1000015, -1000015, 2000015, -2000015}

layer_map = [
    ("VB", 0, 30), ("VB", 1, 30), ("VB", 2, 51), ("VB", 3, 51),
    ("VB", 4, 74), ("VB", 5, 74), ("VB", 6, 102), ("VB", 7, 102),
    ("IB", 0, 127), ("IB", 1, 340), ("IB", 2, 554),
    ("OB", 0, 819), ("OB", 1, 1153), ("OB", 2, 1486),
]

VB_OFFSET = 0.2 
radius_to_idxs = {}
for i,(det,lay,r) in enumerate(layer_map):
    radius_to_idxs.setdefault((det,r), []).append(i)

layer_offsets = []
for i,(det,lay,r) in enumerate(layer_map):
    group = radius_to_idxs[(det,r)]
    if len(group) == 1:
        layer_offsets.append(0.0)
    else:
        k = group.index(i)
        centered = k - (len(group)-1)/2.0
        off = centered * VB_OFFSET
        layer_offsets.append(off if det == "VB" else 0.0)

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

def hit_is_stau(hit, rel_nav):
    encoding = rel_nav["_ENCODING"]
    decoder = UTIL.BitField64(encoding)
    decoder.setValue(int(hit.getCellID0()))
    system = decoder["system"].value()
    det_name = system_to_det.get(system)
    if det_name is None:
        return False
    related = rel_nav[det_name].getRelatedToObjects(hit)
    if not related:
        return False
    for sim in related:
        mcp = sim.getMCParticle()
        if mcp and abs(mcp.getPDG()) in stau_ids:
            return True
    return False

def get_xyz(hit):
    pos = hit.getPosition()
    try:
        return float(pos[0]), float(pos[1]), float(pos[2])
    except Exception:
        p = np.asarray(pos).ravel()
        return float(p[0]), float(p[1]), float(p[2])

ITBARREL_SYS = 3

def decode_system_layer(hit, encoding_str):
    dec = UTIL.BitField64(encoding_str)
    dec.setValue(int(hit.getCellID0()))
    system = dec["system"].value()
    layer  = dec["layer"].value()
    return system, layer

def collect_hits_pre_refit(track, rel_nav):
    zs, rs, ts, mask, sys_arr, lay_arr = [], [], [], [], [], []
    enc = rel_nav["_ENCODING"]

    for hit in track.getTrackerHits():
        x, y, z = get_xyz(hit)
        r = np.hypot(x, y)
        zs.append(z); rs.append(r); ts.append(hit.getTime())
        # stau / non-stau classification via relations (unchanged)
        mask.append(hit_is_stau(hit, rel_nav))
        s, l = decode_system_layer(hit, enc)
        sys_arr.append(s); lay_arr.append(l)

    if not zs:
        emptyB = np.array([], dtype=bool)
        return np.array([]), np.array([]), np.array([]), emptyB, np.array([]), np.array([])

    return (np.asarray(zs), np.asarray(rs), np.asarray(ts),
            np.asarray(mask, dtype=bool), np.asarray(sys_arr), np.asarray(lay_arr))

def annotate_times_with_crowding(ax, x, y, t, drawn_px=None, min_sep_px=12, fontsize=15):
    if drawn_px is None:
        drawn_px = []
    to_disp = ax.transData.transform
    for xi, yi, tt in zip(x, y, t):
        px, py = to_disp((xi, yi))
        if all((px - dx)**2 + (py - dy)**2 > (min_sep_px**2) for dx, dy in drawn_px):
            ax.annotate(f"{tt:.2f} ns", (xi, yi), xytext=(3, 3),
                        textcoords="offset points", fontsize=fontsize)
            drawn_px.append((px, py))
    return drawn_px

VB_radii = [r for det,lay,r in layer_map if det == "VB"]
IB_rads  = [(lay, r) for det,lay,r in layer_map if det == "IB"]
OB_rads  = [(lay, r) for det,lay,r in layer_map if det == "OB"]

def plot_rz_preref(ax, z, r, t, is_stau, sys_arr, lay_arr):
    if r.size == 0:
        ax.set_axis_off()
        ax.set_title("No hits", fontsize=18)
        return

    ax.scatter(z[ is_stau], r[ is_stau], s=70, color=STAU_COLOR, label="Stau hits")
    ax.scatter(z[~is_stau], r[~is_stau], s=70, color=BIB_COLOR,  label="BIB hits")

    for _, _, rr in layer_map:
        ax.axhline(rr, ls="--", lw=0.9, color="lightgray", zorder=0)

    x_right = ax.get_xlim()[1] if ax.get_xlim()[0] < ax.get_xlim()[1] else ax.get_xlim()[0]
    for lay, rr in IB_rads:
        ax.text(x_right, rr, f"IB L{lay}", va="bottom", ha="right", fontsize=10, color="dimgray", fontweight="bold")
    for lay, rr in OB_rads:
        ax.text(x_right, rr, f"OB L{lay}", va="bottom", ha="right", fontsize=10, color="dimgray", fontweight="bold")
    if VB_radii:
        ax.text(x_right, float(np.median(VB_radii)), "VB", va="center", ha="right",
                fontsize=11, color="dimgray", fontweight="bold")

    drawn = []
    bib_mask = ~is_stau
    drawn = annotate_times_with_crowding(ax, z[bib_mask], r[bib_mask], t[bib_mask],
                                         drawn_px=drawn, min_sep_px=12, fontsize=17)
    ib0_stau = (sys_arr == ITBARREL_SYS) & (lay_arr == 0) & (is_stau == True)
    drawn = annotate_times_with_crowding(ax, z[ib0_stau], r[ib0_stau], t[ib0_stau],
                                         drawn_px=drawn, min_sep_px=12, fontsize=17)

    ax.set_xlabel("z [mm]", fontsize=20)
    ax.set_ylabel("r [mm]", fontsize=20)
    ax.tick_params(labelsize=17)
    ax.legend(loc="best",fontsize=11,frameon=True, facecolor="white", edgecolor="black", fancybox=True)
    
reader = pyLCIO.IOIMPL.LCFactory.getInstance().createLCReader()

with PdfPages(outfile) as pdf:
    for i in range(n_files):
        fname = os.path.join(in_dir, f"4000_10_reco{i}.slcio")
        if not os.path.exists(fname):
            continue

        reader.open(fname)
        for event in reader:
            rel_nav = build_rel_nav(event)

            rel_collection   = event.getCollection("MCParticle_SiTracks")
            track_collection = event.getCollection("SiTracks")

            for idx in range(rel_collection.getNumberOfElements()):
                rel   = rel_collection.getElementAt(idx)
                mcp   = rel.getFrom()
                track = rel.getTo()
                if not mcp or abs(mcp.getPDG()) not in stau_ids:
                    continue

                fig, ax = plt.subplots(figsize=(8.5, 6.5))
                z, r, t, is_stau, sys_arr, lay_arr = collect_hits_pre_refit(track, rel_nav)
                plot_rz_preref(ax, z, r, t, is_stau, sys_arr, lay_arr)

                fig.tight_layout()
                pdf.savefig(fig)
                plt.close(fig)

        reader.close()

print(f"Saved plots to: {outfile}")
