import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from collections import OrderedDict
from scipy import optimize
import math

import pyLCIO
from pyLCIO import UTIL
import ROOT

from tqdm import tqdm
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--rebuild", action="store_true")
args = parser.parse_args()

# plot of loose window, most stringent hit req, only 4tev staus, 
# track purity (% of hits in "stau" track actually caused by stau).
# x axis increasing bib percentage, y axis track purity

bib_percent = np.arange(0, 110, 10)
bibdir = "/ospool/uc-shared/project/futurecolliders/wandriscok/reco/stau/loose/bib/"

pdf_path = "/scratch/miralittmann/analysis/mira_analysis_code/n-hits-study/bib_purity.pdf"

n_files = 2500
stau_ids = {1000015, -1000015, 2000015, -2000015}

system_to_det = {
    1: "VB", 2: "VE",
    3: "IB",  4: "IE",
    5: "OB",  6: "OE",
}

system_to_relname = {
    1: "VXDBarrel", 2: "VXDEndcap",
    3: "ITBarrel",  4: "ITEndcap",
    5: "OTBarrel", 6: "OTEndcap",
}

rebuild = args.rebuild
cache_path = "/scratch/miralittmann/analysis/mira_analysis_code/n-hits-study/track_n_hits_cache_bibpurity.pkl"

# Track-level selection cuts
ETA_MAX = 0.8
CHI2NDF_MAX = 3.0
PT_MAX = 10000.0
WRMS_MAX = 1.4

Bfield = 3.57
speedoflight = 299792458 / 1000000  # mm/ns
guess_velo = 180


def build_rel_nav(event, coned=True):
    if coned == False:
        nav = {
            "VXDBarrel": UTIL.LCRelationNavigator(event.getCollection("VXDBarrelHitsRelations")),
            "VXDEndcap": UTIL.LCRelationNavigator(event.getCollection("VXDEndcapHitsRelations")),
            "ITBarrel": UTIL.LCRelationNavigator(event.getCollection("ITBarrelHitsRelations")),
            "ITEndcap": UTIL.LCRelationNavigator(event.getCollection("ITEndcapHitsRelations")),
            "OTBarrel": UTIL.LCRelationNavigator(event.getCollection("OTBarrelHitsRelations")),
            "OTEndcap": UTIL.LCRelationNavigator(event.getCollection("OTEndcapHitsRelations")),
        }
        enc = event.getCollection("ITBarrelHits").getParameters().getStringVal(
            pyLCIO.EVENT.LCIO.CellIDEncoding
        )
    else:
        nav = {
            "VXDBarrel": UTIL.LCRelationNavigator(event.getCollection("VXDBarrelHitsRelationsConed")),
            "VXDEndcap": UTIL.LCRelationNavigator(event.getCollection("VXDEndcapHitsRelationsConed")),
            "ITBarrel": UTIL.LCRelationNavigator(event.getCollection("ITBarrelHitsRelationsConed")),
            "ITEndcap": UTIL.LCRelationNavigator(event.getCollection("ITEndcapHitsRelationsConed")),
            "OTBarrel": UTIL.LCRelationNavigator(event.getCollection("OTBarrelHitsRelationsConed")),
            "OTEndcap": UTIL.LCRelationNavigator(event.getCollection("OTEndcapHitsRelationsConed")),
        }
        enc = event.getCollection("ITBarrelHitsConed").getParameters().getStringVal(
            pyLCIO.EVENT.LCIO.CellIDEncoding
        )

    nav["_ENCODING"] = enc
    nav["_DECODER"] = UTIL.BitField64(enc)
    return nav


def decode_system_layer(hit, decoder):
    decoder.setValue(int(hit.getCellID0()))
    system = decoder["system"].value()
    layer = decoder["layer"].value()
    return system, layer


def time_rms_from_fit(v, t, r, time_unc, b=0.0):
    t = np.asarray(t, float)
    r = np.asarray(r, float)
    st = np.asarray(time_unc, float)

    m = np.isfinite(t) & np.isfinite(r) & np.isfinite(st) & (st > 0)
    t, r, st = t[m], r[m], st[m]

    if t.size < 3 or (not np.isfinite(v)) or abs(v) < 1e-12:
        return np.nan, np.nan

    t_pred = (r - b) / v
    dt = t - t_pred

    rms_t = float(np.sqrt(np.mean(dt * dt)))
    rms_pull = float(np.sqrt(np.mean((dt / st) ** 2)))

    return rms_t, rms_pull


def reco_velo_no_intercept(times, pos, spatial_unc, time_unc):
    x = np.asarray(times, dtype=float)
    y = np.asarray(pos, dtype=float)
    s = np.asarray(spatial_unc, dtype=float)
    st = np.asarray(time_unc, dtype=float)

    m = (
        np.isfinite(x) & np.isfinite(y) &
        np.isfinite(s) & (s > 0) &
        np.isfinite(st) & (st > 0)
    )
    x, y, s, st = x[m], y[m], s[m], st[m]

    if x.size < 3 or np.allclose(x, x.mean()):
        return np.nan, np.nan, np.nan, np.nan

    v0 = np.array([guess_velo])

    def residual0(v, times, pos, spatial_unc):
        vv = float(np.atleast_1d(v)[0])
        return (vv * times - pos) / spatial_unc

    fit = optimize.least_squares(
        residual0,
        v0,
        args=(x, y, s),
        jac="2-point"
    )

    v = float(fit.x[0])

    rms_t, rms_pull = time_rms_from_fit(v, x, y, st, b=0.0)

    try:
        r = v * x - y
        J = fit.jac
        dof = max(1, x.size - 1)
        chi2 = np.sum((r / s) ** 2)
        sigma2 = chi2 / dof
        cov = np.linalg.inv(J.T @ J) * sigma2
        v_err = float(np.sqrt(cov[0, 0]))
    except Exception:
        v_err = np.nan

    return v, v_err, rms_t, rms_pull


def pass_track_level(t):
    return (
        (t["vb_hits"] >= 3) and
        (t["ib_hits"] >= 2) and
        (t["ob_hits"] >= 2) and
        np.isfinite(t["eta"]) and (abs(t["eta"]) < ETA_MAX) and
        np.isfinite(t["chi2ndf"]) and (t["chi2ndf"] < CHI2NDF_MAX) and
        np.isfinite(t["pT"]) and (t["pT"] < PT_MAX) and
        np.isfinite(t["beta"]) and (0 < t["beta"] < 1) and
        np.isfinite(t["mass"]) and
        np.isfinite(t["vrmsw"]) and (t["vrmsw"] < WRMS_MAX)
    )


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

    # Check whether this is an old cache without the track-level quantities
    if not all(
        "vrmsw" in track_n_hits[percent]
        for percent in bib_percent
    ):
        print("Old cache does not contain track-level quantities. Rebuilding...")
        track_n_hits = None


if track_n_hits is None:
    print("rebuilding analysis arrays...")

    track_n_hits = {
        percent: {
            "n_vb": [],
            "n_ib": [],
            "n_ob": [],
            "n_tot": [],
            "pcent_stau": [],
            "chi2": [],
            "eta": [],
            "pT": [],
            "beta": [],
            "mass": [],
            "vrmsw": [],
        }
        for percent in bib_percent
    }

    for percent in bib_percent:
        for i in tqdm(range(n_files)):

            if percent == 0:
                nobibdir = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/loose4/seeding_10GeV/nobib/"
                fname = os.path.join(
                    nobibdir,
                    f"4000_30/4000_30_reco{i}.slcio"
                )
                coned_yn = False
            else:
                fname = os.path.join(
                    bibdir,
                    f"{percent}/4000_30/4000_30_reco{i}.slcio"
                )
                coned_yn = True

            if not os.path.exists(fname):
                print(f"cannot find {fname}")
                continue

            reader.open(fname)

            for event in reader:
                rel_nav = build_rel_nav(event, coned=coned_yn)

                all_collections = event.getCollectionNames()
                rel_collection = (
                    event.getCollection("MCParticle_SiTracks")
                    if "MCParticle_SiTracks" in all_collections
                    else None
                )

                if rel_collection is None:
                    continue

                for idx in range(rel_collection.getNumberOfElements()):
                    rel = rel_collection.getElementAt(idx)
                    mcp = rel.getFrom()
                    track = rel.getTo()

                    # Only looking at stau tracks
                    if not mcp or abs(mcp.getPDG()) not in stau_ids:
                        continue

                    n_vb_hits = 0
                    n_ib_hits = 0
                    n_ob_hits = 0

                    stau_hit_count = 0
                    hit_count = 0

                    track_hits = track.getTrackerHits()

                    track_times = []
                    track_pos = []
                    spatial_unc = []
                    time_unc = []

                    for hit in track_hits:
                        system, layer = decode_system_layer(
                            hit, rel_nav["_DECODER"]
                        )
                        det_name = system_to_det.get(system)

                        if det_name == "VB":
                            n_vb_hits += 1
                            spatial_unc.append(0.005)
                            time_unc.append(0.03)

                        elif det_name == "IB":
                            n_ib_hits += 1
                            spatial_unc.append(0.007)
                            time_unc.append(0.06)

                        elif det_name == "OB":
                            n_ob_hits += 1
                            spatial_unc.append(0.007)
                            time_unc.append(0.06)

                        elif det_name == "VE":
                            spatial_unc.append(0.005)
                            time_unc.append(0.03)

                        elif det_name in ("IE", "OE"):
                            spatial_unc.append(0.007)
                            time_unc.append(0.06)

                        for sim in rel_nav[system_to_relname[system]].getRelatedToObjects(hit):
                            sim_mcp = sim.getMCParticle()
                            if sim_mcp and abs(sim_mcp.getPDG()) in stau_ids:
                                stau_hit_count += 1

                        hit_count += 1

                        x = hit.getPosition()[0]
                        y = hit.getPosition()[1]
                        z = hit.getPosition()[2]

                        hit_pos = np.sqrt(x**2 + y**2 + z**2)
                        tof = hit_pos / speedoflight

                        hit_time = hit.getTime()
                        corrected_t = hit_time - tof
                        corrected_corrected_t = 2 * hit_time - corrected_t

                        track_times.append(corrected_corrected_t)
                        track_pos.append(hit_pos)

                    total_track_hits = n_vb_hits + n_ib_hits + n_ob_hits

                    track_pcent_stau = (
                        stau_hit_count / len(track_hits)
                        if len(track_hits) else np.nan
                    )

                    red_chi2 = (
                        track.getChi2() / track.getNdf()
                        if track.getNdf() > 0 else np.nan
                    )

                    v_fit, v_err, rms_t, rms_w = reco_velo_no_intercept(
                        track_times,
                        track_pos,
                        spatial_unc,
                        time_unc
                    )

                    beta = (
                        v_fit / speedoflight
                        if np.isfinite(v_fit) else np.nan
                    )

                    try:
                        reco_pT = 0.3 * Bfield / abs(track.getOmega() * 1000.)
                    except Exception:
                        reco_pT = np.nan

                    try:
                        tanL = track.getTanLambda()
                    except Exception:
                        tanL = np.nan

                    if np.isfinite(tanL):
                        eta = float(np.arcsinh(tanL))
                        reco_pz = reco_pT * tanL
                        reco_p = math.sqrt(reco_pT**2 + reco_pz**2)
                    else:
                        eta = np.nan
                        reco_p = np.nan

                    if np.isfinite(reco_p) and np.isfinite(beta) and (0 < beta < 1):
                        mass = reco_p * math.sqrt(1.0 / (beta * beta) - 1.0)
                    else:
                        mass = np.nan

                    track_n_hits[percent]["n_tot"].append(total_track_hits)
                    track_n_hits[percent]["n_vb"].append(n_vb_hits)
                    track_n_hits[percent]["n_ib"].append(n_ib_hits)
                    track_n_hits[percent]["n_ob"].append(n_ob_hits)
                    track_n_hits[percent]["pcent_stau"].append(track_pcent_stau)
                    track_n_hits[percent]["chi2"].append(red_chi2)

                    track_n_hits[percent]["eta"].append(eta)
                    track_n_hits[percent]["pT"].append(reco_pT)
                    track_n_hits[percent]["beta"].append(beta)
                    track_n_hits[percent]["mass"].append(mass)
                    track_n_hits[percent]["vrmsw"].append(rms_w)

            reader.close()

    os.makedirs(os.path.dirname(cache_path), exist_ok=True)

    with open(cache_path, "wb") as f:
        print("saving to cache...")
        pickle.dump(track_n_hits, f, protocol=pickle.HIGHEST_PROTOCOL)
        print("cache saved.")


markers = {"VB": "v", "IB": "o", "OB": "^"}

with PdfPages(pdf_path) as pdf:
    plt.style.use("seaborn-v0_8-colorblind")

    dets = ["VB", "IB", "OB", "Total"]
    x = np.arange(len(dets))
    width = 0.24

    fig, ax = plt.subplots(figsize=(8, 6))

    purities = []
    purities_track_level = []

    for percent in bib_percent:
        vb = np.asarray(track_n_hits[percent]["n_vb"], dtype=float)
        ib = np.asarray(track_n_hits[percent]["n_ib"], dtype=float)
        ob = np.asarray(track_n_hits[percent]["n_ob"], dtype=float)
        chi2 = np.asarray(track_n_hits[percent]["chi2"], dtype=float)

        pcent = (
            np.asarray(track_n_hits[percent]["pcent_stau"], dtype=float) * 100
        )

        eta = np.asarray(track_n_hits[percent]["eta"], dtype=float)
        pT = np.asarray(track_n_hits[percent]["pT"], dtype=float)
        beta = np.asarray(track_n_hits[percent]["beta"], dtype=float)
        mass = np.asarray(track_n_hits[percent]["mass"], dtype=float)
        vrmsw = np.asarray(track_n_hits[percent]["vrmsw"], dtype=float)

        ob_mask = (vb >= 1.5) & (ib >= 2) & (ob >= 2)

        ob_values = pcent[ob_mask]
        purities.append(
            np.mean(ob_values) if len(ob_values) > 0 else np.nan
        )

        track_level_mask = (
            (vb >= 3) &
            (ib >= 2) &
            (ob >= 2) &
            np.isfinite(eta) & (np.abs(eta) < ETA_MAX) &
            np.isfinite(chi2) & (chi2 < CHI2NDF_MAX) &
            np.isfinite(pT) & (pT < PT_MAX) &
            np.isfinite(beta) & (beta > 0) & (beta < 1) &
            np.isfinite(mass) &
            np.isfinite(vrmsw) & (vrmsw < WRMS_MAX)
        )

        track_level_values = pcent[track_level_mask]

        purities_track_level.append(
            np.mean(track_level_values)
            if len(track_level_values) > 0 else np.nan
        )

    ax.plot(
        bib_percent,
        purities,
        marker="o",
        linewidth=3,
        ms=7,
        label="Hit requirements only"
    )

    ax.plot(
        bib_percent,
        purities_track_level,
        marker="s",
        linewidth=3,
        ms=7,
        label="With track-level cuts"
    )

    ax.set_xlabel("Percent BIB inclusion in simulation", fontsize=20)
    ax.set_ylim(94, 105)
    ax.grid(True, alpha=0.2)
    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

    ax.set_ylabel("Average % of stau hits in stau tracks", fontsize=20)
    ax.legend(fontsize=13, loc="upper left")

    ax.text(
        0.98, 0.98,
        "Muon Collider",
        ha="right", va="top",
        transform=ax.transAxes,
        fontsize=24,
        fontweight="bold",
        style="italic",
    )

    ax.text(
        0.98, 0.91,
        "Simulation 10% BIB, $\sqrt{s}$=10 TeV",
        ha="right", va="top",
        transform=ax.transAxes,
        fontsize=18
    )

    ax.text(
        0.98, 0.85,
        "4 TeV Staus, Loose window",
        ha="right", va="top",
        transform=ax.transAxes,
        fontsize=18
    )

    ax.text(
        0.98, 0.79,
        "MuColl_v1",
        ha="right", va="top",
        transform=ax.transAxes,
        fontsize=18
    )

    fig.tight_layout(rect=[0, 0, 1, 0.86])
    pdf.savefig(fig)
    plt.close(fig)


print(f"saved plots to {pdf_path}")