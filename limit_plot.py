import os
import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.colors as colors
from scipy.interpolate import RBFInterpolator

# ============================================================
# INPUTS / CONFIG
# ============================================================

stau_cache_loose = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_with_track_level_loose.pkl"
stau_cache_nominal = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_with_track_level_nominal.pkl"
stau_cache_medium = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_with_track_level_medium.pkl"

stau_cache_01_loose = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_with_track_level_loose_01only.pkl"
stau_cache_01_medium = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_with_track_level_medium_01only.pkl"
stau_cache_01_nominal = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_with_track_level_nominal_01only.pkl"

stau_cache_03_1_loose = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_with_track_level_loose_03_1.pkl"
stau_cache_03_1_medium = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_with_track_level_medium_03_1.pkl"
stau_cache_03_1_nominal = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_with_track_level_nominal_03_1.pkl"

stau_base_dir_loose = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/loose4/seeding_10GeV/nobib"
stau_base_dir_nominal = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/nominal/seeding_10GeV/nobib"
stau_base_dir_medium = "/ospool/uc-shared/project/futurecolliders/miralittmann/reco/reder_timing/p26_p5_p6/seeding_10GeV/nobib"

bib_cache_loose = "/scratch/wandriscok/kate_mucoll_scripts/bib_analysis/cache/wrms_100bib_event_plot_leadsub_loose.pkl"
bib_cache_medium = "/scratch/wandriscok/kate_mucoll_scripts/bib_analysis/cache/wrms_100bib_event_plot_leadsub_medium.pkl"
bib_cache_nominal = "/scratch/wandriscok/kate_mucoll_scripts/bib_analysis/cache/wrms_100bib_event_plot_leadsub_tight.pkl"

muon_cache = "/scratch/miralittmann/analysis/mira_analysis_code/cache/mumu_bkg_stats_nominal_nobib_byevent.pkl"

configs = ["nominal", "medium", "loose"]
signal_masses = ["1000", "1500", "2000", "2500", "3000", "3500", "4000", "4500"]
extra_lifetimes_ns = []

SIGNAL_CUTOFF = 3.0
l_abinv = 10.0
l_pb_inv = l_abinv * 1e6
w_rms_cut = 1.6
xs_pb = {
    "muon_bkg": 0.4312,
    "bib_bkg": 150000,
    "1000": 0.0005108,
    "1500": 0.0004715,
    "2000": 0.0004184,
    "2500": 0.0003528,
    "3000": 0.0002780,
    "3500": 0.0001980,
    "4000": 0.0001173,
    "4500": 0.0000450,
}

event_cuts = {
    "1000": {"pT": {"lead": 1500, "sublead": 1500}, "mass": {"lead": (500, 2000), "sublead": (500, 2000)}, "beta": {"lead": (0.95, 0.983), "sublead": (0.95, 0.983)}},
    "1500": {"pT": {"lead": 2000, "sublead": 2000}, "mass": {"lead": (1000, 2000), "sublead": (950, 2000)}, "beta": {"lead": (0.90, 0.98), "sublead": (0.90, 0.98)}},
    "2000": {"pT": {"lead": 2000, "sublead": 2000}, "mass": {"lead": (1500, 3000), "sublead": (1500, 2500)}, "beta": {"lead": (0.85, 0.95), "sublead": (0.85, 0.95)}},
    "2500": {"pT": {"lead": 2000, "sublead": 2000}, "mass": {"lead": (2000, 4000), "sublead": (1900, 3100)}, "beta": {"lead": (0.80, 0.90), "sublead": (0.80, 0.90)}},
    "3000": {"pT": {"lead": 2000, "sublead": 2000}, "mass": {"lead": (2000, 6000), "sublead": (2000, 4500)}, "beta": {"lead": (0.70, 0.85), "sublead": (0.70, 0.85)}},
    "3500": {"pT": {"lead": 2000, "sublead": 2000}, "mass": {"lead": (2900, 5000), "sublead": (2900, 4200)}, "beta": {"lead": (0.65, 0.80), "sublead": (0.65, 0.80)}},
    "4000": {"pT": {"lead": 1900, "sublead": 1900}, "mass": {"lead": (3200, 5500), "sublead": (3200, 4500)}, "beta": {"lead": (0.55, 0.70), "sublead": (0.55, 0.70)}},
    "4500": {"pT": {"lead": 1400, "sublead": 1400}, "mass": {"lead": (3500, 6000), "sublead": (3500, 5500)}, "beta": {"lead": (0.40, 0.55), "sublead": (0.40, 0.55)}},
}

plt.style.use("seaborn-v0_8-colorblind")

# ============================================================
# HELPERS
# ============================================================

def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)

def lifetime_to_dirlabel(lifetime):
    lf = float(lifetime)
    if np.isclose(lf, 0.1):
        return "01"
    if np.isclose(lf, 0.3):
        return "03"
    if lf.is_integer():
        return str(int(lf))
    return str(lf)

def count_files_in_dir(base_dir, mass, lifetime):
    tau_label = lifetime_to_dirlabel(lifetime)
    path = os.path.join(base_dir, f"{mass}_{tau_label}")
    try:
        return len([f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))])
    except FileNotFoundError:
        print(f"WARNING: directory not found: {path}")
        return 0

def passing_cuts(cache, cut_set, bib=False, muons=False):
    lead_beta = np.asarray(cache["leading_beta"])
    sub_beta  = np.asarray(cache["subleading_beta"])
    lead_pT   = np.asarray(cache["leading_pT"])
    sub_pT    = np.asarray(cache["subleading_pT"])
    lead_mass = np.asarray(cache["leading_mass"])
    sub_mass  = np.asarray(cache["subleading_mass"])

    if bib:
        lead_wrms = np.asarray(cache["leading_w_rms"])
        sub_wrms  = np.asarray(cache["subleading_w_rms"])
    elif muons:
        lead_wrms = np.asarray(cache["leading_wvrms"])
        sub_wrms  = np.asarray(cache["subleading_wvrms"])
    else:
        lead_wrms = np.asarray(cache["leading_vrmsw"])
        sub_wrms  = np.asarray(cache["subleading_vrmsw"])

    lead_beta_lo, lead_beta_hi = event_cuts[cut_set]["beta"]["lead"]
    sub_beta_lo,  sub_beta_hi  = event_cuts[cut_set]["beta"]["sublead"]

    lead_mass_lo, lead_mass_hi = event_cuts[cut_set]["mass"]["lead"]
    sub_mass_lo,  sub_mass_hi  = event_cuts[cut_set]["mass"]["sublead"]

    lead_pT_min = event_cuts[cut_set]["pT"]["lead"]
    sub_pT_min  = event_cuts[cut_set]["pT"]["sublead"]

    num_passing = 0
    num_total = len(lead_pT)

    for i in range(num_total):
        passing_lead_pT   = np.isfinite(lead_pT[i])   and (lead_pT[i] > lead_pT_min)
        passing_lead_beta = np.isfinite(lead_beta[i]) and (lead_beta_lo < lead_beta[i] < lead_beta_hi)
        passing_lead_mass = np.isfinite(lead_mass[i]) and (lead_mass_lo < lead_mass[i] < lead_mass_hi)
        passing_lead_wrms = np.isfinite(lead_wrms[i]) and (lead_wrms[i] < w_rms_cut)

        if not (passing_lead_pT and passing_lead_beta and passing_lead_wrms):
            continue

        passing_sub_pT   = np.isfinite(sub_pT[i])   and (sub_pT[i] > sub_pT_min)
        passing_sub_beta = np.isfinite(sub_beta[i]) and (sub_beta_lo < sub_beta[i] < sub_beta_hi)
        passing_sub_mass = np.isfinite(sub_mass[i]) and (sub_mass_lo < sub_mass[i] < sub_mass_hi)
        passing_sub_wrms = np.isfinite(sub_wrms[i]) and (sub_wrms[i] < w_rms_cut)

        if passing_sub_pT and passing_sub_beta and passing_sub_wrms:
            num_passing += 1

    return num_passing, num_total

def get_loose_cutflow_counts(cache, cut_set):
    lead_pT   = np.asarray(cache["leading_pT"])
    sub_pT    = np.asarray(cache["subleading_pT"])
    lead_beta = np.asarray(cache["leading_beta"])
    sub_beta  = np.asarray(cache["subleading_beta"])

    lead_pT_min = event_cuts[cut_set]["pT"]["lead"]
    sub_pT_min  = event_cuts[cut_set]["pT"]["sublead"]
    lead_beta_lo, lead_beta_hi = event_cuts[cut_set]["beta"]["lead"]
    sub_beta_lo,  sub_beta_hi  = event_cuts[cut_set]["beta"]["sublead"]

    has_two_tracks = np.isfinite(lead_pT) & np.isfinite(sub_pT)
    pass_pt = has_two_tracks & (lead_pT > lead_pT_min) & (sub_pT > sub_pT_min)
    pass_pt_beta = (
        pass_pt
        & np.isfinite(lead_beta) & np.isfinite(sub_beta)
        & (lead_beta_lo < lead_beta) & (lead_beta < lead_beta_hi)
        & (sub_beta_lo  < sub_beta)  & (sub_beta  < sub_beta_hi)
    )

    return {
        "n_two_tracks": int(np.count_nonzero(has_two_tracks)),
        "n_pt": int(np.count_nonzero(pass_pt)),
        "n_final": int(np.count_nonzero(pass_pt_beta)),
    }

def stage_fractions_from_counts(counts, total):
    n_two = counts["n_two_tracks"]
    n_pt = counts["n_pt"]
    n_final = counts["n_final"]
    return {
        "fail_before_two": (total - n_two) / total if total > 0 else 0.0,
        "fail_pt_after_two": (n_two - n_pt) / total if total > 0 else 0.0,
        "fail_beta_after_pt": (n_pt - n_final) / total if total > 0 else 0.0,
        "pass_final": n_final / total if total > 0 else 0.0,
    }

def interp_extrap_1d(x_known, y_known, x_new):
    x = np.asarray(x_known, dtype=float)
    y = np.asarray(y_known, dtype=float)
    order = np.argsort(x)
    x = x[order]
    y = y[order]

    if np.any(np.isclose(x, x_new)):
        return float(y[np.argmin(np.abs(x - x_new))])

    if x[0] < x_new < x[-1]:
        return float(np.interp(x_new, x, y))

    if x_new <= x[0]:
        x0, x1 = x[0], x[1]
        y0, y1 = y[0], y[1]
    else:
        x0, x1 = x[-2], x[-1]
        y0, y1 = y[-2], y[-1]

    slope = (y1 - y0) / (x1 - x0)
    return float(y0 + slope * (x_new - x0))

def interp_extrap_vs_lifetime(tau_known, y_known, tau_new):
    return interp_extrap_1d(np.log(np.asarray(tau_known, dtype=float)), y_known, np.log(float(tau_new)))

def find_crossing_mass(masses, yields, cutoff):
    masses = np.asarray(masses, dtype=float)
    yields = np.asarray(yields, dtype=float)

    for i in range(len(masses) - 1):
        m1, m2 = masses[i], masses[i + 1]
        y1, y2 = yields[i], yields[i + 1]

        if y1 == cutoff:
            return m1

        if (y1 - cutoff) * (y2 - cutoff) < 0:
            t = (cutoff - y1) / (y2 - y1)
            return m1 + t * (m2 - m1)

    return np.nan

def get_contour_xy(contour_obj):
    if len(contour_obj.collections) == 0:
        return None, None
    paths = contour_obj.collections[0].get_paths()
    if len(paths) == 0:
        return None, None
    v = paths[0].vertices
    return v[:, 0], v[:, 1]

def centers_to_edges_linear(x):
    x = np.asarray(x, dtype=float)
    edges = np.empty(len(x) + 1)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - 0.5 * (x[1] - x[0])
    edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
    return edges

def centers_to_edges_log(y):
    y = np.asarray(y, dtype=float)
    logy = np.log10(y)
    edges_log = np.empty(len(y) + 1)
    edges_log[1:-1] = 0.5 * (logy[:-1] + logy[1:])
    edges_log[0] = logy[0] - 0.5 * (logy[1] - logy[0])
    edges_log[-1] = logy[-1] + 0.5 * (logy[-1] - logy[-2])
    return 10**edges_log

# ============================================================
# LOAD DATA
# ============================================================

signal_loose = load_pickle(stau_cache_loose)
signal_nominal = load_pickle(stau_cache_nominal)
signal_medium = load_pickle(stau_cache_medium)

signal_01_loose = load_pickle(stau_cache_01_loose)
signal_01_medium = load_pickle(stau_cache_01_medium)
signal_01_nominal = load_pickle(stau_cache_01_nominal)

signal_03_1_loose = load_pickle(stau_cache_03_1_loose)
signal_03_1_medium = load_pickle(stau_cache_03_1_medium)
signal_03_1_nominal = load_pickle(stau_cache_03_1_nominal)

bib_loose = load_pickle(bib_cache_loose)
bib_medium = load_pickle(bib_cache_medium)
bib_nominal = load_pickle(bib_cache_nominal)
muons = load_pickle(muon_cache)

bib_slice_loose = bib_loose["loose"]["ob"]["bib"]
bib_slice_medium = bib_medium["medium"]["ob"]["bib"]
bib_slice_nominal = bib_nominal["tight"]["ob"]["bib"]

signal_loose["0.1"] = signal_01_loose["01"]
signal_medium["0.1"] = signal_01_medium["01"]
signal_nominal["0.1"] = signal_01_nominal["01"]

signal_loose["0.3"] = signal_03_1_loose["03"]
signal_medium["0.3"] = signal_03_1_medium["03"]
signal_nominal["0.3"] = signal_03_1_nominal["03"]

signal_loose["1"] = signal_03_1_loose["1"]
signal_medium["1"] = signal_03_1_medium["1"]
signal_nominal["1"] = signal_03_1_nominal["1"]

print("Loose lifetimes:  ", sorted(signal_loose.keys(), key=float))
print("Medium lifetimes: ", sorted(signal_medium.keys(), key=float))
print("Nominal lifetimes:", sorted(signal_nominal.keys(), key=float))

assert set(signal_loose.keys()) == set(signal_medium.keys()) == set(signal_nominal.keys())

signal_map = {
    "loose": signal_loose,
    "nominal": signal_nominal,
    "medium": signal_medium,
}

denominator_dirs = {
    "loose": stau_base_dir_loose,
    "nominal": stau_base_dir_nominal,
    "medium": stau_base_dir_medium,
}

# ============================================================
# DENOMINATORS
# ============================================================

lifetimes_known = sorted([float(x) for x in signal_loose.keys()])
lifetimes_known_str = [str(int(x)) if float(x).is_integer() else str(x) for x in lifetimes_known]

denominator_map = {cfg: {} for cfg in configs}
for cfg in configs:
    for tau_str in signal_loose.keys():
        denominator_map[cfg][tau_str] = {}
        for m in signal_masses:
            denominator_map[cfg][tau_str][m] = count_files_in_dir(denominator_dirs[cfg], m, tau_str)

print(f"Assuming L = {l_abinv} ab^-1 = {l_pb_inv:.3e} pb^-1")

# ============================================================
# STACKED CUTFLOW PLOTS (KNOWN LIFETIMES ONLY)
# ============================================================

n_bib_tot = len(np.asarray(bib_slice_loose["leading_pT"]))
n_mu_tot = len(np.asarray(muons["leading_pT"]))
masses_numeric = np.array([int(m) for m in signal_masses])

with PdfPages("loose_window_stacked_cutflow.pdf") as pdf:
    # signal pages
    for tau_str in sorted(signal_loose.keys(), key=float):
        fail_before_two, fail_pt_after_two, fail_beta_after_pt, pass_final = [], [], [], []

        for m in signal_masses:
            total = denominator_map["loose"][tau_str][m]
            counts = get_loose_cutflow_counts(signal_loose[tau_str][m], m)
            fracs = stage_fractions_from_counts(counts, total)

            fail_before_two.append(100 * fracs["fail_before_two"])
            fail_pt_after_two.append(100 * fracs["fail_pt_after_two"])
            fail_beta_after_pt.append(100 * fracs["fail_beta_after_pt"])
            pass_final.append(100 * fracs["pass_final"])

        plt.figure(figsize=(10, 6))
        plt.bar(masses_numeric, fail_before_two, width=300, label="Fail before 2 tracks")
        plt.bar(masses_numeric, fail_pt_after_two, width=300, bottom=fail_before_two, label="Fail pT after 2 tracks")
        plt.bar(masses_numeric, fail_beta_after_pt, width=300,
                bottom=np.array(fail_before_two) + np.array(fail_pt_after_two),
                label="Fail beta after pT")
        plt.bar(masses_numeric, pass_final, width=300,
                bottom=np.array(fail_before_two) + np.array(fail_pt_after_two) + np.array(fail_beta_after_pt),
                label="Pass final")
        plt.xlabel("Stau mass [GeV]", fontsize=18)
        plt.ylabel("Percent of total events", fontsize=18)
        plt.title(f"Loose-window stacked cutflow: signal ({tau_str} ns)", fontsize=15)
        plt.ylim(0, 100)
        plt.xticks(masses_numeric, [str(m) for m in masses_numeric], fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=11)
        plt.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    # BIB page
    x = np.arange(len(signal_masses))
    x_labels = signal_masses

    bib_fail_before_two, bib_fail_pt_after_two, bib_fail_beta_after_pt, bib_pass_final = [], [], [], []
    for m in signal_masses:
        counts = get_loose_cutflow_counts(bib_slice_loose, m)
        fracs = stage_fractions_from_counts(counts, n_bib_tot)
        bib_fail_before_two.append(100 * fracs["fail_before_two"])
        bib_fail_pt_after_two.append(100 * fracs["fail_pt_after_two"])
        bib_fail_beta_after_pt.append(100 * fracs["fail_beta_after_pt"])
        bib_pass_final.append(100 * fracs["pass_final"])

    plt.figure(figsize=(10, 6))
    plt.bar(x, bib_fail_before_two, label="Fail before 2 tracks")
    plt.bar(x, bib_fail_pt_after_two, bottom=bib_fail_before_two, label="Fail pT after 2 tracks")
    plt.bar(x, bib_fail_beta_after_pt,
            bottom=np.array(bib_fail_before_two) + np.array(bib_fail_pt_after_two),
            label="Fail beta after pT")
    plt.bar(x, bib_pass_final,
            bottom=np.array(bib_fail_before_two) + np.array(bib_fail_pt_after_two) + np.array(bib_fail_beta_after_pt),
            label="Pass final")
    plt.xlabel("Signal-tuned cut set [GeV]", fontsize=18)
    plt.ylabel("Percent of total events", fontsize=18)
    plt.title("Loose-window stacked cutflow: BIB", fontsize=15)
    plt.ylim(0, 100)
    plt.xticks(x, x_labels, fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    # muons page
    mu_fail_before_two, mu_fail_pt_after_two, mu_fail_beta_after_pt, mu_pass_final = [], [], [], []
    for m in signal_masses:
        counts = get_loose_cutflow_counts(muons, m)
        fracs = stage_fractions_from_counts(counts, n_mu_tot)
        mu_fail_before_two.append(100 * fracs["fail_before_two"])
        mu_fail_pt_after_two.append(100 * fracs["fail_pt_after_two"])
        mu_fail_beta_after_pt.append(100 * fracs["fail_beta_after_pt"])
        mu_pass_final.append(100 * fracs["pass_final"])

    plt.figure(figsize=(10, 6))
    plt.bar(x, mu_fail_before_two, label="Fail before 2 tracks")
    plt.bar(x, mu_fail_pt_after_two, bottom=mu_fail_before_two, label="Fail pT after 2 tracks")
    plt.bar(x, mu_fail_beta_after_pt,
            bottom=np.array(mu_fail_before_two) + np.array(mu_fail_pt_after_two),
            label="Fail beta after pT")
    plt.bar(x, mu_pass_final,
            bottom=np.array(mu_fail_before_two) + np.array(mu_fail_pt_after_two) + np.array(mu_fail_beta_after_pt),
            label="Pass final")
    plt.xlabel("Signal-tuned cut set [GeV]", fontsize=18)
    plt.ylabel("Percent of total events", fontsize=18)
    plt.title("Loose-window stacked cutflow: muons", fontsize=15)
    plt.ylim(95, 100)
    plt.xticks(x, x_labels, fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

print("Saved stacked cutflow plots to loose_window_stacked_cutflow.pdf")

# ============================================================
# BUILD BASE LIMIT GRID (KNOWN LIFETIMES)
# ============================================================

rows = []

for tau_str in sorted(signal_loose.keys(), key=float):
    for m in signal_masses:
        row = {
            "lifetime_ns": float(tau_str),
            "mass_GeV": int(m),
            "theory_xs_pb": xs_pb[m],
        }

        for cfg in configs:
            cache = signal_map[cfg][tau_str][m]
            n_pass, _ = passing_cuts(cache, m, bib=False)
            denom = denominator_map[cfg][tau_str][m]

            eff = n_pass / denom if denom > 0 else 0.0
            s = eff * xs_pb[m] * l_pb_inv
            sigma_eff = np.sqrt(eff * (1.0 - eff) / denom) if denom > 0 else 0.0
            sigma_s = xs_pb[m] * l_pb_inv * sigma_eff

            row[f"eff_{cfg}"] = eff
            row[f"npass_{cfg}"] = n_pass
            row[f"denom_{cfg}"] = denom
            row[f"s_{cfg}"] = s
            row[f"sigma_s_{cfg}"] = sigma_s
            row[f"s_{cfg}_up"] = s + sigma_s
            row[f"s_{cfg}_down"] = max(0.0, s - sigma_s)
            row[f"excluded_{cfg}"] = int(s > SIGNAL_CUTOFF)
            row[f"excluded_{cfg}_up"] = int((s + sigma_s) > SIGNAL_CUTOFF)
            row[f"excluded_{cfg}_down"] = int(max(0.0, s - sigma_s) > SIGNAL_CUTOFF)

            sigma_limit = SIGNAL_CUTOFF / (eff * l_pb_inv) if eff > 0 else np.inf
            row[f"sigma_limit_pb_{cfg}"] = sigma_limit

        row["is_interpolated_lifetime"] = 0
        rows.append(row)

limit_df = pd.DataFrame(rows)
limit_df["mass_GeV"] = limit_df["mass_GeV"].astype(int)
limit_df["lifetime_ns"] = limit_df["lifetime_ns"].astype(float)

# ============================================================
# AUGMENT WITH 1 ns AND 20 ns
# ============================================================

new_rows = []

for tau_new in extra_lifetimes_ns:
    if np.any(np.isclose(limit_df["lifetime_ns"].unique(), tau_new)):
        continue

    for mass in sorted(limit_df["mass_GeV"].unique()):
        sub = limit_df[limit_df["mass_GeV"] == mass].sort_values("lifetime_ns")

        row = {
            "lifetime_ns": float(tau_new),
            "mass_GeV": int(mass),
            "theory_xs_pb": float(sub["theory_xs_pb"].iloc[0]),
            "is_interpolated_lifetime": 1,
        }

        xs_here = row["theory_xs_pb"]

        for cfg in configs:
            tau_known = sub["lifetime_ns"].to_numpy(dtype=float)
            s_known = sub[f"s_{cfg}"].to_numpy(dtype=float)
            s_up_known = sub[f"s_{cfg}_up"].to_numpy(dtype=float)
            s_down_known = sub[f"s_{cfg}_down"].to_numpy(dtype=float)

            s_new = max(0.0, interp_extrap_vs_lifetime(tau_known, s_known, tau_new))
            s_up_new = max(0.0, interp_extrap_vs_lifetime(tau_known, s_up_known, tau_new))
            s_down_new = max(0.0, interp_extrap_vs_lifetime(tau_known, s_down_known, tau_new))

            eff_new = s_new / (xs_here * l_pb_inv) if xs_here > 0 else 0.0
            sigma_limit = SIGNAL_CUTOFF / (eff_new * l_pb_inv) if eff_new > 0 else np.inf

            row[f"eff_{cfg}"] = eff_new
            row[f"npass_{cfg}"] = np.nan
            row[f"denom_{cfg}"] = np.nan
            row[f"s_{cfg}"] = s_new
            row[f"sigma_s_{cfg}"] = 0.5 * max(0.0, s_up_new - s_down_new)
            row[f"s_{cfg}_up"] = s_up_new
            row[f"s_{cfg}_down"] = s_down_new
            row[f"excluded_{cfg}"] = int(s_new > SIGNAL_CUTOFF)
            row[f"excluded_{cfg}_up"] = int(s_up_new > SIGNAL_CUTOFF)
            row[f"excluded_{cfg}_down"] = int(s_down_new > SIGNAL_CUTOFF)
            row[f"sigma_limit_pb_{cfg}"] = sigma_limit
            row[f"limit_ratio_{cfg}"] = xs_here / sigma_limit if np.isfinite(sigma_limit) else 0.0

        new_rows.append(row)

limit_df = pd.DataFrame(rows + new_rows)
limit_df["mass_GeV"] = limit_df["mass_GeV"].astype(int)
limit_df["lifetime_ns"] = limit_df["lifetime_ns"].astype(float)

# ------------------------------------------------------------
# add interpolated / extrapolated lifetimes
# ------------------------------------------------------------
limit_df = pd.DataFrame(rows + new_rows)
limit_df["mass_GeV"] = limit_df["mass_GeV"].astype(int)
limit_df["lifetime_ns"] = limit_df["lifetime_ns"].astype(float)
limit_df = limit_df.sort_values(["lifetime_ns", "mass_GeV"]).reset_index(drop=True)

limit_df.to_csv("limit_results_full_grid_with_extra_lifetimes.csv", index=False)
print("Saved limit grid to limit_results_full_grid_with_extra_lifetimes.csv")
print("Lifetimes in grid:", sorted(limit_df["lifetime_ns"].unique()))

# ============================================================
# REAL (SIMULATED) POINTS ONLY
# ============================================================

real_df = limit_df[limit_df["is_interpolated_lifetime"] == 0].copy()

# ============================================================
# PLOTTING LIMITS WITH SMALL MARGINS
# ============================================================

XMIN_DATA = real_df["mass_GeV"].min()
XMAX_DATA = real_df["mass_GeV"].max()
YMIN_DATA = real_df["lifetime_ns"].min()
YMAX_DATA = real_df["lifetime_ns"].max()

X_PAD = 120.0

Y_PAD_LOW = 0.85
Y_PAD_HIGH = 1.10

X_LIM = (XMIN_DATA - X_PAD, XMAX_DATA + X_PAD)
Y_LIM = (YMIN_DATA * Y_PAD_LOW, YMAX_DATA * Y_PAD_HIGH)


# ============================================================
# LABEL-ONLY OVERLAY FOR REAL SIMULATED POINTS
# ============================================================

def overlay_real_labels_only(ax, cfg, text_size=8):
    sub = real_df[["mass_GeV", "lifetime_ns", f"s_{cfg}"]].copy()

    xmin = XMIN_DATA
    xmax = XMAX_DATA
    ymin = YMIN_DATA
    ymax = YMAX_DATA

    # Small inward shifts for labels on the boundaries
    dx = 35.0
    low_y_factor = 1.10
    high_y_factor = 1.0 / 1.10

    for _, r in sub.iterrows():
        m = float(r["mass_GeV"])
        tau = float(r["lifetime_ns"])
        val = float(r[f"s_{cfg}"])

        x_text = m
        y_text = tau
        ha = "center"
        va = "center"

        # Move edge labels inward horizontally
        if np.isclose(m, xmin):
            x_text = m + dx
            ha = "left"
        elif np.isclose(m, xmax):
            x_text = m - dx
            ha = "right"

        # Move edge labels inward vertically (important on log axis)
        if np.isclose(tau, ymin):
            y_text = tau * low_y_factor
            va = "bottom"
        elif np.isclose(tau, ymax):
            y_text = tau * high_y_factor
            va = "top"

        txt_color = "black" if val < SIGNAL_CUTOFF else "white"

        ax.text(
            x_text,
            y_text,
            f"{val:.0f}",
            color=txt_color,
            ha=ha,
            va=va,
            fontsize=text_size,
            fontweight="bold",
            zorder=6,
        )


def overlay_real_points_and_labels(ax, cfg, text_size=7, point_size=10):
    """
    Overlay only the REAL simulated grid points and their raw yield labels.
    """
    sub = real_df[["mass_GeV", "lifetime_ns", f"s_{cfg}"]].copy()

    for _, r in sub.iterrows():
        m = float(r["mass_GeV"])
        tau = float(r["lifetime_ns"])
        val = float(r[f"s_{cfg}"])

        txt_color = "black" if val < SIGNAL_CUTOFF else "white"

        ax.scatter(
            m, tau,
            color="black",
            s=point_size,
            zorder=5
        )

        ax.text(
            m, tau,
            f"{val:.0f}",
            color=txt_color,
            ha="center",
            va="center",
            fontsize=text_size,
            fontweight="bold",
            zorder=6,
        )

# ============================================================
# HELPERS FOR PLOTS / CONTOURS
# ============================================================

def find_crossing_mass(masses, yields, cutoff):
    masses = np.asarray(masses, dtype=float)
    yields = np.asarray(yields, dtype=float)

    for i in range(len(masses) - 1):
        y1, y2 = yields[i], yields[i + 1]
        m1, m2 = masses[i], masses[i + 1]

        if y1 == cutoff:
            return m1

        if (y1 - cutoff) * (y2 - cutoff) < 0:
            t = (cutoff - y1) / (y2 - y1)
            return m1 + t * (m2 - m1)

    return np.nan


def centers_to_edges_linear(x):
    x = np.asarray(x, dtype=float)
    edges = np.empty(len(x) + 1)
    edges[1:-1] = 0.5 * (x[:-1] + x[1:])
    edges[0] = x[0] - 0.5 * (x[1] - x[0])
    edges[-1] = x[-1] + 0.5 * (x[-1] - x[-2])
    return edges


def centers_to_edges_log(y):
    y = np.asarray(y, dtype=float)
    logy = np.log10(y)
    edges_log = np.empty(len(y) + 1)
    edges_log[1:-1] = 0.5 * (logy[:-1] + logy[1:])
    edges_log[0] = logy[0] - 0.5 * (logy[1] - logy[0])
    edges_log[-1] = logy[-1] + 0.5 * (logy[-1] - logy[-2])
    return 10 ** edges_log


def get_contour_xy(contour_obj):
    if len(contour_obj.collections) == 0:
        return None, None
    paths = contour_obj.collections[0].get_paths()
    if len(paths) == 0:
        return None, None
    v = paths[0].vertices
    return v[:, 0], v[:, 1]


# ============================================================
# STACKED CUTFLOW PLOTS
# ============================================================

def get_loose_cutflow_counts_final(cache, cut_set):
    lead_pT   = np.asarray(cache["leading_pT"])
    sub_pT    = np.asarray(cache["subleading_pT"])
    lead_beta = np.asarray(cache["leading_beta"])
    sub_beta  = np.asarray(cache["subleading_beta"])

    lead_pT_min = event_cuts[cut_set]["pT"]["lead"]
    sub_pT_min  = event_cuts[cut_set]["pT"]["sublead"]

    lead_beta_lo, lead_beta_hi = event_cuts[cut_set]["beta"]["lead"]
    sub_beta_lo,  sub_beta_hi  = event_cuts[cut_set]["beta"]["sublead"]

    has_two_tracks = np.isfinite(lead_pT) & np.isfinite(sub_pT)

    pass_pt = has_two_tracks & (lead_pT > lead_pT_min) & (sub_pT > sub_pT_min)

    pass_pt_beta = (
        pass_pt &
        np.isfinite(lead_beta) & np.isfinite(sub_beta) &
        (lead_beta_lo < lead_beta) & (lead_beta < lead_beta_hi) &
        (sub_beta_lo < sub_beta) & (sub_beta < sub_beta_hi)
    )

    return {
        "n_two_tracks": int(np.count_nonzero(has_two_tracks)),
        "n_pt": int(np.count_nonzero(pass_pt)),
        "n_final": int(np.count_nonzero(pass_pt_beta)),
    }


def stage_fractions_from_counts(counts, total):
    n_two = counts["n_two_tracks"]
    n_pt = counts["n_pt"]
    n_final = counts["n_final"]

    return {
        "fail_before_two": (total - n_two) / total if total > 0 else 0.0,
        "fail_pt_after_two": (n_two - n_pt) / total if total > 0 else 0.0,
        "fail_beta_after_pt": (n_pt - n_final) / total if total > 0 else 0.0,
        "pass_final": n_final / total if total > 0 else 0.0,
    }


masses_numeric = np.array([int(m) for m in signal_masses])
n_bib_tot = len(np.asarray(bib_slice_loose["leading_pT"]))
n_mu_tot = len(np.asarray(muons["leading_pT"]))

with PdfPages("loose_window_stacked_cutflow.pdf") as pdf:
    # signal pages: only real simulated lifetimes
    for lifetime in sorted(signal_loose.keys(), key=float):
        fail_before_two = []
        fail_pt_after_two = []
        fail_beta_after_pt = []
        pass_final = []

        for m in signal_masses:
            total = denominator_map["loose"][lifetime][m] 
            counts = get_loose_cutflow_counts_final(signal_loose[lifetime][m], m)
            fracs = stage_fractions_from_counts(counts, total)

            fail_before_two.append(100 * fracs["fail_before_two"])
            fail_pt_after_two.append(100 * fracs["fail_pt_after_two"])
            fail_beta_after_pt.append(100 * fracs["fail_beta_after_pt"])
            pass_final.append(100 * fracs["pass_final"])

        plt.figure(figsize=(10, 6))
        plt.bar(masses_numeric, fail_before_two, width=300, label="Fail before 2 tracks")
        plt.bar(masses_numeric, fail_pt_after_two, width=300, bottom=fail_before_two, label="Fail pT after 2 tracks")
        plt.bar(
            masses_numeric,
            fail_beta_after_pt,
            width=300,
            bottom=np.array(fail_before_two) + np.array(fail_pt_after_two),
            label="Fail beta after pT",
        )
        plt.bar(
            masses_numeric,
            pass_final,
            width=300,
            bottom=np.array(fail_before_two) + np.array(fail_pt_after_two) + np.array(fail_beta_after_pt),
            label="Pass final",
        )

        plt.xlabel("Stau mass [GeV]", fontsize=18)
        plt.ylabel("Percent of total events", fontsize=18)
        plt.title(f"Loose-window stacked cutflow: signal ({float(lifetime):g} ns)", fontsize=15)
        plt.ylim(0, 100)
        plt.xticks(masses_numeric, [str(m) for m in masses_numeric], fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=11)
        plt.grid(axis="y", alpha=0.3)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

    # BIB
    x_labels = [str(m) for m in signal_masses]
    x = np.arange(len(signal_masses))

    bib_fail_before_two = []
    bib_fail_pt_after_two = []
    bib_fail_beta_after_pt = []
    bib_pass_final = []

    for m in signal_masses:
        counts = get_loose_cutflow_counts_final(bib_slice_loose, m)
        fracs = stage_fractions_from_counts(counts, n_bib_tot)
        bib_fail_before_two.append(100 * fracs["fail_before_two"])
        bib_fail_pt_after_two.append(100 * fracs["fail_pt_after_two"])
        bib_fail_beta_after_pt.append(100 * fracs["fail_beta_after_pt"])
        bib_pass_final.append(100 * fracs["pass_final"])

    plt.figure(figsize=(10, 6))
    plt.bar(x, bib_fail_before_two, label="Fail before 2 tracks")
    plt.bar(x, bib_fail_pt_after_two, bottom=bib_fail_before_two, label="Fail pT after 2 tracks")
    plt.bar(
        x,
        bib_fail_beta_after_pt,
        bottom=np.array(bib_fail_before_two) + np.array(bib_fail_pt_after_two),
        label="Fail beta after pT",
    )
    plt.bar(
        x,
        bib_pass_final,
        bottom=np.array(bib_fail_before_two) + np.array(bib_fail_pt_after_two) + np.array(bib_fail_beta_after_pt),
        label="Pass final",
    )
    plt.xlabel("Signal-tuned cut set [GeV]", fontsize=18)
    plt.ylabel("Percent of total events", fontsize=18)
    plt.title("Loose-window stacked cutflow: BIB", fontsize=15)
    plt.ylim(0, 100)
    plt.xticks(x, x_labels, fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

    # Muons
    mu_fail_before_two = []
    mu_fail_pt_after_two = []
    mu_fail_beta_after_pt = []
    mu_pass_final = []

    for m in signal_masses:
        counts = get_loose_cutflow_counts_final(muons, m)
        fracs = stage_fractions_from_counts(counts, n_mu_tot)
        mu_fail_before_two.append(100 * fracs["fail_before_two"])
        mu_fail_pt_after_two.append(100 * fracs["fail_pt_after_two"])
        mu_fail_beta_after_pt.append(100 * fracs["fail_beta_after_pt"])
        mu_pass_final.append(100 * fracs["pass_final"])

    plt.figure(figsize=(10, 6))
    plt.bar(x, mu_fail_before_two, label="Fail before 2 tracks")
    plt.bar(x, mu_fail_pt_after_two, bottom=mu_fail_before_two, label="Fail pT after 2 tracks")
    plt.bar(
        x,
        mu_fail_beta_after_pt,
        bottom=np.array(mu_fail_before_two) + np.array(mu_fail_pt_after_two),
        label="Fail beta after pT",
    )
    plt.bar(
        x,
        mu_pass_final,
        bottom=np.array(mu_fail_before_two) + np.array(mu_fail_pt_after_two) + np.array(mu_fail_beta_after_pt),
        label="Pass final",
    )
    plt.xlabel("Signal-tuned cut set [GeV]", fontsize=18)
    plt.ylabel("Percent of total events", fontsize=18)
    plt.title("Loose-window stacked cutflow: muons", fontsize=15)
    plt.ylim(95, 100)
    plt.xticks(x, x_labels, fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend(fontsize=11)
    plt.grid(axis="y", alpha=0.3)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

print("Saved stacked cutflow plots to loose_window_stacked_cutflow.pdf")

# ============================================================
# MAPS FOR HEATMAPS / DEBUG / CONTOURS
# ============================================================

mass_vals = np.array(sorted(limit_df["mass_GeV"].unique()), dtype=float)
lifetime_vals = np.array(sorted(limit_df["lifetime_ns"].unique()), dtype=float)

M, T = np.meshgrid(mass_vals, lifetime_vals)

mass_to_idx = {m: i for i, m in enumerate(mass_vals)}
life_to_idx = {t: i for i, t in enumerate(lifetime_vals)}

S_maps = {}
S_up_maps = {}
S_down_maps = {}
Z_maps = {}
Z_up_maps = {}
Z_down_maps = {}

for cfg in configs:
    S_maps[cfg] = np.full((len(lifetime_vals), len(mass_vals)), np.nan)
    S_up_maps[cfg] = np.full((len(lifetime_vals), len(mass_vals)), np.nan)
    S_down_maps[cfg] = np.full((len(lifetime_vals), len(mass_vals)), np.nan)

    Z_maps[cfg] = np.zeros((len(lifetime_vals), len(mass_vals)), dtype=int)
    Z_up_maps[cfg] = np.zeros((len(lifetime_vals), len(mass_vals)), dtype=int)
    Z_down_maps[cfg] = np.zeros((len(lifetime_vals), len(mass_vals)), dtype=int)

for _, r in limit_df.iterrows():
    i = life_to_idx[float(r["lifetime_ns"])]
    j = mass_to_idx[float(r["mass_GeV"])]

    for cfg in configs:
        S_maps[cfg][i, j] = r[f"s_{cfg}"]
        S_up_maps[cfg][i, j] = r[f"s_{cfg}_up"]
        S_down_maps[cfg][i, j] = r[f"s_{cfg}_down"]
        Z_maps[cfg][i, j] = int(r[f"excluded_{cfg}"])
        Z_up_maps[cfg][i, j] = int(r[f"excluded_{cfg}_up"])
        Z_down_maps[cfg][i, j] = int(r[f"excluded_{cfg}_down"])

# ============================================================
# RBF SMOOTHED MAPS
# ============================================================

RBF_N_MASS = 400
RBF_N_TAU = 400

RBF_SMOOTHING = 0.05

RBF_KERNEL = "thin_plate_spline"

RBF_EPS_YIELD = 1e-6
RBF_YIELD_SCALE = 1.0

USE_ONLY_REAL_POINTS_FOR_RBF = True

RBF_NEIGHBORS = None


mass_dense = np.linspace(mass_vals.min(), mass_vals.max(), RBF_N_MASS)

tau_dense_rbf = np.logspace(
    np.log10(lifetime_vals.min()),
    np.log10(lifetime_vals.max()),
    RBF_N_TAU
)

M_dense, T_dense = np.meshgrid(mass_dense, tau_dense_rbf)


def scaled_rbf_coordinates(mass, lifetime, mass_min, mass_max, logtau_min, logtau_max):
    """
    RBF distance is Euclidean, so rescale mass and log10(lifetime)
    to comparable ranges. Otherwise the 1000-4500 GeV mass scale dominates.
    """
    mass = np.asarray(mass, dtype=float)
    lifetime = np.asarray(lifetime, dtype=float)

    logtau = np.log10(lifetime)

    x = (mass - mass_min) / (mass_max - mass_min)
    y = (logtau - logtau_min) / (logtau_max - logtau_min)

    return np.column_stack([x, y])

def yield_to_rbf_space(y):
    """
    Transform expected yield before RBF fitting.

    This is a soft-log transform. It compresses large yields but does not
    send zero-yield points to an extreme negative value like log10(y) does.
    """
    y = np.asarray(y, dtype=float)
    y = np.clip(y, 0.0, None)
    return np.log1p(y / RBF_YIELD_SCALE)


def yield_from_rbf_space(z):
    """
    Inverse transform back to expected-yield space.
    """
    z = np.asarray(z, dtype=float)
    y = RBF_YIELD_SCALE * np.expm1(z)
    return np.clip(y, 0.0, None)


def build_rbf_yield_map(limit_df, cfg, value_col, mass_dense, tau_dense):
    """
    Build a smooth 2D RBF surface for one quantity, e.g.
    s_nominal, s_nominal_up, or s_nominal_down.

    Uses a soft-log transform instead of log10(yield). This keeps the surface
    positive and smooth but makes the s = 3 contour much more faithful to the
    actual discrete yields.
    """
    sub = limit_df.copy()

    if USE_ONLY_REAL_POINTS_FOR_RBF and "is_interpolated_lifetime" in sub.columns:
        sub = sub[sub["is_interpolated_lifetime"] == 0].copy()

    train_mass = sub["mass_GeV"].to_numpy(dtype=float)
    train_tau = sub["lifetime_ns"].to_numpy(dtype=float)
    train_y = sub[value_col].to_numpy(dtype=float)

    good = (
        np.isfinite(train_mass)
        & np.isfinite(train_tau)
        & np.isfinite(train_y)
        & (train_tau > 0)
    )

    train_mass = train_mass[good]
    train_tau = train_tau[good]
    train_y = train_y[good]

    train_z = yield_to_rbf_space(train_y)

    mass_min = train_mass.min()
    mass_max = train_mass.max()
    logtau_min = np.log10(train_tau.min())
    logtau_max = np.log10(train_tau.max())

    train_xy = scaled_rbf_coordinates(
        train_mass,
        train_tau,
        mass_min,
        mass_max,
        logtau_min,
        logtau_max,
    )

    rbf = RBFInterpolator(
        train_xy,
        train_z,
        kernel=RBF_KERNEL,
        smoothing=RBF_SMOOTHING,
        degree=1,
        neighbors=RBF_NEIGHBORS,
    )

    eval_mass = M_dense.ravel()
    eval_tau = T_dense.ravel()

    eval_xy = scaled_rbf_coordinates(
        eval_mass,
        eval_tau,
        mass_min,
        mass_max,
        logtau_min,
        logtau_max,
    )

    pred_z = rbf(eval_xy).reshape(len(tau_dense), len(mass_dense))
    pred_y = yield_from_rbf_space(pred_z)

    return pred_y


S_smooth_maps = {}
S_up_smooth_maps = {}
S_down_smooth_maps = {}

for cfg in configs:
    S_smooth_maps[cfg] = build_rbf_yield_map(
        limit_df, cfg, f"s_{cfg}", mass_dense, tau_dense_rbf
    )

    S_up_smooth_maps[cfg] = build_rbf_yield_map(
        limit_df, cfg, f"s_{cfg}_up", mass_dense, tau_dense_rbf
    )

    S_down_smooth_maps[cfg] = build_rbf_yield_map(
        limit_df, cfg, f"s_{cfg}_down", mass_dense, tau_dense_rbf
    )

    S_up_smooth_maps[cfg] = np.maximum(S_up_smooth_maps[cfg], S_smooth_maps[cfg])
    S_down_smooth_maps[cfg] = np.minimum(S_down_smooth_maps[cfg], S_smooth_maps[cfg])

# ============================================================
# PRINT YIELD TABLES
# ============================================================

for cfg in configs:
    print(f"\n===== EXPECTED SIGNAL YIELDS: {cfg.upper()} =====")
    pivot = limit_df.pivot(index="lifetime_ns", columns="mass_GeV", values=f"s_{cfg}")
    print(pivot.to_string(float_format=lambda x: f"{x:8.3f}"))

# ============================================================
# RBF-SMOOTHED HEATMAP OF SIGNAL YIELD (using contourf)
# ============================================================

with PdfPages("expected_signal_yield_heatmaps_rbf.pdf") as pdf:
    for cfg in configs:
        fig, ax = plt.subplots(figsize=(8, 6))

        data = np.array(S_smooth_maps[cfg], dtype=float)

        finite_vals = data[np.isfinite(data)]
        vmin = np.min(finite_vals)
        vmax = np.max(finite_vals)
        levels = np.linspace(vmin, vmax, 120)

        cf = ax.contourf(
            M_dense,
            T_dense,
            data,
            levels=levels,
            cmap="viridis"
        )

        cbar = fig.colorbar(cf, ax=ax)
        cbar.set_label("Expected signal yield", fontsize=14)
        cbar.ax.tick_params(labelsize=12)

        ax.contour(
            M_dense,
            T_dense,
            data,
            levels=[SIGNAL_CUTOFF],
            colors="black",
            linewidths=2.0,
        )

        overlay_real_labels_only(ax, cfg, text_size=8)


        ax.set_xlabel("Stau mass [GeV]", fontsize=20)
        ax.set_ylabel("Lifetime [ns]", fontsize=20)
        ax.set_yscale("log")
        ax.set_title(f"RBF-smoothed expected signal yield: {cfg}", fontsize=14)

        ax.set_xlim(*X_LIM)
        ax.set_ylim(*Y_LIM)

        ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
        ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
        ax.grid(True, alpha=0.25)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

print("Saved RBF-smoothed heatmaps to expected_signal_yield_heatmaps_rbf.pdf")

# ============================================================
# RBF-SMOOTHED SIDE-BY-SIDE HEATMAPS WITH < 3 EVENTS SHOWN IN GRAY
# ============================================================

cmap = plt.cm.viridis.copy()
cmap.set_bad(color="lightgray")

all_valid_vals = []
for cfg in configs:
    vals = np.array(S_smooth_maps[cfg], dtype=float)
    vals = vals[np.isfinite(vals) & (vals >= SIGNAL_CUTOFF)]
    if vals.size > 0:
        all_valid_vals.append(vals)

if len(all_valid_vals) > 0:
    all_valid_vals = np.concatenate(all_valid_vals)
    vmin = np.min(all_valid_vals)
else:
    vmin = SIGNAL_CUTOFF

vmax_display = 1000.0

norm = colors.PowerNorm(gamma=0.65, vmin=vmin, vmax=vmax_display)

levels = np.linspace(vmin, vmax_display, 120)

window_titles = {
    "loose": "Loose",
    "medium": "Medium",
    "nominal": "Nominal",
}

with PdfPages("expected_signal_yield_heatmaps_side_by_side_rbf.pdf") as pdf:
    fig, axes = plt.subplots(1, 3, figsize=(19, 6), sharey=True)

    contourf_obj = None

    for ax, cfg in zip(axes, ["nominal", "medium", "loose"]):
        data = np.array(S_smooth_maps[cfg], dtype=float)

        masked_data = np.ma.masked_where(
            (~np.isfinite(data)) | (data < SIGNAL_CUTOFF),
            data
        )

        contourf_obj = ax.contourf(
            M_dense,
            T_dense,
            masked_data,
            levels=levels,
            cmap=cmap,
            norm=norm,
            extend="max"
        )

        ax.contour(
            M_dense,
            T_dense,
            data,
            levels=[SIGNAL_CUTOFF],
            colors="black",
            linewidths=2.0,
        )

        overlay_real_labels_only(ax, cfg, text_size=10)

        ax.set_xlim(*X_LIM)
        ax.set_ylim(*Y_LIM)

        ax.set_title(window_titles[cfg], fontsize=20)
        ax.set_xlabel("Stau mass [GeV]", fontsize=20)
        ax.set_yscale("log")
        ax.tick_params(axis="x", labelsize=14)
        ax.tick_params(axis="y", labelsize=14)
        ax.grid(True, alpha=0.25)

    axes[0].set_ylabel("Lifetime [ns]", fontsize=20)

    panel = axes[2]

    text_kw = dict(
        ha="right",
        va="bottom",
        transform=panel.transAxes,
        bbox=dict(
            facecolor="white",
            alpha=0.8,     
            edgecolor="none",
            boxstyle="round,pad=0.25"
        )
    )

    panel.text(
        0.98, 0.11, "Muon Collider",
        fontsize=16, fontweight="bold", style="italic",
        **text_kw
    )

    panel.text(
        0.98, 0.07, "Simulation, √s = 10 TeV",
        fontsize=12,
        **text_kw
    )

    panel.text(
        0.98, 0.01, "MuColl_v1",
        fontsize=12,
        **text_kw
    )

    fig.subplots_adjust(right=0.88, wspace=0.10)

    cax = fig.add_axes([0.90, 0.15, 0.015, 0.70])
    cbar = fig.colorbar(contourf_obj, cax=cax)
    cbar.set_label("Expected signal yield", fontsize=14)
    cbar.ax.tick_params(labelsize=11)
    cbar.set_ticks([3, 250, 500, 750, 1000])

    pdf.savefig(fig, bbox_inches="tight")
    plt.close(fig)

print("Saved RBF-smoothed side-by-side heatmaps to expected_signal_yield_heatmaps_side_by_side_rbf.pdf")


# ============================================================
# DEBUG GRID OF EXCLUDED POINTS
# ============================================================

with PdfPages("excluded_points_debug.pdf") as pdf:
    plt.figure(figsize=(8, 6))

    markers = {"loose": "o", "nominal": "s", "medium": "^"}
    colors = {"loose": "tab:orange", "nominal": "tab:blue", "medium": "tab:green"}

    for cfg in configs:
        sub = limit_df[limit_df[f"excluded_{cfg}"] == 1]
        plt.scatter(
            sub["mass_GeV"],
            sub["lifetime_ns"],
            marker=markers[cfg],
            color=colors[cfg],
            s=80,
            label=f"{cfg.capitalize()} excluded",
        )

    plt.xlabel("Stau mass [GeV]", fontsize=20)
    plt.ylabel("Lifetime [ns]", fontsize=20)
    plt.yscale("log")
    plt.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    plt.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    plt.title(f"Excluded grid points (cutoff = {SIGNAL_CUTOFF:.3f})")
    plt.legend()
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    pdf.savefig()
    plt.close()

print("Saved debug grid to excluded_points_debug.pdf")

# ============================================================
# 1D CROSS-SECTION LIMITS BY LIFETIME
# ============================================================

plt.style.use("seaborn-v0_8-colorblind")
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

style_map = {
    "nominal": {"color": colors[0], "linestyle": "-",  "linewidth": 2.5},
    "medium":  {"color": colors[1], "linestyle": ":",  "linewidth": 2.5},
    "loose":   {"color": colors[2], "linestyle": "--", "linewidth": 2.5},
}

with PdfPages("expected_cross_section_limits_by_lifetime.pdf") as pdf:
    for lifetime in sorted(limit_df["lifetime_ns"].unique()):
        plt.figure(figsize=(8, 6))
        sub = limit_df[limit_df["lifetime_ns"] == lifetime].sort_values("mass_GeV")

        plt.plot(
            sub["mass_GeV"],
            sub["theory_xs_pb"],
            marker="o",
            color="black",
            linestyle=":",
            linewidth=2.5,
            label=r"$\sigma_{\mathrm{theory}}$"
        )

        for cfg in configs:
            plt.plot(
                sub["mass_GeV"],
                sub[f"sigma_limit_pb_{cfg}"],
                marker="s",
                linewidth=2.0,
                color=style_map[cfg]["color"],
                label=fr"$\sigma_{{95}}$ ({cfg})"
            )

        plt.yscale("log")
        plt.xlabel("Stau mass [GeV]", fontsize=20)
        plt.ylabel("Cross section [pb]", fontsize=20)
        plt.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
        plt.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
        interp_tag = ""
        if "is_interpolated_lifetime" in sub.columns and int(sub["is_interpolated_lifetime"].iloc[0]) == 1:
            interp_tag = " [interpolated]"
        plt.title(f"Expected 95% CL cross-section limits ({lifetime:g} ns){interp_tag}", fontsize=14)
        plt.legend(fontsize=10)
        plt.grid(True, which="both", alpha=0.3)
        plt.tight_layout()
        pdf.savefig()
        plt.close()

print("Saved 1D cross-section limit plots to expected_cross_section_limits_by_lifetime.pdf")

# ============================================================
# RBF-SMOOTHED EXCLUSION CONTOURS WITH FILLED BAND
# ============================================================

def find_rightmost_crossing_mass_with_edge(masses, yields, cutoff):
    """
    Find the largest-mass crossing of yield = cutoff.

    Returns:
        x_cross:
            Actual crossing location. np.nan if no crossing occurs inside
            the plotted mass range.
        x_display:
            Boundary location to use for filled uncertainty bands.
            If the contour has moved outside the plotted range, this is
            clamped to the left or right mass edge.
        status:
            'crossing', 'all_above', 'all_below', or 'bad'.
    """
    masses = np.asarray(masses, dtype=float)
    yields = np.asarray(yields, dtype=float)

    good = np.isfinite(masses) & np.isfinite(yields)
    masses = masses[good]
    yields = yields[good]

    if len(masses) < 2:
        return np.nan, np.nan, "bad"

    crossings = []

    for i in range(len(masses) - 1):
        m1, m2 = masses[i], masses[i + 1]
        y1, y2 = yields[i], yields[i + 1]

        if y1 == cutoff:
            crossings.append(m1)

        if (y1 - cutoff) * (y2 - cutoff) < 0:
            t = (cutoff - y1) / (y2 - y1)
            crossings.append(m1 + t * (m2 - m1))

    if len(crossings) > 0:
        x = np.max(crossings)
        return x, x, "crossing"

        if np.all(yields >= cutoff):
        return np.nan, masses.max(), "all_above"

    if np.all(yields < cutoff):
        return np.nan, masses.min(), "all_below"

    return np.nan, np.nan, "bad"


contour_rows = []

for cfg in configs:
    for i, tau in enumerate(tau_dense_rbf):
        s_tau = S_smooth_maps[cfg][i, :]
        s_up_tau = S_up_smooth_maps[cfg][i, :]
        s_down_tau = S_down_smooth_maps[cfg][i, :]

        m_c, m_c_disp, stat_c = find_rightmost_crossing_mass_with_edge(
            mass_dense, s_tau, SIGNAL_CUTOFF
        )
        m_up, m_up_disp, stat_up = find_rightmost_crossing_mass_with_edge(
            mass_dense, s_up_tau, SIGNAL_CUTOFF
        )
        m_down, m_down_disp, stat_down = find_rightmost_crossing_mass_with_edge(
            mass_dense, s_down_tau, SIGNAL_CUTOFF
        )

        contour_rows.append({
            "config": cfg,
            "lifetime_ns": tau,

            "m_boundary_central": m_c,
            "m_boundary_up": m_up,
            "m_boundary_down": m_down,

            "m_boundary_central_display": m_c_disp,
            "m_boundary_up_display": m_up_disp,
            "m_boundary_down_display": m_down_disp,
            "status_central": stat_c,
            "status_up": stat_up,
            "status_down": stat_down,
        })

contour_df = pd.DataFrame(contour_rows)
contour_df.to_csv("contour_boundaries_with_uncertainty_rbf.csv", index=False)
print("Saved RBF-smoothed contour boundaries to contour_boundaries_with_uncertainty_rbf.csv")


with PdfPages("exclusion_contours_filled_band_rbf.pdf") as pdf:
    fig, ax = plt.subplots(figsize=(8, 6))
    
    real_points = real_df[["mass_GeV", "lifetime_ns"]].drop_duplicates()

    ax.scatter(
        real_points["mass_GeV"],
        real_points["lifetime_ns"],
        facecolors="none",
        edgecolors="0.35",
        linewidths=1.2,
        s=42,
        marker="o",
        alpha=0.9,
        label="Simulated points",
        zorder=1,
    )

    for cfg in configs:
        sub = contour_df[contour_df["config"] == cfg].sort_values("lifetime_ns")

        y = sub["lifetime_ns"].to_numpy(dtype=float)

        x_c = sub["m_boundary_central"].to_numpy(dtype=float)

        x_up_fill = sub["m_boundary_up_display"].to_numpy(dtype=float)
        x_down_fill = sub["m_boundary_down_display"].to_numpy(dtype=float)

        stat_up = sub["status_up"].to_numpy()
        stat_down = sub["status_down"].to_numpy()

        valid_c = np.isfinite(x_c)

        valid_band = (
            np.isfinite(x_up_fill)
            & np.isfinite(x_down_fill)
        )

        both_all_below = (stat_up == "all_below") & (stat_down == "all_below")
        both_all_above = (stat_up == "all_above") & (stat_down == "all_above")
        valid_band = valid_band & (~both_all_below) & (~both_all_above)

        color = style_map[cfg]["color"]

        if np.any(valid_band):
            x_low = np.minimum(x_up_fill[valid_band], x_down_fill[valid_band])
            x_high = np.maximum(x_up_fill[valid_band], x_down_fill[valid_band])

            ax.fill_betweenx(
                y[valid_band],
                x_low,
                x_high,
                color=color,
                alpha=0.20,
                linewidth=0,
            )

        if np.any(valid_c):
            ax.plot(
                x_c[valid_c],
                y[valid_c],
                color=color,
                linestyle=style_map[cfg]["linestyle"],
                linewidth=style_map[cfg]["linewidth"],
                label=cfg.capitalize(),
            )

    for cfg in configs:
        sub_raw = real_df[real_df[f"excluded_{cfg}"] == 1]
        ax.scatter(
            sub_raw["mass_GeV"],
            sub_raw["lifetime_ns"],
            color=style_map[cfg]["color"],
            s=18,
            alpha=0.45,
            marker="o",
            zorder=2,
        )

    ax.text(
        0.02, 0.97, "Muon Collider",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=16, fontweight="bold", style="italic"
    )

    ax.text(
        0.02, 0.93, "Simulation, √s = 10 TeV",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=12
    )
    
    ax.text(
        0.02, 0.89, "MuColl_v1",
        ha="left", va="top",
        transform=ax.transAxes,
        fontsize=12
    )

    ax.set_xlabel("Stau mass [GeV]", fontsize=20)
    ax.set_ylabel("Lifetime [ns]", fontsize=20)

    ax.set_xlim(*X_LIM)
    ax.set_ylim(*Y_LIM)
    ax.set_yscale("log")

    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)

    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)

print("Saved RBF-smoothed exclusion contour plot to exclusion_contours_filled_band_rbf.pdf")


# ============================================================
# COMBINED LOOSE STACKED CUTFLOW:
# ============================================================

selected_masses = ["1000", "2000", "3000", "4000"]
selected_labels = ["1000", "2000", "3000", "4000", "BIB", r"$\mu^+\mu^-$"]

# choose which signal lifetime to plot
plot_lifetime = sorted(signal_loose.keys(), key=float)[1]   # change to "3" or "30" if you want

fail_before_two = []
fail_pt_after_two = []
fail_beta_after_pt = []
pass_final = []

# -------------------------
# signal bars
# -------------------------
for m in selected_masses:
    total = denominator_map["loose"][plot_lifetime][m]
    counts = get_loose_cutflow_counts_final(signal_loose[plot_lifetime][m], m)
    fracs = stage_fractions_from_counts(counts, total)

    fail_before_two.append(100 * fracs["fail_before_two"])
    fail_pt_after_two.append(100 * fracs["fail_pt_after_two"])
    fail_beta_after_pt.append(100 * fracs["fail_beta_after_pt"])
    pass_final.append(100 * fracs["pass_final"])

# -------------------------
# BIB bar
# -------------------------
bib_cut_set = "3000"
counts = get_loose_cutflow_counts_final(bib_slice_loose, bib_cut_set)
fracs = stage_fractions_from_counts(counts, n_bib_tot)

fail_before_two.append(100 * fracs["fail_before_two"])
fail_pt_after_two.append(100 * fracs["fail_pt_after_two"])
fail_beta_after_pt.append(100 * fracs["fail_beta_after_pt"])
pass_final.append(100 * fracs["pass_final"])

# -------------------------
# muon bar
# -------------------------
mu_cut_set = "3000"
counts = get_loose_cutflow_counts_final(muons, mu_cut_set)
fracs = stage_fractions_from_counts(counts, n_mu_tot)

fail_before_two.append(100 * fracs["fail_before_two"])
fail_pt_after_two.append(100 * fracs["fail_pt_after_two"])
fail_beta_after_pt.append(100 * fracs["fail_beta_after_pt"])
pass_final.append(100 * fracs["pass_final"])

# -------------------------
# plot
# -------------------------
x = np.arange(len(selected_labels))
width = 0.75

plt.figure(figsize=(11, 6))

plt.bar(x, fail_before_two, width=width, label="Fail before 2 tracks")
plt.bar(x, fail_pt_after_two, width=width,
        bottom=fail_before_two, label="Fail pT after 2 tracks")
plt.bar(x, fail_beta_after_pt, width=width,
        bottom=np.array(fail_before_two) + np.array(fail_pt_after_two),
        label="Fail beta after pT")
plt.bar(x, pass_final, width=width,
        bottom=np.array(fail_before_two) + np.array(fail_pt_after_two) + np.array(fail_beta_after_pt),
        label="Pass final")

plt.xlabel("Sample / cut set", fontsize=20)
plt.ylabel("Percent of total events", fontsize=20)
plt.title(f"Loose window stacked cutflow ({plot_lifetime} ns signal)", fontsize=15)
plt.ylim(0, 100)
plt.xticks(x, selected_labels, fontsize=16)
plt.yticks(fontsize=16)
plt.legend(fontsize=12)
plt.grid(axis="y", alpha=0.3)
plt.tight_layout()
plt.savefig("loose_window_stacked_cutflow_combined.pdf")
plt.close()

print("Saved plot to loose_window_stacked_cutflow_combined.pdf")


# ============================================================
# 1D CROSS-SECTION LIMITS: ALL REAL LIFETIMES ON ONE PLOT
# ============================================================

plt.style.use("seaborn-v0_8-colorblind")
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]

style_map = {
    "nominal": {"color": colors[0], "linestyle": "-"},
    "medium":  {"color": colors[1], "linestyle": ":"},
    "loose":   {"color": colors[2], "linestyle": "--"},
}

real_lifetimes = [t for t in sorted(limit_df["lifetime_ns"].unique()) if t in [3.0, 10.0, 30.0]]

alpha_map = {
    30.0: 1.00,
    10.0: 0.70,
    3.0: 0.40,
}

with PdfPages("expected_cross_section_limits_all_real_lifetimes.pdf") as pdf:
    fig, ax = plt.subplots(figsize=(9, 7))

    theory_sub = (
        limit_df[limit_df["lifetime_ns"] == real_lifetimes[0]]
        .sort_values("mass_GeV")
    )
    ax.plot(
        theory_sub["mass_GeV"],
        theory_sub["theory_xs_pb"],
        color="black",
        linestyle="-",
        linewidth=2.5,
        marker="o",
        label=r"$\sigma_{\mathrm{theory}}$"
    )

    for tau in sorted(real_lifetimes, reverse=True):   # 30, 10, 3
        sub_tau = limit_df[limit_df["lifetime_ns"] == tau].sort_values("mass_GeV")

        for cfg in ["nominal", "medium", "loose"]:
            ax.plot(
                sub_tau["mass_GeV"],
                sub_tau[f"sigma_limit_pb_{cfg}"],
                color=style_map[cfg]["color"],
                linestyle=style_map[cfg]["linestyle"],
                linewidth=2.5,
                marker="s",
                alpha=alpha_map[tau],
            )
    from matplotlib.lines import Line2D

    window_handles = [
        Line2D([0], [0],
            color=style_map[cfg]["color"],
            linestyle=style_map[cfg]["linestyle"],
            linewidth=2.5,
            label=cfg.capitalize())
        for cfg in ["nominal", "medium", "loose"]
    ]

    theory_handle = Line2D(
        [0], [0],
        color="black",
        linestyle="-",
        linewidth=2.5,
        marker="o",
        label=r"$\sigma_{\mathrm{theory}}$"
    )

    lifetime_handles = [
        Line2D([0], [0],
            color="black",
            linestyle="-",
            linewidth=3,
            alpha=alpha_map[tau],
            label=fr"${int(tau)}$ ns")
        for tau in [30.0, 10.0, 3.0]
    ]

    legend1 = ax.legend(
        handles=[theory_handle] + window_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.16),
        ncol=4,
        fontsize=14,
        frameon=False
    )

    legend2 = ax.legend(
        handles=lifetime_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, 1.10),
        ncol=3,
        fontsize=14,
        frameon=False,
        # title="Stau lifetime:"
    )

    ax.add_artist(legend1)

    ax.set_yscale("log")
    ax.set_xlabel("Stau mass [GeV]", fontsize=20)
    ax.set_ylabel("Cross section [pb]", fontsize=20)
    ax.tick_params(axis="both", which="major", labelsize=16, length=6, width=1.5)
    ax.tick_params(axis="both", which="minor", labelsize=14, length=4, width=1.0)
    # ax.set_title("Expected 95% CL cross-section limits", fontsize=15)
    ax.grid(True, which="both", alpha=0.3)

    ax.text(0.98, 0.98, "Muon Collider", ha="right", va="top",
            transform=ax.transAxes, fontsize=16, fontweight="bold", style="italic")
    ax.text(0.98, 0.93, r"Expected 95$\%$ CL cross-section limits",
            ha="right", va="top", transform=ax.transAxes, fontsize=12)
    # ax.text(0.02, 0.90, f"N_events={N_events}",
    #         ha="left", va="top", transform=ax.transAxes, fontsize=12)
    
    plt.subplots_adjust(top=0.87)
    pdf.savefig(fig)
    plt.close(fig)

print("Saved combined cross-section limit plot to expected_cross_section_limits_all_real_lifetimes.pdf")