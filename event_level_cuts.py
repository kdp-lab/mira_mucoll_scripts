import pickle
import numpy as np

stau_cache = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_with_track_level.pkl"
bib_cache = "/scratch/wandriscok/kate_mucoll_scripts/bib_analysis/cache/bib_event_plot_lead_sub.pkl"
# muon_cache = "


xs_pb = {
    "mumu_bkg": 0.4312,
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
l_abinv = 10.0
l_pb_inv = l_abinv * 1e6  # 1 ab^-1 = 1e6 pb^-1

total_events = 2500  

event_cuts = {
    "1000": {
        "pT": {
            "lead": 2500, # as a minimum
            "sublead": 1500
        },
        "mass": {
            "lead": (500, 2000),  # (low, high)
            "sublead": (500, 2000)
        },
        "beta": {
            "lead": (0.95, 0.973), 
            "sublead": (0.95, 0.973)
        }
    },
    "1500": {
        "pT": {
            "lead": 2500, 
            "sublead": 2000
        },
        "mass": {
            "lead": (1000, 2000), 
            "sublead": (950, 2000)
        },
        "beta": {
            "lead": (0.91, 0.97), 
            "sublead": (0.91, 0.97)
        }
    },
    "2000": {
        "pT": {
            "lead": 3000, 
            "sublead": 2500
        },
        "mass": {
            "lead": (1500, 3000), 
            "sublead": (1500, 2500)
        },
        "beta": {
            "lead": (0.89, 0.93), 
            "sublead": (0.89, 0.93)
        }
    },
    "2500": {
        "pT": {
            "lead": 2800, 
            "sublead": 2000
        },
        "mass": {
            "lead": (2000, 4000), 
            "sublead": (1900, 3100)
        },
        "beta": {
            "lead": (0.85, 0.89), 
            "sublead": (0.85, 0.89)
        }
    },
    "3000": {
        "pT": {
            "lead": 2700, 
            "sublead": 2000
        },
        "mass": {
            "lead": (2000, 6000), 
            "sublead": (2000, 4500)
        },
        "beta": {
            "lead": (0.70, 0.82), 
            "sublead": (0.70, 0.82)
        }
    },
    "3500": {
        "pT": {
            "lead": 2500, 
            "sublead": 2000
        },
        "mass": {
            "lead": (2900, 5000), 
            "sublead": (2900, 4200)
        },
        "beta": {
            "lead": (0.67, 0.74), 
            "sublead": (0.67, 0.74)
        }
    },
    "4000": {
        "pT": {
            "lead": 2000, 
            "sublead": 1900
        },
        "mass": {
            "lead": (3200, 5500), 
            "sublead": (3200, 4500)
        },
        "beta": {
            "lead": (0.57, 0.61), 
            "sublead": (0.57, 0.61)
        }
    },
    "4500": {
        "pT": {
            "lead": 1500, 
            "sublead": 1400
        },
        "mass": {
            "lead": (3500, 6000), 
            "sublead": (3500, 5500)
        },
        "beta": {
            "lead": (0.41, 0.45), 
            "sublead": (0.41, 0.45)
        }
    }
}
z0_cut = 0.75 # req abs(value) < this
w_rms_cut = 1.6

def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)

def passing_cuts(cache, cut_set, bib=False):
    lead_beta = np.asarray(cache["leading_beta"])
    sub_beta  = np.asarray(cache["subleading_beta"])

    lead_pT = np.asarray(cache["leading_pT"])
    sub_pT  = np.asarray(cache["subleading_pT"])

    lead_mass = np.asarray(cache["leading_mass"])
    sub_mass  = np.asarray(cache["subleading_mass"])

    lead_z0 = np.asarray(cache["leading_z0"])
    sub_z0  = np.asarray(cache["subleading_z0"])

    if bib:
        lead_wrms = np.asarray(cache["leading_w_rms"])
        sub_wrms  = np.asarray(cache["subleading_w_rms"])
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

    for i in range(len(lead_pT)):
        # lead requirements
        passing_lead_pT   = np.isfinite(lead_pT[i])   and (lead_pT[i] > lead_pT_min)
        passing_lead_beta = np.isfinite(lead_beta[i]) and (lead_beta_lo < lead_beta[i] < lead_beta_hi)
        passing_lead_mass = np.isfinite(lead_mass[i]) and (lead_mass_lo < lead_mass[i] < lead_mass_hi)
        passing_lead_z0   = np.isfinite(lead_z0[i])   and (abs(lead_z0[i]) < z0_cut)
        passing_lead_wrms = np.isfinite(lead_wrms[i]) and (lead_wrms[i] < w_rms_cut)

        if not (passing_lead_pT and passing_lead_beta and passing_lead_mass and passing_lead_z0 and passing_lead_wrms):
            continue

        # sublead requirements
        passing_sub_pT   = np.isfinite(sub_pT[i])   and (sub_pT[i] > sub_pT_min)
        passing_sub_beta = np.isfinite(sub_beta[i]) and (sub_beta_lo < sub_beta[i] < sub_beta_hi)
        passing_sub_mass = np.isfinite(sub_mass[i]) and (sub_mass_lo < sub_mass[i] < sub_mass_hi)
        passing_sub_z0   = np.isfinite(sub_z0[i])   and (abs(sub_z0[i]) < z0_cut)
        passing_sub_wrms = np.isfinite(sub_wrms[i]) and (sub_wrms[i] < w_rms_cut)

        if passing_sub_pT and passing_sub_beta and passing_sub_mass and passing_sub_z0 and passing_sub_wrms:
            num_passing += 1

    return num_passing, num_total 

signal = load_pickle(stau_cache)   # signal[lifetime][mass]
bib    = load_pickle(bib_cache)    # bib["loose"]["ob"]["10_bib"]
bib_slice = bib["loose"]["ob"]["10_bib"]
signal_masses = ["1000","1500","2000","2500","3000","3500","4000","4500"]
print(f"Assuming L = {l_abinv} ab^-1 = {l_pb_inv:.3e} pb^-1")
print(f"Raw denominator = {total_events}\n")

for lifetime in sorted(signal.keys()):
    print(f"=== lifetime {lifetime} ns ===")
    for m in signal_masses:
        # raw 
        n_sig_pass, n_sig_tot = passing_cuts(signal[lifetime][m], m, bib=False)
        n_bib_pass, n_bib_tot = passing_cuts(bib_slice, m, bib=True)

        frac_sig = n_sig_pass / n_sig_tot 
        frac_bib = n_bib_pass / n_bib_tot 

        # scaled expected yields at 10 ab^-1
        exp_sig = frac_sig * xs_pb[m] * l_pb_inv 
        exp_bib = frac_bib * xs_pb["bib_bkg"] * l_pb_inv

        print(
            f"m={m:>4} | "
            f"sig raw {n_sig_pass:4d}/{n_sig_tot} ({frac_sig:7.4f}) -> expected {exp_sig:12.4g} | "
            f"bib raw {n_bib_pass:4d}/{n_bib_tot} ({frac_bib:7.4f}) -> expected {exp_bib:12.4g}"
        )
    print()