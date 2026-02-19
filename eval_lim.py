import numpy as np
import pickle
import pandas as pd

SIG_CACHE = "/scratch/miralittmann/analysis/mira_analysis_code/cache/sig_by_event.pkl"
BKG_CACHE = "/scratch/miralittmann/analysis/mira_analysis_code/cache/mumu_bkg_stats_nominal_nobib_byevent.pkl"

WINDOW  = "nominal"
REQ     = "ob"
OPTION  = "nobib"
LIFETIME = 30  

L_ABINV = 10.0
L_PBINV = L_ABINV * 1e6  # 10 ab^-1 = 1e7 pb^-1

# Cross sections in pb, taken from MadGraph
XS_PB = {
    "mumu_bkg": 0.4312,
    "1000": 0.0005108,
    "1500": 0.0004715,
    "2000": 0.0004184,
    "2500": 0.0003528,
    "3000": 0.0002780,
    "3500": 0.0001980,
    "4000": 0.0001173,
    "4500": 0.0000450,
}

SIGNAL_SAMPLES = ["1000","1500","2000","2500","3000","3500","4000","4500"]

M_CUTS = np.arange(0, 5000 + 50, 50)          # GeV, require mass > mass_cut
BETA_CUTS = np.arange(0.90, 1.0001, 0.001)    # require beta < beta_cut

B_TARGET = 1.0  # expected background events allowed after cuts. starting with 1 

MIN_BKG_MC_PASS = 0  # set to 0 to disable, just because don't have enough events yet

# Track requirement mode:
# "ge1" = at least one candidate passes (use leading only)
# "ge2" = require leading AND subleading both pass
MODE = "ge1"


def load_pickle(path):
    with open(path, "rb") as f:
        return pickle.load(f)

def expected_events(xs_pb, eff, L_pbinv=L_PBINV):
    return float(xs_pb) * float(L_pbinv) * float(eff)

def align_and_mask(mL, bL, mS, bS, denom_mask):
    mL = np.asarray(mL, dtype=float)
    bL = np.asarray(bL, dtype=float)
    denom_mask = np.asarray(denom_mask, dtype=bool)

    n = min(len(mL), len(bL), len(denom_mask))
    mL, bL, denom_mask = mL[:n], bL[:n], denom_mask[:n]

    if mS is not None and bS is not None:
        mS = np.asarray(mS, dtype=float)[:n]
        bS = np.asarray(bS, dtype=float)[:n]
    else:
        mS, bS = None, None

    mL = mL[denom_mask]
    bL = bL[denom_mask]
    if mS is not None:
        mS = mS[denom_mask]
        bS = bS[denom_mask]

    return mL, bL, mS, bS

def pass_mask_grid(mL, bL, mS, bS, mc, bc, mode="ge1"):
    finL = np.isfinite(mL) & np.isfinite(bL)
    if mode == "ge1":
        return finL & (mL > mc) & (bL < bc)

    if mode == "ge2":
        if mS is None or bS is None:
            raise ValueError("mode='ge2' requires subleading arrays.")
        finS = np.isfinite(mS) & np.isfinite(bS)
        return (finL & finS &
                (mL > mc) & (bL < bc) &
                (mS > mc) & (bS < bc))

    raise ValueError("mode must be 'ge1' or 'ge2'")

def compute_efficiency_grid(mass_lead, beta_lead,
                            mass_sub, beta_sub,
                            denom_mask,
                            m_cuts, beta_cuts,
                            mode="ge1"):
    """
    returns
      eff_grid[j,i] = efficiency for beta_cuts[j], m_cuts[i]
      n_pass_grid[j,i] = number of passing events in denominator sample
      N = denominator size
    """
    mL, bL, mS, bS = align_and_mask(mass_lead, beta_lead, mass_sub, beta_sub, denom_mask)
    N = len(mL)
    if N == 0:
        eff = np.full((len(beta_cuts), len(m_cuts)), np.nan)
        n_pass = np.zeros_like(eff, dtype=int)
        return eff, n_pass, 0

    eff = np.zeros((len(beta_cuts), len(m_cuts)), dtype=float)
    n_pass = np.zeros((len(beta_cuts), len(m_cuts)), dtype=int)

    for j, bc in enumerate(beta_cuts):
        for i, mc in enumerate(m_cuts):
            pm = pass_mask_grid(mL, bL, mS, bS, mc, bc, mode=mode)
            k = int(np.sum(pm))
            n_pass[j, i] = k
            eff[j, i] = k / N

    return eff, n_pass, N

def pick_working_point_expected_space(
    Zs, Zb, npass_bkg_mc, npass_sig_mc,
    m_cuts, beta_cuts,
    xs_sig_pb, xs_bkg_pb, L_pbinv,
    B_target=1.0,
    min_bkg_mc_pass=0,
):
    """
      choose cuts with Bexp leq B_target
      optionally require at least min_bkg_mc_pass passing bkg MC events
      among allowed cuts, maximize signal (Sexp)
    """
    Bexp = xs_bkg_pb * L_pbinv * Zb
    Sexp = xs_sig_pb * L_pbinv * Zs

    ok = np.isfinite(Bexp) & np.isfinite(Sexp)
    ok &= (Bexp <= B_target)
    if min_bkg_mc_pass > 0:
        ok &= (npass_bkg_mc >= min_bkg_mc_pass)

    if not np.any(ok):
        print("none above min bkg pass")
        return None

    # maximize expected signal
    Sexp_ok = np.where(ok, Sexp, -np.inf)
    flat = int(np.argmax(Sexp_ok))
    j, i = np.unravel_index(flat, Sexp.shape)

    return {
        "m_cut": float(m_cuts[i]),
        "beta_cut": float(beta_cuts[j]),
        "sig_eff": float(Zs[j, i]),
        "bkg_eff": float(Zb[j, i]),
        "Sexp": float(Sexp[j, i]),
        "Bexp": float(Bexp[j, i]),
        "bkg_mc_pass": int(npass_bkg_mc[j, i]),
        "sig_mc_pass": int(npass_sig_mc[j, i])
    }



signal = load_pickle(SIG_CACHE)
muons  = load_pickle(BKG_CACHE)

d_bkg = muons[WINDOW][REQ][OPTION]
bkg_mL = d_bkg.get("leading_mass", [])
bkg_bL = d_bkg.get("leading_beta", [])
bkg_mS = d_bkg.get("subleading_mass", None)
bkg_bS = d_bkg.get("subleading_beta", None)

# background denominator: all events in cache
denom_bkg = np.ones(len(bkg_mL), dtype=bool)

Zb, nb_pass_bkg, Nb_den = compute_efficiency_grid(
    bkg_mL, bkg_bL, bkg_mS, bkg_bS,
    denom_bkg,
    M_CUTS, BETA_CUTS,
    mode=MODE
)

print(f"[bkg] denominator events = {Nb_den}")

rows = []
for sample in SIGNAL_SAMPLES:
    d_sig = signal[LIFETIME][sample][REQ]

    sig_mL = d_sig.get("leading_mass", [])
    sig_bL = d_sig.get("leading_beta", [])
    sig_mS = d_sig.get("subleading_mass", None)
    sig_bS = d_sig.get("subleading_beta", None)

    # signal denominator:
    # denominator = "events with at least 1 accepted stau" --> has1_acc_stau. TODO: need to change this???
    denom_sig = np.asarray(d_sig.get("has1_acc_stau", np.ones(len(sig_mL), dtype=bool)), dtype=bool)

    Zs, ns_pass_sig, Ns_den = compute_efficiency_grid(
        sig_mL, sig_bL, sig_mS, sig_bS,
        denom_sig,
        M_CUTS, BETA_CUTS,
        mode=MODE
    )

    best = pick_working_point_expected_space(
        Zs=Zs, Zb=Zb, npass_bkg_mc=nb_pass_bkg, npass_sig_mc=ns_pass_sig,
        m_cuts=M_CUTS, beta_cuts=BETA_CUTS,
        xs_sig_pb=XS_PB[sample],
        xs_bkg_pb=XS_PB["mumu_bkg"],
        L_pbinv=L_PBINV,
        B_target=B_TARGET,
        min_bkg_mc_pass=MIN_BKG_MC_PASS,
    )

    if best is None:
        rows.append({
            "mass_GeV": int(sample),
            "note": f"No cut found with Bexp <= {B_TARGET} and bkg_mc_pass >= {MIN_BKG_MC_PASS}",
        })
        continue

    rows.append({
        "mass_GeV": int(sample),
        "m_cut_GeV": best["m_cut"],
        "beta_cut": best["beta_cut"],
        "sig_eff": best["sig_eff"],
        "bkg_eff": best["bkg_eff"],
        "Sexp_10ab": best["Sexp"],
        # "Bexp_10ab": best["Bexp"],
        "bkg_pass": best["bkg_mc_pass"],
        # "sig_pass": best["sig_mc_pass"],
        # "exclude_S>=3": (best["Sexp"] >= 3.0),
        # "discover_S>=5": (best["Sexp"] >= 5.0),
        "Ns_den": int(Ns_den),
        "Nb_den": int(Nb_den),
    })

df = pd.DataFrame(rows).sort_values("mass_GeV")
pd.set_option("display.max_columns", None)
print(df.to_string(index=False, float_format=lambda x: f"{x:.4g}"))
