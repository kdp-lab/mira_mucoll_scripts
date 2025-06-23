import numpy as np
import pandas as pd

stau_masses = np.array([1, 2.5, 4, 4.5]) # in TeV. just using what we have rn in the tbl files, could do more
com_energy = 10 # also in TeV, CoM energy in muon collider
stau_energy = com_energy / 2 # assuming it's split evenly (is this fair?)
c = 299792458/1000000  # mm/ns 

# E^2 = m^2c^4 + p^2c^2

# locations of subdetector layers (in mm) (just barrel)
subdet_locs = {
    "vb_0": 30,
    "vb_1": 51,
    "vb_3": 74,
    "vb_4": 102,

    "ib_0": 127,
    "ib_1": 340,
    "ib_2": 554,

    "ob_0": 819,
    "ob_1": 1153,
    "ob_2": 1486
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
print("Corrected Time-of-Flight (ns) for Stau particles:")
print(stau_tof_df)
print(" ")


     