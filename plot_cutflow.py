import matplotlib.pyplot as plt
import numpy as np

loose = [3499110, 551402, 549865, 113591, 113544, 6]
medium = [2278715, 315934, 314930, 65408, 65389, 78]
tight = [1070336, 165640, 165191, 27081, 27076, 857]

loose_arr = np.array(loose, dtype=float)
medium_arr = np.array(medium, dtype=float)
tight_arr = np.array(tight, dtype=float)

windows = ["Nominal", "Medium", "Loose"]
cuts = ["Total",  
    r"$|\eta| \leq 0.8$",
    r"$\chi^2/ndf < 3$",
    "Hit req",
    r"$p_T < 10 TeV$",
    "W.RMS < 1.6"]

plt.style.use("seaborn-v0_8-colorblind")
def make_reg_cutflow():
    x = np.arange(len(cuts))

    plt.figure(figsize=(10, 6))

    plt.plot(x, tight,  marker='o', linewidth=2, label='Nominal window')
    plt.plot(x, medium, marker='o', linewidth=2, label='Medium window')
    plt.plot(x, loose,  marker='o', linewidth=2, label='Loose window')

    plt.xticks(x, cuts, rotation=20, ha='right', fontsize=20)
    plt.yticks(fontsize=16)
    plt.yscale('log')
    plt.ylabel("Surviving tracks", fontsize=20)
    plt.title('100% BIB Track-Level Cutflow', fontsize=20)
    plt.grid(True, which='both', linestyle='--', alpha=0.4)
    plt.legend(fontsize=13)
    plt.tight_layout()
    plt.savefig('/scratch/miralittmann/analysis/mira_analysis_code/cutflow.pdf')

make_reg_cutflow()