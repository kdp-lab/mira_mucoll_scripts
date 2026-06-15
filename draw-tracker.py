import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 

pdf_path = "/scratch/miralittmann/analysis/mira_analysis_code/mucoll_tracker.pdf"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

pdf_path = "/scratch/miralittmann/analysis/mira_analysis_code/mucoll_tracker.pdf"


vb_r = [30, 51, 74, 102]
vb_l = 130

ib_r = [127, 340, 554]
ib_l = [963.2, 963.2, 1384.6]

ob_r = [819, 1153, 1486]
ob_l = 2528.4

# Vertex endcap disks
ve_z    = [80, 82, 120, 122, 200, 202, 280, 282]
ve_rmin = [25, 25,  31,  31,  38,  38,  53,  53]
ve_rmax = [112]*len(ve_z)

# Inner endcap disks
ie_z    = [524, 808, 1093, 1377, 1661, 1946, 2190]
ie_rmin = [95,  147, 190,  212,  237,  264,  284]
ie_rmax = [405, 555, 555,  555,  555,  555,  555]

# Outer endcap disks
oe_z    = [1310, 1617, 1883, 2190]
oe_rmin = [617.5]*len(oe_z)
oe_rmax = [1430.2]*len(oe_z)

def get_eta_lines(ax, etas, zmax, rmax, color='gray', ls='--', alpha=0.8, lw=1.2):
    z_vals = np.linspace(0, zmax, 400)
    for eta in etas:
        theta = 2*np.arctan(np.exp(-eta))
        r_vals = np.tan(theta) * z_vals

        ax.plot( z_vals, r_vals, ls=ls, color=color, alpha=alpha, lw=lw)
        ax.plot(-z_vals, r_vals, ls=ls, color=color, alpha=alpha, lw=lw)

        z_top = rmax / np.tan(theta)
        if z_top <= zmax:
            x_edge, y_edge = z_top, rmax
            ax.annotate(f"$\\eta$={eta:.1f}", xy=(x_edge, y_edge),
                        xytext=(0, 4), textcoords='offset points',
                        rotation=0, color=color, fontsize=14)
        else:
            r_right = np.tan(theta) * zmax
            r_right = min(r_right, rmax)
            x_edge, y_edge = zmax, r_right
            ax.annotate(f"$\\eta$={eta:.1f}", xy=(x_edge, y_edge),
                        xytext=(4, 0), textcoords='offset points',
                        ha='left', va='center', rotation=0, color=color, fontsize=14)
            
def get_theta_lines(ax, thetas_deg, zmax, rmax, color='gray', ls='--', alpha=0.8, lw=1.2):
    z_vals = np.linspace(0, zmax, 400)

    for theta_deg in thetas_deg:
        theta = np.radians(theta_deg)

        r_vals = np.tan(theta) * z_vals

        ax.plot(z_vals, r_vals, ls=ls, color=color, alpha=alpha, lw=lw)
        ax.plot(-z_vals, r_vals, ls=ls, color=color, alpha=alpha, lw=lw)

        z_top = rmax / np.tan(theta)

        if z_top <= zmax:
            x_edge, y_edge = z_top, rmax
            ax.annotate(f"$\\theta$={theta_deg:.0f}°", xy=(x_edge, y_edge),
                        xytext=(0,4), textcoords='offset points',
                        fontsize=14, color=color)

        else:
            r_right = np.tan(theta) * zmax
            r_right = min(r_right, rmax)

            ax.annotate(f"$\\theta$={theta_deg:.0f}°",
                        xy=(zmax, r_right),
                        xytext=(4,0), textcoords='offset points',
                        ha='left', va='center',
                        fontsize=14, color=color)

with PdfPages(pdf_path) as pdf:
    fig, ax = plt.subplots(figsize=(10, 6))

    for r in vb_r:
        ax.plot([-vb_l/2, vb_l/2], [r, r], color="#348ec2")
        ax.plot([-vb_l/2, vb_l/2], [r+2, r+2], color="#348ec2")

    for z, rmin, rmax in zip(ve_z, ve_rmin, ve_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#348ec2")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#348ec2")

    for i, r in enumerate(ib_r):
        ax.plot([-ib_l[i]/2, ib_l[i]/2], [r, r], color="#34b18f")

    for z, rmin, rmax in zip(ie_z, ie_rmin, ie_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#34b18f")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#34b18f")

    for r in ob_r:
        ax.plot([-ob_l/2, ob_l/2], [r, r], color="#dd7e33")

    for z, rmin, rmax in zip(oe_z, oe_rmin, oe_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#dd7e33")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#dd7e33")

    ax.set_xlabel("z [mm]", fontsize=20)
    ax.set_ylabel("r [mm]", fontsize=20)

    ax.plot([], [], marker='o', label="Vertex Detector", linestyle=None, color="#348ec2")
    ax.plot([], [], marker='o', label="Inner Tracker",   linestyle=None, color="#34b18f")
    ax.plot([], [], marker='o', label="Outer Tracker",   linestyle=None, color="#dd7e33")
    ax.legend(loc="upper left", fontsize=14)

    ax.minorticks_on()
    ax.tick_params(labelsize=17)
    # ax.set_xlim(-2500, 2500)
    ax.set_xlim(0, 2500)
    ax.set_ylim(0, 1600)

    etas = [0.5, 1.0, 1.5, 2.0, 2.4, 3.0]
    get_eta_lines(ax, etas, zmax=ax.get_xlim()[1], rmax=ax.get_ylim()[1])

    pdf.savefig(fig)
    plt.close(fig)

    
    fig, ax = plt.subplots(figsize=(10, 6))

    for r in vb_r:
        ax.plot([-vb_l/2, vb_l/2], [r, r], color="#348ec2")
        ax.plot([-vb_l/2, vb_l/2], [r+2, r+2], color="#348ec2")

    for z, rmin, rmax in zip(ve_z, ve_rmin, ve_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#348ec2")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#348ec2")

    for i, r in enumerate(ib_r):
        ax.plot([-ib_l[i]/2, ib_l[i]/2], [r, r], color="#34b18f")

    for z, rmin, rmax in zip(ie_z, ie_rmin, ie_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#34b18f")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#34b18f")

    for r in ob_r:
        ax.plot([-ob_l/2, ob_l/2], [r, r], color="#dd7e33")

    for z, rmin, rmax in zip(oe_z, oe_rmin, oe_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#dd7e33")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#dd7e33")

    ax.set_xlabel("z [mm]", fontsize=20)
    ax.set_ylabel("r [mm]", fontsize=20)

    ax.plot([], [], marker='o', label="Vertex Detector", linestyle=None, color="#348ec2")
    ax.plot([], [], marker='o', label="Inner Tracker",   linestyle=None, color="#34b18f")
    ax.plot([], [], marker='o', label="Outer Tracker",   linestyle=None, color="#dd7e33")
    ax.legend(loc="upper left", fontsize=14)

    ax.minorticks_on()
    ax.tick_params(labelsize=17)
    # ax.set_xlim(-2500, 2500)
    ax.set_xlim(0, 250)
    ax.set_ylim(0, 115)

    etas = [1.6, 1.8, 2.0, 2.2, 2.4, 2.6]
    get_eta_lines(ax, etas, zmax=ax.get_xlim()[1], rmax=ax.get_ylim()[1])

    pdf.savefig(fig)
    plt.close(fig)


    fig, ax = plt.subplots(figsize=(10,6))

    for r in vb_r:
        ax.plot([-vb_l/2, vb_l/2], [r, r], color="#348ec2")
        ax.plot([-vb_l/2, vb_l/2], [r+2, r+2], color="#348ec2")

    for z, rmin, rmax in zip(ve_z, ve_rmin, ve_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#348ec2")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#348ec2")

    for i, r in enumerate(ib_r):
        ax.plot([-ib_l[i]/2, ib_l[i]/2], [r, r], color="#34b18f")

    for z, rmin, rmax in zip(ie_z, ie_rmin, ie_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#34b18f")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#34b18f")

    for r in ob_r:
        ax.plot([-ob_l/2, ob_l/2], [r, r], color="#dd7e33")

    for z, rmin, rmax in zip(oe_z, oe_rmin, oe_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#dd7e33")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#dd7e33")

    ax.set_xlabel("z [mm]", fontsize=20)
    ax.set_ylabel("r [mm]", fontsize=20)

    ax.minorticks_on()
    ax.tick_params(labelsize=17)

    ax.set_xlim(0,2500)
    ax.set_ylim(0,1600)

    # THETA lines
    thetas = [5, 10, 20, 30, 45, 60, 75, 85]
    get_theta_lines(ax, thetas, zmax=ax.get_xlim()[1], rmax=ax.get_ylim()[1])

    pdf.savefig(fig)
    plt.close(fig)


    fig, ax = plt.subplots(figsize=(10,6))

    for r in vb_r:
        ax.plot([-vb_l/2, vb_l/2], [r, r], color="#348ec2")
        ax.plot([-vb_l/2, vb_l/2], [r+2, r+2], color="#348ec2")

    for z, rmin, rmax in zip(ve_z, ve_rmin, ve_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#348ec2")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#348ec2")

    for i, r in enumerate(ib_r):
        ax.plot([-ib_l[i]/2, ib_l[i]/2], [r, r], color="#34b18f")

    for z, rmin, rmax in zip(ie_z, ie_rmin, ie_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#34b18f")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#34b18f")

    for r in ob_r:
        ax.plot([-ob_l/2, ob_l/2], [r, r], color="#dd7e33")

    for z, rmin, rmax in zip(oe_z, oe_rmin, oe_rmax):
        ax.plot([ z,  z], [rmin, rmax], alpha=0.5, color="#dd7e33")
        ax.plot([-z, -z], [rmin, rmax], alpha=0.5, color="#dd7e33")

    ax.set_xlabel("z [mm]", fontsize=20)
    ax.set_ylabel("r [mm]", fontsize=20)

    ax.minorticks_on()
    ax.tick_params(labelsize=17)

    ax.set_xlim(0,250)
    ax.set_ylim(0,115)

    # THETA lines
    thetas = [6, 8, 10, 12, 14, 16, 18, 20, 30]
    get_theta_lines(ax, thetas, zmax=ax.get_xlim()[1], rmax=ax.get_ylim()[1])

    pdf.savefig(fig)
    plt.close(fig)
print(f"saved plots to {pdf_path}")