import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages 

pdf_path = "/scratch/miralittmann/analysis/mira_analysis_code/mucoll_tracker.pdf"

vb_r = [30, 51, 74, 102]
vb_l = 130
ve_z = [80, 120, 200, 280]
ve_h = [0, 15, 35, 50]

ib_r = [127, 340, 554]
ib_l = [963.2, 963.2, 1384.6]
ie_z = [524, 808, 1093, 1377, 1661, 1946, 2190]
ie_h = np.arange(100,300, (300-100)/7) 

ob_r = [819, 1153, 1486]
ob_l = 2528.4
oe_z = [1310, 1617, 1883, 2190]

etas = [0.5, 1.0, 2.0, 3.0]

def get_eta_lines(ax, etas, zmax, rmax, color='gray', ls='--', alpha=0.8, lw=1.2, pad=60):
    z_vals = np.linspace(0, zmax, 400)
    for eta in etas:
        theta = 2*np.arctan(np.exp(-eta))
        r_vals = np.tan(theta) * z_vals

        ax.plot( z_vals, r_vals, ls=ls, color=color, alpha=alpha, lw=lw)
        ax.plot(-z_vals, r_vals, ls=ls, color=color, alpha=alpha, lw=lw)
        z_top = rmax / np.tan(theta)
        if z_top <= zmax:
            # hits top edge
            x_edge, y_edge = z_top, rmax
            # label just above the top edge (in points)
            ax.annotate(f"$\eta$={eta:.1f}", xy=(x_edge, y_edge),
                        xytext=(0, 4), textcoords='offset points',
                        rotation=0,color="gray", fontsize=14)
        else:
            # hits right edge: z = zmax => r_right = m * zmax
            r_right = np.tan(theta) * zmax
            # guard against tiny numeric overshoot
            r_right = min(r_right, rmax)
            x_edge, y_edge = zmax, r_right
            # label just to the right (in points)
            ax.annotate(f"$\eta$={eta:.1f}", xy=(x_edge, y_edge),
                        xytext=(4, 0), textcoords='offset points',
                        ha='left', va='center', rotation=0, color="gray", fontsize=14
                        ) 


with PdfPages(pdf_path) as pdf:
    fig,ax = plt.subplots(figsize=(10,6))
    for layer in vb_r:
        ax.plot([-vb_l/2, vb_l/2], [layer, layer], color="#348ec2")
        ax.plot([-vb_l/2, vb_l/2], [layer+2, layer+2], color="#348ec2")

    for i,layer in enumerate(ve_z):
        ax.plot([layer, layer], [ve_h[i], 125], alpha=0.5, color="#348ec2")
        ax.plot([-layer, -layer], [ve_h[i], 125], alpha=0.5, color="#348ec2")
        ax.plot([layer+2, layer+2], [ve_h[i], 125], alpha=0.5, color="#348ec2")
        ax.plot([-layer+2, -layer+2], [ve_h[i], 125], alpha=0.5, color="#348ec2")

    for i,layer in enumerate(ib_r):
        ax.plot([-ib_l[i]/2, ib_l[i]/2], [layer, layer], color="#34b18f")

    for i,layer in enumerate(ie_z):
        if i == 0 :
            ax.plot([layer, layer], [ie_h[i], 410],  alpha=0.5, color="#34b18f")
            ax.plot([-layer, -layer], [ie_h[i], 410], alpha=0.5, color="#34b18f")
        else:
            ax.plot([layer, layer], [ie_h[i], 570], alpha=0.5, color="#34b18f")
            ax.plot([-layer, -layer], [ie_h[i], 570], alpha=0.5, color="#34b18f")
   
    for layer in ob_r:
        ax.plot([-ob_l/2, ob_l/2], [layer, layer], color="#dd7e33")
    for i,layer in enumerate(oe_z):
        ax.plot([layer,layer], [600, 1430], alpha= 0.5, color="#dd7e33")    
        ax.plot([-layer,-layer], [600, 1430], alpha= 0.5, color="#dd7e33")

       
    ax.set_xlabel("z [mm]", fontsize=20)
    ax.set_ylabel("r [mm]", fontsize=20)

    ax.plot([], [], marker='o',label="Vertex Detector", linestyle=None, color="#348ec2")
    ax.plot([], [], marker='o', label="Inner Tracker", linestyle=None, color="#34b18f")
    ax.plot([], [], marker='o', label="Outer Tracker", linestyle=None, color="#dd7e33")
    ax.legend(loc="upper left", fontsize=14)
    ax.minorticks_on()
    ax.tick_params(labelsize=17)
    ax.set_xlim(-2500, 2500)
    ax.set_ylim(0,1600)

    etas = [0.5, 1.0, 1.5, 2.0, 3.0]
    get_eta_lines(ax, etas, zmax=ax.get_xlim()[1], rmax=ax.get_ylim()[1])


    pdf.savefig(fig)
    plt.close(fig)

print(f"saved plots to {pdf_path}") 

     