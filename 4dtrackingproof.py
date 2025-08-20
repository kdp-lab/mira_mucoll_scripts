import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

layer_map = [
    ("VB", 0, 30), ("VB", 1, 30), ("VB", 2, 51), ("VB", 3, 51), ("VB", 4, 74), ("VB", 5, 74), ("VB", 6, 102), ("VB", 7, 102),
    ("IB", 0, 127), ("IB", 1, 340), ("IB", 2, 554), 
    ("OB", 0, 819), ("OB", 1, 1153), ("OB", 2, 1486)
]

layer_names = [f"{det}{lay}" for det, lay, _ in layer_map]
layer_radii = np.array([r for _,_,r in layer_map])
n_layers_tot = len(layer_map)
timing_windows = {"tight": (0.15, 0.3, 0.3), "medium": (0.32, 1.25, 4.0), "loose": (0.32, 3.30, 10.0)}

bib_color = "tab:red"
sig_color = "tab:green"

