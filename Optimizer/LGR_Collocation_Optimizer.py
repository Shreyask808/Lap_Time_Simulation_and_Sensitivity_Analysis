#====================================================================================================================================================================================
# Legendre-Gauss-Radau (LGR) Pseudospectral Optimal Control
# Reference: Limebeer & Perantoni (2015), JDSMC 137, DOI: 10.1115/1.4029466
#====================================================================================================================================================================================

import os
import sys
import numpy as np
import pandas as pd
import tkinter as tk
import scipy.ndimage as ndimage
from scipy.ndimage import gaussian_filter1d
from scipy.special import roots_jacobi
from tkinter import filedialog
import matplotlib.pyplot as plt
from types import SimpleNamespace
import json
import pickle
import casadi as ca
from casadi import *
import matplotlib
matplotlib.use('Qt5Agg')
import plotly.graph_objects as go
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Tire Models'))
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')

#=======================================================================================================================================================================================================================
# User Inputs
h = 60                                                                                                          # Number of Segments the Track is divided into
p = 5                                                                                                          # Degree of the Polynomial approximating the state in segments

print(f"#================================================================================================================================================================================================================")
print(f"")
# Import Track Data
root = tk.Tk()
root.withdraw()
root.attributes('-topmost',True)

file_path1 = filedialog.askopenfilename(title="Select the Track Data",filetypes=[("pickle files","*.pkl"),])
if file_path1:
    with open(file_path1, 'rb') as f:  # 'rb' is mandatory for loading pickles
        track_data = pickle.load(f)
    print("Loaded Track Data")
    print(f"Number of Data Points on the Track: {len(track_data.arc_s)}")
else:
    raise FileNotFoundError("Track Data file not found")

# Track Segmentation
segment_length = track_data.arc_s[-1]/h
print(f"track length in m is: {track_data.arc_s[-1]}")
print(f"Segment Length in m is: {segment_length}")
print(f"")
print(f"#================================================================================================================================================================================================================")
print(f"#================================================================================================================================================================================================================")
print(f"")
# Import Vehicle Data
root = tk.Tk()
root.withdraw()
root.attributes('-topmost',True)

file_path2 = filedialog.askopenfilename(title="Select the Vehicle Data",filetypes=[("json files","*.json")])
if file_path2:
    with open(file_path2, 'r') as f:
        vehicledata = json.load(f)
        vehicle = SimpleNamespace(**vehicledata)
    print("Loaded Vehicle Data")
    print(f"Mass of the Vehicle: {vehicle.m}")
else:
    raise FileNotFoundError("Vehicle Data not found")
print(f"")
print(f"#================================================================================================================================================================================================================")

#=======================================================================================================================================================================================================================
#=======================================================================================================================================================================================================================
# Functions
## Flipped Legendre-Gauss-Radau Points
def flipped_LGR_points(N):

    nodes, weights = roots_jacobi(N,1,0)
    col_points = np.sort(np.concatenate([nodes,[1.0]]))
    return col_points

col_p = flipped_LGR_points(p)

print(f"#================================================================================================================================================================================================================")
print(f"")
print(f"Flipped LGR Points used in Analysis are: {col_p}")
print(f"")
print(f"#================================================================================================================================================================================================================")
