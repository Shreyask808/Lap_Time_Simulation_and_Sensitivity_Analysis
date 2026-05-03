#====================================================================================================================================================================================
# This project utilizes CasADi for symbolic differentiation and optimization.
# Andersson, J. A. E., et al. (2019). "CasADi: a software framework for 
# nonlinear optimization and optimal control." Mathematical Programming Computation.
#====================================================================================================================================================================================

import os
import sys
import numpy as np
import pandas as pd
import tkinter as tk
import scipy.ndimage as ndimage
from scipy.ndimage import gaussian_filter1d
from tkinter import filedialog
import matplotlib.pyplot as plt
from types import SimpleNamespace
import json
import math
import casadi as ca
from casadi import *
import matplotlib
import pickle
matplotlib.use('Qt5Agg') # Or 'TkAgg' if you don't have Qt installed
import plotly.graph_objects as go
import matplotlib.pyplot as plt
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')

#=====================================================================================================================================================================================================================
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

#=====================================================================================================================================================================================================================
# Import Vehicle Data
root = tk.Tk()
root.withdraw()
root.attributes('-topmost',True)

file_path2 = filedialog.askopenfilename(title="Select the Vehicle Data",filetypes=[("json files","*.json")])
if file_path2:
    with open(file_path2, 'r') as f:
        vehicledata = json.load(f)
        vehicle = SimpleNamespace(vehicledata)
    print("Loaded Vehicle Data")
    print(f"Mass of the Vehicle: {vehicle.m}")
else:
    raise FileNotFoundError("Vehicle Data not found")

