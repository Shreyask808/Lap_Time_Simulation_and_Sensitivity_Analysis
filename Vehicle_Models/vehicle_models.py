#=========================================================================================================================================================================================================================
# This project utilizes CasADi for symbolic differentiation and optimization.
# Andersson, J. A. E., et al. (2019). "CasADi: a software framework for 
# nonlinear optimization and optimal control." Mathematical Programming Computation.
#=========================================================================================================================================================================================================================

import os
import sys
import matplotlib
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
matplotlib.use('Qt5Agg')
import plotly.graph_objects as go
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Tire Models'))
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')

#=======================================================================================================================================================================================================================
# User Inputs
h = 250                                                                              # Number of Segments the Track is divided into
p = 3                                                                              # Degree of the Polynomial approximating the state in segments

print(f"#================================================================================================================================================================================================================")
print(f"")

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

# Track Segmentation
segment_length = track_data.arc_s[-1]/h
print(f"track length in m is: {track_data.arc_s[-1]} [m]")
print(f"Segment Length in m is: {segment_length} [m]")
print(f"")
print(f"#================================================================================================================================================================================================================")
print(f"#================================================================================================================================================================================================================")
print(f"")

#=====================================================================================================================================================================================================================
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
    print(f"Mass of the Vehicle: {vehicle.m} [kg]")
else:
    raise FileNotFoundError("Vehicle Data not found")
print(f"")
print(f"#================================================================================================================================================================================================================")

#=====================================================================================================================================================================================================================
# Import Point Mass Optimization
root = tk.Tk()
root.withdraw()
root.attributes('-topmost',True)

file_path3 = filedialog.askopenfilename(title="Select the Point Mass Optimization input",filetypes=[("pickle files","*.pkl"),])
if file_path3:
    with open(file_path3, 'rb') as F:  # 'rb' is mandatory for loading pickles
        x0_ini = pickle.load(F)
    print("Loaded Point Mass Optimization data")
else:
    raise FileNotFoundError("Point Mass Optimization data not found")


#=======================================================================================================================================================================================================================
#=======================================================================================================================================================================================================================
# Functions
## Flipped Legendre-Gauss-Radau Points
def flipped_LGR_points(N):

    nodes, weights = roots_jacobi(N,0,1)
    col_points = np.sort(np.concatenate([[-1],nodes]))
    return col_points,weights

col_p, weights = flipped_LGR_points(p)

## Lagrange Polynomial and its Derivative
def lagrange_polynomials(col_p,I,N):
    tau = SX.sym('tau')
    L_ki = SX(1)
    D_ki = SX(1)
    for j in range(N+1):
        if j != I:
            L_ki = L_ki*(tau - col_p[j])/(col_p[I] - col_p[j])
    
    D_ki = jacobian(L_ki,tau)

    return tau,L_ki,D_ki

print(f"#================================================================================================================================================================================================================")
print(f"")
print(f"Flipped LGR Points used in Analysis are: {col_p}")
print(f"")
print(f"#================================================================================================================================================================================================================")

#====================================================================================================================================================================================================================
# Results
#====================================================================================================================================================================================================================
# Plots 
## Figure 1 - Racing Line
fig = go.Figure()

# Optimized Right Boundary Data
fig.add_trace(go.Scatter3d(
    x=track_data.br[0,:], y=track_data.br[1,:], z=track_data.br[2,:],
    mode='lines',
    name='Right Boundary',
    line=dict(color='blue', width=4)
))

# Optimized Left Boundary Data
fig.add_trace(go.Scatter3d(
    x=track_data.bl[0,:], y=track_data.bl[1,:], z=track_data.bl[2,:], 
    mode='lines',
    name='Left Boundary',
    line=dict(color='blue', width=4)
))

# Optimized Centerline
fig.add_trace(go.Scatter3d(
    x=track_data.xc, y=track_data.yc, z=track_data.zc, 
    mode='lines',
    name='Centerline',
    line=dict(color='red', width=4, dash='dashdot')
))

# Optimal Racing Line Data
fig.add_trace(go.Scatter3d(
    x=x0_ini.x_opt, y=x0_ini.y_opt, z=x0_ini.z_opt,
    mode='lines',
    name='Optimized Racing Line',
    line=dict(color='black', width=7)
))

fig.update_layout(
    title="Autodromo Nazionale Monza - Optimal Racing Line",
    scene=dict(
        xaxis_title='x-coordinate (m)',
        yaxis_title='y-coordinate (m)',
        zaxis_title='Elevation (m)',
        aspectmode='data'),
    margin=dict(l=0, r=0, b=0, t=50)
)
fig.show()