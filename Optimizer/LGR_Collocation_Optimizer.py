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
h = 100                                                                              # Number of Segments the Track is divided into
p = 5                                                                               # Degree of the Polynomial approximating the state in segments

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
print(f"track length in m is: {track_data.arc_s[-1]} [m]")
print(f"Segment Length in m is: {segment_length} [m]")
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
    print(f"Mass of the Vehicle: {vehicle.m} [kg]")
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
    col_points = np.sort(np.concatenate([[-1],nodes]))
    return col_points

col_p = flipped_LGR_points(p)

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

#=======================================================================================================================================================================================================================
# Point mass optimization for intial conditions
#=======================================================================================================================================================================================================================
# Point Mass Dynamics Definition
## Symbolic Variable Definition
X_sym = SX.sym('X_sym',1,4)
U_sym = SX.sym('U_sym',1,2)
s = SX.sym('s')

## Interpolation Functions 
f_kappa = ca.interpolant('f_kappa','bspline',[track_data.arc_s.tolist()],track_data.omega_z.tolist())
f_nl = ca.interpolant('f_nl','linear',[track_data.arc_s.tolist()],track_data.nl.tolist())
f_nr = ca.interpolant('f_nr','linear',[track_data.arc_s.tolist()],track_data.nr.tolist())
f_theta = ca.interpolant('f_theta','linear',[track_data.arc_s.tolist()],track_data.theta.tolist())

## Scaling Factors
force_scale = 1/(vehicle.m*vehicle.g)                                                                               # Force Scale [N^-1]
length_scale = 1/vehicle.l                                                                                          # Length Scale [m^-1]
time_scale = np.sqrt(vehicle.g/vehicle.l)                                                                           # Time Scale [sec^-1]
speed_scale = length_scale/time_scale                                                                               # Speed Scale [sec/m]
angle_scale = 1                                                                                                     # Angle Scale [rad^-1]

### States Definition
t_N = X_sym[0]                                                                                                      # Normalized Time
n_N = X_sym[1]                                                                                                      # Normalized Lateral offset from track centerline
u_N = X_sym[2]                                                                                                      # Normalized Longitudinal velocity of the point mass
v_N = X_sym[3]                                                                                                      # Normalized Lateral velocity of the point mass

### Control Definition
F_d_N = U_sym[0]                                                                                                    # Normalized Drive Force of the point mass
F_n_N = U_sym[1]                                                                                                    # Normalized Centripital Force of the point mass

### States
t_sym = t_N/time_scale
n_sym = n_N/length_scale
u_sym = u_N/speed_scale
v_sym = v_N/speed_scale

### Controls
F_d_sym = F_d_N/force_scale
F_n_sym = F_n_N/force_scale

### Aerodynamic Forces
Drag = (1/2)*vehicle.Cd*vehicle.rho*vehicle.A*(u_sym**2)                                                            # Drag Force Definition [N]
Lift = (1/2)*vehicle.Cl*vehicle.rho*vehicle.A*(u_sym**2)
f_Drag = Function('f_Drag',[X_sym[2]],[Drag])
f_Lift = Function('f_Lift',[X_sym[2]],[Lift])
curv = f_kappa(s)

#### Dynamics
s_dot = (1/time_scale)*(u_sym)/(1 - n_sym*curv)
n_dot_N = (length_scale/time_scale)*v_sym
u_dot_N = (length_scale/time_scale**2)*((F_d_sym - Drag)/vehicle.m + v_sym*u_sym*curv/(1 - n_sym*curv))
v_dot_N = (length_scale/time_scale**2)*(F_n_sym/vehicle.m - (u_sym**2)*curv/(1 - n_sym*curv))

ODE = horzcat(1/s_dot,n_dot_N/s_dot,u_dot_N/s_dot,v_dot_N/s_dot)
f_dynamics = Function('f_dynamics',[X_sym,U_sym,s],[ODE])

# Global State Definition
X = SX.sym('X',h*(p+1),4)
U = SX.sym('U',h*p,2)
states = vertcat(reshape(X, -1, 1), reshape(U, -1, 1))

# Gauss Pseudospectral Differentiation Matrix Definition
D = np.zeros((p,p+1))
for i in range(p+1):
    tau,L_ki, D_ki = lagrange_polynomials(col_p,i,p)
    f_Dki = Function('f_Dki',[tau],[D_ki])
    for j in range(p):
        D[j,i] = f_Dki(col_p[j+1])

print(f"#================================================================================================================================================================================================================")
print(f"")
print(f"D shape: {D.shape}")
print(f"D:\n{D}")
print(f"")
print(f"#================================================================================================================================================================================================================")

# Nonlinear Program Solver (NLP) Definition for Point Mass
## Dynamics Constraints
g_dynamics = []
lbg_dynamics = []
ubg_dynamics = []

## Continuity Constraints
g_continuity = []
lbg_continuity = []
ubg_continuity = []

## Time Constraints
g_time = []
lbg_time = []
ubg_time = []

## Power Constraints
g_power = []
lbg_power = []
ubg_power = []

## Tire Force Limits
g_tire = []
lbg_tire = []
ubg_tire = []

## Main Loop
for k in range(h):
    F_dynamics = []
    X_segment = []
    U_segment = []
    s_segment = []
    
    s0 = k*segment_length
    sf = (k+1)*segment_length
    s_segment = ((sf - s0)/2)*col_p + ((s0+sf)/2)
    X_segment = X[k*(p+1):((k+1)*(p+1)),:]
    U_segment = U[k*p:(k+1)*p,:]
    
    if k < (h-1):
        g_continuity.append(X[(k+1)*(p+1),:] - X[(k+1)*(p+1)-1,:])
        lbg_continuity.append(np.zeros((1,4)))
        ubg_continuity.append(np.zeros((1,4)))
    
    for z in range(p):
        state = X_segment[z+1,:]
        ctrl = U_segment[z,:]

        # Dynamics Derivation
        F_dynamics = vertcat(F_dynamics,f_dynamics(state,ctrl,s_segment[z+1]))

        # Forward Time Constraint
        g_time.append(X_segment[z+1,0] - X_segment[z,0])
        lbg_time.append(1e-6)
        ubg_time.append(np.inf)

        # Vehicle Power Constraint
        g_power.append(state[2]*ctrl[0]/(force_scale*speed_scale))
        lbg_power.append(vehicle.peakbrakingpower)
        ubg_power.append(vehicle.peakdrivingpower)

        # Tire Force Constraint
        g_tire.append((ctrl[0]/force_scale)**2 + (ctrl[1]/force_scale)**2 - (vehicle.mu0*(vehicle.m*vehicle.g - f_Lift[state[2]]))**2)
        lbg_tire.append(-np.inf)
        ubg_tire.append(0)
           
    g_dynamics.append(ca.mtimes(ca.DM(D), X_segment) - F_dynamics)
    lbg_dynamics.append(np.zeros(F_dynamics.shape))
    ubg_dynamics.append(np.zeros(F_dynamics.shape))
    

print(f"Shape of g_dynamics:{len(g_dynamics)}")
print(f"Shape of lbg_dynamics:{len(lbg_dynamics)}")
print(f"Shape of ubg_dynamics:{len(ubg_dynamics)}")