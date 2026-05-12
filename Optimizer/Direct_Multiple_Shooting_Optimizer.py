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
# User Inputs
N = 4000                                                                                                    # Number of Segments on the track          

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

##====================================================================================================================================================================================================================
# Point Mass Model Optimization [to generate initial guess for the vehicle model]
## System Model Definition
lap_length = track_data.arc_s[-1]                                                                           # Total Lap Length along the track centerline in m                         
s_grid = np.linspace(0,lap_length,N+1)                                                                      # Discretization of centerline track length in m
ds = lap_length/N                                                                                           # Step length in m

## Optimization States Definition
X = SX.sym('X',4,N+1)
U = SX.sym('U',2,N)
states = vertcat(reshape(X, -1, 1), reshape(U, -1, 1))
print(f"Length of Optimal Control State Vector is :{states.numel()}")

## Interpolation Functions 
f_kappa = ca.interpolant('f_kappa','linear',[track_data.arc_s.tolist()],track_data.omega_z.tolist())
f_nl = ca.interpolant('f_nl','linear',[track_data.arc_s.tolist()],track_data.nl.tolist())
f_nr = ca.interpolant('f_nr','linear',[track_data.arc_s.tolist()],track_data.nr.tolist())
f_theta = ca.interpolant('f_theta','linear',[track_data.arc_s.tolist()],track_data.theta.tolist())

## System Definition
X_sym = SX.sym('X_sym',4)
U_sym = SX.sym('U_sym',2)
s = SX.sym('s')

### Scaling X[.] = X_scale [1/units]*X [Units]
force_scale = 1/(vehicle.m*vehicle.g)                                                                       # Force Scale in N^-1
length_scale = 1/vehicle.l                                                                                  # Length Scale in m^-1
time_scale = np.sqrt(vehicle.g/vehicle.l)                                                                   # Time Scale in sec^-1
speed_scale = 1/np.sqrt(vehicle.g*vehicle.l)                                                                # Speed Scale in sec/m
angle_scale = 1                                                                                             # Angle Scale in rad^-1

#### States
t_N = X_sym[0]  
n_N = X_sym[1]
n_dot_N = X_sym[2]
sigma_N = X_sym[3]

#### Controls
n_ddot_N = U_sym[0]
F_d_N = U_sym[1]

### Dynamics Definition
t = t_N/time_scale
n = n_N/length_scale
n_dot = n_dot_N*time_scale/length_scale
sigma = sigma_N/speed_scale
n_ddot = n_ddot_N*time_scale/speed_scale
F_d = F_d_N/force_scale

#### Aerodynamic Forces
Drag = (1/2)*vehicle.Cd*vehicle.rho*vehicle.A*(sigma**2)
Lift = (1/2)*vehicle.Cl*vehicle.rho*vehicle.A*(sigma**2)
f_Drag = Function('f_Drag',[X_sym],[Drag])
f_Lift = Function('f_Lift',[X_sym],[Lift])
curv = f_kappa(s)

#### Dynamics
s_dot = (1/time_scale)*sigma/(1 - n*curv)
n_ddot = (speed_scale/time_scale)*curv*(sigma**2)/(1 - n*curv)
sigma_dot = (speed_scale/time_scale)*(F_d - f_Drag(s))/vehicle.m

ODE = vertcat(1/s_dot, n_dot_N/s_dot, n_ddot/s_dot, sigma_dot/s_dot)
f_dynamics = Function('f_dynamics',[X_sym,U_sym,s],[ODE])

### Cost Function & Constraints
#### Cost Function
cost = 0
e1 = 1e-4

#### Constraints
g_dynamics = []
lbg_dynamics = []
ubg_dynamics = []

g_time = []
lbg_time = []
ubg_time = []

g_power = []
lbg_power = []
ubg_power = []

g_tire = []
lbg_tire = []
ubg_tire = []

g_end = []
lbg_end = []
ubg_end = []

### Nonlinear Program Solver (NLP) Definition for Point Mass
for k in range(N):
    s_current = s_grid[k]
    state = X[:,k]
    ctrl = U[:,k]

    #### 4th Order Ranga Kutta Method for Discretization
    k1 = f_dynamics(state,ctrl,s_current)
    k2 = f_dynamics(state + (ds/2)*k1,ctrl,s_current)
    k3 = f_dynamics(state + (ds/2)*k2, ctrl,s_current)
    k4 = f_dynamics(state + (ds)*k3, ctrl,s_current)

    #### Dynamics Constraints
    state_next = state + (ds/6)*(k1 + 2*k2 + 2*k3 + k4)
    g_dynamics.append(X[:,k+1] - state_next)
    lbg_dynamics.append(np.zeros((4,1)))
    ubg_dynamics.append(np.zeros((4,1)))

    #### Time Constraints
    g_time.append(state_next[0] - state[0])
    lbg_time.append(1e-6)
    ubg_time.append(np.inf)

    #### Vehicle Power Constraints
    g_power.append(ctrl[1]*state[3])
    lbg_power.append(vehicle.peakbrakingpower)
    ubg_power.append(vehicle.peakdrivingpower)

    #### Tire Force Constraints
    g_tire.append(((ctrl[1])**2 + (vehicle.m*ctrl[0])**2) - (vehicle.mu0*(vehicle.m*vehicle.g - f_Lift(state[3])))**2)
    lbg_tire.append(-np.inf)
    ubg_tire.append(0)

    #### Cost Function
    dt = (state_next[0] - state[0])/time_scale
    cost = cost + dt*e1*(ctrl[1]**2)

cost = cost + X[0,-1]/time_scale

g_end.append(X[[1,3],-1] - X[[1,3],0])
lbg_end.append(np.zeros((2,1)))
ubg_end.append(np.zeros((2,1)))

g = vertcat(*g_dynamics, *g_time, *g_power, *g_tire, *g_end)
lbg = vertcat(*lbg_dynamics, *lbg_time, *lbg_power, *lbg_tire, *lbg_end)
ubg = vertcat(*ubg_dynamics, *ubg_time, *ubg_power, *ubg_tire, *ubg_end)

### States and Control Input Limits
lbx = np.full(states.shape, -np.inf)
ubx = np.full(states.shape, np.inf)
x0 = np.zeros(states.numel())

for r in range(N+1):
    idx_t = r*4
    idx_n = r*4 + 1
    idx_n_dot = r*4 + 2
    idx_sigma = r*4 + 3

    if r == 0:
        lbx[idx_t] = 0
        ubx[idx_t] = 0
    else:
        lbx[idx_t] = 0
        ubx[idx_t] = 500*time_scale

    lbx[idx_n] = f_nr(s_grid[r])*length_scale
    ubx[idx_n] = f_nl(s_grid[r])*length_scale

    lbx[idx_sigma] = 2*speed_scale
    ubx[idx_sigma] = vehicle.umax*speed_scale

    x0[idx_t] = time_scale*(s_grid[r]/vehicle.umax)
    x0[idx_sigma] = vehicle.umax*speed_scale

    if r < N:
        idx_n_ddot = 4*(N+1) + 2*r
        idx_F_d = 4*(N+1) + 2*r + 1
        lbx[idx_F_d] = (vehicle.peakbrakingtorque/vehicle.R)*force_scale
        ubx[idx_F_d] = (vehicle.peakdrivingtorque/vehicle.R)*force_scale
        x0[idx_F_d] = (vehicle.peakdrivingtorque/vehicle.R)*force_scale

### Nonlinear Solver for Point Mass Model
nlp = {'x': states,'f': cost,'g': g}
opts = {'ipopt': {'max_iter': 2000,'print_level': 5}}

solver = nlpsol('solver', 'ipopt', nlp, opts)
sol = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
full_sol = np.array(sol['x']).flatten()
n_x_total =  4*(N + 1)
X_opt = full_sol[:n_x_total].reshape((4, N + 1), order='F')
U_opt = full_sol[n_x_total:].reshape((2, N), order='F')
print(f"Shape of X_opt is: {X_opt.shape}")
print(f"Shape of U_opt is: {U_opt.shape}")