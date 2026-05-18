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
N = 2000                                                                                                    # Number of Segments on the track          

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
        vehicle = SimpleNamespace(**vehicledata)
    print("Loaded Vehicle Data")
    print(f"Mass of the Vehicle: {vehicle.m}")
else:
    raise FileNotFoundError("Vehicle Data not found")

#====================================================================================================================================================================================================================
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
f_kappa = ca.interpolant('f_kappa','bspline',[track_data.arc_s.tolist()],track_data.omega_z.tolist())
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
u_N = X_sym[2]
v_N = X_sym[3]

#### Controls
F_d_N = U_sym[0]
F_n_N = U_sym[1]

### Dynamics Definition
t = t_N/time_scale
n = n_N/length_scale
u = u_N/speed_scale
v = v_N/speed_scale

F_d = F_d_N/force_scale
F_n = F_n_N/force_scale

#### Aerodynamic Forces
Drag = (1/2)*vehicle.Cd*vehicle.rho*vehicle.A*(u**2)
Lift = (1/2)*vehicle.Cl*vehicle.rho*vehicle.A*(u**2)
f_Drag = Function('f_Drag',[X_sym],[Drag])
f_Lift = Function('f_Lift',[X_sym],[Lift])
curv = f_kappa(s)

#### Dynamics
s_dot = (1/time_scale)*(u)/(1 - n*curv)
n_dot_N = (length_scale/time_scale)*v
u_dot_N = (length_scale/time_scale**2)*((F_d - Drag)/vehicle.m + v*u*curv/(1 - n*curv))
v_dot_N = (length_scale/time_scale**2)*(F_n/vehicle.m - (u**2)*curv/(1 - n*curv))

ODE = vertcat(1/s_dot,n_dot_N/s_dot,u_dot_N/s_dot,v_dot_N/s_dot)
f_dynamics = Function('f_dynamics',[X_sym,U_sym,s],[ODE])

### Cost Function & Constraints
#### Cost Function
cost = 0
e1 = 0
e2 = 0
e3 = 0
e4 = 0

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
### Nonlinear Program Solver (NLP) Definition for Point Mass
for k in range(N):
    s_current = s_grid[k]
    state = X[:,k]
    ctrl = U[:,k]

    #### Hermite-Simpson Collocation for Discretization
    f0 = f_dynamics(state,        ctrl, s_current)
    f1 = f_dynamics(X[:,k+1],     ctrl, s_grid[k+1])

    # Hermite interpolation for midpoint state
    X_mid = 0.5*(state + X[:,k+1]) + (ds/8)*(f0 - f1)

    # Dynamics at midpoint
    f_m = f_dynamics(X_mid, ctrl, s_current + ds/2)

    #### Dynamics Constraints — Simpson quadrature
    state_next = state + (ds/6)*(f0 + 4*f_m + f1)
    g_dynamics.append(X[:,k+1] - state_next)
    lbg_dynamics.append(np.zeros((4,1)))
    ubg_dynamics.append(np.zeros((4,1)))

    #### Time Constraints
    g_time.append(X[0,k+1] - X[0,k])
    lbg_time.append(1e-6)
    ubg_time.append(np.inf)

    #### Vehicle Power Constraints
    g_power.append(ctrl[0]*state[2]/(force_scale*speed_scale))
    lbg_power.append(vehicle.peakbrakingpower)
    ubg_power.append(vehicle.peakdrivingpower)

    #### Tire Force Constraints
    g_tire.append((ctrl[0]/force_scale)**2 + (ctrl[1]/force_scale)**2 - (vehicle.mu0*(vehicle.m*vehicle.g - f_Lift(state)))**2)
    lbg_tire.append(-np.inf)
    ubg_tire.append(0)

    #### Cost Function
    dt = (X[0,k+1] - X[0,k])/time_scale
    cost = cost + dt*(e1*((ctrl[0])**2) + e2*((ctrl[1])**2) + e3*((state[1])**2) + e4*((state[3])**2))
cost = cost + X[0,-1]/time_scale

g_end.append(X[[1,2,3],-1] - X[[1,2,3],0])
lbg_end.append(np.zeros((3,1)))
ubg_end.append(np.zeros((3,1)))

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
    idx_u = r*4 + 2
    idx_v = r*4 + 3

    if r == 0:
        lbx[idx_t] = 0
        ubx[idx_t] = 0
    else:
        lbx[idx_t] = 0
        ubx[idx_t] = 500*time_scale

    lbx[idx_n] = f_nr(s_grid[r])*length_scale
    ubx[idx_n] = f_nl(s_grid[r])*length_scale

    lbx[idx_u] = 2*speed_scale
    ubx[idx_u] = vehicle.umax*speed_scale

    lbx[idx_v] = -vehicle.vmax*speed_scale
    ubx[idx_v] = vehicle.vmax*speed_scale

    x0[idx_t] = time_scale*(s_grid[r]/vehicle.umax)
    x0[idx_u] = vehicle.umax*speed_scale

    if r < N:
        idx_F_d = 4*(N+1) + 2*r
        idx_F_n = 4*(N+1) + 2*r + 1

        lbx[idx_F_d] = (4*vehicle.peakbrakingtorque/vehicle.R)*force_scale
        ubx[idx_F_d] = (2*vehicle.peakdrivingtorque/vehicle.R)*force_scale
        x0[idx_F_d] = (2*vehicle.peakdrivingtorque/vehicle.R)*force_scale

        lbx[idx_F_n] = -(vehicle.mu0*vehicle.m*vehicle.g - (0.5*vehicle.Cl)*vehicle.rho*vehicle.A*vehicle.umax**2)*force_scale
        ubx[idx_F_n] = (vehicle.mu0*vehicle.m*vehicle.g - (0.5*vehicle.Cl)*vehicle.rho*vehicle.A*vehicle.umax**2)*force_scale
        x0[idx_F_n]  = 0

### Nonlinear Solver for Point Mass Model
nlp = {'x': states,'f': cost,'g': g}
opts = {'ipopt': {'max_iter': 4000,'print_level': 5,'mu_strategy': 'adaptive','acceptable_obj_change_tol': 5e-2}}
solver = nlpsol('solver', 'ipopt', nlp, opts)
sol = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
full_sol = np.array(sol['x']).flatten()
n_x_total =  4*(N + 1)
X_opt = full_sol[:n_x_total].reshape((4, N + 1), order='F')
U_opt = full_sol[n_x_total:].reshape((2, N), order='F')
print(f"Shape of X_opt is: {X_opt.shape}")
print(f"Shape of U_opt is: {U_opt.shape}")

time_opt = X_opt[0,:]/time_scale
n_opt = X_opt[1,:]/length_scale
u_opt = X_opt[2,:]/speed_scale
v_opt = X_opt[3,:]/speed_scale

F_d_opt = U_opt[0,:]/force_scale
F_n_opt = U_opt[1,:]/force_scale

Power_opt = F_d_opt*u_opt[:-1]

print(f"Optimized Lap Time is: {time_opt[-1]}")

#====================================================================================================================================================================================================================
# Frenet Coordinates to Cartesian Coordinates
f_normal1 = ca.interpolant('f_normal1','linear',[track_data.arc_s.tolist()],track_data.n[0,:].tolist())
f_normal2 = ca.interpolant('f_normal2','linear',[track_data.arc_s.tolist()],track_data.n[1,:].tolist())
f_normal3 = ca.interpolant('f_normal3','linear',[track_data.arc_s.tolist()],track_data.n[2,:].tolist())

f_xc = ca.interpolant('f_xc','linear',[track_data.arc_s.tolist()],track_data.xc.tolist())
f_yc = ca.interpolant('f_yc','linear',[track_data.arc_s.tolist()],track_data.yc.tolist())
f_zc = ca.interpolant('f_zc','linear',[track_data.arc_s.tolist()],track_data.zc.tolist())

x_opt = np.array(f_xc(s_grid)).flatten() + n_opt * np.array(f_normal1(s_grid)).flatten()
y_opt = np.array(f_yc(s_grid)).flatten() + n_opt * np.array(f_normal2(s_grid)).flatten()
z_opt = np.array(f_zc(s_grid)).flatten() + n_opt * np.array(f_normal3(s_grid)).flatten()

#%%
#====================================================================================================================================================================================================================
# Results and Plots 
## Figure 1 - Racing Line
fig = go.Figure()

# Optimized Right Boundary Data
fig.add_trace(go.Scatter3d(
    x=track_data.br[0,:], y=track_data.br[1,:], z=track_data.br[2,:],
    mode='lines',
    name='Right Boundary - Optimized',
    line=dict(color='blue', width=4)
))

# Optimized Left Boundary Data
fig.add_trace(go.Scatter3d(
    x=track_data.bl[0,:], y=track_data.bl[1,:], z=track_data.bl[2,:], 
    mode='lines',
    name='Left Boundary - Optimized',
    line=dict(color='blue', width=4)
))

# Optimized Centerline
fig.add_trace(go.Scatter3d(
    x=track_data.xc, y=track_data.yc, z=track_data.zc, 
    mode='lines',
    name='Centerline Coordinates - Optimized',
    line=dict(color='red', width=4, dash='dashdot')
))

# Optimal Racing Line Data
fig.add_trace(go.Scatter3d(
    x=x_opt, y=y_opt, z=z_opt,
    mode='lines',
    name='Optimized Racing Line',
    line=dict(color='black', width=4)
))

fig.update_layout(
    title="Circuit de Barcelona-Catalunya - 3D Track Geometry",
    scene=dict(
        xaxis_title='x-coordinate (m)',
        yaxis_title='y-coordinate (m)',
        zaxis_title='Elevation (m)',
        aspectmode='data'),
    margin=dict(l=0, r=0, b=0, t=50)
)
fig.show()

# Figure 2 - Lateral and Drive Forces
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Drive Force
ax1.plot(s_grid[0:-1],F_d_opt,color='black')
ax1.set_xlabel('Track Centerline Arc Length [m]')
ax1.set_ylabel('Drive Force [N]')
ax1.grid(True, alpha=0.3)
ax1.axhline(0,color='red')

# Lateral Force
ax2.plot(s_grid[0:-1],F_n_opt,color='black')
ax2.set_xlabel('Track Centerline Arc Length [m]')
ax2.set_ylabel('Lateral Force [N]')
ax2.grid(True, alpha=0.3)
ax2.axhline(0,color='red')

plt.tight_layout()
plt.show()

# Figure 3 - Lateral Offset, Longitudinal and Lateral Velocities
fig, (ax1,ax2,ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# Lateral Offset
ax1.plot(s_grid,n_opt,color='black')
ax1.plot(track_data.arc_s,track_data.nl,color='blue',label='Left Boundary')
ax1.plot(track_data.arc_s,track_data.nr,color='red',label='Right Boundary')
ax1.set_xlabel('Track Centerline Arc Length [m]')
ax1.set_ylabel('Lateral Offset [m]')
ax1.grid(True, alpha=0.3)
ax1.axhline(0,color='red')
ax1.legend()

# Longitudinal Speed
ax2.plot(s_grid,u_opt*(18/5),color='black')
ax2.set_xlabel('Track Centerline Arc Length [m]')
ax2.set_ylabel('Longitudinal Speed [kmph]')
ax2.grid(True, alpha=0.3)
ax2.axhline(0,color='red')
ax2.axhline(vehicle.umax*(18/5),color='red')

# Lateral Speed
ax3.plot(s_grid,v_opt*(18/5),color='black')
ax3.set_xlabel('Track Centerline Arc Length [m]')
ax3.set_ylabel('Lateral Speed [kmph]')
ax3.grid(True, alpha=0.3)
ax3.axhline(0,color='red')
ax3.axhline(vehicle.vmax*(18/5),color='red')
ax3.axhline(-vehicle.vmax*(18/5),color='red')

plt.tight_layout()
plt.show()


fig, (ax1) = plt.subplots(1, 1, figsize=(10, 12), sharex=True)

# Lateral Offset
ax1.plot(s_grid[:-1],Power_opt/1000,color='black')
ax1.set_xlabel('Track Centerline Arc Length [m]')
ax1.set_ylabel('Vehicle Power [kW]')
ax1.grid(True, alpha=0.3)
ax1.axhline(0,color='red')
ax1.legend()

plt.tight_layout()
plt.show()
