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

#====================================================================================================================================================================================
# Track Direction (Clockwise / Anti-clockwise)
direction = input("Is the track clockwise ? (Y/N): ").strip().upper()

if direction == 'Y':
    direction_angle = -2*math.pi
elif direction == 'N':
    direction_angle = 2*math.pi
else:
    raise ValueError("Invalid Input. Please Enter Y/N")

#====================================================================================================================================================================================
# Import Track GPS Data
# Right Boundary Data
root = tk.Tk()
root.withdraw()
root.attributes('-topmost',True)

file_path1 = filedialog.askopenfilename(title="Select the Right Boundary Data",filetypes=[("text files","*.txt"),("CSV files","*.csv")])
if file_path1:
    right_data = np.genfromtxt(file_path1,skip_header=2)
    right = SimpleNamespace()
    right.lat = ndimage.gaussian_filter1d(right_data[:,1]*np.pi/180, sigma=2.5)
    right.lon = ndimage.gaussian_filter1d(right_data[:,2]*np.pi/180, sigma=2.5)
    right.height = ndimage.gaussian_filter1d(right_data[:,3], sigma=10)
    num_points_r = len(right.lat)
    Xr = np.zeros(num_points_r)
    Yr = np.zeros(num_points_r)
    Zr = np.zeros(num_points_r)
    print("Loaded Right Boundary Data")
    print(f"Number of Right Boundary Data Points: {len(right_data)}")
else:
    raise FileNotFoundError("Right boundary file not found")

# Left Boundary Data
root = tk.Tk()
root.withdraw()
root.attributes('-topmost',True)

file_path2 = filedialog.askopenfilename(title="Select the Left Boundary Data",filetypes=[("text files","*.txt"),("CSV files","*.csv")])
if file_path2:
    left_data = np.genfromtxt(file_path2,skip_header=2)
    left = SimpleNamespace()
    left.lat = ndimage.gaussian_filter1d(left_data[:,1]*np.pi/180, sigma=2.5)
    left.lon = ndimage.gaussian_filter1d(left_data[:,2]*np.pi/180, sigma=2.5)
    left.height = ndimage.gaussian_filter1d(left_data[:,3], sigma=2.5)
    num_points_l = len(left.lat)
    Xl = np.zeros(num_points_l)
    Yl = np.zeros(num_points_l)
    Zl = np.zeros(num_points_l)
    print("Loaded Left Boundary Data")
    print(f"Number of Left Boundary Data Points: {len(left_data)}")
else:
    raise FileNotFoundError("Left boundary file not found")

#===================================================================================================================================================================================
# Conversion of GPS data to Cartisian Coordinates
origin = SimpleNamespace()
origin.lat = ((left.lat[0] + right.lat[0])/2)
origin.lon = ((left.lon[0] + right.lon[0])/2)
origin.height = ((left.height[0] + right.height[0])/2)

# Earth Dimensions as an Ellipse
a = 6378137
f = 1/298.257223563
e_sq = f*(2-f)

N_o = a/np.sqrt(1 - e_sq*(np.sin(origin.lat)**2))
Xo = (N_o + origin.height)*np.cos(origin.lat)*np.cos(origin.lon)
Yo = (N_o + origin.height)*np.cos(origin.lat)*np.sin(origin.lon)
Zo = (N_o*(1 - e_sq) + origin.height)*np.sin(origin.lat)
sin_lat = np.sin(origin.lat)
cos_lat = np.cos(origin.lat)
sin_lon = np.sin(origin.lon)
cos_lon = np.cos(origin.lon)

# GPS data conversion to Cartisian coordinates
for i in range(len(right.lat)):
    N = a/np.sqrt(1 - e_sq*(np.sin(right.lat[i]))**2)
    xr = (N + right.height[i])*np.cos(right.lat[i])*np.cos(right.lon[i]) - Xo
    yr = (N + right.height[i])*np.cos(right.lat[i])*np.sin(right.lon[i]) - Yo
    zr = (N*(1 - e_sq) + right.height[i])*np.sin(right.lat[i]) - Zo
    Xr[i] = -sin_lon*xr + cos_lon*yr
    Yr[i] = -sin_lat*cos_lon*xr + -sin_lat*sin_lon*yr + cos_lat*zr
    Zr[i] = cos_lat*cos_lon*xr + cos_lat*sin_lon*yr + sin_lat*zr

Xr = np.append(Xr,Xr[0])
Yr = np.append(Yr,Yr[0])
Zr = np.append(Zr,Zr[0])

for j in range(len(left.lat)):
    N = a/np.sqrt(1 - e_sq*(np.sin(left.lat[j]))**2)
    xl = (N + left.height[j])*np.cos(left.lat[j])*np.cos(left.lon[j]) - Xo
    yl = (N + left.height[j])*np.cos(left.lat[j])*np.sin(left.lon[j]) - Yo
    zl = (N*(1 - e_sq) + left.height[j])*np.sin(left.lat[j]) - Zo
    Xl[j] = -sin_lon*xl + cos_lon*yl
    Yl[j] = -sin_lat*cos_lon*xl + -sin_lat*sin_lon*yl + cos_lat*zl
    Zl[j] = cos_lat*cos_lon*xl + cos_lat*sin_lon*yl + sin_lat*zl
Xl = np.append(Xl,Xl[0])
Yl = np.append(Yl,Yl[0])
Zl = np.append(Zl,Zl[0])

#=======================================================================================================================================================================================================================
## Frenet Coordinate Estimation
# Centerline Coordinate Estimation
if len(right.lat) >= len(left.lat):
    coord1 = np.array([Xr, Yr, Zr])
    coord2 = np.array([Xl, Yl, Zl])
else:
    coord1 = np.array([Xl, Yl, Zl])
    coord2 = np.array([Xr, Yr, Zr])
dist = np.zeros(coord1.shape[1])

centerline = SimpleNamespace()
centerline.coord = np.zeros(coord2.shape)
new_coord = np.zeros(coord2.shape)

for i in range(coord2.shape[1]):
    dist = np.linalg.norm(coord1 - coord2[:,i].reshape(-1,1), axis=0)
    idx_min = np.argmin(dist)
    centerline.coord[:,i] = (coord2[:,i] + coord1[:,idx_min])/2
    new_coord[:,i] = coord1[:,idx_min]

print(f"Shape of centerline.coord: {centerline.coord.shape}")

#=====================================================================================================================================================================================================================
# Centerline Arc Length Estimation
arc_s = np.zeros(coord2.shape[1])
for k in range(1,coord2.shape[1]):
    arc_dist = centerline.coord[:,k] - centerline.coord[:,k-1]
    arc_s[k] = arc_s[k-1] + np.linalg.norm(arc_dist, axis=0)
lap_length = arc_s[-1]
print(f"Total length of the lap in km is:{lap_length/1000}")
print(f"Length of arc_s: {len(arc_s)}")

#=====================================================================================================================================================================================================================
# Direct Multiple Shooting Optimal Control Track Optimization [CASADI]
N = 4000                                                                                            # Step Length

#=====================================================================================================================================================================================================================
# Weights
w_c = 1e-3
w_l = w_c
w_r = w_c

# Track Euler Angle Weights
w_theta = 3e3
w_mu = 1e10
w_phi = w_mu

# Track Width Weights
w_nl = 1e-2
w_nr = w_nl

# Initial Track Width Guess
nl_0 = 5.3
nr_0 = -5.3

#======================================================================================================================================================================================================================
# Limits
# Track Pitch Limits [rad]
mu_min = -np.pi/2 + 0.1
mu_max = np.pi/2 - 0.1

# Left Half Track Limits [m]
nl_min = 3
nl_max = 10

# Right Half Track Limits [m]
nr_min = -10
nr_max = -3

# Track Yaw about z-axis [rad]
theta_min = -np.inf
theta_max = np.inf

# Track Relative Twist about x-axis [rad]
phi_min = -np.inf
phi_max = np.inf

s_grid = np.linspace(arc_s[0],lap_length,N+1)
ds = (lap_length - arc_s[0])/N
print(f"Step Length in m is:{ds}")

# States
X = SX.sym('X',11,N+1)
U = SX.sym('U',5,N)
states = vertcat(reshape(X, -1, 1), reshape(U, -1, 1))
print(f"Length of State Vector is :{states.numel()}")

# Objective Function
cost = 0

# Constraints
g = []
lbg = []
ubg = []
opts = {}

# Interpolant Functions for Optimal Control
cline_coord = centerline.coord.flatten(order='F').tolist()

if len(right.lat) >= len(left.lat):
    br_coord = new_coord.flatten(order='F').tolist()
    bl_coord = np.array([Xl, Yl, Zl]).flatten(order='F').tolist()
else:
    br_coord = np.array([Xr, Yr, Zr]).flatten(order='F').tolist()
    bl_coord = new_coord.flatten(order='F').tolist()

cline_fun = ca.interpolant('cline_fun','linear',[arc_s.tolist()],cline_coord,opts)
bl_fun = ca.interpolant('bl_fun','linear',[arc_s.tolist()],bl_coord,opts)
br_fun = ca.interpolant('br_fun','linear',[arc_s.tolist()],br_coord,opts)

# Symbolic Track Dynamics Definition
X_sym = SX.sym('X_sym',11,1)
U_sym = SX.sym('U_sym',5,1)

# Symbolic States
x = X_sym[0]
y = X_sym[1]
z = X_sym[2]
theta_f = X_sym[3]
mu_f = X_sym[4]
phi_f = X_sym[5]
theta_dot_f = X_sym[6]
mu_dot_f = X_sym[7]
phi_dot_f = X_sym[8]

# Control States
theta_ddot_f = U_sym[0]
mu_ddot_f = U_sym[1]
phi_ddot_f = U_sym[2]
nl_dot_f = U_sym[3]
nr_dot_f = U_sym[4]

ode = vertcat(cos(theta_f)*cos(mu_f), sin(theta_f)*cos(mu_f), -sin(mu_f), theta_dot_f, mu_dot_f, phi_dot_f, theta_ddot_f, mu_ddot_f, phi_ddot_f, nl_dot_f, nr_dot_f)
f_dynamics = Function('f_dynamics',[X_sym,U_sym],[ode])

#======================================================================================================================================================================================================================
# Nonlinear Program Solver (NLP) Definition
for k in range(N):
    s = s_grid[k]
    state = X[:,k]
    ctrl = U[:,k]

    # State Definition
    cline = vertcat(state[0],state[1],state[2])
    theta = state[3]
    mu = state[4]
    phi = state[5]
    theta_dot = state[6]
    mu_dot = state[7]
    phi_dot = state[8]
    nl = state[9]
    nr = state[10]

    # Control Definition
    theta_ddot = ctrl[0]
    mu_ddot = ctrl[1]
    phi_ddot = ctrl[2]
    nl_dot = ctrl[3]
    nr_dot = ctrl[4]

    # Normal Defintion
    normal = vertcat(cos(theta)*sin(mu)*sin(phi) - sin(theta)*cos(phi),sin(theta)*sin(mu)*sin(phi) + cos(theta)*cos(phi), cos(mu)*sin(phi))
    Boundary_left = cline + nl*normal
    Boundary_right = cline + nr*normal

    # Objective Function Construction
    e_cline = w_c*(ca.sumsqr(cline - cline_fun(s)))
    e_bl = w_l*(ca.sumsqr(Boundary_left - bl_fun(s)))
    e_br = w_r*(ca.sumsqr(Boundary_right - br_fun(s)))
    r_c = w_theta*(theta_ddot**2) + w_mu*(mu_ddot**2) + w_phi*(phi_ddot**2)
    r_w = w_nl*(nl_dot**2) + w_nr*(nr_dot**2)
    cost = cost + ds*(e_cline + e_bl + e_br + r_c + r_w)

    # Dynamics Definition
    # 4th Order Ranga Kutta Method for Discretization
    k1 = f_dynamics(state,ctrl)
    k2 = f_dynamics(state + (ds/2)*k1, ctrl)
    k3 = f_dynamics(state + (ds/2)*k2, ctrl)
    k4 = f_dynamics(state + (ds)*k3, ctrl)

    state_next = state + (ds/6)*(k1 + 2*k2 + 2*k3 + k4)
    g.append(X[:,k+1] - state_next)

# Bounds and Limits
lbx = np.full(states.shape, -np.inf)
ubx = np.full(states.shape, np.inf)

# Indices in states
for r in range(N+1):
    idx_theta = r*11 + 3
    idx_mu = r*11 + 4
    idx_phi = r*11 + 5
    idx_nl = r*11 + 9
    idx_nr = r*11 + 10

    # Limits on mu
    lbx[idx_mu] = mu_min
    ubx[idx_mu] = mu_max

    lbx[idx_theta] = theta_min
    ubx[idx_theta] = theta_max

    lbx[idx_phi] = phi_min
    ubx[idx_phi] = phi_max    

    lbx[idx_nl] = nl_min
    ubx[idx_nl] = nl_max

    lbx[idx_nr] = nr_min
    ubx[idx_nr] = nr_max

#======================================================================================================================================================================================================================
# End Constraints
for w in range(11):
    if w == 3:
        g.append(X[w,N] - X[w,0] - direction_angle)
    else:
       g.append(X[w, N] - X[w,0])

#====================================================================================================================================================================================================================
# Initial Guess
x0 = np.zeros(states.shape)
x_coords = np.interp(s_grid, arc_s, centerline.coord[0,:])
y_coords = np.interp(s_grid, arc_s, centerline.coord[1,:])
z_coords = np.zeros(len(s_grid))
x_guess = np.vstack([x_coords, y_coords, z_coords]).T
dx = np.diff(x_coords)
dy = np.diff(y_coords)
theta_guess = np.unwrap(np.arctan2(dy, dx))
theta_guess = ndimage.gaussian_filter1d(theta_guess, sigma=5) 
theta_guess = np.append(theta_guess, theta_guess[-1])

for o in range(N+1):
    idx = o*11
    x0[idx:idx+3] = x_guess[o,:].reshape(-1,1)
    x0[idx+3] = theta_guess[o]
    x0[idx+9] = nl_0
    x0[idx+10] = nr_0

#=======================================================================================================================================================================================================================
## Nonlinear Solver
g_vector = vertcat(*g)
lbg = np.zeros(g_vector.shape[0])
ubg = np.zeros(g_vector.shape[0])
nlp = {'x': states,'f': cost,'g': g_vector}
opts = {'ipopt': {'max_iter': 2000,'print_level': 5}}

solver = nlpsol('solver', 'ipopt', nlp, opts)
sol = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
full_sol = np.array(sol['x']).flatten()
n_x_total = 11 * (N + 1)
X_opt = full_sol[:n_x_total].reshape((11, N + 1), order='F')
U_opt = full_sol[n_x_total:].reshape((5, N), order='F')
print(f"Shape of X_opt is: {X_opt.shape}")
print(f"Shape of U_opt is: {U_opt.shape}")

##=====================================================================================================================================================================================================================
# Store Track Data
track_data = SimpleNamespace()
track_data.arc_s = s_grid
# Centerline Coordinates
track_data.xc = X_opt[0,:]
track_data.yc = X_opt[1,:]
track_data.zc = X_opt[2,:]
track_data.cline = np.vstack([track_data.xc, track_data.yc, track_data.zc])

# Track Euler Angle
track_data.theta = X_opt[3,:]
track_data.mu = X_opt[4,:]
track_data.phi = X_opt[5,:]

# Track Euler Angle Rates
track_data.theta_dot = X_opt[6,:]
track_data.mu_dot = X_opt[7,:]
track_data.phi_dot = X_opt[8,:]

# Left and Right Track Width
track_data.nl = X_opt[9,:]
track_data.nr = X_opt[10,:]

# Tangent [t], Normal [n], Unit Normal [m], Left & Right Boundary Coordinates
# Track Tangent
track_data.t = np.zeros(track_data.cline.shape)
track_data.t[0,:] = np.cos(track_data.theta)*np.cos(track_data.mu)
track_data.t[1,:] = np.sin(track_data.theta)*np.cos(track_data.mu)
track_data.t[2,:] = -np.sin(track_data.mu)

# Track Normal
track_data.n = np.zeros(track_data.cline.shape)
track_data.n[0,:] = np.cos(track_data.theta)*np.sin(track_data.mu)*np.sin(track_data.phi) - np.sin(track_data.theta)*np.cos(track_data.phi)
track_data.n[1,:] = np.sin(track_data.theta)*np.sin(track_data.mu)*np.sin(track_data.phi) + np.cos(track_data.theta)*np.cos(track_data.phi)
track_data.n[2,:] = np.cos(track_data.mu)*np.sin(track_data.phi)

# Track Unit Normal
track_data.m = np.cross(track_data.t,track_data.n,axis=0)

# Left Boundary of the Track
track_data.bl = track_data.cline + track_data.nl*track_data.n

# Right Boundary of the Track
track_data.br = track_data.cline + track_data.nr*track_data.n

# Track Euler Angle Angular Velocities
# Track Roll Rate
track_data.omega_x = np.zeros(track_data.theta.shape)
track_data.omega_x = track_data.phi_dot -track_data.theta_dot*np.sin(track_data.mu)

# Track Pitch Rate
track_data.omega_y = np.zeros(track_data.theta.shape)
track_data.omega_y = track_data.theta_dot*np.cos(track_data.mu)*np.sin(track_data.phi) + track_data.mu_dot*np.cos(track_data.phi)

# Track Yaw Rate
track_data.omega_z = np.zeros(track_data.theta.shape)
track_data.omega_z = track_data.theta_dot*np.cos(track_data.mu)*np.cos(track_data.phi) - track_data.mu_dot*np.sin(track_data.phi)

# Save Track Data
root = tk.Tk()
root.withdraw()
files = [('Pickle Files','*.pkl'), ('All Files', '*.*')]
full_save_path = filedialog.asksaveasfilename(
    defaultextension=".pkl",
    filetypes=files,
    title="Save Track Data"
)
with open(full_save_path, 'wb') as f:
    pickle.dump(track_data, f, protocol=pickle.HIGHEST_PROTOCOL)
print(f"Track data successfully saved to: {full_save_path}")

#=====================================================================================================================================================================================================================
# Results and Plots
# Track Map
fig = go.Figure()

# Right Boundary raw GPS Data
fig.add_trace(go.Scatter3d(
    x=Xr, y=Yr, z=Zr,
    mode='lines',
    name='Right Boundary - Raw GPS Data',
    line=dict(color='blue', width=4,dash='dash')
))

# Optimized Right Boundary Data
fig.add_trace(go.Scatter3d(
    x=track_data.br[0,:], y=track_data.br[1,:], z=track_data.br[2,:],
    mode='lines',
    name='Right Boundary - Optimized',
    line=dict(color='blue', width=4)
))

# Left Boundary raw GPS data
fig.add_trace(go.Scatter3d(
    x=Xl, y=Yl, z=Zl, 
    mode='lines',
    name='Left Boundary - Raw GPS Data',
    line=dict(color='red', width=4,dash='dash')
))

# Optimized Left Boundary Data
fig.add_trace(go.Scatter3d(
    x=track_data.bl[0,:], y=track_data.bl[1,:], z=track_data.bl[2,:], 
    mode='lines',
    name='Left Boundary - Optimized',
    line=dict(color='red', width=4)
))

# Raw GPS data centerline
fig.add_trace(go.Scatter3d(
    x=centerline.coord[0,:], y=centerline.coord[1,:], z=centerline.coord[2,:], 
    mode='lines',
    name='Centerline Coordinates - Raw GPS Data',
    line=dict(color='black',width=4, dash='dash')
))

# Optimized Centerline
fig.add_trace(go.Scatter3d(
    x=track_data.xc, y=track_data.yc, z=track_data.zc, 
    mode='lines',
    name='Centerline Coordinates - Optimized',
    line=dict(color='black',width=4)
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

plt.figure(1)
# Euler Angles
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# Yaw Angle
ax1.plot(track_data.arc_s,track_data.theta,color='black')
ax1.set_xlabel('Track Centerline Arc Length [m]')
ax1.set_ylabel('Track Yaw Angle [rad]')
ax1.set_title('Track Yaw Angle vs Centerline Arc Length')
ax1.axhline(0,color='red')
ax1.grid(True,alpha=0.3)

# Pitch Angle
ax2.plot(track_data.arc_s,track_data.mu,color='black')
ax2.set_xlabel('Track Centerline Arc Length [m]')
ax2.set_ylabel('Track Pitch Angle [rad]')
ax2.axhline(0,color='red')
ax2.set_title('Track Pitch Angle vs Centerline Arc Length')
ax2.grid(True,alpha=0.3)

# Roll Angle
ax3.plot(track_data.arc_s,track_data.phi,color='black')
ax3.set_xlabel('Track Centerline Arc Length [m]')
ax3.set_ylabel('Track Roll Angle [rad]')
ax3.axhline(0,color='red')
ax3.set_title('Track Roll Angle vs Centerline Arc Length')
ax3.grid(True,alpha=0.3)
plt.tight_layout()
plt.show()

plt.figure(2)
# Track Half Widths
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Left Half Track Width
ax1.plot(track_data.arc_s,track_data.nl,color='black')
ax1.set_xlabel('Track Centerline Arc Length [m]')
ax1.set_ylabel('Track Left Half Width (n_l) [m]')
ax1.grid(True,alpha=0.3)
ax1.axhline(nl_0,color='red')

# Right Half Track Width
ax2.plot(track_data.arc_s,track_data.nr,color='black')
ax2.set_xlabel('Track Centerline Arc Length [m]')
ax2.set_ylabel('Track Right Half Width (n_r) [m]')
ax2.grid(True,alpha=0.3)
ax2.axhline(nr_0,color='red')
plt.tight_layout()
plt.show()

plt.figure(3)
# Euler Angle Angular Velocities
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 12), sharex=True)

# Yaw Angle Rate
ax1.plot(track_data.arc_s, track_data.omega_z, color='black')
ax1.set_ylabel('Track Yaw Rate [rad/sec]')
ax1.set_title('Track Yaw Rate vs Centerline Arc Length')
ax1.grid(True, alpha=0.3)
ax1.axhline(0,color='red')

# Pitch Angle Rate
ax2.plot(track_data.arc_s, track_data.omega_y, color='black')
ax2.set_ylabel('Track Pitch Rate [rad/sec]')
ax2.set_title('Track Pitch Rate vs Centerline Arc Length')
ax2.grid(True, alpha=0.3)
ax2.axhline(0,color='red')

# Roll Angle Rate
ax3.plot(track_data.arc_s, track_data.omega_x, color='black')
ax3.set_ylabel('Track Roll Rate [rad/sec]')
ax3.set_xlabel('Track Centerline Arc Length [m]')
ax3.set_title('Track Roll Rate vs Centerline Arc Length')
ax3.grid(True, alpha=0.3)
ax3.axhline(0,color='red')

plt.tight_layout()
plt.show()

## End
#===================================================================================================================================================================================================================