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
h = 150                                                                              # Number of Segments the Track is divided into
p = 5                                                                              # Degree of the Polynomial approximating the state in segments

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
print(f"Shape of the Gauss Pseudospectral Differentiation Matrix (D): {D.shape}")
print(f"Gauss Pseudospectral Differentiation Matrix (D):\n{D}")
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

## Tire Force Constraints
g_tire = []
lbg_tire = []
ubg_tire = []

## Closed Loop Racing Line Constraints
g_end = []
lbg_end = []
ubg_end = []

## Cost Function Definition
cost = 0
e0 = 1e-7
e1 = 1e-4
e2 = 1e-4

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
        g_tire.append((ctrl[0]/force_scale)**2 + (ctrl[1]/force_scale)**2 - (vehicle.mu0*(vehicle.m*vehicle.g - f_Lift(state[2])))**2)
        lbg_tire.append(-np.inf)
        ubg_tire.append(0)

        cost = cost + (segment_length/2)*weights[z]*(e0*(state[1]**2) + e1*(ctrl[0]**2) + e2*(ctrl[1]**2))
           
    g_dynamics.append(ca.mtimes(ca.DM(D), X_segment) - (segment_length/2)*F_dynamics)
    lbg_dynamics.append(np.zeros(F_dynamics.shape))
    ubg_dynamics.append(np.zeros(F_dynamics.shape))

### Complete Cost Function with Terminal and Integral Cost
cost = cost + X[-1,0]/time_scale

print(f"#================================================================================================================================================================================================================")
print(f"")    
print(f"Shape of g_dynamics:{len(g_dynamics)}")
print(f"Shape of lbg_dynamics:{len(lbg_dynamics)}")
print(f"Shape of ubg_dynamics:{len(ubg_dynamics)}")
print(f"")
print(f"#================================================================================================================================================================================================================")

## Closed Loop Constraint
g_end.append(X[-1,[1,2,3]] - X[0,[1,2,3]])
lbg_end.append(np.zeros((1,3)))
ubg_end.append(np.zeros((1,3)))

## Net Constraints
g = vertcat(*[ca.vec(g) for g in g_dynamics],*[ca.vec(g) for g in g_continuity],*[ca.vec(g) for g in g_time],*[ca.vec(g) for g in g_power],*[ca.vec(g) for g in g_tire],*[ca.vec(g) for g in g_end])
# Convert every single list element to a flattened column vector (n, 1)
# You MUST be consistent: if you use ca.vec() for 'g', you must use it for 'lbg' and 'ubg'
lbg = vertcat(
    *[vec(b) for b in lbg_dynamics],
    *[vec(b) for b in lbg_continuity],
    *[vec(b) for b in lbg_time],
    *[vec(b) for b in lbg_power],
    *[vec(b) for b in lbg_tire],
    *[vec(b) for b in lbg_end]
)

ubg = vertcat(
    *[vec(b) for b in ubg_dynamics],
    *[vec(b) for b in ubg_continuity],
    *[vec(b) for b in ubg_time],
    *[vec(b) for b in ubg_power],
    *[vec(b) for b in ubg_tire],
    *[vec(b) for b in ubg_end]
)

## Limits and Initial Guess Definition
lbx = np.full(states.shape, -np.inf)
ubx = np.full(states.shape, np.inf)
x0 = np.zeros(states.numel())
arc_s = []
states_num = h*(p+1)*4

## Indexing Loop
for node in range(h*(p+1)):
    r = node//(p + 1)   # segment index
    q = node%(p + 1)   # collocation point index within segment
    s_0 = r*segment_length
    s_f = (r+1)*segment_length
    s_idx = ((s_f - s_0)/2)*col_p + ((s_0 + s_f)/2)
    s_current = s_idx[q]
    arc_s.append(s_current)

    idx_t = node
    idx_n = h*(p+1) + node
    idx_u = 2*h*(p+1) + node
    idx_v = 3*h*(p+1) + node
    idx_Fd = 4*h*(p+1) + r*p + q
    idx_Fn = 4*h*(p+1) + h*p + r*p + q

    # Time Limits and Initial Guess
    if r == 0 and q == 0:
        lbx[idx_t] = 0
        ubx[idx_t] = 0
    else:
        lbx[idx_t] = 0
        ubx[idx_t] = 500*time_scale
    x0[idx_t] = (s_current/vehicle.umax)*time_scale

    # Lateral Offset Limits and Initial Guess
    s_for_bounds = float(s_f) if q == p else float(s_current)
    lbx[idx_n] = float(f_nr(s_for_bounds))*length_scale
    ubx[idx_n] = float(f_nl(s_for_bounds))*length_scale
    x0[idx_n] = 0

    # Longitudinal Velocity Limits and Initial Guess
    lbx[idx_u] = 2*speed_scale
    ubx[idx_u] = vehicle.umax*speed_scale
    x0[idx_u] = vehicle.umax*speed_scale

    # Lateral Velocity Limits and Initial Guess
    lbx[idx_v] = -vehicle.vmax*speed_scale
    ubx[idx_v] = vehicle.vmax*speed_scale
    x0[idx_v] = 0

    if q < p:
        # Drive / Brake Force Limits and Initial Guess
        lbx[idx_Fd] = (4*vehicle.peakbrakingtorque/vehicle.R)*force_scale
        ubx[idx_Fd] = (2*vehicle.peakdrivingtorque/vehicle.R)*force_scale
        x0[idx_Fd] = (2*vehicle.peakdrivingtorque/vehicle.R)*force_scale

        # Lateral Force Limits and Initial Guess
        lbx[idx_Fn] = -(vehicle.mu0*vehicle.m*vehicle.g - (0.5*vehicle.Cl)*vehicle.rho*vehicle.A*vehicle.umax**2)*force_scale
        ubx[idx_Fn] = (vehicle.mu0*vehicle.m*vehicle.g - (0.5*vehicle.Cl)*vehicle.rho*vehicle.A*vehicle.umax**2)*force_scale
        x0[idx_Fn] = 0

# Nonlinear Solver for Point Mass Model
print(f"#================================================================================================================================================================================================================")
print(f"")
print(f"IPOPT Solver Running - ")
print(f"")
nlp = {'x': states,'f': cost,'g': g}
opts = {'ipopt': {'max_iter': 4000,'print_level': 5,'mu_strategy': 'adaptive'}}
solver = nlpsol('solver', 'ipopt', nlp, opts)
sol = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
full_sol = np.array(sol['x']).flatten()

arc_s = np.array(arc_s)
time_opt = full_sol[:h*(p+1)]/time_scale
full_sol = full_sol[h*(p+1):]

n_opt = full_sol[:h*(p+1)]/length_scale
full_sol = full_sol[h*(p+1):]

u_opt = full_sol[:h*(p+1)]/speed_scale
full_sol = full_sol[h*(p+1):]

v_opt = full_sol[:h*(p+1)]/speed_scale
full_sol = full_sol[h*(p+1):]

F_d_opt = full_sol[:h*p]/force_scale
full_sol = full_sol[h*p:]

F_n_opt = full_sol[:h*p]/force_scale
print(f"")
print(f"#================================================================================================================================================================================================================")

# Results
print(f"#================================================================================================================================================================================================================")
print(f"")
print(f"Optimal Lap Time in Seconds is: {time_opt[-1]} [sec]")
print(f"Shape of arc_s : {arc_s.shape}")
print(f"Shape of Time Vector: {time_opt.shape}")
print(f"Shape of N Vector: {n_opt.shape}")
print(f"Shape of U Vector: {u_opt.shape}")
print(f"Shape of V Vector: {v_opt.shape}")
print(f"Shape of F_d Vector: {F_d_opt.shape}")
print(f"Shape of F_n Vector: {F_n_opt.shape}")
print(f"")
print(f"#================================================================================================================================================================================================================")

# Frenet Coordinates to Cartesian Coordinates
f_normal1 = ca.interpolant('f_normal1','linear',[track_data.arc_s.tolist()],track_data.n[0,:].tolist())
f_normal2 = ca.interpolant('f_normal2','linear',[track_data.arc_s.tolist()],track_data.n[1,:].tolist())
f_normal3 = ca.interpolant('f_normal3','linear',[track_data.arc_s.tolist()],track_data.n[2,:].tolist())

f_xc = ca.interpolant('f_xc','linear',[track_data.arc_s.tolist()],track_data.xc.tolist())
f_yc = ca.interpolant('f_yc','linear',[track_data.arc_s.tolist()],track_data.yc.tolist())
f_zc = ca.interpolant('f_zc','linear',[track_data.arc_s.tolist()],track_data.zc.tolist())

x_opt = np.array(f_xc(arc_s)).flatten() + n_opt * np.array(f_normal1(arc_s)).flatten()
y_opt = np.array(f_yc(arc_s)).flatten() + n_opt * np.array(f_normal2(arc_s)).flatten()
z_opt = np.array(f_zc(arc_s)).flatten() + n_opt * np.array(f_normal3(arc_s)).flatten()

#====================================================================================================================================================================================================================
# Plots 
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


# Vehicle Speed
fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Longitudinal Velocity
ax1.plot(arc_s,u_opt*18/5,color='black')
ax1.set_xlabel('Track Centerline Arc Length [m]')
ax1.set_ylabel('Longitudinal Vehicle Velocity [kmph]')
ax1.grid(True, alpha=0.3)
ax1.axhline(0,color='green')

# Lateral Velocity
ax2.plot(arc_s,v_opt*18/5,color='black')
ax2.set_xlabel('Track Centerline Arc Length [m]')
ax2.set_ylabel('Lateral Vehicle Velocity [kmph]')
ax2.grid(True, alpha=0.3)
ax2.axhline(0,color='green')

plt.tight_layout()
plt.show()

# Lap Time and Lateral Offset Plot
fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Lap Time Plot
ax1.plot(arc_s,time_opt,color='black')
ax1.set_xlabel('Track Centerline Arc Length [m]')
ax1.set_ylabel('Lap Time [sec]')
ax1.grid(True, alpha=0.3)
ax1.axhline(0,color='green')

# Lateral Offset Plot
ax2.plot(arc_s,n_opt,color='black',label = 'Optimal Racing Line')
ax2.plot(track_data.arc_s,track_data.nl,color='red',label = 'Left Track Boundary')
ax2.plot(track_data.arc_s,track_data.nr,color='blue',label = 'Right Track Boundary')
ax2.set_xlabel('Track Centerline Arc Length [m]')
ax2.set_ylabel('Lateral Offset [m]')
ax2.grid(True, alpha=0.3)
ax2.axhline(0,color='green')


plt.tight_layout()
plt.show()