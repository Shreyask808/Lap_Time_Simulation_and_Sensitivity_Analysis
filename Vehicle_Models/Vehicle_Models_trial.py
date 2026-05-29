#=========================================================================================================================================================================================================================
# This project utilizes CasADi for symbolic differentiation and optimization.
# Andersson, J. A. E., et al. (2019). "CasADi: a software framework for 
# nonlinear optimization and optimal control." Mathematical Programming Computation.
#=========================================================================================================================================================================================================================

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
import pickle
import casadi as ca
from casadi import *
import matplotlib
matplotlib.use('Qt5Agg')
import plotly.graph_objects as go
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Tire Models'))
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')
 

 #=====================================================================================================================================================================================================================
# User Inputs
N = 700                                                                                                    # Number of Segments on the track          
name = 'brush'

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

n_states = 6
u_states = 13

lap_length = track_data.arc_s[-1]
s_grid = np.linspace(0,lap_length,N+1)
ds = lap_length/N

#==========================================================================
# 3. INTERPOLANTS
#==========================================================================
f_kappa = ca.interpolant('f_kappa', 'bspline', [track_data.arc_s.tolist()], track_data.omega_z.tolist())
f_theta = ca.interpolant('f_theta', 'bspline', [track_data.arc_s.tolist()], track_data.theta.tolist())
f_nl    = ca.interpolant('f_nl',    'bspline', [track_data.arc_s.tolist()], track_data.nl.tolist())
f_nr    = ca.interpolant('f_nr',    'bspline', [track_data.arc_s.tolist()], track_data.nr.tolist())

#==========================================================================
# 4. DECISION VARIABLES
#==========================================================================
X      = SX.sym('X', n_states, N+1)    # states
U      = SX.sym('U', u_states, N)      # controls
states = vertcat(reshape(X, -1, 1),reshape(U, -1, 1))

#==========================================================================
# 5. SYMBOLIC VARIABLES FOR FUNCTION DEFINITIONS
#==========================================================================
X_sym  = SX.sym('X_sym',  n_states)    # symbolic state
U_sym  = SX.sym('U_sym',  u_states)    # symbolic control
s      = SX.sym('s')                   # arc length
F = SX.sym('F')

#==========================================================================
# 6. SCALING
#==========================================================================
force_scale  = 1/(vehicle.m*vehicle.g)
length_scale = 1/vehicle.l
time_scale   = np.sqrt(vehicle.g/vehicle.l)
speed_scale  = 1/np.sqrt(vehicle.g*vehicle.l)
angle_scale  = 1

#==========================================================================
# 7. UNPACK STATES AND CONTROLS
#==========================================================================
# Normalized states
t_N       = X_sym[0]
n_N       = X_sym[1]
psi_N     = X_sym[2]
psi_dot_N = X_sym[3]
u_N       = X_sym[4]
v_N       = X_sym[5]

# Normalized controls
delta_N  = U_sym[0]
Fxfl_N = U_sym[1]
Fyfl_N = U_sym[2]

Fxfr_N = U_sym[3]
Fyfr_N = U_sym[4]

Fxrl_N = U_sym[5]
Fyrl_N = U_sym[6]

Fxrr_N = U_sym[7]
Fyrr_N = U_sym[8]

Fzfl_tN = U_sym[9]
Fzfr_tN = U_sym[10]
Fzrl_tN = U_sym[11]
Fzrr_tN = U_sym[12]

# Descaled States
t       = t_N/time_scale
n       = n_N/length_scale
psi     = psi_N/angle_scale
psi_dot = psi_dot_N*time_scale/angle_scale
u       = u_N/speed_scale
v       = v_N/speed_scale

# Descaled Control Inputs
delta  = delta_N/angle_scale
Fxfl = Fxfl_N/force_scale
Fyfl = Fyfl_N/force_scale

Fxfr = Fxfr_N/force_scale
Fyfr = Fyfr_N/force_scale

Fxrl = Fxrl_N/force_scale
Fyrl = Fyrl_N/force_scale

Fxrr = Fxrr_N/force_scale
Fyrr = Fyrr_N/force_scale

Fzfl_t = Fzfl_tN/force_scale
Fzfr_t = Fzfr_tN/force_scale
Fzrl_t = Fzrl_tN/force_scale
Fzrr_t = Fzrr_tN/force_scale

#==========================================================================
# 9. TRACK DATA
#==========================================================================
curv  = f_kappa(s)
theta = f_theta(s)
xi    = psi - theta

#==========================================================================
# 10. AERODYNAMICS
#==========================================================================
Drag = (1/2)*vehicle.Cd*vehicle.rho*vehicle.A*(u**2)
Lift =  (1/2)*vehicle.Cl*vehicle.rho*vehicle.A*(u**2)
f_Lift = Function('f_Lift',[X_sym[4]],[Lift])
f_Drag = Function('f_Drag',[X_sym[4]],[Drag])

Fr = vehicle.m*vehicle.g*(vehicle.l - vehicle.d)/(2*vehicle.l) + Lift*(vehicle.d - vehicle.l - vehicle.a)/(2*vehicle.l) + Drag*(vehicle.hcp - vehicle.hcg)/(2*vehicle.l)
Ff = 0.5*(vehicle.m*vehicle.g - Lift) - Fr
f_Fr = Function('f_Fr',[X_sym],[Fr])
f_Ff = Function('f_Ff',[X_sym],[Ff])

A = ca.vertcat(ca.horzcat(1,1,1,1),ca.horzcat((vehicle.Wf/2), -(vehicle.Wf/2), (vehicle.Wr/2), (-vehicle.Wr/2)), ca.horzcat((vehicle.d - vehicle.l), (vehicle.d - vehicle.l), vehicle.d, vehicle.d), ca.horzcat((1 - vehicle.Droll), (vehicle.Droll - 1), -vehicle.Droll, vehicle.Droll))
B = ca.vertcat(vehicle.m*vehicle.g - Lift, -(Fyrl + Fyrr + (Fyfl + Fyfr)*cos(delta) + (Fxfl + Fxfr)*sin(delta)), Drag*(vehicle.hcg - vehicle.hcp) - vehicle.a*Lift - vehicle.hcg*((Fxfr + Fxfl)*cos(delta) - (Fyfr + Fyfl)*sin(delta) - (Fxrl + Fxrr)),0)
Fz_t = vertcat(Fzfl_t,Fzfr_t,Fzrl_t,Fzrr_t)

eps = 1e-4
#Fz_phys = 0.5*(Fz_t + ca.sqrt(Fz_t**2 + eps))
#Fzfl_p = Fz_phys[0]
#Fzfr_p = Fz_phys[1]
#Fzrl_p = Fz_phys[2]
#Fzrr_p = Fz_phys[3]

mechanical_residual = ca.mtimes(A[0:3,:],Fz_t) - B[0:3]
suspension_residual = ca.mtimes(A[3,:],Fz_t) - B[3]
f_algebraic = Function('f_algebraic',[X_sym,U_sym],[mechanical_residual*force_scale,suspension_residual*force_scale])


mu = vehicle.constant*(F**vehicle.n)
f_mu = Function('f_mu',[F],[mu])

#==========================================================================
# 15. DYNAMICS ODE — using Fz_sym
#==========================================================================
s_dot = (1/time_scale)*(u*cos(xi) - v*sin(xi))/(1 - n*curv)
n_N_dot = (length_scale/time_scale)*(u*sin(xi) + v*cos(xi))
psi_N_dot = (angle_scale/time_scale)*psi_dot
psi_ddot_N = (angle_scale/(time_scale**2))*(1/vehicle.Iz)*(((vehicle.Wf/2)*Fxrr - vehicle.d*Fyrr) - ((vehicle.Wr/2)*Fxrl + vehicle.d*Fyrl) + ((Fyfl*sin(delta) - Fxfl*cos(delta))*(vehicle.Wf/2) + (Fyfl*cos(delta) + Fxfl*sin(delta))*(vehicle.l - vehicle.d)) + ((Fxfr*cos(delta) - Fyfr*sin(delta))*(vehicle.Wf/2) + (Fxfr*sin(delta) + Fyfr*cos(delta))*(vehicle.l - vehicle.d)))
u_N_dot = (length_scale/(time_scale**2))*(v*psi_dot + (Fxrr + Fxrl + (Fxfl*cos(delta) - Fyfl*sin(delta)) + (Fxfr*cos(delta) - Fyfr*sin(delta)) - Drag)/(vehicle.m))
v_N_dot = (length_scale/(time_scale**2))*(-u*psi_dot + (Fyrr + Fyrl + (Fxfl*sin(delta) + Fyfl*cos(delta)) + (Fxfr*sin(delta) + Fyfr*cos(delta)))/(vehicle.m))

ODE = vertcat(1/s_dot, n_N_dot/s_dot, psi_N_dot/s_dot,psi_ddot_N/s_dot, u_N_dot/s_dot, v_N_dot/s_dot)
f_dynamics = Function('f_dynamics',[X_sym,U_sym,s],[ODE])

#==========================================================================
# 11. CONTACT PATCH KINEMATICS
#==========================================================================
# Rear right
u_rr = u + (vehicle.Wr/2)*psi_dot + eps
v_rr = v - vehicle.d*psi_dot
f_rr = Function('f_rr',[X_sym,U_sym],[u_rr,v_rr])

# Rear left
u_rl = u - (vehicle.Wr/2)*psi_dot + eps
v_rl = v - vehicle.d*psi_dot
f_rl = Function('f_rl',[X_sym,U_sym],[u_rl,v_rl])

# Front left
u_fl = u*cos(delta) + v*sin(delta) - (vehicle.Wf/2)*cos(delta)*psi_dot + (vehicle.l-vehicle.d)*sin(delta)*psi_dot + eps
v_fl = -u*sin(delta) + v*cos(delta) + (vehicle.Wf/2)*sin(delta)*psi_dot + (vehicle.l-vehicle.d)*cos(delta)*psi_dot
f_fl = Function('f_fl',[X_sym,U_sym],[u_fl,v_fl])

# Front right
u_fr = u*cos(delta) + v*sin(delta) + (vehicle.Wf/2)*cos(delta)*psi_dot + (vehicle.l-vehicle.d)*sin(delta)*psi_dot + eps
v_fr = -u*sin(delta) + v*cos(delta) - (vehicle.Wf/2)*sin(delta)*psi_dot + (vehicle.l-vehicle.d)*cos(delta)*psi_dot
f_fr = Function('f_fr',[X_sym,U_sym],[u_fr,v_fr])

#==========================================================================
# 16. NLP LOOP
#==========================================================================
g_dynamics  = []
lbg_dynamics = []
ubg_dynamics = []

g_time      = []
lbg_time     = []
ubg_time     = []

g_power = []
lbg_power = []
ubg_power = []

cost = 0
e0 = 1e-6
e1 = 0
e2 = 1e-6

for k in range(N):
    s_current = s_grid[k]
    state = X[:,k]
    ctrl = U[:,k]
    k1 = f_dynamics(state, ctrl, s_current)
    k2 = f_dynamics(state + (ds/2)*k1, ctrl, s_current + ds/2)
    k3 = f_dynamics(state + (ds/2)*k2, ctrl, s_current + ds/2)
    k4 = f_dynamics(state + ds*k3, ctrl, s_grid[k+1])

    state_next = state + (ds/6)*(k1 + 2*k2 + 2*k3 + k4)
    mech_eval,susp_eval = f_algebraic(state,ctrl)
    
    g_dynamics.append(X[:,k+1] - state_next)
    lbg_dynamics.append(np.zeros((n_states,1)))
    ubg_dynamics.append(np.zeros((n_states,1)))

    g_dynamics.append(ctrl[1]**2 + ctrl[2]**2 - (f_mu(ctrl[9]/force_scale)*(ctrl[9]))**2)
    lbg_dynamics.append(-np.inf)
    ubg_dynamics.append(0)

    g_dynamics.append(ctrl[3]**2 + ctrl[4]**2 - (f_mu(ctrl[10]/force_scale)*(ctrl[10]))**2)
    lbg_dynamics.append(-np.inf)
    ubg_dynamics.append(0)

    g_dynamics.append(ctrl[5]**2 + ctrl[6]**2 - (f_mu(ctrl[11]/force_scale)*(ctrl[11]))**2)
    lbg_dynamics.append(-np.inf)
    ubg_dynamics.append(0)

    g_dynamics.append(ctrl[7]**2 + ctrl[8]**2 - (f_mu(ctrl[12]/force_scale)*(ctrl[12]))**2)
    lbg_dynamics.append(-np.inf)
    ubg_dynamics.append(0)

    g_dynamics.append(mech_eval)
    lbg_dynamics.append(np.zeros((3,1)))
    ubg_dynamics.append(np.zeros((3,1)))

    g_dynamics.append(susp_eval)
    lbg_dynamics.append(-np.inf)
    ubg_dynamics.append(0)

    U_FR,V_FR = f_fr(state,ctrl)
    U_FL,V_FL = f_fl(state,ctrl)
    U_RR,V_RR = f_rr(state,ctrl)
    U_RL,V_RL = f_rl(state,ctrl)

    g_dynamics.append(arctan(V_RL/U_RL))
    lbg_dynamics.append(-vehicle.alpha_max)
    ubg_dynamics.append(vehicle.alpha_max)

    g_dynamics.append(arctan(V_FL/U_FL))
    lbg_dynamics.append(-vehicle.alpha_max)
    ubg_dynamics.append(vehicle.alpha_max)

    g_dynamics.append(arctan(V_FR/U_FR))
    lbg_dynamics.append(-vehicle.alpha_max)
    ubg_dynamics.append(vehicle.alpha_max)

    g_dynamics.append(arctan(V_RR/U_RR))
    lbg_dynamics.append(-vehicle.alpha_max)
    ubg_dynamics.append(vehicle.alpha_max)

    g_time.append(state_next[0] - state[0])
    lbg_time.append((ds/vehicle.umax)*time_scale)
    ubg_time.append(np.inf)

    g_power.append((ctrl[1] + ctrl[3] + ctrl[5] + ctrl[7])*state[4])
    lbg_power.append(vehicle.peakbrakingpower*force_scale*speed_scale)
    ubg_power.append(vehicle.peakdrivingpower*force_scale*speed_scale)

    dt = (state_next[0] - state[0])/time_scale
    cost = cost + dt*(e1*state[1]**2 + e2*state[5]**2 + e0*(ctrl[0]**2 + ctrl[1]**2 +ctrl[2]**2 + ctrl[3]**2 + ctrl[4]**2 + ctrl[5]**2 + ctrl[6]**2 + ctrl[7]**2 + ctrl[8]**2 + ctrl[9]**2 + ctrl[10]**2 + ctrl[11]**2 + ctrl[12]**2))

cost = cost + X[0,-1]/time_scale

g_end = []
lbg_end = []
ubg_end = []

g_end.append(X[[1,3,4,5],-1] - X[[1,3,4,5],0])
lbg_end.append(np.zeros((4,1)))
ubg_end.append(np.zeros((4,1)))

g_end.append(X[2,-1] - X[2,0] - (track_data.theta[-1] - track_data.theta[0]))
lbg_end.append(0)
ubg_end.append(0)

#==========================================================================
# 18. ASSEMBLE
#==========================================================================
g   = vertcat(*g_dynamics, *g_time, *g_end)
lbg = vertcat(*lbg_dynamics, *lbg_time, *lbg_end)
ubg = vertcat(*ubg_dynamics, *ubg_time, *ubg_end)

#==========================================================================
# 19. BOUNDS AND INITIAL GUESS
#==========================================================================
lbx = np.full(states.numel(), -np.inf)
ubx = np.full(states.numel(),  np.inf)
x0  = np.zeros(states.numel())

for r in range(N+1):
    idx_t       = r*n_states + 0
    idx_n       = r*n_states + 1
    idx_psi     = r*n_states + 2
    idx_psi_dot = r*n_states + 3
    idx_u       = r*n_states + 4
    idx_v       = r*n_states + 5

    idx_delta = n_states*(N+1) + u_states*r + 0
    idx_Fxfl = n_states*(N+1) + u_states*r + 1
    idx_Fyfl = n_states*(N+1) + u_states*r + 2

    idx_Fxfr = n_states*(N+1) + u_states*r + 3
    idx_Fyfr = n_states*(N+1) + u_states*r + 4

    idx_Fxrl = n_states*(N+1) + u_states*r + 5
    idx_Fyrl = n_states*(N+1) + u_states*r + 6

    idx_Fxrr = n_states*(N+1) + u_states*r + 7
    idx_Fyrr = n_states*(N+1) + u_states*r + 8

    idx_Fzfl = n_states*(N+1) + u_states*r + 9
    idx_Fzfr = n_states*(N+1) + u_states*r + 10
    idx_Fzrl = n_states*(N+1) + u_states*r + 11
    idx_Fzrr = n_states*(N+1) + u_states*r + 12

    # State bounds
    if r == 0:
        lbx[idx_t]   = 0;  ubx[idx_t]   = 0
        lbx[idx_psi] = float(f_theta(s_grid[0]))*angle_scale
        ubx[idx_psi] = float(f_theta(s_grid[0]))*angle_scale
    else:
        lbx[idx_t]   = 0;  ubx[idx_t]   = 200*time_scale
        lbx[idx_psi] = float((f_theta(s_grid[r]) - np.pi/3))*angle_scale
        ubx[idx_psi] = float((f_theta(s_grid[r]) + np.pi/3))*angle_scale

    lbx[idx_n] = float(f_nr(s_grid[r]))*length_scale
    ubx[idx_n] = float(f_nl(s_grid[r]))*length_scale
    lbx[idx_u]    = 2*speed_scale
    ubx[idx_u]    = vehicle.umax*speed_scale

    #lbx[idx_v]    = -vehicle.vmax*speed_scale
    #ubx[idx_v]    = vehicle.vmax*speed_scale

    # State initial guess
    x0[idx_t]       = x0_ini.time_opt[r]*time_scale
    x0[idx_n]       = x0_ini.n_opt[r]*length_scale
    x0[idx_psi]     = float(f_theta(s_grid[r]))*angle_scale
    x0[idx_u]       = x0_ini.u_opt[r]*speed_scale
    x0[idx_v]       = x0_ini.v_opt[r]*speed_scale
    x0[idx_psi_dot] = x0_ini.u_opt[r]*float(f_kappa(s_grid[r]))*(angle_scale/time_scale)

    if r < N:
        # Control bounds
        lbx[idx_delta] = vehicle.delta_min*angle_scale
        ubx[idx_delta] = vehicle.delta_max*angle_scale

        lbx[idx_Fxfl] = vehicle.peakbrakingtorque*vehicle.R*force_scale
        ubx[idx_Fxfl] = 0

        lbx[idx_Fxfr] = vehicle.peakbrakingtorque*vehicle.R*force_scale
        ubx[idx_Fxfr] = 0

        lbx[idx_Fxrl] = vehicle.peakbrakingtorque*vehicle.R*force_scale
        ubx[idx_Fxrl] = vehicle.peakdrivingtorque*vehicle.R*force_scale

        lbx[idx_Fxrr] = vehicle.peakbrakingtorque*vehicle.R*force_scale
        ubx[idx_Fxrr] = vehicle.peakdrivingtorque*vehicle.R*force_scale
        
        lbx[idx_Fzrr] = 0
        lbx[idx_Fzrl] = 0
        lbx[idx_Fzfr] = 0
        lbx[idx_Fzfl] = 0

        # Control initial guess
        x0[idx_delta] = 0

        x0[idx_Fxfl] = 0
        x0[idx_Fyfl] = x0_ini.F_n_opt[r]*force_scale/4
        
        x0[idx_Fxfr] = 0
        x0[idx_Fyfr] = x0_ini.F_n_opt[r]*force_scale/4

        x0[idx_Fxrl] = x0_ini.F_d_opt[r]*force_scale/2
        x0[idx_Fyrl] = x0_ini.F_n_opt[r]*force_scale/4
        
        x0[idx_Fxrr] = x0_ini.F_d_opt[r]*force_scale/2
        x0[idx_Fyrr] = x0_ini.F_n_opt[r]*force_scale/4

        x0[idx_Fzrr] = vehicle.m*vehicle.g*(vehicle.l - vehicle.d)*force_scale/(2*vehicle.l)
        x0[idx_Fzrl] = vehicle.m*vehicle.g*(vehicle.l - vehicle.d)*force_scale/(2*vehicle.l)

        x0[idx_Fzfl] = vehicle.m*vehicle.g*vehicle.d*force_scale/(2*vehicle.l)
        x0[idx_Fzfr] = vehicle.m*vehicle.g*vehicle.d*force_scale/(2*vehicle.l)
        

nlp = {'x': states,'f': cost,'g': g}
opts = {'ipopt': {'max_iter': 10000,'print_level': 5,'acceptable_tol': 1e-1}}
solver = nlpsol('solver', 'ipopt', nlp, opts)
sol = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
full_sol = np.array(sol['x']).flatten()
n_x_total =  6*(N + 1)
X_opt = full_sol[:n_x_total].reshape((6, N + 1), order='F')
U_opt = full_sol[n_x_total:].reshape((13, N), order='F')
print(f"Shape of X_opt is: {X_opt.shape}")
print(f"Shape of U_opt is: {U_opt.shape}")

time_opt = X_opt[0,:]/time_scale
n_opt = X_opt[1,:]/length_scale
psi_opt = X_opt[2,:]/angle_scale
psi_dot_opt = X_opt[3,:]*time_scale/angle_scale
u_opt = X_opt[4,:]/speed_scale
v_opt = X_opt[5,:]/speed_scale

delta_opt = U_opt[0,:]/angle_scale
Fxfl_opt = U_opt[1,:]/force_scale
Fyfl_opt = U_opt[2,:]/force_scale

Fxfr_opt = U_opt[3,:]/force_scale
Fyfr_opt = U_opt[4,:]/force_scale

Fxrl_opt = U_opt[5,:]/force_scale
Fyrl_opt = U_opt[6,:]/force_scale

Fxrr_opt = U_opt[7,:]/force_scale
Fyrr_opt = U_opt[8,:]/force_scale

Fzfl_opt = U_opt[9,:]/force_scale
Fzfr_opt = U_opt[10,:]/force_scale
Fzrl_opt = U_opt[11,:]/force_scale
Fzrr_opt = U_opt[12,:]/force_scale

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

# Tire Forces
fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Longitudinal Forces
ax1.plot(s_grid[:-1],Fxfl_opt,color='black',label='Front Left')
ax1.plot(s_grid[:-1],Fxfr_opt,color='blue',label='Front Right')
ax1.plot(s_grid[:-1],Fxrl_opt,color='red',label='Rear Left')
ax1.plot(s_grid[:-1],Fxrr_opt,color='orange',label='Rear Right')
ax1.set_xlabel('Track Centerline Arc Length [m]')
ax1.set_ylabel('Longitudinal Forces [N]')
ax1.grid(True, alpha=0.3)
ax1.axhline(0,color='green')
ax1.legend()

# Lateral Forces
ax2.plot(s_grid[:-1],Fyfl_opt,color='black',label='Front Left')
ax2.plot(s_grid[:-1],Fyfr_opt,color='blue',label='Front Right')
ax2.plot(s_grid[:-1],Fyrl_opt,color='red',label='Rear Left')
ax2.plot(s_grid[:-1],Fyrr_opt,color='orange',label='Rear Right')
ax2.set_xlabel('Track Centerline Arc Length [m]')
ax2.set_ylabel('Lateral Forces [N]')
ax2.grid(True, alpha=0.3)
ax2.axhline(0,color='green')
ax2.legend()

plt.tight_layout()
plt.show()


# Vehicle Velocities
fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Longitudinal Speed
ax1.plot(s_grid,u_opt*(18/5),color='black')
ax1.set_xlabel('Track Centerline Arc Length [m]')
ax1.set_ylabel('Longitudinal Speed [kmph]')
ax1.grid(True, alpha=0.3)
ax1.axhline(0,color='green')

# Lateral Speed
ax2.plot(s_grid,v_opt*(18/5),color='black')
ax2.set_xlabel('Track Centerline Arc Length [m]')
ax2.set_ylabel('Lateral Speed [kmph]')
ax2.grid(True, alpha=0.3)
ax2.axhline(0,color='green')

plt.tight_layout()
plt.show()

# Vehicle Normal Forces
fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)

# Longitudinal Speed
ax1.plot(s_grid[:-1],Fzfl_opt,color='black',label='Front Left')
ax1.plot(s_grid[:-1],Fzfr_opt,color='blue',label='Front Right')
ax1.plot(s_grid[:-1],Fzrl_opt,color='red',label='Rear Left')
ax1.plot(s_grid[:-1],Fzrr_opt,color='orange',label='Rear Right')
ax1.set_xlabel('Track Centerline Arc Length [m]')
ax1.set_ylabel('Normal Forces')
ax1.grid(True, alpha=0.3)
ax1.axhline(0,color='green')

# Lateral Speed
ax2.plot(s_grid[:-1],Fzfl_opt + Fzfr_opt + Fzrl_opt + Fzrr_opt,color='black',label='Front Left')
ax2.set_xlabel('Track Centerline Arc Length [m]')
ax2.set_ylabel('Total Normal Force')
ax2.grid(True, alpha=0.3)
ax2.axhline(0,color='green')

plt.tight_layout()
plt.show()
