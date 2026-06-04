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
from tire_models import get_tire_model
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')
 

 #=====================================================================================================================================================================================================================
# User Inputs
N = 1500                                                                                                    # Number of Segments on the track          
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

n_states = 10
u_states = 9

lap_length = track_data.arc_s[-1]
s_grid = np.linspace(0,lap_length,N+1)
ds = lap_length/N

#==========================================================================
# 3. INTERPOLANTS
#==========================================================================
f_kappa = ca.interpolant('f_kappa', 'linear', [track_data.arc_s.tolist()], track_data.omega_z.tolist())
f_theta = ca.interpolant('f_theta', 'linear', [track_data.arc_s.tolist()], track_data.theta.tolist())
f_nl    = ca.interpolant('f_nl',    'linear', [track_data.arc_s.tolist()], track_data.nl.tolist())
f_nr    = ca.interpolant('f_nr',    'linear', [track_data.arc_s.tolist()], track_data.nr.tolist())

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
O_fl_N    = X_sym[6]
O_fr_N    = X_sym[7]
O_rl_N    = X_sym[8]
O_rr_N    = X_sym[9]

# Normalized controls
delta_N  = U_sym[0]
Md_fl_N  = U_sym[1]
Md_fr_N  = U_sym[2]
Md_rl_N  = U_sym[3]
Md_rr_N  = U_sym[4]
Fzfl_t_N = U_sym[5]
Fzfr_t_N = U_sym[6]
Fzrl_t_N = U_sym[7]
Fzrr_t_N = U_sym[8]

Fz_t = (1/force_scale)*ca.vertcat(Fzfl_t_N,Fzfr_t_N,Fzrl_t_N,Fzrr_t_N)

eps = 1e-4
Fz_phys = 0.5*(Fz_t + ca.sqrt(Fz_t**2 + eps))
Fzfl_p = Fz_phys[0]
Fzfr_p = Fz_phys[1]
Fzrl_p = Fz_phys[2]
Fzrr_p = Fz_phys[3]

#==========================================================================
# 8. DESCALE
#==========================================================================
t       = t_N/time_scale
n       = n_N/length_scale
psi     = psi_N/angle_scale
psi_dot = psi_dot_N*time_scale/angle_scale
u       = u_N/speed_scale
v       = v_N/speed_scale
O_fl    = O_fl_N*time_scale/angle_scale
O_fr    = O_fr_N*time_scale/angle_scale
O_rl    = O_rl_N*time_scale/angle_scale
O_rr    = O_rr_N*time_scale/angle_scale

delta  = delta_N/angle_scale
Md_fl  = Md_fl_N/(force_scale*length_scale)
Md_fr  = Md_fr_N/(force_scale*length_scale)
Md_rl  = Md_rl_N/(force_scale*length_scale)
Md_rr  = Md_rr_N/(force_scale*length_scale)
#==========================================================================
# 9. TRACK DATA
#==========================================================================
curv  = f_kappa(s)
theta = f_theta(s)
xi    = psi - theta

#==========================================================================
# 10. AERODYNAMICS
#==========================================================================
Drag = -(1/2)*vehicle.Cd*vehicle.rho*vehicle.A*(u**2)
Lift =  (1/2)*vehicle.Cl*vehicle.rho*vehicle.A*(u**2)
f_Lift = Function('f_Lift',[X_sym[4]],[Lift])
f_Drag = Function('f_Drag',[X_sym[4]],[Drag])

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
# 12. SLIP ANGLES
#==========================================================================
alpha_rr = ca.atan2(v_rr, u_rr)
alpha_rrt = 0.5*(alpha_rr + vehicle.alpha_max - sqrt((alpha_rr - vehicle.alpha_max)**2 + 1e-4)) - 0.5*(alpha_rr - vehicle.alpha_max - sqrt((alpha_rr + vehicle.alpha_max)**2 + 1e-4)) - vehicle.alpha_max

alpha_rl = ca.atan2(v_rl, u_rl)
alpha_rlt = 0.5*(alpha_rl + vehicle.alpha_max - sqrt((alpha_rl - vehicle.alpha_max)**2 + 1e-4)) - 0.5*(alpha_rl - vehicle.alpha_max - sqrt((alpha_rl + vehicle.alpha_max)**2 + 1e-4)) - vehicle.alpha_max

alpha_fl = ca.atan2(v_fl, u_fl)
alpha_flt = 0.5*(alpha_fl + vehicle.alpha_max - sqrt((alpha_fl - vehicle.alpha_max)**2 + 1e-4)) - 0.5*(alpha_fl - vehicle.alpha_max - sqrt((alpha_fl + vehicle.alpha_max)**2 + 1e-4)) - vehicle.alpha_max

alpha_fr = ca.atan2(v_fr, u_fr)
alpha_frt = 0.5*(alpha_fr + vehicle.alpha_max - sqrt((alpha_fr - vehicle.alpha_max)**2 + 1e-4)) - 0.5*(alpha_fr - vehicle.alpha_max - sqrt((alpha_fr + vehicle.alpha_max)**2 + 1e-4)) - vehicle.alpha_max

#==========================================================================
# 13. TIRE FORCES — using Fz_sym
#==========================================================================
tire_model = get_tire_model(name)

Fxrr, Fyrr, Mzrr, rrr = tire_model(u_rr, O_rr, alpha_rrt, Fzrr_p, vehicle)
Fxrl, Fyrl, Mzrl, rrl = tire_model(u_rl, O_rl, alpha_rlt, Fzrl_p, vehicle)
Fxfr, Fyfr, Mzfr, rfr = tire_model(u_fr, O_fr, alpha_frt, Fzfr_p, vehicle)
Fxfl, Fyfl, Mzfl, rfl = tire_model(u_fl, O_fl, alpha_flt, Fzfl_p, vehicle)

#==========================================================================
# 14. FZ RESIDUAL
#==========================================================================

A = ca.vertcat(ca.horzcat(1,1,1,1),ca.horzcat((vehicle.Wf/2), -(vehicle.Wf/2), (vehicle.Wr/2), (-vehicle.Wr/2)), ca.horzcat((vehicle.d - vehicle.l), (vehicle.d - vehicle.l), vehicle.d, vehicle.d), ca.horzcat((1 - vehicle.Droll), (vehicle.Droll - 1), -vehicle.Droll, vehicle.Droll))
B = ca.vertcat(vehicle.m*vehicle.g - Lift, -(Fyrl + Fyrr + (Fyfl + Fyfr)*cos(delta) + (Fxfl + Fxfr)*sin(delta)), Drag*(vehicle.hcg - vehicle.hcp) - vehicle.a*Lift - vehicle.hcg*((Fxfr + Fxfl)*cos(delta) - (Fyfr + Fyfl)*sin(delta) - (Fxrl + Fxrr)),0)

mechanical_residual = ca.mtimes(A[0:3,:], Fz_phys) - B[0:3]
suspension_residual = ca.mtimes(A[3,:],Fz_t) - B[3]
f_algebraic = Function('f_algebraic',[X_sym,U_sym],[mechanical_residual*force_scale,suspension_residual*force_scale])

#==========================================================================
# 15. DYNAMICS ODE — using Fz_sym
#==========================================================================
s_dot = (1/time_scale)*(u*cos(xi) - v*sin(xi))/(1 - n*curv)
n_N_dot = (length_scale/time_scale)*(u*sin(xi) + v*cos(xi))
psi_N_dot = (angle_scale/time_scale)*psi_dot
psi_ddot_N = (angle_scale/(time_scale**2))*(1/vehicle.Iz)*(((vehicle.Wf/2)*Fxrr - vehicle.d*Fyrr + Mzrr) - ((vehicle.Wr/2)*Fxrl + vehicle.d*Fyrl - Mzrl) + ((Fyfl*sin(delta) - Fxfl*cos(delta))*(vehicle.Wf/2) + (Fyfl*cos(delta) + Fxfl*sin(delta))*(vehicle.l - vehicle.d) + Mzfl) + ((Fxfr*cos(delta) - Fyfr*sin(delta))*(vehicle.Wf/2) + (Fxfr*sin(delta) + Fyfr*cos(delta))*(vehicle.l - vehicle.d) + Mzfr))
u_N_dot = (length_scale/(time_scale**2))*(v*psi_dot + (Fxrr + Fxrl + (Fxfl*cos(delta) - Fyfl*sin(delta)) + (Fxfr*cos(delta) - Fyfr*sin(delta)) - Drag)/(vehicle.m))
v_N_dot = (length_scale/(time_scale**2))*(-u*psi_dot + (Fyrr + Fyrl + (Fxfl*sin(delta) + Fyfl*cos(delta)) + (Fxfr*sin(delta) + Fyfr*cos(delta)))/(vehicle.m))
O_fl_N_dot = (angle_scale/(time_scale**2))*(Md_fl - Fxfl*rfl)*(1/vehicle.Ifl)
O_fr_N_dot = (angle_scale/(time_scale**2))*(Md_fr - Fxfr*rfr)*(1/vehicle.Ifr)
O_rl_N_dot = (angle_scale/(time_scale**2))*(Md_rl - Fxrl*rrl)*(1/vehicle.Irl)
O_rr_N_dot = (angle_scale/(time_scale**2))*(Md_rr - Fxrr*rrr)*(1/vehicle.Irr)

ODE = vertcat(1/s_dot, n_N_dot/s_dot, psi_N_dot/s_dot,psi_ddot_N/s_dot, u_N_dot/s_dot, v_N_dot/s_dot, O_fl_N_dot/s_dot, O_fr_N_dot/s_dot, O_rl_N_dot/s_dot, O_rr_N_dot/s_dot)
f_dynamics = Function('f_dynamics',[X_sym,U_sym,s],[ODE])

penalty_weight = 1.0  # Start with a high penalty
#==========================================================================
# 16. NLP LOOP
#==========================================================================
g_dynamics  = []
lbg_dynamics = []
ubg_dynamics = []

g_time      = []
lbg_time     = []
ubg_time     = []

g_power     = []
lbg_power    = []
ubg_power    = []

cost        = 0

e0 = 1e-6
e1 = 1e-6
e2 = 1e-6
e3 = 1e-6
e4 = 1e-6

for k in range(N):
    s_current = s_grid[k]
    state = X[:,k]
    ctrl = U[:,k]
    k1 = f_dynamics(state, ctrl, s_current)
    k2 = f_dynamics(state + (ds/2)*k1, ctrl, s_current + ds/2)
    k3 = f_dynamics(state + (ds/2)*k2, ctrl, s_current + ds/2)
    k4 = f_dynamics(state + ds*k3, ctrl, s_grid[k+1])

    state_next = state + (ds/6)*(k1 + 2*k2 + 2*k3 + k4)

    g_dynamics.append(X[:,k+1] - state_next)
    lbg_dynamics.append(np.zeros((n_states,1)))
    ubg_dynamics.append(np.zeros((n_states,1)))

    mech_eval,susp_eval = f_algebraic(state,ctrl)

    #g_dynamics.append(mech_eval*force_scale)
    #lbg_dynamics.append(np.zeros((3,1)))
    #ubg_dynamics.append(np.zeros((3,1)))

    #g_dynamics.append(susp_eval*force_scale)
    #lbg_dynamics.append(0)
    #ubg_dynamics.append(0)

    g_time.append(state_next[0] - state[0])
    lbg_time.append((ds/vehicle.umax)*time_scale)
    ubg_time.append(np.inf)

    u_rr_eval,v_rr_eval = f_rr(state,ctrl)
    g_dynamics.append(ca.atan2(v_rr_eval,u_rr_eval))
    lbg_dynamics.append(-vehicle.alpha_max)
    ubg_dynamics.append(vehicle.alpha_max)

    u_rl_eval,v_rl_eval = f_rl(state,ctrl)
    g_dynamics.append(ca.atan2(v_rl_eval,u_rl_eval))
    lbg_dynamics.append(-vehicle.alpha_max)
    ubg_dynamics.append(vehicle.alpha_max)

    u_fr_eval,v_fr_eval = f_fr(state,ctrl)
    g_dynamics.append(ca.atan2(v_fr_eval,u_fr_eval))
    lbg_dynamics.append(-vehicle.alpha_max)
    ubg_dynamics.append(vehicle.alpha_max)

    u_fl_eval,v_fl_eval = f_fl(state,ctrl)
    g_dynamics.append(ca.atan2(v_fl_eval,u_fl_eval))
    lbg_dynamics.append(-vehicle.alpha_max)
    ubg_dynamics.append(vehicle.alpha_max)

    g_power.append((state[6]*ctrl[1] + state[7]*ctrl[2] + state[8]*ctrl[3] + state[9]*ctrl[4]))
    lbg_power.append(vehicle.peakbrakingpower*force_scale*speed_scale)
    ubg_power.append(vehicle.peakdrivingpower*force_scale*speed_scale)

    dt = (state_next[0] - state[0])/time_scale
    cost = cost + dt*(e0*ctrl[0]**2 + e1*ctrl[1]**2 + e2*ctrl[2]**2 + e3*ctrl[3]**2 + e4*ctrl[4]**2) + penalty_weight * ca.sumsqr(susp_eval) + penalty_weight*ca.sumsqr(mech_eval)

cost = cost + X[0,-1]/time_scale
    


#==========================================================================
# 17. END CONSTRAINTS
#==========================================================================

g_end = []
lbg_end = []
ubg_end = []

g_end.append(X[[1,3,4,5,6,7,8,9],-1] - X[[1,3,4,5,6,7,8,9],0])
lbg_end.append(np.zeros((8,1)))
ubg_end.append(np.zeros((8,1)))

g_end.append(X[2,-1] - X[2,0] - (track_data.theta[-1] - track_data.theta[0]))
lbg_end.append(0)
ubg_end.append(0)

#==========================================================================
# 18. ASSEMBLE
#==========================================================================
g   = vertcat(*g_dynamics, *g_time, *g_power, *g_end)
lbg = vertcat(*lbg_dynamics, *lbg_time, *lbg_power, *lbg_end)
ubg = vertcat(*ubg_dynamics, *ubg_time, *ubg_power, *ubg_end)

#==========================================================================
# 19. BOUNDS AND INITIAL GUESS
#==========================================================================
lbx = np.full(states.numel(), -np.inf)
ubx = np.full(states.numel(),  np.inf)
x0  = np.zeros(states.numel())
Force_mag = 1.2
Weq = 0.5*(vehicle.Wf + vehicle.Wr)
Fr_static = (vehicle.l - vehicle.d)*(vehicle.m*vehicle.g)/(vehicle.l)
Ff_static = vehicle.m*vehicle.g*vehicle.d/vehicle.l

for r in range(N+1):
    idx_t       = r*n_states + 0
    idx_n       = r*n_states + 1
    idx_psi     = r*n_states + 2
    idx_psi_dot = r*n_states + 3
    idx_u       = r*n_states + 4
    idx_v       = r*n_states + 5
    idx_O_fl    = r*n_states + 6
    idx_O_fr    = r*n_states + 7
    idx_O_rl    = r*n_states + 8
    idx_O_rr    = r*n_states + 9

    idx_delta = n_states*(N+1) + u_states*r + 0
    idx_Mdfl  = n_states*(N+1) + u_states*r + 1
    idx_Mdfr  = n_states*(N+1) + u_states*r + 2
    idx_Mdrl  = n_states*(N+1) + u_states*r + 3
    idx_Mdrr  = n_states*(N+1) + u_states*r + 4

    idx_Fzt_fl = n_states*(N+1) + u_states*r + 5
    idx_Fzt_fr = n_states*(N+1) + u_states*r + 6
    idx_Fzt_rl = n_states*(N+1) + u_states*r + 7
    idx_Fzt_rr = n_states*(N+1) + u_states*r + 8

    # State bounds
    if r == 0:
        lbx[idx_t]   = 0;  ubx[idx_t]   = 0
        lbx[idx_psi] = float(f_theta(s_grid[0]))*angle_scale
        ubx[idx_psi] = float(f_theta(s_grid[0]))*angle_scale
    else:
        lbx[idx_t]   = 0;  ubx[idx_t]   = 200*time_scale
        lbx[idx_psi] = (f_theta(s_grid[r]) - np.pi/3)*angle_scale
        ubx[idx_psi] = (f_theta(s_grid[r]) + np.pi/3)*angle_scale

    lbx[idx_n]    = f_nr(s_grid[r])*length_scale
    ubx[idx_n]    = f_nl(s_grid[r])*length_scale
    lbx[idx_u]    = 2*speed_scale
    ubx[idx_u]    = vehicle.umax*speed_scale
    lbx[idx_O_fl] = (2/vehicle.R)*(angle_scale/time_scale)
    ubx[idx_O_fl] = (vehicle.umax/vehicle.R)*(angle_scale/time_scale)
    lbx[idx_O_fr] = (5/vehicle.R)*(angle_scale/time_scale)
    ubx[idx_O_fr] = (vehicle.umax/vehicle.R)*(angle_scale/time_scale)
    lbx[idx_O_rl] = (2/vehicle.R)*(angle_scale/time_scale)
    ubx[idx_O_rl] = (vehicle.umax/vehicle.R)*(angle_scale/time_scale)
    lbx[idx_O_rr] = (2/vehicle.R)*(angle_scale/time_scale)
    ubx[idx_O_rr] = (vehicle.umax/vehicle.R)*(angle_scale/time_scale)

    # State initial guess
    x0[idx_t]       = x0_ini.time_opt[r]*time_scale
    x0[idx_n]       = x0_ini.n_opt[r]*length_scale
    x0[idx_psi]     = float(f_theta(s_grid[r]))*angle_scale
    x0[idx_u]       = x0_ini.u_opt[r]*speed_scale
    x0[idx_v]       = 0
    x0[idx_O_fl]    = (x0_ini.u_opt[r]/vehicle.R)*(angle_scale/time_scale)
    x0[idx_O_fr]    = (x0_ini.u_opt[r]/vehicle.R)*(angle_scale/time_scale)
    x0[idx_O_rl]    = (x0_ini.u_opt[r]/vehicle.R)*(angle_scale/time_scale)
    x0[idx_O_rr]    = (x0_ini.u_opt[r]/vehicle.R)*(angle_scale/time_scale)
    x0[idx_psi_dot] = x0_ini.u_opt[r]*float(f_kappa(s_grid[r]))*(angle_scale/time_scale)

    if r < N:
        # Control bounds
        lbx[idx_delta] = vehicle.delta_min*angle_scale
        ubx[idx_delta] = vehicle.delta_max*angle_scale
        lbx[idx_Mdfl]  = vehicle.peakbrakingtorque*force_scale*length_scale
        ubx[idx_Mdfl]  = vehicle.peakdrivingtorque*force_scale*length_scale
        lbx[idx_Mdfr]  = vehicle.peakbrakingtorque*force_scale*length_scale
        ubx[idx_Mdfr]  = vehicle.peakdrivingtorque*force_scale*length_scale
        lbx[idx_Mdrl]  = vehicle.peakbrakingtorque*force_scale*length_scale
        ubx[idx_Mdrl]  = vehicle.peakdrivingtorque*force_scale*length_scale
        lbx[idx_Mdrr]  = vehicle.peakbrakingtorque*force_scale*length_scale
        ubx[idx_Mdrr]  = vehicle.peakdrivingtorque*force_scale*length_scale
        
        lbx[idx_Fzt_fl] = 0
        ubx[idx_Fzt_fl] = max(float((force_scale*Force_mag)*((f_Drag(x0_ini.u_opt[r]*speed_scale)*(vehicle.hcg - vehicle.hcp) + (vehicle.a - vehicle.d)*f_Lift(x0_ini.u_opt[r]*speed_scale) + vehicle.m*vehicle.g*vehicle.d)/(2*vehicle.l) - vehicle.m*vehicle.hcg*x0_ini.ax_opt[r]/(2*vehicle.l) - vehicle.Droll*vehicle.m*x0_ini.ay_opt[r]*vehicle.hcg/Weq)),Ff_static*force_scale)

        lbx[idx_Fzt_fr] = 0
        ubx[idx_Fzt_fr] = max(float((force_scale*Force_mag)*((f_Drag(x0_ini.u_opt[r]*speed_scale)*(vehicle.hcg - vehicle.hcp) + (vehicle.a - vehicle.d)*f_Lift(x0_ini.u_opt[r]*speed_scale) + vehicle.m*vehicle.g*vehicle.d)/(2*vehicle.l) - vehicle.m*vehicle.hcg*x0_ini.ax_opt[r]/(2*vehicle.l) + vehicle.Droll*vehicle.m*x0_ini.ay_opt[r]*vehicle.hcg/Weq)),Ff_static*force_scale)

        lbx[idx_Fzt_rl] = 0
        ubx[idx_Fzt_rl] = max(float((force_scale*Force_mag)*(((vehicle.m*vehicle.g - f_Lift(x0_ini.u_opt[r]*speed_scale))*(vehicle.l - vehicle.d) - vehicle.a*f_Lift(x0_ini.u_opt[r]*speed_scale) + f_Drag(x0_ini.u_opt[r]*speed_scale)*(vehicle.hcp - vehicle.hcg))/(2*vehicle.l) + vehicle.m*vehicle.hcg*x0_ini.ax_opt[r]/(2*vehicle.l) - (1 - vehicle.Droll)*(vehicle.m*x0_ini.ay_opt[r]*vehicle.hcg)/(Weq))),Fr_static*force_scale)

        lbx[idx_Fzt_rr] = 0
        ubx[idx_Fzt_rr] = max(float((force_scale*Force_mag)*(((vehicle.m*vehicle.g - f_Lift(x0_ini.u_opt[r]*speed_scale))*(vehicle.l - vehicle.d) - vehicle.a*f_Lift(x0_ini.u_opt[r]*speed_scale) + f_Drag(x0_ini.u_opt[r]*speed_scale)*(vehicle.hcp - vehicle.hcg))/(2*vehicle.l) + vehicle.m*vehicle.hcg*x0_ini.ax_opt[r]/(2*vehicle.l) + (1 - vehicle.Droll)*(vehicle.m*x0_ini.ay_opt[r]*vehicle.hcg)/(Weq))),Fr_static*force_scale)

        # Control initial guess
        x0[idx_delta] = 0
        x0[idx_Mdfl]  = (x0_ini.F_d_opt[r]*vehicle.R/4)*force_scale*length_scale
        x0[idx_Mdfr]  = (x0_ini.F_d_opt[r]*vehicle.R/4)*force_scale*length_scale
        x0[idx_Mdrl]  = (x0_ini.F_d_opt[r]*vehicle.R/4)*force_scale*length_scale
        x0[idx_Mdrr]  = (x0_ini.F_d_opt[r]*vehicle.R/4)*force_scale*length_scale
        
        x0[idx_Fzt_fl] = ((vehicle.m*vehicle.g*(vehicle.l - vehicle.d))/(2*vehicle.l))*force_scale
        x0[idx_Fzt_fr] = ((vehicle.m*vehicle.g*(vehicle.l - vehicle.d))/(2*vehicle.l))*force_scale
        x0[idx_Fzt_rl] = ((vehicle.m*vehicle.g*vehicle.d)/(2*vehicle.l))*force_scale
        x0[idx_Fzt_rr] = ((vehicle.m*vehicle.g*vehicle.d)/(2*vehicle.l))*force_scale

nlp = {'x': states,'f': cost,'g': g}
opts = {'ipopt': {'max_iter': 4000,'print_level': 5,'mu_strategy': 'adaptive','acceptable_tol': 1e-3}}
solver = nlpsol('solver', 'ipopt', nlp, opts)
sol = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
full_sol = np.array(sol['x']).flatten()