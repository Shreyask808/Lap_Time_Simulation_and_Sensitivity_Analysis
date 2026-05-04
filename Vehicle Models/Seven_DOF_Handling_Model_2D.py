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
import casadi as ca
from casadi import *
import matplotlib
matplotlib.use('Qt5Agg')
import plotly.graph_objects as go
import matplotlib.pyplot as plt
from tire_models import get_tire_model
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')
 
def Seven_DOF_Handling_Model_2D(vehicle,track_data,N,name):
#=========================================================================================================================================================================================================================
# System Dimensions    
    n_states = 10
    u_states = 9

#==========================================================================================================================================================================================================================
# System Model Definition
    lap_length = track_data.arc_s[-1]                                                                           # Total Lap Length along the track centerline in m                         
    s_grid = np.linspace(0,lap_length,N+1)                                                                      # Discretization of centerline track length in m
    ds = lap_length/N                                                                                           # Step length in m

# CASADI Interpolation Function
    f_kappa = ca.interpolant('f_kappa','linear',[track_data.arc_s.tolist()],track_data.omega_z.tolist())
    f_nl = ca.interpolant('f_nl','linear',[track_data.arc_s.tolist()],track_data.nl.tolist())
    f_nr = ca.interpolant('f_nr','linear',[track_data.arc_s.tolist()],track_data.nr.tolist())
    f_theta = ca.interpolant('f_theta','linear',[track_data.arc_s.tolist()],track_data.theta.tolist())

# States Definition
    X = SX.sym('X',n_states,N+1)
    U = SX.sym('U',u_states,N)
    states = vertcat(reshape(X, -1, 1), reshape(U, -1, 1))
    print(f"Length of State Vector is :{states.numel()}")

# System Definition
    X_sym = SX.sym('X_sym',n_states)                                                                            # Symbolic State Definition
    U_sym = SX.sym('U_sym',u_states)                                                                            # Symbolic Control Input Definition
    s = SX.sym('s')                                                                                             # Symbolic Centerline Arc Length Definition
    # Scaling - X_N [.] = X_scale [1/units]*X [units]
    force_scale = 1/(vehicle.m*vehicle.g)                                                                       # Force Scale in N^-1
    length_scale = 1/vehicle.l                                                                                  # Length Scale in m^-1
    time_scale = np.sqrt(vehicle.g/vehicle.l)                                                                   # Time Scale in sec^-1
    speed_scale = 1/np.sqrt(vehicle.g*vehicle.l)                                                                # Speed Scale in sec/m
    angle_scale = 1                                                                                             # Angle Scale in rad^-1

    # Normalized States
    t_N = X_sym[0]                                                                                              # Normalized Time [.]
    n_N = X_sym[1]
    psi_N = X_sym[2]
    psi_dot_N = X_sym[3]
    u_N = X_sym[4]
    v_N = X_sym[5]
    O_fl_N = X_sym[6]
    O_fr_N = X_sym[7]
    O_rl_N = X_sym[8]
    O_rr_N = X_sym[9]

    # Normalized Control Inputs
    delta_N = U_sym[0]
    Md_fl_N = U_sym[1]
    Md_fr_N = U_sym[2]
    Md_rl_N = U_sym[3]
    Md_rr_N = U_sym[4]
    Fzfl_N = U_sym[5]
    Fzfr_N = U_sym[6]
    Fzrl_N = U_sym[7]
    Fzrr_N = U_sym[8]

    # Dynamics Definition
    # States
    t = t_N/time_scale
    n = n_N/length_scale
    psi = psi_N/angle_scale
    psi_dot = psi_dot_N*time_scale/angle_scale
    u = u_N/speed_scale
    v = v_N/speed_scale
    O_fl = O_fl_N*time_scale/angle_scale
    O_fr = O_fr_N*time_scale/angle_scale
    O_rl = O_rl_N*time_scale/angle_scale
    O_rr = O_rr_N*time_scale/angle_scale

    # Control Inputs
    delta = delta_N/angle_scale
    Md_fl = Md_fl_N/(force_scale*length_scale)
    Md_fr = Md_fr_N/(force_scale*length_scale)
    Md_rl = Md_rl_N/(force_scale*length_scale)
    Md_rr = Md_rr_N/(force_scale*length_scale)
    Fzfl = Fzfl_N/force_scale
    Fzfr = Fzfr_N/force_scale
    Fzrl = Fzrl_N/force_scale
    Fzrr = Fzrr_N/force_scale

    # Symbolic Track Data
    curv = f_kappa(s)
    theta = f_theta(s)
    xi = psi - theta

    # Aerodynamic Forces
    Drag = -(1/2)*vehicle.Cd*vehicle.rho*vehicle.A*(u**2)
    f_Drag = Function('f_Drag',{X_sym},{Drag})
    Downforce = -(1/2)*vehicle.Cl*vehicle.rho*vehicle.A*(u**2)
    f_Downforce = Function('f_Downforce',{X_sym},{Downforce})

    # Tire Contact Patch Velocities in Vehicle Frame
    # Rear Right Tire
    u_rr = u + (vehicle.Wr/2)*psi_dot
    v_rr = v - (vehicle.d)*psi_dot

    # Rear Left Tire
    u_rl = u - (vehicle.Wr/2)*psi_dot
    v_rl = v - (vehicle.d)*psi_dot

    # Front Left Tire
    u_fl = u - (vehicle.Wf/2)*psi_dot
    v_fl = v - (vehicle.l - vehicle.d)*psi_dot

    # Front Right Tire
    u_fr = u + (vehicle.Wf/2)*psi_dot
    v_fr = v - (vehicle.l - vehicle.d)*psi_dot

    # Tire Slip Angle Definition
    alpha_rr = - ca.arctan(v_rr/u_rr)                                                                           # Rear Right Slip Angle in rad
    alpha_rl = -ca.arctan(v_rl/u_rl)                                                                            # Rear Left Slip Angle in rad
    alpha_fr = -ca.arctan((-u_fr*sin(delta) + v_fr*cos(delta))/(u_fr*cos(delta) + v_fr*sin(delta)))             # Front Right Slip Angle in rad
    alpha_fl = -ca.arctan((-u_fl*sin(delta) + v_fl*cos(delta))/(u_fl*cos(delta) + v_fl*sin(delta)))             # Front Left Slip Angle in rad

    # Tire Forces and Moments
    tire_model = get_tire_model(name)

    # Rear Right Tire Forces
    Fxrr,Fyrr,Mzrr,rrr = tire_model(u_rr,O_rr,alpha_rr,Fzrr,vehicle)

    # Rear Left Tire Forces
    Fxrl,Fyrl,Mzrl,rrl = tire_model(u_rl,O_rl,alpha_rl,Fzrl,vehicle)

    # Front Right Tire Forces
    Fxfr,Fyfr,Mzfr,rfr = tire_model(u_fr,O_fr,alpha_fr,Fzfr,vehicle)
    
    # Front Left Tire Forces
    Fxfl,Fyfl,Mzfl,rfl = tire_model(u_fl,O_fl,alpha_fl,Fzfl,vehicle)

    # Dynamics ODE Definition
    s_dot = (1/length_scale)*(u*cos(xi) - v*cos(xi))/(1 - n*curv)
    n_N_dot = (length_scale/time_scale)*(u*sin(xi) + v*cos(xi))
    psi_N_dot = psi_dot_N
    psi_ddot_N = (angle_scale/(time_scale**2))*(1/2)*(-2*vehicle.d*(Fyrl + Fyrr) + 2*(Mzfl + Mzfr + Mzrl + Mzrr) + (Fxrr - Fxrl)*vehicle.Wr + (2*(vehicle.d - vehicle.l)*(Fyfl + Fyfr) + vehicle.Wf*(Fxfr - Fxfl))*cos(delta) + (2*(vehicle.d - vehicle.l)*(Fxfl+Fxfr) + vehicle.Wf*(Fyfl - Fyfr))*sin(delta))*(1/vehicle.Iz)
    u_N_dot = (speed_scale/time_scale)*(v*psi_dot + (Drag + Fxrl + Fxrr + Fxfl*cos(delta) + Fxfr*cos(delta) - Fyfl*sin(delta) - Fyfr*sin(delta))/vehicle.m)
    v_N_dot = (speed_scale/time_scale)*(-u*psi_dot + (Fyrl + Fyrr + Fyfl*cos(delta) + Fyfr*cos(delta) + Fxfl*sin(delta) + Fxfr*sin(delta))/vehicle.m)
    O_fl_N_dot = (angle_scale/(time_scale**2))*(Md_fl - Fxfl*rfl)*(1/vehicle.Ifl)
    O_fr_N_dot = (angle_scale/(time_scale**2))*(Md_fr - Fxfr*rfr)*(1/vehicle.Ifr)
    O_rl_N_dot = (angle_scale/(time_scale**2))*(Md_rl - Fxrl*rrl)*(1/vehicle.Irl)
    O_rr_N_dot = (angle_scale/(time_scale**2))*(Md_rr - Fxrr*rrr)*(1/vehicle.Irr)

    ODE = vertcat(1/s_dot, n_N_dot/s_dot, psi_N_dot/s_dot,psi_ddot_N/s_dot, u_N_dot/s_dot, v_N_dot/s_dot, O_fl_N_dot/s_dot, O_fr_N_dot/s_dot, O_rl_N_dot/s_dot, O_rr_N_dot/s_dot)
    f_dynamics = Function('f_dynamics',[X_sym,U_sym],[ODE])

# Cost Function & Constraint Definition
    # Cost Function
    cost = 0

    # Dynamics Constraints
    g_dynamics = []
    lbg_dynamics = []
    ubg_dynamics = []

    # Time Constraints
    g_time = []
    lbg_time = []
    ubg_time = []

    # Vehicle Power Constraints
    g_power = []
    lbg_power = []
    ubg_power = []

    # Nonlinear Programing Definition
    for k in range(N):
        state = X[:,k]
        ctrl = U[:,k]
        # Dynamics Definition
        # 4th Order Ranga Kutta Method for Discretization
        k1 = f_dynamics(state,ctrl)
        k2 = f_dynamics(state + (ds/2)*k1, ctrl)
        k3 = f_dynamics(state + (ds/2)*k2, ctrl)
        k4 = f_dynamics(state + (ds)*k3, ctrl)
        state_next = state + (ds/6)*(k1 + 2*k2 + 2*k3 + k4)

        # Dynamics Constraint Definition
        g_dynamics.append(X[:,k+1] - state_next)
        lbg_dynamics.append(np.zeros((n_states),1))
        ubg_dynamics.append(np.zeros((n_states),1))

        # Time Constraint Definition (Time is forced to move ahead)
        g_time.append(state_next[1] - state[1])
        lbg_time.append(1e-6)
        ubg_time.append(np.inf)

        