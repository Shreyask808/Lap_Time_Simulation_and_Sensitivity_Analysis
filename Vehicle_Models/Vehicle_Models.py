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
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Tire Models'))
from tire_models import get_tire_model
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')
 
def Seven_DOF_Handling_Model_2D(vehicle,track_data,N,name,x0_ini):
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
    f_Drag = Function('f_Drag',[X_sym],[Drag])
    Lift = (1/2)*vehicle.Cl*vehicle.rho*vehicle.A*(u**2)
    f_Lift = Function('f_Lift',[X_sym],[Lift])

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

    # Normal Force 
    Ffl = -vehicle.d*(Lift - vehicle.m*vehicle.g)/(2*vehicle.l) - (vehicle.Droll*vehicle.hcg*(Fyrl + Fyrr + (Fyfl + Fyfr)*cos(delta) + (Fxfl + Fxfr)*sin(delta)))/(vehicle.Droll*(vehicle.Wf - vehicle.Wr) + vehicle.Wr) - ((Fxrl + Fxrr)*vehicle.hcg - vehicle.a*Lift + Drag*(vehicle.hcp - vehicle.hcg) - (Fxfl + Fxfr)*vehicle.hcg*cos(delta) + (Fyfl + Fyfr)*vehicle.hcg*sin(delta))/(2*vehicle.l)
    f_Ffl = Function('f_Ffl',[X_sym,U_sym,s],[Ffl])

    Ffr = (vehicle.d*(vehicle.m*vehicle.g - Lift))/(2*vehicle.l) + (2*vehicle.Droll*vehicle.hcg*(Fyrl + Fyrr + (Fyfl + Fyfr)*cos(delta) + (Fxfl + Fxfr)*sin(delta)))/(vehicle.Droll*(vehicle.Wf - vehicle.Wr) + vehicle.Wr) + ((vehicle.hcg - vehicle.hcp)*Drag - (Fxrl + Fxrr)*vehicle.hcg + vehicle.a*Lift + (Fxfl + Fxfr)*vehicle.hcg*cos(delta) - (Fyfl + Fyfr)*vehicle.hcg*sin(delta))/(vehicle.l)
    f_Ffr = Function('f_Ffr',[X_sym,U_sym,s],[Ffr])

    Frl = ((vehicle.d - vehicle.l)*(Lift - vehicle.m*vehicle.g))/(2*vehicle.l) + (2*vehicle.hcg*(vehicle.Droll - 1)*(Fyrl + Fyrr + (Fyfl + Fyfr)*cos(delta) + (Fxfl + Fxfr)*sin(delta)))/(vehicle.Droll*(vehicle.Wf - vehicle.Wr) + vehicle.Wr) + (Drag*(vehicle.hcp - vehicle.hcg) + (Fxrl + Fxrr)*vehicle.hcg - Lift*vehicle.a - (Fxfl + Fxfr)*vehicle.hcg*cos(delta) + (Fyfl + Fyfr)*vehicle.hcp*sin(delta))/(vehicle.l)
    f_Frl = Function('f_Frl',[X_sym,U_sym,s],[Frl])

    Frr = ((vehicle.d - vehicle.l)*(Lift - vehicle.m*vehicle.g))/(2*vehicle.l) - (2*vehicle.hcg*(vehicle.Droll - 1)*(Fyrl + Fyrr + (Fyfl + Fyfr)*cos(delta) + (Fxfl + Fxfr)*sin(delta)))/(vehicle.Droll*(vehicle.Wf - vehicle.Wr) + vehicle.Wr) + (Drag*(vehicle.hcp - vehicle.hcg) + (Fxrl + Fxrr)*vehicle.hcg -  Lift*vehicle.a - (Fxfl + Fxfr)*vehicle.hcg*cos(delta) + (Fyfl + Fyfr)*vehicle.hcg*sin(delta))/(vehicle.l)
    f_Frr = Function('f_Frr',[X_sym,U_sym,s],[Frr])
    
    # Dynamics ODE Definition
    s_dot = (1/length_scale)*(u*cos(xi) - v*sin(xi))/(1 - n*curv)
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
    f_dynamics = Function('f_dynamics',[X_sym,U_sym,s],[ODE])

#==========================================================================================================================================================================================================================
# Cost Function & Constraint Definition
    # Cost Function
    cost = 0
    
    # Cost Function Weights
    e0 = 1e-6
    e1 = 1e-6
    e2 = 1e-6
    e3 = 1e-6
    e4 = 1e-6

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

    # Vehicle Normal Force Constraints
    g_force = []
    lbg_force = []
    ubg_force = []

    # Nonlinear Programing Definition
    for k in range(N):
        s_current = s_grid[k]
        state = X[:,k]
        ctrl = U[:,k]

        # Dynamics Definition
        # 4th Order Ranga Kutta Method for Discretization
        k1 = f_dynamics(state,ctrl,s_current)
        k2 = f_dynamics(state + (ds/2)*k1,ctrl,s_current + ds/2)
        k3 = f_dynamics(state + (ds/2)*k2, ctrl,s_current + ds/2)
        k4 = f_dynamics(state + (ds)*k3, ctrl,s_grid[k+1])
        state_next = state + (ds/6)*(k1 + 2*k2 + 2*k3 + k4)

        # Dynamics Constraint Definition
        g_dynamics.append(X[:,k+1] - state_next)
        lbg_dynamics.append(np.zeros((n_states,1)))
        ubg_dynamics.append(np.zeros((n_states,1)))

        # Time Constraint Definition (Time is forced to move ahead)
        g_time.append(state_next[0] - state[0])
        lbg_time.append(1e-6)
        ubg_time.append(np.inf)

        # Driving Power Constraint Definition
        g_power.append(vehicle.peakdrivingpower - (1/(force_scale*speed_scale))*(state[6]*ctrl[1] + state[7]*ctrl[2] + state[8]*ctrl[3] + state[9]*ctrl[4]))
        lbg_power.append(-np.inf)
        ubg_power.append(0)

        # Normal Force Constraints
        g_force.append(ctrl[5]/force_scale - f_Ffl(state,ctrl,s_current))
        lbg_force.append(0)
        ubg_force.append(0)

        g_force.append(ctrl[6]/force_scale - f_Ffr(state,ctrl,s_current))
        lbg_force.append(0)
        ubg_force.append(0)

        g_force.append(ctrl[7]/force_scale - f_Frl(state,ctrl,s_current))
        lbg_force.append(0)
        ubg_force.append(0)

        g_force.append(ctrl[8]/force_scale - f_Frr(state,ctrl,s_current))
        lbg_force.append(0)
        ubg_force.append(0)

        # Cost Function
        dt = (state_next[0] - state[0])/time_scale
        cost = cost + dt*(e0*ctrl[0]**2 + e1*ctrl[1]**2 + e2*ctrl[2]**2 + e3*ctrl[3]**2 + e4*ctrl[4]**2)

    cost = cost + X[0,-1]/time_scale
    
    # Closed Loop Constraints
    g_end = []
    lbg_end = []
    ubg_end = []

    g_end.append(X[[1,3,4,5,6,7,8,9],-1] - X[[1,3,4,5,6,7,8,9],0])
    lbg_end.append(np.zeros((8,1)))
    ubg_end.append(np.zeros((8,1)))

    g_end.append(X[2,-1] - X[2,0] - (track_data.theta[-1] - track_data.theta[0]))
    lbg_end.append(0)
    ubg_end.append(0)

    g = vertcat(*g_dynamics, *g_time, *g_power, *g_force, *g_end)
    lbg = vertcat(*lbg_dynamics, *lbg_time, *lbg_power, *lbg_force, *lbg_end)
    ubg = vertcat(*ubg_dynamics, *ubg_time, *ubg_power, *ubg_force, *ubg_end)

    #==========================================================================================================================================================================================================================
    # States and Control Input Limits
    lbx = np.full(states.shape, -np.inf)
    ubx = np.full(states.shape, np.inf)
    x0 = np.zeros(states.shape)

    for r in range(N+1):
        idx_t = r*n_states
        idx_n = r*n_states + 1
        idx_psi = r*n_states + 2
        idx_u = r*n_states + 4
        idx_v = r*n_states + 5
        idx_O_fl = r*n_states + 6
        idx_O_fr = r*n_states + 7
        idx_O_rl = r*n_states + 8
        idx_O_rr = r*n_states + 9

        idx_delta = n_states*(N+1) + u_states*r
        idx_Mdfl = n_states*(N+1) + u_states*r + 1
        idx_Mdfr = n_states*(N+1) + u_states*r + 2
        idx_Mdrl = n_states*(N+1) + u_states*r + 3
        idx_Mdrr = n_states*(N+1) + u_states*r + 4
        idx_Fzfl = n_states*(N+1) + u_states*r + 5
        idx_Fzfr = n_states*(N+1) + u_states*r + 6
        idx_Fzrl = n_states*(N+1) + u_states*r + 7
        idx_Fzrr = n_states*(N+1) + u_states*r + 8
        
        if r == 0:
            lbx[idx_t] = 0
            ubx[idx_t] = 0

            lbx[idx_psi] = float(f_theta(s_grid[r]))*angle_scale
            ubx[idx_psi] = float(f_theta(s_grid[r]))*angle_scale
        else:
            lbx[idx_t] = 0
            ubx[idx_t] = 200*time_scale

            lbx[idx_psi] = (f_theta(s_grid[r]) - np.pi/2)*angle_scale
            ubx[idx_psi] = (f_theta(s_grid[r]) + np.pi/2)*angle_scale
            
        lbx[idx_n] = f_nr(s_grid[r])*length_scale
        ubx[idx_n] = f_nl(s_grid[r])*length_scale
        
        lbx[idx_u] = 5*speed_scale
        ubx[idx_u] = vehicle.umax*speed_scale

        lbx[idx_O_fl] = (5/vehicle.R)*(angle_scale/time_scale)
        ubx[idx_O_fl] = (vehicle.umax/vehicle.R)*(angle_scale/time_scale)

        lbx[idx_O_fr] = (5/vehicle.R)*(angle_scale/time_scale)
        ubx[idx_O_fr] = (vehicle.umax/vehicle.R)*(angle_scale/time_scale)

        lbx[idx_O_rl] = (5/vehicle.R)*(angle_scale/time_scale)
        ubx[idx_O_rl] = (vehicle.umax/vehicle.R)*(angle_scale/time_scale)

        lbx[idx_O_rr] = (5/vehicle.R)*(angle_scale/time_scale)
        ubx[idx_O_rr] = (vehicle.umax/vehicle.R)*(angle_scale/time_scale)

        # States
        x0[idx_t] = x0_ini.time_opt[r]*time_scale
        x0[idx_n] = x0_ini.n_opt[r]*length_scale
        x0[idx_psi]  = float(f_theta(s_grid[r]))*angle_scale
        x0[idx_u] = x0_ini.u_opt[r]*speed_scale
        x0[idx_v] = x0_ini.v_opt[r]*speed_scale
        x0[idx_O_fl] = (x0_ini.u_opt[r]/vehicle.R)*(angle_scale/time_scale)
        x0[idx_O_fr] = (x0_ini.u_opt[r]/vehicle.R)*(angle_scale/time_scale)
        x0[idx_O_rl] = (x0_ini.u_opt[r]/vehicle.R)*(angle_scale/time_scale)
        x0[idx_O_rr] = (x0_ini.u_opt[r]/vehicle.R)*(angle_scale/time_scale)

        if r < N:
            lbx[idx_delta] = vehicle.delta_min*angle_scale
            ubx[idx_delta] = vehicle.delta_max*angle_scale

            lbx[idx_Mdfl] = vehicle.peakbrakingtorque*force_scale*length_scale
            ubx[idx_Mdfl] = vehicle.peakdrivingtorque*force_scale*length_scale

            lbx[idx_Mdfr] = vehicle.peakbrakingtorque*force_scale*length_scale
            ubx[idx_Mdfr] = vehicle.peakdrivingtorque*force_scale*length_scale

            lbx[idx_Mdrl] = vehicle.peakbrakingtorque*force_scale*length_scale
            ubx[idx_Mdrl] = vehicle.peakdrivingtorque*force_scale*length_scale

            lbx[idx_Mdrr] = vehicle.peakbrakingtorque*force_scale*length_scale
            ubx[idx_Mdrr] = vehicle.peakdrivingtorque*force_scale*length_scale

            lbx[idx_Fzfl] = 0
            ubx[idx_Fzfl] = np.inf

            lbx[idx_Fzfr] = 0
            ubx[idx_Fzfr] = np.inf
            
            lbx[idx_Fzrl] = 0
            ubx[idx_Fzrl] = np.inf
            
            lbx[idx_Fzrr] = 0
            ubx[idx_Fzrr] = np.inf

           # Controls
            x0[idx_delta] = 0

            # Step 1 - rough Fz guess
            x0[idx_Fzfl] = (vehicle.m*vehicle.g/4)*force_scale
            x0[idx_Fzfr] = (vehicle.m*vehicle.g/4)*force_scale
            x0[idx_Fzrl] = (vehicle.m*vehicle.g/4)*force_scale
            x0[idx_Fzrr] = (vehicle.m*vehicle.g/4)*force_scale

    return states, cost, g, lbg, ubg, lbx, ubx, x0