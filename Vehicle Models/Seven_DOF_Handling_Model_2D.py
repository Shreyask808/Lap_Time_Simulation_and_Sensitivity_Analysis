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
matplotlib.use('Qt5Agg') # Or 'TkAgg' if you don't have Qt installed
import plotly.graph_objects as go
import matplotlib.pyplot as plt
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')
 
def Seven_DOF_Handling_Model_2D(vehicle,track_data,N,loc4):
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
