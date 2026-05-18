#====================================================================================================================================================================================
# This project utilizes CasADi for symbolic differentiation and optimization.
# Andersson, J. A. E., et al. (2019). "CasADi: a software framework for
# nonlinear optimization and optimal control." Mathematical Programming Computation.
#
# Collocation method: Legendre-Gauss-Radau (LGR) Pseudospectral Method
# Reference: Limebeer & Perantoni (2015), JDSMC Vol. 137, DOI: 10.1115/1.4029466
#            Patterson & Rao (2013), GPOPS-II, University of Florida
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
matplotlib.use('Qt5Agg')  # Or 'TkAgg' if you don't have Qt installed
import plotly.graph_objects as go
import matplotlib.pyplot as plt
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')

#=====================================================================================================================================================================================================================
# LGR Pseudospectral Collocation Utilities
#=====================================================================================================================================================================================================================

def lgr_nodes_weights(N):
    """
    Compute N Legendre-Gauss-Radau collocation nodes and weights on [-1, 1).
    The LGR nodes are roots of P_N(s) + P_{N-1}(s), where P_N is the Nth-degree
    Legendre polynomial. One node is always fixed at s = -1.

    Returns
    -------
    tau : np.ndarray (N,)   nodes in [-1, 1), ascending, tau[0] = -1
    w   : np.ndarray (N,)   quadrature weights (sum to 2)
    """
    # Start from Gauss-Legendre nodes as initial guess for the interior nodes
    if N == 1:
        return np.array([-1.0]), np.array([2.0])

    # Build the companion matrix for P_N + P_{N-1} = 0
    # Equivalently, the LGR nodes satisfy the generalised eigenvalue problem.
    # We use the well-known recurrence: construct (N-1) interior nodes via
    # the eigenvalues of the Jacobi matrix for P_{N-1}, then prepend -1.

    # Jacobi matrix for Legendre polynomials (symmetric tridiagonal)
    n_int = N - 1  # number of interior nodes
    beta = np.array([k / np.sqrt(4*k*k - 1) for k in range(1, n_int + 1)])
    J = np.diag(beta[:-1], -1) + np.diag(beta[:-1], 1)

    # The LGR interior nodes are NOT simply GL nodes; we use the direct
    # approach: roots of (1+s)*P_{N-1}(s) excluding s=-1 via companion matrix.
    # A numerically stable alternative is the Golub-Welsch-style computation
    # for LGR nodes specifically.

    # Build the (N x N) LGR eigenvalue matrix following Trefethen (2000)
    # and Hesthaven et al. — differentiation matrix approach
    # tau[0] = -1 always (the "left" Radau node)
    tau = np.zeros(N)
    tau[0] = -1.0

    if N > 1:
        # Interior LGR nodes: eigenvalues of the (N-1)-order Jacobi-like matrix
        # for the LGR distribution. We compute them as roots of
        # (1 + s) * P_{N-1}(s) shifted; a stable route is via the fact that
        # the N-1 interior nodes are the GL nodes of order N-1 mapped with the
        # s = -1 pin. We use numpy.polynomial.legendre to evaluate and find roots.
        from numpy.polynomial.legendre import leggauss
        gl_nodes, _ = leggauss(N - 1)  # (N-1) GL nodes on (-1,1)
        # Refine: the true LGR nodes (excl. -1) satisfy P_N(tau) + P_{N-1}(tau) = 0
        # Use Newton iteration from GL nodes as starting guess.
        from numpy.polynomial.legendre import legval
        for i, s0 in enumerate(sorted(gl_nodes)):
            s = s0
            for _ in range(100):
                # Evaluate P_{N-1} and P_N via recurrence
                Pm1, P0 = 1.0, s
                for k in range(1, N - 1):
                    Pm1, P0 = P0, ((2*k+1)*s*P0 - k*Pm1) / (k+1)
                PN_1 = P0  # P_{N-1}
                Pm1, P0 = 1.0, s
                for k in range(1, N):
                    Pm1, P0 = P0, ((2*k+1)*s*P0 - k*Pm1) / (k+1)
                PN = P0    # P_N
                f   =  PN + PN_1
                # Derivative: d/ds [P_N + P_{N-1}]
                # Using: P_N'(s) = N*(P_{N-1}(s) - s*P_N(s)) / (1 - s^2)
                if abs(1 - s**2) > 1e-14:
                    dPN   = N   * (PN_1 - s*PN)   / (1 - s**2)
                    # P_{N-2}
                    Pm2, Pm1_ = 1.0, s
                    for k in range(1, N - 2):
                        Pm2, Pm1_ = Pm1_, ((2*k+1)*s*Pm1_ - k*Pm2) / (k+1)
                    PN_2 = Pm2 if N == 2 else Pm1_
                    dPN_1 = (N-1) * (PN_2 - s*PN_1) / (1 - s**2)
                    df = dPN + dPN_1
                else:
                    break
                if abs(df) < 1e-14:
                    break
                delta = f / df
                s -= delta
                if abs(delta) < 1e-14:
                    break
            tau[i + 1] = s
        tau[1:] = np.sort(tau[1:])

    # LGR weights
    # w_i = 1/N^2 * (1 - tau_i) / [P_{N-1}(tau_i)]^2  for i >= 1
    # w_0 = 2/N^2  (the -1 node)
    w = np.zeros(N)
    w[0] = 2.0 / (N * N)
    for i in range(1, N):
        s = tau[i]
        Pm1, P0 = 1.0, s
        for k in range(1, N - 1):
            Pm1, P0 = P0, ((2*k+1)*s*P0 - k*Pm1) / (k+1)
        PN_1 = P0
        w[i] = (1.0 - s) / (N * PN_1)**2

    return tau, w


def lgr_differentiation_matrix(tau):
    """
    Build the N x N LGR differentiation matrix D such that
    (D @ x_col)[i] approximates dx/dtau at tau[i].

    For LGR collocation the state derivative at the collocation points is
    D @ X ≈ dX/dtau, matching the scaled ODE rhs F.
    """
    N = len(tau)
    D = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            if i != j:
                # Barycentric Lagrange differentiation
                num = 1.0
                den = 1.0
                for k in range(N):
                    if k != j:
                        num *= (tau[i] - tau[k])
                    if k != i:
                        den *= (tau[i] - tau[k])
                D[i, j] = num / (den * (tau[i] - tau[j]))
        # Diagonal via sum condition: each row of D sums to 0
        D[i, i] = -np.sum([D[i, j] for j in range(N) if j != i])
    return D


#=====================================================================================================================================================================================================================
# User Inputs
#=====================================================================================================================================================================================================================

# Number of mesh intervals (segments). With LGR we also choose the
# polynomial degree K within each interval (K=3 gives a good balance).
N_intervals = 500       # Number of mesh intervals
K = 4                   # LGR polynomial degree per interval (collocation points per segment)
                        # Total collocation points = N_intervals * K

#=====================================================================================================================================================================================================================
# Import Track Data
#=====================================================================================================================================================================================================================
root = tk.Tk()
root.withdraw()
root.attributes('-topmost', True)

file_path1 = filedialog.askopenfilename(title="Select the Track Data",
                                        filetypes=[("pickle files", "*.pkl")])
if file_path1:
    with open(file_path1, 'rb') as f:
        track_data = pickle.load(f)
    print("Loaded Track Data")
    print(f"Number of Data Points on the Track: {len(track_data.arc_s)}")
else:
    raise FileNotFoundError("Track Data file not found")

#=====================================================================================================================================================================================================================
# Import Vehicle Data
#=====================================================================================================================================================================================================================
root = tk.Tk()
root.withdraw()
root.attributes('-topmost', True)

file_path2 = filedialog.askopenfilename(title="Select the Vehicle Data",
                                        filetypes=[("json files", "*.json")])
if file_path2:
    with open(file_path2, 'r') as f:
        vehicledata = json.load(f)
        vehicle = SimpleNamespace(**vehicledata)
    print("Loaded Vehicle Data")
    print(f"Mass of the Vehicle: {vehicle.m}")
else:
    raise FileNotFoundError("Vehicle Data not found")

#=====================================================================================================================================================================================================================
# LGR Pseudospectral Optimal Control Problem
# Replaces the Hermite-Simpson single-step transcription with multi-interval
# LGR collocation following Limebeer & Perantoni (2015) / GPOPS-II.
#
# Problem structure:
#   States  X = [t, n, u, v]   (time, lateral offset, long. speed, lat. speed)
#   Controls U = [F_d, F_n]    (longitudinal force, lateral force)
#
# The track arc length s ∈ [0, L] is divided into N_intervals equal intervals.
# Within each interval m, the state is approximated by a degree-K Lagrange
# polynomial at the K LGR nodes tau ∈ [-1,1). The ODE constraint is enforced
# at every collocation point via the LGR differentiation matrix D.
#
# Continuity between intervals is enforced by matching the extrapolated
# end-of-interval state (using the LGR integration weights) to the initial
# state of the next interval.
#=====================================================================================================================================================================================================================

lap_length = track_data.arc_s[-1]
h = lap_length / N_intervals           # Length of each mesh interval [m]

# --- Compute LGR nodes, weights, differentiation matrix (local, on [-1,1)) ---
tau_lgr, w_lgr = lgr_nodes_weights(K)
D_lgr = lgr_differentiation_matrix(tau_lgr)    # K x K

# Map local LGR nodes to arc-length positions within interval m:
#   s_{m,k} = s_m + h/2 * (1 + tau_k)   where s_m = m * h  (left endpoint)
# The derivative transform: d/ds = (2/h) * d/dtau

# Scaling
force_scale = 1 / (vehicle.m * vehicle.g)
length_scale = 1 / vehicle.l
time_scale   = np.sqrt(vehicle.g / vehicle.l)
speed_scale  = 1 / np.sqrt(vehicle.g * vehicle.l)

# --- Interpolation functions for track geometry ---
f_kappa = ca.interpolant('f_kappa', 'bspline', [track_data.arc_s.tolist()], track_data.omega_z.tolist())
f_nl    = ca.interpolant('f_nl',    'linear',  [track_data.arc_s.tolist()], track_data.nl.tolist())
f_nr    = ca.interpolant('f_nr',    'linear',  [track_data.arc_s.tolist()], track_data.nr.tolist())

# --- Continuous dynamics (identical to original, reused symbolically) ---
X_sym = SX.sym('X_sym', 4)
U_sym = SX.sym('U_sym', 2)
s_sym = SX.sym('s')

t_N  = X_sym[0];  n_N = X_sym[1];  u_N = X_sym[2];  v_N = X_sym[3]
F_d_N = U_sym[0]; F_n_N = U_sym[1]

t = t_N / time_scale;   n = n_N / length_scale
u = u_N / speed_scale;  v = v_N / speed_scale
F_d = F_d_N / force_scale;  F_n = F_n_N / force_scale

Drag = (1/2) * vehicle.Cd * vehicle.rho * vehicle.A * (u**2)
Lift = (1/2) * vehicle.Cl * vehicle.rho * vehicle.A * (u**2)
curv = f_kappa(s_sym)

s_dot   = (1 / time_scale) * u / (1 - n * curv)
n_dot_N = (length_scale / time_scale) * v
u_dot_N = (length_scale / time_scale**2) * ((F_d - Drag) / vehicle.m + v * u * curv / (1 - n * curv))
v_dot_N = (length_scale / time_scale**2) * (F_n / vehicle.m - (u**2) * curv / (1 - n * curv))

# ODE expressed as dX/ds  (scaled state, arc-length independent variable)
ODE_s = vertcat(1 / s_dot, n_dot_N / s_dot, u_dot_N / s_dot, v_dot_N / s_dot)
f_ode = Function('f_ode', [X_sym, U_sym, s_sym], [ODE_s])

f_Lift = Function('f_Lift', [X_sym], [Lift])

# ==============================================================================
# Decision variable layout
# ------------------------------------------------------------------------------
# For each interval m = 0 … N_intervals-1:
#   State at K LGR collocation points:  X_col[m, k, :]  (4 values each)
# Initial state of each interval (= final state propagated from previous):
#   X_edge[m, :]  m = 0 … N_intervals   (4 values each)
# Controls (piecewise constant per interval, evaluated at each collocation pt):
#   U_col[m, k, :] (2 values each) — or one control per interval if desired.
#   Here we keep K controls per interval (one per collocation point) for
#   maximum fidelity, matching the state polynomial degree.
# ==============================================================================

n_states   = 4
n_controls = 2
n_col_pts  = N_intervals * K          # total collocation points

# Symbolic variables
# X_edges: state at interval boundaries  (N_intervals+1) x n_states
X_edges = SX.sym('X_edges', n_states, N_intervals + 1)

# X_col: state at internal collocation points  n_states x (N_intervals * K)
X_col = SX.sym('X_col', n_states, n_col_pts)

# U_col: controls at collocation points  n_controls x n_col_pts
U_col = SX.sym('U_col', n_controls, n_col_pts)

all_vars = vertcat(
    reshape(X_edges, -1, 1),
    reshape(X_col,   -1, 1),
    reshape(U_col,   -1, 1)
)
print(f"Total NLP decision variables: {all_vars.numel()}")

# --- Precompute the LGR integration row vector w_int for end-of-interval
#     extrapolation:  X_end = X_edge_m + (h/2) * sum_k w_lgr[k] * F(X_col_mk, ...)
#     This is exact for the polynomial of degree K representing the state.
# --- Also precompute D as a CasADi DM for fast matrix-vector products ---
D_ca  = ca.DM(D_lgr)     # K x K differentiation matrix
w_ca  = ca.DM(w_lgr)     # K x 1 quadrature weights (column)

# --- Assemble constraints and cost ---
g_all  = [];  lbg_all = [];  ubg_all = []
cost   = 0.0

# Regularisation weights (small, following eq. 56 in the paper)
e1 = 1e-3;  e2 = 1e-3

for m in range(N_intervals):
    s_left  = m * h                       # arc-length at left edge of interval
    col_idx = slice(m * K, (m + 1) * K)  # index range in the flat collocation arrays

    X_m  = X_col[:, col_idx]   # 4 x K  states at collocation pts in interval m
    U_m  = U_col[:, col_idx]   # 2 x K  controls at collocation pts
    Xm0  = X_edges[:, m]       # 4      state at left boundary of interval m

    # ------------------------------------------------------------------
    # 1. LGR COLLOCATION CONSTRAINT (core of the pseudospectral method)
    #    D_ca @ X_m^T ≈ (h/2) * F(X_m, U_m, s)
    #    i.e., for each collocation point k:
    #       sum_j D[k,j] * X_m[:,j]  =  (h/2) * f_ode(X_m[:,k], U_m[:,k], s_mk)
    # ------------------------------------------------------------------
    F_col = SX.zeros(n_states, K)
    for k in range(K):
        tau_k = tau_lgr[k]
        s_mk  = s_left + (h / 2) * (1.0 + tau_k)   # arc-length at collocation pt
        F_col[:, k] = f_ode(X_m[:, k], U_m[:, k], s_mk)

    # LGR collocation residual:  D @ X_m.T = (h/2) * F_col.T
    # (K x K) @ (K x 4) = K x 4
    DX  = mtimes(D_ca, X_m.T)        # K x 4
    hF  = (h / 2) * F_col.T          # K x 4

    collocation_residual = DX - hF    # K x 4  (should be zero)
    g_all.append(reshape(collocation_residual, -1, 1))
    lbg_all.append(np.zeros(K * n_states))
    ubg_all.append(np.zeros(K * n_states))

    # ------------------------------------------------------------------
    # 2. INITIAL CONDITION / CONTINUITY CONSTRAINT
    #    The first collocation point (tau = -1) equals the left-edge state.
    #    X_m[:, 0] = X_edges[:, m]
    # ------------------------------------------------------------------
    g_all.append(X_m[:, 0] - Xm0)
    lbg_all.append(np.zeros(n_states))
    ubg_all.append(np.zeros(n_states))

    # ------------------------------------------------------------------
    # 3. INTERVAL END-STATE via LGR quadrature
    #    X_end = Xm0 + (h/2) * sum_k w[k] * F_col[:,k]
    #    This is also the left-edge state of the next interval:
    #    X_edges[:, m+1] = X_end
    # ------------------------------------------------------------------
    X_end = Xm0 + (h / 2) * mtimes(F_col, w_ca)   # 4 x 1
    g_all.append(X_edges[:, m + 1] - X_end)
    lbg_all.append(np.zeros(n_states))
    ubg_all.append(np.zeros(n_states))

    # ------------------------------------------------------------------
    # 4. TIME MONOTONICITY (t must increase along the track)
    # ------------------------------------------------------------------
    g_all.append(X_edges[0, m + 1] - X_edges[0, m])
    lbg_all.append(np.array([1e-6]))
    ubg_all.append(np.array([np.inf]))

    # ------------------------------------------------------------------
    # 5. POWER CONSTRAINT  (at each collocation point)
    # ------------------------------------------------------------------
    for k in range(K):
        power_k = U_m[0, k] * X_m[2, k] / (force_scale * speed_scale)
        g_all.append(power_k)
        lbg_all.append(np.array([vehicle.peakbrakingpower]))
        ubg_all.append(np.array([vehicle.peakdrivingpower]))

    # ------------------------------------------------------------------
    # 6. TIRE FRICTION ELLIPSE CONSTRAINT  (at each collocation point)
    # ------------------------------------------------------------------
    for k in range(K):
        tau_k = tau_lgr[k]
        s_mk  = s_left + (h / 2) * (1.0 + tau_k)
        lift_k = f_Lift(X_m[:, k])
        tire_k = ((U_m[0, k] / force_scale)**2
                + (U_m[1, k] / force_scale)**2
                - (vehicle.mu0 * (vehicle.m * vehicle.g - lift_k))**2)
        g_all.append(tire_k)
        lbg_all.append(np.array([-np.inf]))
        ubg_all.append(np.array([0.0]))

    # ------------------------------------------------------------------
    # 7. COST ACCUMULATION  (regularised lap time, eq. 56 in the paper)
    #    J = integral_0^L [1 + e1*F_d^2 + e2*F_n^2] ds
    #      ≈ sum_m (h/2) * sum_k w[k] * integrand(X_mk, U_mk)
    # ------------------------------------------------------------------
    dt_m = (X_edges[0, m + 1] - X_edges[0, m]) / time_scale
    reg_m = 0
    for k in range(K):
        reg_m += w_lgr[k] * (e1 * U_m[0, k]**2 + e2 * U_m[1, k]**2)
    cost += dt_m + (h / 2) * reg_m

# ------------------------------------------------------------------
# 8. PERIODICITY / CYCLIC CONSTRAINT  (states[n, u, v] match at start/finish)
# ------------------------------------------------------------------
g_all.append(X_edges[[1, 2, 3], -1] - X_edges[[1, 2, 3], 0])
lbg_all.append(np.zeros(3))
ubg_all.append(np.zeros(3))

# Final lap time is an explicit state — add it to the cost
cost += X_edges[0, -1] / time_scale

# Concatenate all constraints
g   = vertcat(*g_all)
lbg = np.concatenate(lbg_all)
ubg = np.concatenate(ubg_all)

# ==============================================================================
# Variable bounds and initial guess
# ==============================================================================
n_edge_vars = n_states * (N_intervals + 1)
n_col_vars  = n_states * n_col_pts
n_ctrl_vars = n_controls * n_col_pts
n_total     = n_edge_vars + n_col_vars + n_ctrl_vars

lbx = np.full(n_total, -np.inf)
ubx = np.full(n_total,  np.inf)
x0  = np.zeros(n_total)

# Helper: flat index for X_edges[state, interval_idx]
def idx_edge(state_i, m):
    return state_i * (N_intervals + 1) + m

# Helper: flat index for X_col[state_i, m*K + k]
def idx_col(state_i, m, k):
    return n_edge_vars + state_i * n_col_pts + m * K + k

# Helper: flat index for U_col[ctrl_i, m*K + k]
def idx_ctrl(ctrl_i, m, k):
    return n_edge_vars + n_col_vars + ctrl_i * n_col_pts + m * K + k

# --- Bounds on edge states ---
for m in range(N_intervals + 1):
    s_m = m * h

    # Time
    if m == 0:
        lbx[idx_edge(0, m)] = 0.0
        ubx[idx_edge(0, m)] = 0.0
    else:
        lbx[idx_edge(0, m)] = 0.0
        ubx[idx_edge(0, m)] = 500 * time_scale

    # Lateral offset
    lbx[idx_edge(1, m)] = float(f_nr(s_m)) * length_scale
    ubx[idx_edge(1, m)] = float(f_nl(s_m)) * length_scale

    # Longitudinal speed
    lbx[idx_edge(2, m)] = 2 * speed_scale
    ubx[idx_edge(2, m)] = vehicle.umax * speed_scale

    # Lateral speed
    lbx[idx_edge(3, m)] = -vehicle.vmax * speed_scale
    ubx[idx_edge(3, m)] =  vehicle.vmax * speed_scale

    # Initial guess: constant forward speed at umax, zero lateral
    x0[idx_edge(0, m)] = time_scale * (s_m / vehicle.umax)
    x0[idx_edge(2, m)] = vehicle.umax * speed_scale

# --- Bounds on collocation states and controls ---
for m in range(N_intervals):
    for k in range(K):
        tau_k = tau_lgr[k]
        s_mk  = m * h + (h / 2) * (1.0 + tau_k)

        # State bounds at collocation points
        lbx[idx_col(0, m, k)] = 0.0
        ubx[idx_col(0, m, k)] = 500 * time_scale
        x0[idx_col(0, m, k)]  = time_scale * (s_mk / vehicle.umax)

        lbx[idx_col(1, m, k)] = float(f_nr(s_mk)) * length_scale
        ubx[idx_col(1, m, k)] = float(f_nl(s_mk)) * length_scale

        lbx[idx_col(2, m, k)] = 2 * speed_scale
        ubx[idx_col(2, m, k)] = vehicle.umax * speed_scale
        x0[idx_col(2, m, k)]  = vehicle.umax * speed_scale

        lbx[idx_col(3, m, k)] = -vehicle.vmax * speed_scale
        ubx[idx_col(3, m, k)] =  vehicle.vmax * speed_scale

        # Control bounds at collocation points
        lbx[idx_ctrl(0, m, k)] = (4 * vehicle.peakbrakingtorque / vehicle.R) * force_scale
        ubx[idx_ctrl(0, m, k)] = (2 * vehicle.peakdrivingtorque / vehicle.R) * force_scale
        x0[idx_ctrl(0, m, k)]  = (2 * vehicle.peakdrivingtorque / vehicle.R) * force_scale

        lat_lim = (vehicle.mu0 * vehicle.m * vehicle.g
                   - 0.5 * vehicle.Cl * vehicle.rho * vehicle.A * vehicle.umax**2) * force_scale
        lbx[idx_ctrl(1, m, k)] = -lat_lim
        ubx[idx_ctrl(1, m, k)] =  lat_lim

# ==============================================================================
# NLP solver (IPOPT)
# ==============================================================================
nlp = {'x': all_vars, 'f': cost, 'g': g}
opts = {
    'ipopt': {
        'max_iter':               4000,
        'print_level':            5,
        'mu_strategy':            'adaptive',
        'acceptable_obj_change_tol': 5e-1,
        'linear_solver':          'mumps',   # change to 'ma57' if HSL available
    }
}
solver = nlpsol('solver', 'ipopt', nlp, opts)
sol    = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
full_sol = np.array(sol['x']).flatten()

# ==============================================================================
# Extract solution
# ==============================================================================
# Edge states: one per interval boundary  (N_intervals+1 points)
X_edge_sol = np.zeros((n_states, N_intervals + 1))
for i in range(n_states):
    for m in range(N_intervals + 1):
        X_edge_sol[i, m] = full_sol[idx_edge(i, m)]

# Collocation states and controls: K per interval → dense arrays
X_col_sol = np.zeros((n_states,   n_col_pts))
U_col_sol = np.zeros((n_controls, n_col_pts))
for i in range(n_states):
    for j in range(n_col_pts):
        X_col_sol[i, j] = full_sol[n_edge_vars + i * n_col_pts + j]
for i in range(n_controls):
    for j in range(n_col_pts):
        U_col_sol[i, j] = full_sol[n_edge_vars + n_col_vars + i * n_col_pts + j]

# s-grid corresponding to edge points and collocation points
s_edges = np.linspace(0, lap_length, N_intervals + 1)
s_col   = np.array([m * h + (h / 2) * (1.0 + tau_lgr[k])
                    for m in range(N_intervals) for k in range(K)])

# Unscale
time_opt = X_edge_sol[0, :] / time_scale
n_opt    = X_edge_sol[1, :] / length_scale
u_opt    = X_edge_sol[2, :] / speed_scale
v_opt    = X_edge_sol[3, :] / speed_scale

F_d_opt  = U_col_sol[0, :] / force_scale
F_n_opt  = U_col_sol[1, :] / force_scale
u_col    = X_col_sol[2, :]  / speed_scale

Power_opt = F_d_opt * u_col

print(f"\nOptimized Lap Time (LGR): {time_opt[-1]:.4f} s")

#====================================================================================================================================================================================================================
# Frenet Coordinates to Cartesian Coordinates
#====================================================================================================================================================================================================================
f_normal1 = ca.interpolant('f_normal1', 'linear', [track_data.arc_s.tolist()], track_data.n[0, :].tolist())
f_normal2 = ca.interpolant('f_normal2', 'linear', [track_data.arc_s.tolist()], track_data.n[1, :].tolist())
f_normal3 = ca.interpolant('f_normal3', 'linear', [track_data.arc_s.tolist()], track_data.n[2, :].tolist())
f_xc = ca.interpolant('f_xc', 'linear', [track_data.arc_s.tolist()], track_data.xc.tolist())
f_yc = ca.interpolant('f_yc', 'linear', [track_data.arc_s.tolist()], track_data.yc.tolist())
f_zc = ca.interpolant('f_zc', 'linear', [track_data.arc_s.tolist()], track_data.zc.tolist())

x_opt = np.array(f_xc(s_edges)).flatten() + n_opt * np.array(f_normal1(s_edges)).flatten()
y_opt = np.array(f_yc(s_edges)).flatten() + n_opt * np.array(f_normal2(s_edges)).flatten()
z_opt = np.array(f_zc(s_edges)).flatten() + n_opt * np.array(f_normal3(s_edges)).flatten()

#====================================================================================================================================================================================================================
# Results and Plots (identical to original)
#====================================================================================================================================================================================================================

## Figure 1 - Racing Line
fig = go.Figure()
fig.add_trace(go.Scatter3d(x=track_data.br[0, :], y=track_data.br[1, :], z=track_data.br[2, :],
    mode='lines', name='Right Boundary', line=dict(color='blue', width=4)))
fig.add_trace(go.Scatter3d(x=track_data.bl[0, :], y=track_data.bl[1, :], z=track_data.bl[2, :],
    mode='lines', name='Left Boundary',  line=dict(color='blue', width=4)))
fig.add_trace(go.Scatter3d(x=track_data.xc, y=track_data.yc, z=track_data.zc,
    mode='lines', name='Centreline',     line=dict(color='red', width=4, dash='dashdot')))
fig.add_trace(go.Scatter3d(x=x_opt, y=y_opt, z=z_opt,
    mode='lines', name='Optimal Racing Line (LGR)', line=dict(color='black', width=4)))
fig.update_layout(title="Circuit de Barcelona-Catalunya — LGR Optimal Racing Line",
    scene=dict(xaxis_title='x (m)', yaxis_title='y (m)', zaxis_title='Elevation (m)', aspectmode='data'),
    margin=dict(l=0, r=0, b=0, t=50))
fig.show()

## Figure 2 - Forces
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
ax1.plot(s_col, F_d_opt, color='black')
ax1.set_ylabel('Drive Force [N]');  ax1.grid(True, alpha=0.3);  ax1.axhline(0, color='red')
ax2.plot(s_col, F_n_opt, color='black')
ax2.set_xlabel('Arc Length [m]');   ax2.set_ylabel('Lateral Force [N]')
ax2.grid(True, alpha=0.3);          ax2.axhline(0, color='red')
plt.tight_layout();  plt.show()

## Figure 3 - States
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
ax1.plot(s_edges, n_opt, color='black')
ax1.plot(track_data.arc_s, track_data.nl, color='blue',  label='Left Boundary')
ax1.plot(track_data.arc_s, track_data.nr, color='red',   label='Right Boundary')
ax1.set_ylabel('Lateral Offset [m]');  ax1.grid(True, alpha=0.3);  ax1.legend()
ax2.plot(s_edges, u_opt * (18/5), color='black')
ax2.set_ylabel('Long. Speed [km/h]');  ax2.grid(True, alpha=0.3)
ax2.axhline(vehicle.umax * (18/5), color='red')
ax3.plot(s_edges, v_opt * (18/5), color='black')
ax3.set_xlabel('Arc Length [m]');  ax3.set_ylabel('Lat. Speed [km/h]')
ax3.grid(True, alpha=0.3)
ax3.axhline( vehicle.vmax * (18/5), color='red')
ax3.axhline(-vehicle.vmax * (18/5), color='red')
plt.tight_layout();  plt.show()

## Figure 4 - Power
fig, ax = plt.subplots(1, 1, figsize=(10, 4))
ax.plot(s_col, Power_opt / 1000, color='black')
ax.set_xlabel('Arc Length [m]');  ax.set_ylabel('Power [kW]')
ax.grid(True, alpha=0.3);         ax.axhline(0, color='red')
plt.tight_layout();  plt.show()