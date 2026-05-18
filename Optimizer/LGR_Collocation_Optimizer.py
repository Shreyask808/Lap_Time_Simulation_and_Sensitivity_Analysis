#====================================================================================================================================================================================
# Legendre-Gauss-Radau (LGR) Pseudospectral Optimal Control
# Reference: Limebeer & Perantoni (2015), JDSMC 137, DOI: 10.1115/1.4029466
#            Patterson & Rao (2013), GPOPS-II
#
# Variable layout (following GPOPS-II convention):
#   - X_edge[m]  : state at the LEFT boundary of interval m  (m = 0..N)
#                  X_edge[0] is the initial state (t=0 fixed)
#                  X_edge[N] is the final state (periodic match with X_edge[0])
#   - X_col[m,k] : state at the K-1 INTERIOR LGR collocation points of interval m
#                  (tau[0]=-1 is NOT duplicated here; it is represented by X_edge[m])
#   - U_col[m,k] : control at each of the K nodes (incl. left edge node)
#
# Collocation constraint (augmented differentiation matrix):
#   The degree-K polynomial over interval m is pinned at K+1... wait, K nodes:
#     node 0 (tau=-1) = X_edge[m]          <- boundary, shared with prev interval
#     nodes 1..K-1    = X_col[m, 0..K-2]   <- interior, free variables
#   D_aug (shape (K-1) x K) maps the K nodal values to the K-1 interior derivatives:
#     D_aug @ [X_edge[m]; X_col[m].T]  =  h2 * F_col_interior.T
#
# End-of-interval propagation (LGR quadrature):
#   X_edge[m+1] = X_edge[m] + h2 * F_all @ w_lgr
#   where F_all contains the ODE rhs at all K nodes of interval m.
#====================================================================================================================================================================================

import os
import numpy as np
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from types import SimpleNamespace
import json
import casadi as ca
from casadi import *
import matplotlib
import pickle
matplotlib.use('Qt5Agg')
import plotly.graph_objects as go
plt.close('all')
os.system('cls' if os.name == 'nt' else 'clear')

# =====================================================================================
# LGR nodes, weights, and differentiation matrices
# =====================================================================================

def lgr_nodes_weights(K):
    """
    K-point LGR rule on [-1, 1).
    tau[0] = -1 (fixed left Radau point).
    Interior nodes tau[1..K-1] solve  P_K(s) + P_{K-1}(s) = 0.
    Returns tau (K,) and w (K,) with sum(w) = 2.
    """
    if K == 1:
        return np.array([-1.0]), np.array([2.0])

    from numpy.polynomial.legendre import leggauss

    def legendre(n, x):
        """Return P_n(x) and P_{n-1}(x) via three-term recurrence."""
        if n == 0:
            return 1.0, 0.0
        p0, p1 = 1.0, float(x)
        for k in range(1, n):
            p0, p1 = p1, ((2*k + 1)*x*p1 - k*p0) / (k + 1)
        return p1, p0

    tau = np.empty(K)
    tau[0] = -1.0

    # Newton iteration for the K-1 interior nodes
    gl, _ = leggauss(K - 1)
    for i, s0 in enumerate(sorted(gl)):
        s = float(s0)
        for _ in range(200):
            PN,  PN1 = legendre(K,     s)
            PK1, PK2 = legendre(K - 1, s)
            f = PN + PK1
            denom = 1.0 - s**2
            if abs(denom) < 1e-15:
                break
            dPN  = K       * (PN1  - s * PN)  / denom
            dPK1 = (K - 1) * (PK2  - s * PK1) / denom
            df   = dPN + dPK1
            if abs(df) < 1e-15:
                break
            step = f / df
            s -= step
            if abs(step) < 1e-14:
                break
        tau[i + 1] = s
    tau[1:] = np.sort(tau[1:])

    # LGR weights
    w = np.empty(K)
    w[0] = 2.0 / (K * K)
    for i in range(1, K):
        PK1, _ = legendre(K - 1, tau[i])
        w[i] = (1.0 - tau[i]) / (K * PK1) ** 2

    return tau, w


def lgr_differentiation_matrix(tau):
    """Full K x K differentiation matrix for LGR nodes tau."""
    K = len(tau)
    D = np.zeros((K, K))
    for i in range(K):
        for j in range(K):
            if i != j:
                num = np.prod([tau[i] - tau[m] for m in range(K) if m != j])
                den = np.prod([tau[j] - tau[m] for m in range(K) if m != j])
                D[i, j] = num / (den * (tau[i] - tau[j]))
        D[i, i] = -np.sum(D[i, :])
    return D


# =====================================================================================
# User inputs
# =====================================================================================

N = 200   # mesh intervals
K = 4     # LGR nodes per interval (incl. tau=-1 left boundary)
          # interior collocation points per interval = K-1 = 3

# =====================================================================================
# Load data
# =====================================================================================

root = tk.Tk(); root.withdraw(); root.attributes('-topmost', True)
fp1  = filedialog.askopenfilename(title="Select Track Data",
                                  filetypes=[("pickle files", "*.pkl")])
if fp1:
    with open(fp1, 'rb') as f:
        track_data = pickle.load(f)
    print(f"Loaded track  |  {len(track_data.arc_s)} points")
else:
    raise FileNotFoundError("Track data not found")

root = tk.Tk(); root.withdraw(); root.attributes('-topmost', True)
fp2  = filedialog.askopenfilename(title="Select Vehicle Data",
                                  filetypes=[("json files", "*.json")])
if fp2:
    with open(fp2, 'r') as f:
        vehicle = SimpleNamespace(**json.load(f))
    print(f"Loaded vehicle  |  mass = {vehicle.m} kg")
else:
    raise FileNotFoundError("Vehicle data not found")

# =====================================================================================
# LGR precomputation
# =====================================================================================

tau, w_lgr = lgr_nodes_weights(K)
# tau[0] = -1  (left boundary node, mapped to X_edge[m])
# tau[1:] = K-1 interior collocation nodes

D_full  = lgr_differentiation_matrix(tau)   # K x K
D_aug   = D_full[1:, :]                     # (K-1) x K  — rows for interior nodes only
D_aug_ca = ca.DM(D_aug)
w_ca     = ca.DM(w_lgr.reshape(-1, 1))      # K x 1

Kc = K - 1   # interior collocation points per interval

# =====================================================================================
# Problem geometry and scaling
# =====================================================================================

lap_length = float(track_data.arc_s[-1])
h  = lap_length / N
h2 = h / 2.0

force_scale  = 1.0 / (vehicle.m * vehicle.g)
length_scale = 1.0 / vehicle.l
time_scale   = np.sqrt(vehicle.g / vehicle.l)
speed_scale  = 1.0 / np.sqrt(vehicle.g * vehicle.l)

# Precompute s-coordinates for all K nodes of every interval
# s_nodes_all[m, k] = arc-length at node k of interval m
s_nodes_all = np.array([
    m * h + h2 * (1.0 + tau)   # shape (K,)
    for m in range(N)
])   # shape (N, K)

# =====================================================================================
# Track interpolants
# =====================================================================================

f_kappa = ca.interpolant('f_kappa', 'bspline',
                         [track_data.arc_s.tolist()], track_data.omega_z.tolist())
f_nl    = ca.interpolant('f_nl', 'linear',
                         [track_data.arc_s.tolist()], track_data.nl.tolist())
f_nr    = ca.interpolant('f_nr', 'linear',
                         [track_data.arc_s.tolist()], track_data.nr.tolist())

# =====================================================================================
# Continuous ODE  f(X, U, s)  →  dX/ds  (scaled)
# =====================================================================================

X_sym = SX.sym('X_sym', 4)
U_sym = SX.sym('U_sym', 2)
s_sym = SX.sym('s')

n_s = X_sym[1] / length_scale
u_s = X_sym[2] / speed_scale
v_s = X_sym[3] / speed_scale
Fd_s = U_sym[0] / force_scale
Fn_s = U_sym[1] / force_scale

Drag  = 0.5 * vehicle.Cd * vehicle.rho * vehicle.A * u_s**2
Lift  = 0.5 * vehicle.Cl * vehicle.rho * vehicle.A * u_s**2
kappa = f_kappa(s_sym)
denom = 1.0 - n_s * kappa

s_dot   = (1.0 / time_scale) * u_s / denom
n_dot_N = (length_scale / time_scale) * v_s
u_dot_N = (length_scale / time_scale**2) * ((Fd_s - Drag) / vehicle.m
           + v_s * u_s * kappa / denom)
v_dot_N = (length_scale / time_scale**2) * (Fn_s / vehicle.m
           - u_s**2 * kappa / denom)

ODE = vertcat(1.0 / s_dot, n_dot_N / s_dot, u_dot_N / s_dot, v_dot_N / s_dot)
f_ode  = Function('f_ode',  [X_sym, U_sym, s_sym], [ODE])
f_Lift = Function('f_Lift', [X_sym],               [Lift])

# =====================================================================================
# Decision variables
#   X_edge : 4 x (N+1)         state at interval left boundaries + final boundary
#   X_col  : 4 x (N * Kc)      state at interior collocation points (Kc = K-1 per interval)
#   U_col  : 2 x (N * K)       control at all K nodes per interval (incl. left boundary node)
# =====================================================================================

n_xe = 4 * (N + 1)
n_xc = 4 * N * Kc
n_uc = 2 * N * K

X_edge = SX.sym('X_edge', 4, N + 1)
X_col  = SX.sym('X_col',  4, N * Kc)
U_col  = SX.sym('U_col',  2, N * K)

all_vars = vertcat(reshape(X_edge, -1, 1),
                   reshape(X_col,  -1, 1),
                   reshape(U_col,  -1, 1))
print(f"NLP variables: {all_vars.numel()}")

# Flat index helpers
def ie(si, m):   return si * (N + 1) + m
def ic(si, m, k): return n_xe + si * (N * Kc) + m * Kc + k
def iu(ci, m, k): return n_xe + n_xc + ci * (N * K) + m * K + k

# =====================================================================================
# Constraint and cost assembly
# =====================================================================================

g_all   = []
lbg_all = []
ubg_all = []
cost    = 0.0
e1, e2  = 1e-3, 1e-3

for m in range(N):
    sn    = s_nodes_all[m]       # arc-length at K nodes of this interval
    Xe_m  = X_edge[:, m]         # 4 — left boundary state
    Xe_m1 = X_edge[:, m + 1]     # 4 — right boundary state
    Xc_m  = X_col[:, m*Kc : (m+1)*Kc]   # 4 x Kc  interior collocation states
    Uc_m  = U_col[:, m*K  : (m+1)*K ]   # 2 x K   controls at all K nodes

    # Full K-node state matrix: [Xe_m | Xc_m]  shape 4 x K
    Xall  = horzcat(Xe_m, Xc_m)          # 4 x K

    # ODE rhs at all K nodes
    Fall  = horzcat(*[f_ode(Xall[:, k], Uc_m[:, k], sn[k]) for k in range(K)])  # 4 x K

    # ------------------------------------------------------------------
    # 1. COLLOCATION  (enforces the ODE at K-1 interior nodes)
    #    D_aug @ Xall.T  =  h2 * Fall[:, 1:].T
    #    Shape: (K-1 x K) @ (K x 4)  =  (K-1) x 4
    # ------------------------------------------------------------------
    DX = mtimes(D_aug_ca, Xall.T)              # (K-1) x 4
    hF = h2 * Fall[:, 1:].T                   # (K-1) x 4
    g_all.append(reshape(DX - hF, -1, 1))
    lbg_all.append(np.zeros(Kc * 4))
    ubg_all.append(np.zeros(Kc * 4))

    # ------------------------------------------------------------------
    # 2. CONTINUITY  (quadrature propagation to next edge)
    #    X_edge[m+1] = X_edge[m] + h2 * Fall @ w_lgr
    # ------------------------------------------------------------------
    X_end = Xe_m + h2 * mtimes(Fall, w_ca)    # 4 x 1
    g_all.append(Xe_m1 - X_end)
    lbg_all.append(np.zeros(4))
    ubg_all.append(np.zeros(4))

    # ------------------------------------------------------------------
    # 3. TIME MONOTONICITY
    # ------------------------------------------------------------------
    g_all.append(Xe_m1[0] - Xe_m[0])
    lbg_all.append(np.array([1e-6]))
    ubg_all.append(np.array([np.inf]))

    # ------------------------------------------------------------------
    # 4. POWER CONSTRAINT  (at each of the K nodes)
    # ------------------------------------------------------------------
    for k in range(K):
        pw = Uc_m[0, k] * Xall[2, k] / (force_scale * speed_scale)
        g_all.append(pw)
        lbg_all.append(np.array([vehicle.peakbrakingpower]))
        ubg_all.append(np.array([vehicle.peakdrivingpower]))

    # ------------------------------------------------------------------
    # 5. TIRE FRICTION ELLIPSE  (at each of the K nodes)
    # ------------------------------------------------------------------
    for k in range(K):
        mu_Fz = vehicle.mu0 * (vehicle.m * vehicle.g - f_Lift(Xall[:, k]))
        tc = (Uc_m[0, k] / force_scale)**2 + (Uc_m[1, k] / force_scale)**2 - mu_Fz**2
        g_all.append(tc)
        lbg_all.append(np.array([-np.inf]))
        ubg_all.append(np.array([0.0]))

    # ------------------------------------------------------------------
    # 6. COST  (regularised lap time, eq. 56 in paper)
    # ------------------------------------------------------------------
    dt_m  = (Xe_m1[0] - Xe_m[0]) / time_scale
    reg_m = sum([w_lgr[k] * (e1 * Uc_m[0, k]**2 + e2 * Uc_m[1, k]**2)
                 for k in range(K)])
    cost += dt_m + h2 * reg_m

# ------------------------------------------------------------------
# 7. PERIODICITY  (n, u, v match across start/finish)
# ------------------------------------------------------------------
g_all.append(X_edge[[1, 2, 3], -1] - X_edge[[1, 2, 3], 0])
lbg_all.append(np.zeros(3))
ubg_all.append(np.zeros(3))

g   = vertcat(*g_all)
lbg = np.concatenate(lbg_all)
ubg = np.concatenate(ubg_all)

# =====================================================================================
# Variable bounds and warm initial guess
# =====================================================================================

lbx = np.full(all_vars.numel(), -np.inf)
ubx = np.full(all_vars.numel(),  np.inf)
x0  = np.zeros(all_vars.numel())

# Reference: drive at 60% of umax; compute drag-balancing force as control guess
u_ref    = 0.6 * vehicle.umax
Drag_ref = 0.5 * vehicle.Cd * vehicle.rho * vehicle.A * u_ref**2
Fd_ref   = Drag_ref * force_scale          # scaled control that sustains u_ref
Fd_max   = (2 * vehicle.peakdrivingtorque / vehicle.R) * force_scale
Fd_min   = (4 * vehicle.peakbrakingtorque / vehicle.R) * force_scale
Fn_lim   = (vehicle.mu0 * vehicle.m * vehicle.g
            - 0.5 * vehicle.Cl * vehicle.rho * vehicle.A * u_ref**2) * force_scale

s_edges = np.linspace(0, lap_length, N + 1)

# --- Edge state bounds and guess ---
for m in range(N + 1):
    sm   = s_edges[m]
    t_m  = sm / u_ref
    nl_m = float(f_nl(sm))
    nr_m = float(f_nr(sm))

    # time: fix t=0 at start
    lbx[ie(0, m)] = 0.0 if m == 0 else 0.0
    ubx[ie(0, m)] = 0.0 if m == 0 else 500 * time_scale
    x0 [ie(0, m)] = t_m * time_scale

    # lateral offset: centreline (always feasible)
    lbx[ie(1, m)] = nr_m * length_scale
    ubx[ie(1, m)] = nl_m * length_scale
    x0 [ie(1, m)] = 0.0

    # longitudinal speed: non-trivial guess at u_ref
    lbx[ie(2, m)] = 2.0 * speed_scale
    ubx[ie(2, m)] = vehicle.umax * speed_scale
    x0 [ie(2, m)] = u_ref * speed_scale

    # lateral speed: zero
    lbx[ie(3, m)] = -vehicle.vmax * speed_scale
    ubx[ie(3, m)] =  vehicle.vmax * speed_scale
    x0 [ie(3, m)] = 0.0

# Fix t=0 at start
lbx[ie(0, 0)] = 0.0
ubx[ie(0, 0)] = 0.0

# --- Interior collocation state bounds and guess ---
for m in range(N):
    sn = s_nodes_all[m]
    for k in range(Kc):
        s_mk = sn[k + 1]          # skip tau[0]=-1 (that's the edge)
        t_mk = s_mk / u_ref
        nl_mk = float(f_nl(s_mk))
        nr_mk = float(f_nr(s_mk))

        lbx[ic(0, m, k)] = 0.0
        ubx[ic(0, m, k)] = 500 * time_scale
        x0 [ic(0, m, k)] = t_mk * time_scale

        lbx[ic(1, m, k)] = nr_mk * length_scale
        ubx[ic(1, m, k)] = nl_mk * length_scale
        x0 [ic(1, m, k)] = 0.0

        lbx[ic(2, m, k)] = 2.0 * speed_scale
        ubx[ic(2, m, k)] = vehicle.umax * speed_scale
        x0 [ic(2, m, k)] = u_ref * speed_scale

        lbx[ic(3, m, k)] = -vehicle.vmax * speed_scale
        ubx[ic(3, m, k)] =  vehicle.vmax * speed_scale
        x0 [ic(3, m, k)] = 0.0

# --- Control bounds and guess ---
for m in range(N):
    for k in range(K):
        lbx[iu(0, m, k)] = Fd_min
        ubx[iu(0, m, k)] = Fd_max
        x0 [iu(0, m, k)] = Fd_ref   # sustain speed — non-zero, non-trivial

        lbx[iu(1, m, k)] = -Fn_lim
        ubx[iu(1, m, k)] =  Fn_lim
        x0 [iu(1, m, k)] = 0.0

# =====================================================================================
# Solve
# =====================================================================================

nlp  = {'x': all_vars, 'f': cost, 'g': g}
opts = {
    'ipopt': {
        'max_iter':                  3000,
        'print_level':               5,
        'mu_strategy':               'adaptive',
        'acceptable_obj_change_tol': 1e-2,
        'linear_solver':             'mumps',
    }
}
solver = nlpsol('solver', 'ipopt', nlp, opts)
sol    = solver(x0=x0, lbx=lbx, ubx=ubx, lbg=lbg, ubg=ubg)
sv     = np.array(sol['x']).flatten()

# =====================================================================================
# Extract solution
# =====================================================================================

time_opt = np.array([sv[ie(0, m)] / time_scale  for m in range(N + 1)])
n_opt    = np.array([sv[ie(1, m)] / length_scale for m in range(N + 1)])
u_opt    = np.array([sv[ie(2, m)] / speed_scale  for m in range(N + 1)])
v_opt    = np.array([sv[ie(3, m)] / speed_scale  for m in range(N + 1)])

s_col_all, F_d_all, F_n_all, u_col_all = [], [], [], []
for m in range(N):
    sn = s_nodes_all[m]
    for k in range(K):
        s_col_all.append(sn[k])
        F_d_all.append(sv[iu(0, m, k)] / force_scale)
        F_n_all.append(sv[iu(1, m, k)] / force_scale)
        u_col_all.append(sv[ie(2, m)] / speed_scale if k == 0
                         else sv[ic(2, m, k - 1)] / speed_scale)

s_col_all = np.array(s_col_all)
F_d_all   = np.array(F_d_all)
F_n_all   = np.array(F_n_all)
Power_opt = F_d_all * np.array(u_col_all)

print(f"\nOptimised Lap Time (LGR): {time_opt[-1]:.4f} s")

# =====================================================================================
# Cartesian racing line
# =====================================================================================

f_n1 = ca.interpolant('f_n1','linear',[track_data.arc_s.tolist()],track_data.n[0,:].tolist())
f_n2 = ca.interpolant('f_n2','linear',[track_data.arc_s.tolist()],track_data.n[1,:].tolist())
f_n3 = ca.interpolant('f_n3','linear',[track_data.arc_s.tolist()],track_data.n[2,:].tolist())
f_xc = ca.interpolant('f_xc','linear',[track_data.arc_s.tolist()],track_data.xc.tolist())
f_yc = ca.interpolant('f_yc','linear',[track_data.arc_s.tolist()],track_data.yc.tolist())
f_zc = ca.interpolant('f_zc','linear',[track_data.arc_s.tolist()],track_data.zc.tolist())

x_opt = np.array(f_xc(s_edges)).flatten() + n_opt*np.array(f_n1(s_edges)).flatten()
y_opt = np.array(f_yc(s_edges)).flatten() + n_opt*np.array(f_n2(s_edges)).flatten()
z_opt = np.array(f_zc(s_edges)).flatten() + n_opt*np.array(f_n3(s_edges)).flatten()

# =====================================================================================
# Plots
# =====================================================================================

fig = go.Figure()
fig.add_trace(go.Scatter3d(x=track_data.br[0,:],y=track_data.br[1,:],z=track_data.br[2,:],
    mode='lines', name='Right boundary', line=dict(color='blue', width=4)))
fig.add_trace(go.Scatter3d(x=track_data.bl[0,:],y=track_data.bl[1,:],z=track_data.bl[2,:],
    mode='lines', name='Left boundary',  line=dict(color='blue', width=4)))
fig.add_trace(go.Scatter3d(x=track_data.xc, y=track_data.yc, z=track_data.zc,
    mode='lines', name='Centreline', line=dict(color='red', width=4, dash='dashdot')))
fig.add_trace(go.Scatter3d(x=x_opt, y=y_opt, z=z_opt,
    mode='lines', name='Optimal line (LGR)', line=dict(color='black', width=4)))
fig.update_layout(title="LGR Optimal Racing Line",
    scene=dict(xaxis_title='x (m)', yaxis_title='y (m)',
               zaxis_title='Elevation (m)', aspectmode='data'),
    margin=dict(l=0,r=0,b=0,t=50))
fig.show()

fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
axes[0].plot(s_col_all, F_d_all, 'k'); axes[0].axhline(0,'r')
axes[0].set_ylabel('Drive Force [N]'); axes[0].grid(True, alpha=0.3)
axes[1].plot(s_col_all, F_n_all, 'k'); axes[1].axhline(0,'r')
axes[1].set_ylabel('Lateral Force [N]'); axes[1].set_xlabel('Arc Length [m]')
axes[1].grid(True, alpha=0.3)
plt.tight_layout(); plt.show()

fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
axes[0].plot(s_edges, n_opt, 'k')
axes[0].plot(track_data.arc_s, track_data.nl, 'b', label='Left')
axes[0].plot(track_data.arc_s, track_data.nr, 'r', label='Right')
axes[0].set_ylabel('Lateral Offset [m]'); axes[0].grid(True, alpha=0.3); axes[0].legend()
axes[1].plot(s_edges, u_opt*(18/5), 'k')
axes[1].axhline(vehicle.umax*(18/5), color='red')
axes[1].set_ylabel('Long. Speed [km/h]'); axes[1].grid(True, alpha=0.3)
axes[2].plot(s_edges, v_opt*(18/5), 'k')
axes[2].axhline( vehicle.vmax*(18/5), color='red')
axes[2].axhline(-vehicle.vmax*(18/5), color='red')
axes[2].set_ylabel('Lat. Speed [km/h]'); axes[2].set_xlabel('Arc Length [m]')
axes[2].grid(True, alpha=0.3)
plt.tight_layout(); plt.show()

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(s_col_all, Power_opt/1000, 'k'); ax.axhline(0, color='red')
ax.set_ylabel('Power [kW]'); ax.set_xlabel('Arc Length [m]'); ax.grid(True, alpha=0.3)
plt.tight_layout(); plt.show()