%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This project utilizes CasADi for symbolic differentiation and optimization.
% Andersson, J. A. E., et al. (2019). "CasADi: a software framework for 
% nonlinear optimization and optimal control." Mathematical Programming Computation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%% 7 DOF Plannar Vehicle Handling Model 
function [N,n_states,u_states,f_dynamics] = Handling_Model(Car,track_data,N)
import casadi.*

%% Direct Multiple Shooting Method Inputs
n_states = 10;                                                              % Number of States (Time, )
u_states = 5;                                                               % Number of Control Inputs

% Independent Variable
lap_length = track_data(arc_s(end));                                        % Total length of the lap (m)
s_grid = linspace(0,lap_length,N+1);                                        % Discretized Centerline Arc Length (m)
ds = lap_length/N;                                                          % Track Centerline Arc Length Step (m)
f_kappa = casadi.interpolant('f_kappa', 'bspline', {track_data.arc_s},track_data.Omega_z);
f_nl = casadi.interpolant('f_nl', 'bspline', {track_data.arc_s},track_data.nl);
f_nr = casadi.interpolant('f_nr', 'bspline', {track_data.arc_s},track_data.nr);
f_theta = casadi.interpolant('f_theta', 'bspline', {track_data.arc_s},track_data.theta);

X = SX.sym('X',n_states,N+1);
U = SX.sym('U',u_states,N);
X_reshape = reshape(X,n_states*(N+1),1);                                    % Column vector of all States
U_reshape = reshape(U,u_states*(N),1);                                      % Column vector of all Control Inputs
states = [X_reshape; U_reshape];                                            % Direct Multiple Shooting States

%% System Definition
X_sym = SX.sym('X_sym',n_states);
U_sym = SX.sym('U_sym',u_states);
s_sym = SX.sym('s_sym');

% Scale States
length_scale = 1/Car.l;                                                     % Length Scaling (m^-1) 
time_scale = sqrt(Car.g/Car.l);                                             % Time Scale (s^-1)
speed_scale = 1/sqrt(Car.g*Car.l);                                          % Speed Scale (s/m)
angle_scale = 1;                                                            % angle Scale (rad^-1)
force_scale = 1/(Car.m*Car.g);                                              % Force Scale (N^-1)

% Normalized States
t_N = X_sym(1);
n_N = X_sym(2);
psi_N = X_sym(3);
psi_dot_N = X_sym(4);
u_N = X_sym(5);
v_N = X_sym(6);
O_fl_N = X_sym(7);
O_fr_N = X_sym(8);
O_rl_N = X_sym(9);
O_rr_N = X_sym(10);

% Normalized Control Inputs
delta_N = U_sym(1);
Md_fl_N = U_sym(2);
Md_fr_N = U_sym(3);
Md_rl_N = U_sym(4);
Md_rr_N = U_sym(5);

%% Dynamics Definition
t = t_N/time_scale;
n = n_N/length_scale;

psi = psi_N/angle_scale;
psi_dot = psi_dot_N*time_scale/angle_scale;

u = u_N/speed_scale;
v = v_N/speed_scale;

O_fl = O_fl_N*time_scale/angle_scale;
O_fr = O_fr_N*time_scale/angle_scale;
O_rl = O_rl_N*time_scale/angle_scale;
O_rr = O_rr_N*time_scale/angle_scale;

delta = delta_N/angle_scale;
Md_fl = Md_fl_N/(length_scale*force_scale);
Md_fr = Md_fr_N/(length_scale*force_scale);
Md_rl = Md_rl_N/(length_scale*force_scale);
Md_rr = Md_rr_N/(length_scale*force_scale);

% Symbolic Track Data
curv = f_kappa(s_sym);
theta = f_theta(s_sym);
xi = psi - theta;

%% Tire Velocities in Vehile Frame
% Rear Right tire
u_rr = u + (Car.Wr/2)*psi_dot;
v_rr = v - (Car.d)*psi_dot;

% Rear Left tire
u_rl = u - (Car.Wr/2)*psi_dot;
v_rl = v - (Car.d)*psi_dot;

% Front Left Tire
u_fl = u - (Car.Wf/2)*psi_dot;
v_fl = v - (Car.l - Car.d)*psi_dot;

% Front Right Tire
u_fr = u + (Car.Wf/2)*psi_dot;
v_fr = v - (Car.l - Car.d)*psi_dot;

%% Tire Slip Angle Definition
alpha_rr = atan(v_rr/u_rr);                                                 % Slip Angle of Rear Right Tire
alpha_rl = atan(v_rl/u_rl);                                                 % Slip Angle of Rear Left Tire
alpha_fl = atan(v_fl/u_fl) + delta;                                         % Slip Angle of Front Left Tire
alpha_fr = atan(v_fr/u_fr) + delta;                                         % Slip Angle of Front Right Tire

%% Normal Forces Definition
Drag = -(1/2)*Car.Cd*Car.rho*Car.A*(u^2);
Downforce = -(1/2)*Car.Cl*Car.rho*Car.A*(u^2);
Fzrr = SX.sym('Fzrr');
Fzrl = SX.sym('Fzrl');
Fzfr = SX.sym('Fzfr');
Fzfl = SX.sym('Fzfl');

%% Tire Forces
% Rear Right Tire Forces
[Fxrr,Fyrr,Mzrr,p_trailrr] = combined_brsuh_model(u_rr,v_rr,Car.R,Car.w,Car.C_p,Car.kv,O_rr,alpha_rr,Fzrr,Car.mu0,Car.mu);

% Rear Left Tire Forces
[Fxrl,Fyrl,Mzrl,p_trailrl] = combined_brsuh_model(u_rl,v_rl,Car.R,Car.w,Car.C_p,Car.kv,O_rl,alpha_rl,Fzrl,Car.mu0,Car.mu);

% Front Right Tire Forces
[Fxfr,Fyfr,Mzfr,p_trailfr] = combined_brsuh_model(u_fr,v_fr,Car.R,Car.w,Car.C_p,Car.kv,O_fr,alpha_fr,Fzfr,Car.mu0,Car.mu);

% Front Right Tire Forces
[Fxfl,Fyfl,Mzfl,p_trailfl] = combined_brsuh_model(u_fl,v_fl,Car.R,Car.w,Car.C_p,Car.kv,O_fl,alpha_fl,Fzfl,Car.mu0,Car.mu);

%% ODE Definitio (Dynamics) 
s_dot = (1/t_scale)*(sigma/(1 - n*curv))*(cos(xi) - (d/l)*sin(xi)*tan(delta));
n_N_dot = (n_scale/t_scale)*(sigma)*(sin(xi) + (d/l)*cos(xi)*tan(delta));
psi_N_dot = psi_dot_N;
psi_ddot_N = (1/2)*(-2*Car.d*(Fyrl + Fyrr) + 2*(Md_fl + Md_fr + Md_rl + Md_rr) + (Fxrr - Fxrl)*Car.Wr + (2*(Car.d - Car.l)*(Fyfl + Fyfr) + (Fxfr - Fxfl)*Car.Wf)*cos(delta) + (2*(Car.d - Car.l)*(Fxfl + Fxfr) + (Fyfl - Fyfr)*Car.Wf)*sin(delta))/Car.Iz;
u_N_dot = (Drag + (Fxrl + Fxrr + (Fxfl + Fxfr)*cos(delta)) - (Fyfr*sin(delta) + Fyfl*sin(delta)))/Car.m + v*psi_dot;
v_N_dot = (Fyrl + Fyrr + (Fyfl + Fyfr)*cos(delta) + (Fxfl + Fxfr)*sin(delta))/Car.m - u*psi_dot;
O_fl_N_dot = Md_fl_N/(force_scale*length_scale) - Crr*Fflz - Fflx*rfl;
O_fr_N_dot = Md_fr_N/(force_scale*length_scale) - Crr*Ffrz - Ffrx*rfr;
O_rl_N_dot = Md_rl_N/(force_scale*length_scale) - Crr*Frlz - Frlx*rrl;
O_rr_N_dot = Md_rr_N/(force_scale*length_scale) - Crr*Frrz - Frrx*rrr;

ODE = [1/s_dot; n_N_dot/s_dot; psi_N_dot/s_dot; psi_ddot_N/s_dot; u_N_dot/s_dot; v_N_dot/s_dot; O_fl_N_dot/s_dot; O_fr_N_dot/s_dot; O_rl_N_dot/s_dot; O_rr_N_dot/s_dot];
f_dynamics = Function('f_dynamics',{X_sym,U_sym,s_sym},{ODE});

%% Constraint Definition
% Dynamics Constraints
g_dynamics = {};
lbg_dynamics = [];
ubg_dynamics = [];

% Time Constraints
g_time = {};
lbg_time = [];
ubg_time = [];

% Normal Force Constraints
g_force = {};
lbg_force = {};
ubg_force = {};

% Engine Power Definition
g_power = {};
lbg_power = [];
ubg_power = [];

% Main Loop
for k = 1:N
    k1 = f_dynamics(X(:,k),U(:,k), s_grid(k));
    k2 = f_dynamics(X(:,k) + ds/2*k1, U(:,k), s_grid(k) + ds/2);
    k3 = f_dynamics(X(:,k) + ds/2*k2, U(:,k), s_grid(k) + ds/2);
    k4 = f_dynamics(X(:,k) + ds*k3,   U(:,k), s_grid(k+1));
    x_next = X(:,k) + (ds/6)*(k1 + 2*k2 + 2*k3 + k4);
   
    % Time Constraint Definition (Time is forced to move ahead)
    g_time{end+1} = x_next(1) - X(1,k);
    lbg_time = [lbg_time; 1e-6];
    ubg_time = [ubg_time; inf];

    % Dynamics Constraint Definition
    g_dynamics{end+1} = X(:,k+1) - x_next;
    lbg_dynamics = [lbg_dynamics; zeros(n_states,1)];
    ubg_dynamics = [ubg_dynamics; zeros(n_states,1)];

    % Normal Force Constraint Definition
    % Front Left Wheel
    g_force{end+1} = Fzfl - (Car.m*Car.g*Car.d/(2*Car.l) - Car.m*(k1(5)/k1(1))*Car.h/(2*Car.l) - Car.Droll*Car.m*(k1(6)/k1(1))*Car.h/(Car.Weq) - Downforce*(Car.d - Car.a)/(2*Car.l));
    lbg_force = [lbg_force; 0];
    ubg_force = [ubg_force; 0];

    % Front Right Wheel
    g_force{end+1} = Fzfr - (Car.m*Car.g*Car.d/(2*Car.l) - Car.m*(k1(5)/k1(1))*Car.h/(2*Car.l) + Car.Droll*Car.m*(k1(6)/k1(1))*Car.h/(Car.Weq) - Downforce*(Car.d - Car.a)/(2*Car.l));
    lbg_force = [lbg_force; 0];
    ubg_force = [ubg_force; 0];

    % Rear Left Wheel
    g_force{end+1} = Fzrl - (Car.m*(Car.g/2)*(1 - (Car.d/Car.l)) + Car.m*(k1(5)/k1(1))*Car.h/(2*Car.l) - (1 - Car.Droll)*Car.m*(k1(6)/k1(1))*Car.h/(Car.Weq) - Downforce*(Car.l - Car.d + Car.a)/(2*Car.l));
    lbg_force = [lbg_force; 0];
    ubg_force = [ubg_force; 0];

    % Rear Right Wheel
    g_force{end+1} = Fzrr - (Car.m*(Car.g/2)*(1 - (Car.d/Car.l)) + Car.m*(k1(5)/k1(1))*Car.h/(2*Car.l) + (1 - Car.Droll)*Car.m*(k1(6)/k1(1))*Car.h/(Car.Weq) - Downforce*(Car.l - Car.d + Car.a)/(2*Car.l));
    lbg_force = [lbg_force; 0];
    ubg_force = [ubg_force; 0];

    % Engine Power Constraint Definition
    g_power{end+1} = Car.PeakPower - time_scale/((force_scale*length_scale)*angle_scale)*(U(2,k)*X(7,k) + U(3,k)*X(8,k) + U(4,k)*X(9,k) + U(5,k)*X(10,k));
    lbg_power = [lbg_power; -inf];
    ubg_power = [ubg_power; 0];
end

%% Net Constraint Definition
g = [vertcat(g_dynamics{:}); vertcat(g_time{:}); vertcat(g_force{:}); vertcat(g_power{:})];
lbg = [lbg_dynamics; lbg_time; lbg_force; lbg_power];
ubg = [ubg_dynamics; ubg_time; ubg_force; ubg_power];

% Closed Loop Constraint
g   = [g;X([2:n_states],end) - X([2:n_states],1)]; 
lbg = [lbg; zeros(n_states-1,1)]; 
ubg = [ubg; zeros(n_states-1,1)];

end