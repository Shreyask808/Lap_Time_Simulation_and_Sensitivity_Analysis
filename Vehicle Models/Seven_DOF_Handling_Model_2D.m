%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This project utilizes CasADi for symbolic differentiation and optimization.
% Andersson, J. A. E., et al. (2019). "CasADi: a software framework for 
% nonlinear optimization and optimal control." Mathematical Programming Computation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 7 DOF Plannar Vehicle Handling Model 
function [n_states,u_states,g,lbg,ubg,lbx,ubx,x0,states,cost,time_scale] = Seven_DOF_Handling_Model_2D(Car,track_data,N,loc4)
import casadi.*

%% Direct Multiple Shooting Method Inputs
n_states = 10;                                                              % Number of States
u_states = 9;                                                               % Number of Control Inputs

u_max = 100;
v_max = 10;
v_min = -10;

% Cost Function Weights;
e1 = 1e-6;
e2 = 1e-6;
e3 = 1e-6;
e4 = 1e-6;
e5 = 1e-6;

% Independent Variable
lap_length = track_data.arc_s(end);                                         % Total length of the lap (m)
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
length_states = n_states*(N+1) + u_states*N;

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
Fzfl_N = U_sym(6);
Fzfr_N = U_sym(7);
Fzrl_N = U_sym(8);
Fzrr_N = U_sym(9);

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
Fzfl = Fzfl_N/force_scale;
Fzfr = Fzfr_N/force_scale;
Fzrl = Fzrl_N/force_scale;
Fzrr = Fzrr_N/force_scale;

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
f_Drag = Function('f_Drag',{X_sym},{Drag});
Downforce = -(1/2)*Car.Cl*Car.rho*Car.A*(u^2);
f_Downforce = Function('f_Downforce',{X_sym},{Downforce});

%% Tire Forces
addpath(loc4);
% Rear Right Tire Forces
[Fxrr,Fyrr,Mzrr,rrr] = Combined_Brush_Tire_Model(u_rr,v_rr,Car.R,Car.w,Car.Cp,Car.kv,O_rr,alpha_rr,Fzrr,Car.mu0,Car.mu);

% Rear Left Tire Forces
[Fxrl,Fyrl,Mzrl,rrl] = Combined_Brush_Tire_Model(u_rl,v_rl,Car.R,Car.w,Car.Cp,Car.kv,O_rl,alpha_rl,Fzrl,Car.mu0,Car.mu);

% Front Right Tire Forces
[Fxfr,Fyfr,Mzfr,rfr] = Combined_Brush_Tire_Model(u_fr,v_fr,Car.R,Car.w,Car.Cp,Car.kv,O_fr,alpha_fr,Fzfr,Car.mu0,Car.mu);

% Front Right Tire Forces
[Fxfl,Fyfl,Mzfl,rfl] = Combined_Brush_Tire_Model(u_fl,v_fl,Car.R,Car.w,Car.Cp,Car.kv,O_fl,alpha_fl,Fzfl,Car.mu0,Car.mu);

%% ODE Definitio (Dynamics) 
s_dot = (1/time_scale)*(u/(1 - n*curv))*(cos(xi) - (Car.d/Car.l)*sin(xi)*tan(delta));
n_N_dot = (length_scale/time_scale)*(u)*(sin(xi) + (Car.d/Car.l)*cos(xi)*tan(delta));
psi_N_dot = psi_dot_N;
psi_ddot_N = (1/2)*(-2*Car.d*(Fyrl + Fyrr) + 2*(Md_fl + Md_fr + Md_rl + Md_rr) + (Fxrr - Fxrl)*Car.Wr + (2*(Car.d - Car.l)*(Fyfl + Fyfr) + (Fxfr - Fxfl)*Car.Wf)*cos(delta) + (2*(Car.d - Car.l)*(Fxfl + Fxfr) + (Fyfl - Fyfr)*Car.Wf)*sin(delta))/Car.Iz;
u_N_dot = (Drag + (Fxrl + Fxrr + (Fxfl + Fxfr)*cos(delta)) - (Fyfr*sin(delta) + Fyfl*sin(delta)))/Car.m + v*psi_dot;
v_N_dot = (Fyrl + Fyrr + (Fyfl + Fyfr)*cos(delta) + (Fxfl + Fxfr)*sin(delta))/Car.m - u*psi_dot;
O_fl_N_dot = Md_fl_N/(force_scale*length_scale) - Car.Crr*Fzfl - Fxfl*rfl;
O_fr_N_dot = Md_fr_N/(force_scale*length_scale) - Car.Crr*Fzfr - Fxfr*rfr;
O_rl_N_dot = Md_rl_N/(force_scale*length_scale) - Car.Crr*Fzrl - Fxrl*rrl;
O_rr_N_dot = Md_rr_N/(force_scale*length_scale) - Car.Crr*Fzrr - Fxrr*rrr;

ODE = [1/s_dot; n_N_dot/s_dot; psi_N_dot/s_dot; psi_ddot_N/s_dot; u_N_dot/s_dot; v_N_dot/s_dot; O_fl_N_dot/s_dot; O_fr_N_dot/s_dot; O_rl_N_dot/s_dot; O_rr_N_dot/s_dot];
f_dynamics = Function('f_dynamics',{X_sym,U_sym,s_sym},{ODE});

%% Constraint Definition and Cost Function
% Cost Function
cost = 0;

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
    g_force{end+1} = U(6,k)/force_scale - (Car.m*Car.g*Car.d/(2*Car.l) - Car.m*(k1(5)/k1(1))*Car.h/(2*Car.l) - Car.Droll*Car.m*(k1(6)/k1(1))*Car.h/(Car.Weq) - (f_Downforce(X(5,k)))*(Car.d - Car.a)/(2*Car.l));
    lbg_force = [lbg_force; 0];
    ubg_force = [ubg_force; 0];

    % Front Right Wheel
    g_force{end+1} = U(7,k)/force_scale - (Car.m*Car.g*Car.d/(2*Car.l) - Car.m*(k1(5)/k1(1))*Car.h/(2*Car.l) + Car.Droll*Car.m*(k1(6)/k1(1))*Car.h/(Car.Weq) - (f_Downforce(X(5,k)))*(Car.d - Car.a)/(2*Car.l));
    lbg_force = [lbg_force; 0];
    ubg_force = [ubg_force; 0];

    % Rear Left Wheel
    g_force{end+1} = U(8,k)/force_scale - (Car.m*(Car.g/2)*(1 - (Car.d/Car.l)) + Car.m*(k1(5)/k1(1))*Car.h/(2*Car.l) - (1 - Car.Droll)*Car.m*(k1(6)/k1(1))*Car.h/(Car.Weq) - (f_Downforce(X(5,k)))*(Car.l - Car.d + Car.a)/(2*Car.l));
    lbg_force = [lbg_force; 0];
    ubg_force = [ubg_force; 0];

    % Rear Right Wheel
    g_force{end+1} = U(9,k)/force_scale - (Car.m*(Car.g/2)*(1 - (Car.d/Car.l)) + Car.m*(k1(5)/k1(1))*Car.h/(2*Car.l) + (1 - Car.Droll)*Car.m*(k1(6)/k1(1))*Car.h/(Car.Weq) - (f_Downforce(X(5,k)))*(Car.l - Car.d + Car.a)/(2*Car.l));
    lbg_force = [lbg_force; 0];
    ubg_force = [ubg_force; 0];

    % Engine Power Constraint Definition
    g_power{end+1} = Car.PeakPower - time_scale/((force_scale*length_scale)*angle_scale)*(U(2,k)*X(7,k) + U(3,k)*X(8,k) + U(4,k)*X(9,k) + U(5,k)*X(10,k));
    lbg_power = [lbg_power; -inf];
    ubg_power = [ubg_power; 0];
    
    % Cost Function
    dt = x_next(1) - X(1,k);
    cost = cost + dt*(e1*(U(1,k)^2) + e2*(U(2,k)^2) + e3*(U(3,k)^2) + e4*(U(4,k)^2) + e5*(U(5,k)^2));
end
cost = cost + X(1,end);

%% Net Constraint Definition
g   = vertcat(g_dynamics{:}, g_time{:}, g_force{:}, g_power{:});
lbg = vertcat(lbg_dynamics, lbg_time, lbg_force, lbg_power);
ubg = vertcat(ubg_dynamics, ubg_time, ubg_force, ubg_power);

% Closed Loop Constraint
g   = [g;X([2:n_states],end) - X([2:n_states],1)]; 
lbg = [lbg; zeros(n_states-1,1)]; 
ubg = [ubg; zeros(n_states-1,1)];

% Yaw Rate Constraint
g = [g;X(3,end) - X(3,1) - 2*pi*angle_scale];
lbg = [lbg;0];
ubg = [ubg;0];
%% Limits
lbx = -inf(length_states,1);
ubx = inf(length_states,1);

% Indices for States
idx_t = [];
idx_n = [];
idx_psi = [];
idx_psidot = [];
idx_u = [];
idx_v = [];
idx_Ofl = [];
idx_Ofr = [];
idx_Orl = [];
idx_Orr = [];

% Indices for Controls
idx_delta = [];
idx_Mdfl = [];
idx_Mdfr = [];
idx_Mdrl = [];
idx_Mdrr = [];
idx_Fzfl = [];
idx_Fzfr = [];
idx_Fzrl = [];
idx_Fzrr = [];

for j = 1:N+1
    idx = n_states*(j-1);
    idx_t = [idx_t; idx+1];
    idx_n = [idx_n; idx+2];
    idx_psi = [idx_psi; idx+3];
    idx_psidot = [idx_psidot; idx+4];
    idx_u = [idx_u; idx+5];
    idx_v = [idx_v; idx+6];
    idx_Ofl = [idx_Ofl; idx+7];
    idx_Ofr = [idx_Ofr; idx+8];
    idx_Orl = [idx_Orl; idx+9];
    idx_Orr = [idx_Orr; idx+10];

    if j < N+1
        idx_delta = [idx_delta; n_states*(N+1) + u_states*(j-1) + 1];
        idx_Mdfl = [idx_Mdfl; n_states*(N+1) + u_states*(j-1) + 2];
        idx_Mdfr = [idx_Mdfr; n_states*(N+1) + u_states*(j-1) + 3];
        idx_Mdrl = [idx_Mdrl; n_states*(N+1) + u_states*(j-1) + 4];
        idx_Mdrr = [idx_Mdrr; n_states*(N+1) + u_states*(j-1) + 5];
        idx_Fzfl = [idx_Fzfl; n_states*(N+1) + u_states*(j-1) + 6];
        idx_Fzfr = [idx_Fzfr; n_states*(N+1) + u_states*(j-1) + 7];
        idx_Fzrl = [idx_Fzrl; n_states*(N+1) + u_states*(j-1) + 8];
        idx_Fzrr = [idx_Fzrr; n_states*(N+1) + u_states*(j-1) + 9];
    end
end

% Lower and Upper Bounds
for h = 1:N+1
    track_theta = full(f_theta(s_grid(h)));
    track_nl = full(f_nl(s_grid(h)));
    track_nr = full(f_nr(s_grid(h)));

    if h == 1
        ubx(idx_t(h)) = 0;
        lbx(idx_t(h)) = 0;

        ubx(idx_n(h)) = track_nl*length_scale;
        lbx(idx_n(h)) = track_nr*length_scale;

        ubx(idx_psi(h)) = track_data.theta(1)*angle_scale;
        lbx(idx_psi(h)) = track_data.theta(1)*angle_scale;

        ubx(idx_psidot(h)) = 0;
        lbx(idx_psidot(h)) = 0;

        ubx(idx_u(h)) = u_max*speed_scale;
        lbx(idx_u(h)) = 0.75*u_max*speed_scale;

        ubx(idx_v(h)) = v_max*speed_scale;
        lbx(idx_v(h)) = v_min*speed_scale;

        ubx(idx_Ofl(h)) = (u_max/Car.R)*time_scale/angle_scale;
        lbx(idx_Ofl(h)) = (5/Car.R)*time_scale/angle_scale;

        ubx(idx_Ofr(h)) = (u_max/Car.R)*time_scale/angle_scale;
        lbx(idx_Ofr(h)) = (5/Car.R)*time_scale/angle_scale;

        ubx(idx_Orl(h)) = (u_max/Car.R)*time_scale/angle_scale;
        lbx(idx_Orl(h)) = (5/Car.R)*time_scale/angle_scale;

        ubx(idx_Orr(h)) = (u_max/Car.R)*time_scale/angle_scale;
        lbx(idx_Orr(h)) = (5/Car.R)*time_scale/angle_scale;
    else
        ubx(idx_t(h)) = 300*time_scale;
        lbx(idx_t(h)) = 0;

        ubx(idx_n(h)) = track_nl*length_scale;
        lbx(idx_n(h)) = track_nr*length_scale;

        ubx(idx_psi(h)) = (track_theta + pi/2)*angle_scale;
        lbx(idx_psi(h)) = (track_theta - pi/2)*angle_scale;

        ubx(idx_psidot(h)) = 0;
        lbx(idx_psidot(h)) = 0;

        ubx(idx_u(h)) = u_max*speed_scale;
        lbx(idx_u(h)) = 5*speed_scale;

        ubx(idx_v(h)) = v_max*speed_scale;
        lbx(idx_v(h)) = v_min*speed_scale;

        ubx(idx_Ofl(h)) = (u_max/Car.R)*time_scale/angle_scale;
        lbx(idx_Ofl(h)) = (5/Car.R)*time_scale/angle_scale;

        ubx(idx_Ofr(h)) = (u_max/Car.R)*time_scale/angle_scale;
        lbx(idx_Ofr(h)) = (5/Car.R)*time_scale/angle_scale;

        ubx(idx_Orl(h)) = (u_max/Car.R)*time_scale/angle_scale;
        lbx(idx_Orl(h)) = (5/Car.R)*time_scale/angle_scale;

        ubx(idx_Orr(h)) = (u_max/Car.R)*time_scale/angle_scale;
        lbx(idx_Orr(h)) = (5/Car.R)*time_scale/angle_scale;
    end
    if h < N+1
        ubx(idx_delta(h)) = Car.delta_max*angle_scale;
        lbx(idx_delta(h)) = Car.delta_min*angle_scale;

        ubx(idx_Mdfl(h)) = Car.PeakTorque*(length_scale*force_scale);
        lbx(idx_Mdfl(h)) = -Car.PeakTorque*(length_scale*force_scale);

        ubx(idx_Mdfr(h)) = Car.PeakTorque*(length_scale*force_scale);
        lbx(idx_Mdfr(h)) = -Car.PeakTorque*(length_scale*force_scale);

        ubx(idx_Mdrl(h)) = Car.PeakTorque*(length_scale*force_scale);
        lbx(idx_Mdrl(h)) = -Car.PeakTorque*(length_scale*force_scale);

        ubx(idx_Mdrr(h)) = Car.PeakTorque*(length_scale*force_scale);
        lbx(idx_Mdrr(h)) = -Car.PeakTorque*(length_scale*force_scale);

        ubx(idx_Fzfl(h)) = inf;
        lbx(idx_Fzfl(h)) = 0;

        ubx(idx_Fzfr(h)) = inf;
        lbx(idx_Fzfr(h)) = 0;

        ubx(idx_Fzrl(h)) = inf;
        lbx(idx_Fzrl(h)) = 0;

        ubx(idx_Fzrr(h)) = inf;
        lbx(idx_Fzrr(h)) = 0;
    end
end

%% Initial Conditions
x0 = zeros(length_states,1);

for z = 1:N+1
    kappa_c = full(f_kappa(s_grid(z)));
    if z == 1
        x0(idx_t(z)) = 0;
        x0(idx_n(z)) = 0;
        x0(idx_psi(z)) = track_data.theta(1)*angle_scale;
        x0(idx_u(z)) = u_max*speed_scale;
        x0(idx_v(z)) = v_max*speed_scale;
        x0(idx_Ofl(z)) = (u_max/Car.R)*time_scale/angle_scale;
        x0(idx_Ofr(z)) = (u_max/Car.R)*time_scale/angle_scale;
        x0(idx_Orl(z)) = (u_max/Car.R)*time_scale/angle_scale;
        x0(idx_Orr(z)) = (u_max/Car.R)*time_scale/angle_scale;

        x0(idx_delta(z)) = 0;
        x0(idx_Mdfl(z)) = Car.PeakTorque*(length_scale*force_scale);
        x0(idx_Mdfr(z)) = Car.PeakTorque*(length_scale*force_scale);
        x0(idx_Mdrl(z)) = Car.PeakTorque*(length_scale*force_scale);
        x0(idx_Mdrr(z)) = Car.PeakTorque*(length_scale*force_scale);
    else
        x0(idx_t(z)) = time_scale*(s_grid(z)/u_max);
        x0(idx_n(z)) = 0;
        x0(idx_psi(z)) = full(f_theta(s_grid(z)))*angle_scale;
        x0(idx_u(z)) = u_max*speed_scale;
        x0(idx_v(z)) = v_max*speed_scale;
        x0(idx_Ofl(z)) = (u_max/Car.R)*time_scale/angle_scale;
        x0(idx_Ofr(z)) = (u_max/Car.R)*time_scale/angle_scale;
        x0(idx_Orl(z)) = (u_max/Car.R)*time_scale/angle_scale;
        x0(idx_Orr(z)) = (u_max/Car.R)*time_scale/angle_scale;

        if z < N+1
            x0(idx_delta(z)) = 0;
            x0(idx_Mdfl(z)) = Car.PeakTorque*(length_scale*force_scale);
            x0(idx_Mdfr(z)) = Car.PeakTorque*(length_scale*force_scale);
            x0(idx_Mdrl(z)) = Car.PeakTorque*(length_scale*force_scale);
            x0(idx_Mdrr(z)) = Car.PeakTorque*(length_scale*force_scale);
        
            x0(idx_Fzfl(z)) = Car.m*Car.g*Car.d/(2*Car.l);
            x0(idx_Fzfr(z)) = Car.m*Car.g*Car.d/(2*Car.l);

            x0(idx_Fzrl(z)) = (Car.m*Car.g/2)*(1 - Car.d/(Car.l));
            x0(idx_Fzrr(z)) = (Car.m*Car.g/2)*(1 - Car.d/(Car.l));
        end
    end
end
end
