clc
clear
close all

%% User Controlled Data
% Track Direction - Clockwise / Anti-Clockwise
prompt = "Is the track clockwise? (Y/N): ";
direction = input(prompt, 's');
if strcmpi(direction, 'Y')
    direction_angle = 2*pi;
elseif strcmpi(direction, 'N')
    direction_angle = -2*pi;
else
    error('Invalid input. Please enter Y or N.');
end
N = 1200;                                                                   % Number of individual Segments

% Track Pitch Limits [rad]
mu_min = -pi/2 + 0.1;                                                       % Minimum Track Pitch in rad       
mu_max = pi/2 - 0.1;                                                        % Maximum Track Pitch in rad

% Left Half Track Width [m]
nl_min = 3;                                                                 % Minimum Left Half Track Width in m
nl_max = 10;                                                                % Maximum Left Half Track Width in m

% Right Half Track Width [m]
nr_min = -10;                                                               % Minimum Right Half Track Width in m
nr_max = -3;                                                                % Maximum Right Half Track Width in m

% Track Yae about z-axis [rad]
theta_min = -inf;                                                           % Minimum Track Yaw in rad
theta_max = inf;                                                            % Maximum Track Yaw in rad

% Track Relative Twist about x-axis [rad]
phi_min = -inf;                                                             % Minimum Track Roll in rad
phi_max = inf;                                                              % Maximum Track Roll in rad

% Weights
w_c = 1e-3;
w_l = w_c;
w_r = w_c;

% Track Euler Angle Weights
w_theta = 3e3;
w_mu = 1e10;
w_phi = w_mu;

% Track Width Weights
w_nl = 1e-2;
w_nr = w_nl;

% Initial Track Width Guess
nl_0 = 5.3;
nr_0 = -5.3;

%% Import Track GPS Data
% Right Boundary Data
[file1,loc1] = uigetfile({'*.txt';'*.csv'}, 'Select the Right Boundary Data');
if isequal(file1,0)
    disp('Select Right Boundary Data')
else
    location_right = fullfile(loc1,file1);
    right_boundary = readtable(location_right);
    fprintf('Loaded: %s\n', location_right);
end

% Left Boundary Data
[file2,loc2] = uigetfile({'*.txt';'*.csv'}, 'Select the Left Boundary Data');
if isequal(file2,0)
    disp('Select Left Boundary Data')
else
    location_left = fullfile(loc2,file2);
    left_boundary = readtable(location_left);
    fprintf('Loaded: %s\n', location_left);
end

% Right Coordinates
right.lattitude = table2array(right_boundary(:,2));
right.longitude = table2array(right_boundary(:,3));
right.altitude = table2array(right_boundary(:,4));

% Apply a smoothing filter to the raw data rows
right.lattitude = smoothdata(right.lattitude, 'gaussian', 15);
right.longitude = smoothdata(right.longitude, 'gaussian', 15);
right.altitude = smoothdata(right.altitude, 'gaussian', 25);

% Left Coordinates
left.lattitude = table2array(left_boundary(:,2));
left.longitude = table2array(left_boundary(:,3));
left.altitude = table2array(left_boundary(:,4));

left.lattitude = smoothdata(left.lattitude, 'gaussian', 15);
left.longitude = smoothdata(left.longitude, 'gaussian', 15);
left.altitude = smoothdata(left.altitude, 'gaussian', 40);

% Centerline Origin
origin.lat = (right.lattitude(1) + left.lattitude(1))/2;
origin.lon = (right.longitude(1) + left.longitude(1))/2;
origin.altitude = (right.altitude(1) + left.altitude(1))/2;
wgs84 = wgs84Ellipsoid;

% Convert Right Boundary
[xR, yR, zR] = geodetic2enu(right.lattitude, right.longitude, right.altitude, origin.lat, origin.lon, origin.altitude, wgs84);
bR.coordinates = [xR, yR, zR; xR(1), yR(1), zR(1)]; 
new_bL_coords = zeros(size(bR.coordinates));

% Convert Left Boundary
[xL, yL, zL] = geodetic2enu(left.lattitude, left.longitude, left.altitude, origin.lat, origin.lon, origin.altitude, wgs84);
bL.coordinates = [xL, yL, zL; xL(1), yL(1), zL(1)];
dist = [];

%Centerline Coordinates
for i = 1:length(bR.coordinates)
    for j = 1:length(bL.coordinates)
        distance = bR.coordinates(i,:) - bL.coordinates(j,:);
        dist(j) = norm(distance,2);
    end
    [~,idx_min] = min(dist);
    new_bL_coords(i,:) = bL.coordinates(idx_min,:);
    cline.coordinates(i,:) = (bR.coordinates(i,:) + bL.coordinates(idx_min,:))/2;
end
bL.coordinates = new_bL_coords;

%% Centerline Arc Length
arc_s = zeros(length(cline.coordinates),1);
for a = 2:length(arc_s)
    arc = cline.coordinates(a,:) - cline.coordinates(a-1,:);
    arc_s(a) = arc_s(a-1) + norm(arc,2);
end

%% Track Boundary Optimal Control
import casadi.*
s_grid = linspace(arc_s(1),arc_s(end),N+1);                                 % Coarse Grid for Direct Multiple Shooting 
ds = (arc_s(end) - arc_s(1))/(N);                                           % Step Size for Coarse Grid

%% State and Control Vector
X = SX.sym('X',11,N+1);                                                     % States in the Problem
U = SX.sym('U',5,N);                                                        % Control Inputs in the Problem

X_reshape = reshape(X,11*(N+1),1);                                          % Column vector of all States
U_reshape = reshape(U,5*(N),1);                                             % Column vector of all Control Inputs
states = [X_reshape; U_reshape];                                            % Direct Multiple Shooting States
obj = 0;                                                                    % Cost Function Definition
g = {};                                                                     % Constraints Definition
lbg = [];                                                                   % Lower Bound on Constraints
ubg = [];                                                                   % Upper Bound on Constraints
arc_s = arc_s(:).';
cline_vals = cline.coordinates.';
bL_vals    = bL.coordinates.';
bR_vals    = bR.coordinates.';
opts = struct;
cline_fun = casadi.interpolant('cline_fun','linear',{arc_s},cline_vals(:),opts);
bL_fun    = casadi.interpolant('bL_fun',   'linear',{arc_s},bL_vals(:),opts);
bR_fun    = casadi.interpolant('bR_fun',   'linear',{arc_s},bR_vals(:),opts);

%% Track Dynamics Definition 
s_st = SX.sym('s_st',11);                                                   % State Vector
s_ctrl = SX.sym('s_ctrl',5);                                                % Control Input Vector
theta_f = s_st(4);
mu_f = s_st(5);
phi_f = s_st(6);
theta_dot_f = s_st(7);
mu_dot_f = s_st(8);
phi_dot_f = s_st(9);
theta_ddot_f = s_ctrl(1);
mu_ddot_f = s_ctrl(2);
phi_ddot_f = s_ctrl(3);
nl_dot_f = s_ctrl(4);
nr_dot_f = s_ctrl(5);

ode = [cos(theta_f)*cos(mu_f); sin(theta_f)*cos(mu_f); -sin(mu_f); theta_dot_f; mu_dot_f; phi_dot_f; theta_ddot_f; mu_ddot_f; phi_ddot_f; nl_dot_f; nr_dot_f];
f_dynamics = Function('f_dynamics', {s_st,s_ctrl},{ode});                   % Track Dynamics Definition

%% Problem Definition
for k = 1:N
    s_current = s_grid(k);
    state = X(:,k);
    ctrl = U(:,k);
    % States
    cline_coord = state(1:3);
    nl = state(10);
    nr = state(11);
    theta = state(4);
    mu = state(5);
    phi = state(6);
    theta_dot = state(7);
    mu_dot = state(8);
    phi_dot = state(9); 

    % Controls
    theta_ddot = ctrl(1);
    mu_ddot = ctrl(2);
    phi_ddot = ctrl(3);
    nl_dot = ctrl(4);
    nr_dot = ctrl(5);
    
    normal = [cos(theta)*sin(mu)*sin(phi) - sin(theta)*cos(phi); sin(theta)*sin(mu)*sin(phi) + cos(theta)*cos(phi); cos(mu)*sin(phi)];
    bL_coord = cline_coord + nl*normal;
    bR_coord = cline_coord + nr*normal;

    % Objective Function Construction
    e_cline = w_c*sumsqr(cline_coord - cline_fun(s_current));
    e_bl = w_l*sumsqr(bL_coord - bL_fun(s_current));
    e_br = w_r*sumsqr(bR_coord - bR_fun(s_current));
    r_c = w_theta*(theta_ddot^2) + w_mu*(mu_ddot^2) + w_phi*(phi_ddot^2);
    r_w = w_nl*(nl_dot^2) + w_nr*(nr_dot^2);

    obj = obj + e_cline + e_bl + e_br + r_c + r_w;

    % Dynamics

    % 4th Order Runga Kutta Method for Discretization
    k1 = f_dynamics(state,ctrl);
    k2 = f_dynamics(state + ds/2*k1, ctrl);
    k3 = f_dynamics(state + ds/2*k2, ctrl);
    k4 = f_dynamics(state + ds*k3,   ctrl);
    state_next = state + (ds/6)*(k1 + 2*k2 + 2*k3 + k4);
    g{end+1} = X(:,k+1) - state_next;
end

%% Bounds on States
lbx = - inf(size(states));
ubx = inf(size(states));

for r = 1:N+1
    idx_mu = (r-1)*11 + 5;
    idx_theta = (r-1)*11 + 4;
    idx_phi = (r-1)*11 + 6;
    idx_nl = (r-1)*11 + 10;
    idx_nr = (r-1)*11 + 11;

    lbx(idx_mu) = mu_min;
    ubx(idx_mu) = mu_max;

    lbx(idx_theta) = theta_min;
    ubx(idx_theta) = theta_max;

    lbx(idx_phi) = phi_min;
    ubx(idx_phi) = phi_max;    

    lbx(idx_nl) = nl_min;
    ubx(idx_nl) = nl_max;

    lbx(idx_nr) = nr_min;
    ubx(idx_nr) = nr_max;
end

%% End Constraints
for w = 1:11
    if w == 4
        g{end+1} = X(w,N+1) - X(w,1) + direction_angle;                    
    else
        g{end+1} = X(w,N+1) - X(w,1);                                       % Other States are the same as starting point
    end
end
lbg = zeros(size(vertcat(g{:})));
ubg = zeros(size(vertcat(g{:})));

%% Initial Guess
x0 = zeros(size(states));
x_guess = interp1(arc_s, cline.coordinates, s_grid);
dx = diff(x_guess(:,1));
dy = diff(x_guess(:,2));
theta_guess = unwrap(atan2(dy, dx));
theta_guess(end+1) = theta_guess(end); 

for k = 1:N+1
   idx = (k-1)*11;
   x0(idx+1:idx+3) = x_guess(k, :); 
   x0(idx+4) = theta_guess(k);      
   x0(idx+10) = nl_0; 
   x0(idx+11) = nr_0;
end

%% NLP Sover
nlp = struct('x', states, 'f', obj, 'g', vertcat(g{:}));
opts = struct;
opts.ipopt.max_iter = 2000; % Set a limit for testing
solver = nlpsol('S', 'ipopt', nlp, opts);
sol = solver('x0', x0, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);
full_sol = full(sol.x);
X_opt = reshape(full_sol(1 : 11*(N+1)), 11, N+1);
U_opt = reshape(full_sol(11*(N+1) + 1 : end), 5, N);

%% Store Track Data
% Track Coordinates
x_opt = X_opt(1,:);
y_opt = X_opt(2,:);
z_opt = X_opt(3,:);
cline = [x_opt; y_opt; z_opt];

% Euler Angles
theta_opt = X_opt(4,:);
mu_opt = X_opt(5,:);
phi_opt = X_opt(6,:);

% Euler Angle Rates
theta_dot_opt = X_opt(7,:);
mu_dot_opt = X_opt(8,:);
phi_dot_opt = X_opt(9,:);

% Track Widths
nl_opt = X_opt(10,:);
nr_opt = X_opt(11,:);

% Right and Left Boundaries
bl_opt = zeros(size(cline));
br_opt = zeros(size(cline));

% Direction Vectors
tangent = zeros(size(cline));
normal = zeros(size(cline));
m = zeros(size(cline));

% Angular Velocities
Omega_x = zeros(size(x_opt));
Omega_y = zeros(size(x_opt));
Omega_z = zeros(size(x_opt));

for q = 1:length(x_opt)
    tangent(:,q) = [cos(theta_opt(q))*cos(mu_opt(q)); sin(theta_opt(q))*cos(mu_opt(q)); - sin(mu_opt(q))];
    normal(:,q) = [cos(theta_opt(q))*sin(mu_opt(q))*sin(phi_opt(q)) - sin(theta_opt(q))*cos(phi_opt(q)); sin(theta_opt(q))*sin(mu_opt(q))*sin(phi_opt(q)) + cos(theta_opt(q))*cos(phi_opt(q)); cos(mu_opt(q))*sin(phi_opt(q))];
    m(:,q) = cross(tangent(:,q), normal(:,q));
    bl_opt(:,q) = cline(:,q) + nl_opt(q)*normal(:,q);
    br_opt(:,q) = cline(:,q) + nr_opt(q)*normal(:,q);
    J = [1, 0, -sin(mu_opt(q)); 0, cos(phi_opt(q)), cos(mu_opt(q))*sin(phi_opt(q)); 0, -sin(phi_opt(q)), cos(mu_opt(q))*cos(phi_opt(q))];
    omega_b = J*[phi_dot_opt(q); mu_dot_opt(q); theta_dot_opt(q)];
    Omega_x(q) = omega_b(1,:);
    Omega_y(q) = omega_b(2,:);
    Omega_z(q) = omega_b(3,:);
end

%% Export Track Data
default_name = 'Track_Model.mat';
[file, path] = uiputfile('*.mat', 'Select Folder and Name to Save Track', default_name);

% Check if the user pressed "Cancel"
if isequal(file, 0) || isequal(path, 0)
   disp('Save cancelled by user.');
else
    % Centerline Coordinates
    track_data.arc_s = s_grid;
    track_data.xc = x_opt;
    track_data.yc = y_opt;
    track_data.zc = z_opt;
    
    % Line Data
    track_data.cline = cline;
    track_data.bl = bl_opt;
    track_data.br = br_opt;
    
    % Euler Angles
    track_data.theta = theta_opt;
    track_data.mu = mu_opt;
    track_data.phi = phi_opt;
    
    % Euler Angle Rates
    track_data.theta_dot = theta_dot_opt;
    track_data.mu_dot = mu_dot_opt;
    track_data.phi_dot = phi_dot_opt;
    
    % Left & Right Boundaries
    track_data.nl = nl_opt;
    track_data.nr = nr_opt;
    
    % Directions
    track_data.tangent = tangent;
    track_data.normal = normal;
    track_data.m = m;
    
    % Angular Velocities
    track_data.Omega_x = Omega_x;
    track_data.Omega_y = Omega_y;
    track_data.Omega_z = Omega_z;

    full_save_path = fullfile(path, file); 
    save(full_save_path, 'track_data');    
    fprintf('Successfully saved track model to: %s\n', full_save_path);
end

%% Results and Plots
% Track Map
figure(1)
plot3(track_data.bl(1,:),track_data.bl(2,:),track_data.bl(3,:),'r','LineWidth',2)
hold on
plot3(bL.coordinates(:,1),bL.coordinates(:,2),bL.coordinates(:,3),'ro')
grid minor
plot3(track_data.br(1,:),track_data.br(2,:),track_data.br(3,:),'b','LineWidth',2)
plot3(bR.coordinates(:,1),bR.coordinates(:,2),bR.coordinates(:,3),'bo')
plot3(track_data.cline(1,:),track_data.cline(2,:),track_data.cline(3,:),'k','LineWidth',2)
plot3(cline(1,:),cline(2,:),cline(3,:),'ko')
grid minor
xlabel('$x - coordinate (m) $', 'Interpreter', 'latex', 'FontSize', 16)
ylabel('$y - coordinate (m) $', 'Interpreter', 'latex', 'FontSize', 16)
zlabel('$z - coordinate (m) $', 'Interpreter', 'latex', 'FontSize', 16)
legend('Optimized Left Boundary','Left Boundary Points','Optimized Right Boundary','Right Boundary Points','Optimized Centerline','Centerline Points')
title('Circuit 3D Boundaries')
axis equal
image_save_path1 = fullfile(path,'Track_Map.png');
saveas(figure(1),image_save_path1)

% Euler Angles
figure(2)
subplot(3,1,1)
plot(track_data.arc_s,track_data.theta,'k','LineWidth',2)
hold on
yline(0,'--r','LineWidth',2)
xlabel('$ Centerline\ Arc\ Length\ (m) $', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$ \theta $', 'Interpreter', 'latex', 'FontSize', 14)
grid minor

subplot(3,1,2)
plot(track_data.arc_s,track_data.mu,'k','LineWidth',2)
hold on
yline(0,'--r','LineWidth',2)
xlabel('$ Centerline\ Arc\ Length\ (m) $', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$ \mu $', 'Interpreter', 'latex', 'FontSize', 14)
grid minor

subplot(3,1,3)
plot(track_data.arc_s,track_data.phi,'k','LineWidth',2)
hold on
yline(0,'--r','LineWidth',2)
xlabel('$ Centerline\ Arc\ Length\ (m) $', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$ \phi $', 'Interpreter', 'latex', 'FontSize', 14)
grid minor
image_save_path2 = fullfile(path,'Track_Euler_Angles.png');
saveas(figure(2),image_save_path2)

% Track Half Widths
figure(3)
subplot(2,1,1)
plot(track_data.arc_s,track_data.nl,'k','LineWidth',2)
hold on
yline(nl_max,'--r','LineWidth',2)
xlabel('$ Centerline\ Arc\ Length\ (m) $', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$ n_l $', 'Interpreter', 'latex', 'FontSize', 14)
grid minor

subplot(2,1,2)
plot(track_data.arc_s,track_data.nr,'k','LineWidth',2)
hold on
yline(nr_min,'--r','LineWidth',2)
xlabel('$ Centerline\ Arc\ Length\ (m) $', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$ n_r $', 'Interpreter', 'latex', 'FontSize', 14)
grid minor
image_save_path3 = fullfile(path,'Track_Widths.png');
saveas(figure(3),image_save_path3)

% Track Angular Velocities
figure(4)
subplot(3,1,1)
plot(track_data.arc_s,track_data.Omega_x,'k','LineWidth',2)
hold on
yline(0,'--r','LineWidth',2)
xlabel('$ Centerline\ Arc\ Length\ (m) $', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$ \Omega_x\ (rad/m) $', 'Interpreter', 'latex', 'FontSize', 14)
grid minor

subplot(3,1,2)
plot(track_data.arc_s,track_data.Omega_y,'k','LineWidth',2)
hold on
yline(0,'--r','LineWidth',2)
xlabel('$ Centerline\ Arc\ Length\ (m) $', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$ \Omega_y\ (rad/m) $', 'Interpreter', 'latex', 'FontSize', 14)
grid minor

subplot(3,1,3)
plot(track_data.arc_s,track_data.Omega_z,'k','LineWidth',2)
hold on
yline(0,'--r','LineWidth',2)
xlabel('$ Centerline\ Arc\ Length\ (m) $', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$ \Omega_z\ (rad/m) $', 'Interpreter', 'latex', 'FontSize', 14)
grid minor
image_save_path4 = fullfile(path,'Track_Angular_Velocities.png');
saveas(figure(4),image_save_path4)