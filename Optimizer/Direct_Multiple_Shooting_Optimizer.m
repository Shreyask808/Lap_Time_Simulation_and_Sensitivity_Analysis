%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This project utilizes CasADi for symbolic differentiation and optimization.
% Andersson, J. A. E., et al. (2019). "CasADi: a software framework for 
% nonlinear optimization and optimal control." Mathematical Programming Computation.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc
clear
close all

%% Load Vehicle Model
[file1,loc1] = uigetfile({'*.mat'}, 'Select the Vehicle Data');
if isequal(file1,0)
    disp('Select the Vehicle Data')
else
    location_vehicle = fullfile(loc1,file1);
    load(location_vehicle);
    fprintf('Loaded Vehicle Data: %s\n', location_vehicle);
end

%% Load Track Model
[file2,loc2] = uigetfile({'*.mat'}, 'Select the Track Data');
if isequal(file2,0)
    disp('Select the Track Data')
else
    location_track = fullfile(loc2,file2);
    load(location_track);
    fprintf('Loaded Vehicle Data: %s\n', location_track);
end

N = 500;

%% Vehicle Model Folder Selection
[file4, loc4] = uigetfile({'*.m'}, 'Select the Tire Model.m file');

if isequal(file4, 0)
    error('No tire model file selected.');
else
    fprintf('Tire Model Path Added: %s\n', loc4);
end

%% Vehicle Model Folder Selection
[file3, loc3] = uigetfile({'*.m'}, 'Select the Vehicle Model.m file');

if isequal(file3, 0)
    error('No vehicle model file selected.');
else
    addpath(loc3); 
    [n_states,u_states,f_dynamics,g,lbg,ubg,lbx,ubx,x0,states,cost] = Seven_DOF_Handling_Model_2D(Car,track_data,N,loc4)
    fprintf('Vehicle Model Path Added: %s\n', loc3);
end

%% Optimizer Definition
nlp = struct('x', states, 'f', cost, 'g',g);
opts = struct;
opts.ipopt.max_iter = 5000;
opts.ipopt.obj_scaling_factor = 1;
opts.ipopt.nlp_scaling_method = 'gradient-based';
solver = nlpsol('S', 'ipopt', nlp, opts);
sol = solver('x0', x0, 'lbx', lbx, 'ubx', ubx, 'lbg', lbg, 'ubg', ubg);
full_sol = full(sol.x);
X_opt = reshape(full_sol(1:n_states*(N+1)), n_states, N+1);
U_opt = reshape(full_sol(n_states*(N+1)+1:end), u_states, N);
fprintf('Optimal Lap Time: %.3f seconds\n', X_opt(1, end)/t_scale);
