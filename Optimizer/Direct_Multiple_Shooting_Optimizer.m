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

%% Optimizer Definition