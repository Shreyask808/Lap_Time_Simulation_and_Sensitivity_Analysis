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
    location = fullfile(loc1,file1);
    load(location);
    fprintf('Loaded Vehicle Data: %s\n', location);
end

%% Load Track Model
[file1,loc1] = uigetfile({'*.mat'}, 'Select the Vehicle Data');
if isequal(file1,0)
    disp('Select the Vehicle Data')
else
    location = fullfile(loc1,file1);
    load(location);
    fprintf('Loaded Vehicle Data: %s\n', location);
end