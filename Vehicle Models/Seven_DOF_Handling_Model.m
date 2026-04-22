clc
clear
close all

function [ode] = Handling_Model(Car,track_data,)
import casadi.*

%% Vehicle Parameters
l = Car.l;                                                                  % Vehicle Wheelbase in m
%% Direct Multiple Shooting Method 
N = 2000;
n_states = 10;
u_states = 5;

X = SX.sym('X',n_states,N+1);
U = SX.sym('U',u_states,N);
X_reshape = reshape(X,n_states*(N+1),1);                                    % Column vector of all States
U_reshape = reshape(U,u_states*(N),1);                                      % Column vector of all Control Inputs
states = [X_reshape; U_reshape];                                            % Direct Multiple Shooting States

%% Dynamics Definition
X_sym = SX.sym('X_sym',n_states);
U_sym = SX.sym('U_sym',u_states);
s_sym = SX.sym('s_sym');

% Scale States
lenght_scale = l;
end