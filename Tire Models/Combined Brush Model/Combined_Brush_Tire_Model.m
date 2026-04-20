clc
clear
close all
import casadi.*

%% Tire Data
Fz = SX.sym('Fz');                                                          % Nomral Force at the wheel in N
R = SX.sym('R');                                                            % Radius of the Tire in m
r = SX.sym('r');                                                            % Effective Radius of the Tire in m
a = SX.sym('a');                                                            % Half the length of the contact patch in m
w = SX.sym('w');                                                            % Width of the tire in m
P_in = SX.sym('P_in');                                                      % Tire internal pressure in Pa
s = SX.sym('s');                                                            % Longitudinal Slip Ratio
alpha = SX.sym('alpha');                                                    % Slip angle of the tire in rad

% Effective Radius Estimation
a = Fz/(2*P_in*w);                                                          % Half the contact patch length in m
r = sqrt(R^2 - a^2);                                                        % Effective Radius of the tire in m

% Longitudinal Slip Ratio Definition
s = (r*omega - u*cos(alpha)/(r*omega);


% Adhesion Length of the tire
