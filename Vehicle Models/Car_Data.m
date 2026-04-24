clc
clear
close all

%% Vehicle and Modeling Data
%% Vehicle Dimension Data
Car.l = 3.6;                                                                % Wheelbase of the car (m)
Car.d = 1.6;                                                                % Distance between rear axle and cg (m)
Car.h = 0.3;                                                                % Height of cg above the ground level (m)
Car.a = 0.5;                                                                % Distance between cg and cp along the centerline (m)
Car.Wf = 2;                                                                 % Front Trackwidth (m)
Car.Wr = 2;                                                                 % Rear TrackWidth (m)
Car.Droll = 0.5;                                                            % Rolling Moment Distribution (.)
Car.delta_max = 0.5;                                                        % Max steering angle (rad)
Car.delta_min = -0.5;                                                       % Min steering angle (rad)
Car.A = 1.5;                                                                % Frontal Area (m^2)
Car.Cd = 0.22;                                                              % Vehicle Drag Coefficient (.)
Car.Cl = -0.5;                                                              % Vehicle Drag Coefficient (.)

%% Vehicle Mass and Inertia Properties
Car.m = 660;                                                                % Vehicle Mass (kg)
Car.Iz = 450;                                                               % Vehicle Moment of Interia about cg (kg.m^2)

Car.Ifl = 50;                                                               % Front Left Wheel Moment of Intertia (kg.m^2)
Car.Ifr = 50;                                                               % Front Right Wheel Moment of Intertia (kg.m^2)
Car.Irl = 50;                                                               % Rear Left Wheel Moment of Intertia (kg.m^2)
Car.Irr = 50;                                                               % Rear Right Wheel Moment of Intertia (kg.m^2)

%% Tire Data
Car.R = 0.2;                                                                % Tire Radius (m)
Car.w = 0.225;                                                              % Tire Width (m)
Car.Cp = 4e6;                                                               % Tire Stiffness per unit length (N/m^2)
Car.mu0 = 0.9;                                                              % Static Friction Coefficient (.)
Car.mu = 0.6;                                                               % Sliding Friction Coefficient (.)
Car.kv = 1e6;                                                               % Vertical Stiffness of the Tire (N/m^2)

%% Vehicle Powertrain Limit
Car.PeakPower = 450e3;                                                      % Vehicle Peak Power (W)
Car.PeakTorque = 1000;                                                      % Vehicle Peak Torque (N.m)

%% Miscellaneous Data
Car.rho = 1.2;                                                              % Density of Air (kg/m^3)
Car.g = 9.81;                                                               % Gravitational Constant (m/s^2)

save('VehicleData.mat', 'Car');
disp('Vehicle Data saved successfully');