clc
clear
close all
import casadi.*

function [Fx,Fy] = combined_brsuh_model(u,v,P_in,R,w,C_p,omega,alpha,Fz,mu_0,mu)
% Effective Radius Calculation
a = Fz/(2*P_in*w);                                                          % Half of contact patch length (m)
r = sqrt(R^2 - a^2);                                                        % Effective Radius of the Tire (m)
P_max = 3*Fz/(4*a*w);

% Slip Ratios
kappa = (r*omega - u)/u;
sigma_x = kappa/(1+kappa);
sigma_y = tan(alpha)/(1 + kappa);
sigma = [sigma_x; sigma_y];

% Adhesion Length
lambda = 2*a - (C_p*norm(sigma)*a^2)/(mu_0*w*P_max);

% Forces
Fx = C_p*sigma_x*(lambda^2)/2 + (sigma_x/norm(sigma))*mu*w*P_max*((4*a/3) - (lambda^2/a) + (lambda^3/(3*a^2)));
Fy = C_p*sigma_y*(lambda^2)/2 + (sigma_y/norm(sigma))*mu*w*P_max*((4*a/3) - (lambda^2/a) + (lambda^3/(3*a^2)));
end