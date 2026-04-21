clc
clear
close all

function [Fx,Fy] = combined_brsuh_model(u,v,P_in,R,w,C_p,omega,alpha,Fz,mu_0,mu)
import casadi.*
% Effective Radius Calculation
a = Fz/(2*P_in*w);                                                          % Half of contact patch length (m)
r = sqrt(R^2 - a^2);                                                        % Effective Radius of the Tire (m)
P_max = 3*Fz/(4*a*w);

% Slip Ratios
kappa = (r*omega - u)/(u + 1e-6);
sigma_x = kappa/(1+kappa);
sigma_y = tan(alpha)/(1 + kappa);
sigma = [sigma_x; sigma_y];
sigma_norm = sqrt(sigma_x^2 + sigma_y^2);
% Adhesion Length
lambda = 2*a - (C_p*sigma_norm*a^2)/(mu_0*w*P_max);
lambda = fmax(0,fmin(2*a,lambda));
% Forces
Fx = C_p*sigma_x*(lambda^2)/2 + (sigma_x/(sigma_norm + 1e-6))*mu*w*P_max*((4*a/3) - (lambda^2/a) + (lambda^3/(3*a^2)));
Fy = C_p*sigma_y*(lambda^2)/2 + (sigma_y/(sigma_norm + 1e-6))*mu*w*P_max*((4*a/3) - (lambda^2/a) + (lambda^3/(3*a^2)));
end