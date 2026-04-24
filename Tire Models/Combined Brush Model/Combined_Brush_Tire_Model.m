clc
clear
close all

function [Fx,Fy,Mz] = combined_brsuh_model(u,v,P_in,R,w,C_p,kv,omega,alpha,Fz,mu_0,mu)
import casadi.*
% Effective Radius Calculation
r = R - (Fz/kv);                                                            % Effective Raius of the Tire (m)
a = sqrt(R^2 - r^2);                                                        % Half the length of the contact patch (m)
P_max = 3*Fz/(4*a*w);                                                       % Maximum Contact Pressure (Pa)
    
% Slip Ratios
kappa = (r*omega - u)/(u + 1e-6);
sigma_x = kappa/(1+kappa);                                                  % Longitudinal Theoretical Slip
sigma_y = tan(alpha)/(1 + kappa);                                           % Lateral Theoretical Slip
sigma = [sigma_x; sigma_y];                                                 % Net Slip Vector
sigma_norm = sqrt(sigma_x^2 + sigma_y^2);                                   % Norm of the Slip Vector

% Adhesion Length
lambda = 2*a - (C_p*sigma_norm*a^2)/(mu_0*w*P_max);                         % Adhesion Length estimate (m)
lambda = fmax(0,fmin(2*a,lambda));

% Forces and Moments
Fx = C_p*sigma_x*(lambda^2)/2 + (sigma_x/(sigma_norm + 1e-6))*mu*w*P_max*((4*a/3) - (lambda^2/a) + (lambda^3/(3*a^2)));                         % Lateral Force (N)
Fy = C_p*sigma_y*(lambda^2)/2 + (sigma_y/(sigma_norm + 1e-6))*mu*w*P_max*((4*a/3) - (lambda^2/a) + (lambda^3/(3*a^2)));                         % Longitudinal Force (N)
Mz = C_p*sigma_y*(((a*lambda^2)/2) - (lambda^3)/3) + (sigma_y/(sigma_norm + 1e-6))*mu*w*P_max*((lambda^3/a) - (lambda^4)/(4*a^2) - lambda^2);   % Aligning Moment (N.m)
end