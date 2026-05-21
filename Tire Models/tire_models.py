#=======================================================================================================================================================================================================================
# This project utilizes CasADi for symbolic differentiation and optimization.
# Andersson, J. A. E., et al. (2019). "CasADi: a software framework for 
# nonlinear optimization and optimal control." Mathematical Programming Computation.
#=======================================================================================================================================================================================================================

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
import json
import math
import casadi as ca
from casadi import *

os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')

#=======================================================================================================================================================================================================================
# Combined Brush Tire Model
#=======================================================================================================================================================================================================================
def combined_brush_tire_model(u,omega,alpha,Fz,vehicle):
    mu = vehicle.mu
    mu0 = vehicle.mu0
    R = vehicle.R
    w = vehicle.w
    Cp = vehicle.Cp
    kv = vehicle.kv
    alpha_max = vehicle.alpha_max
    
    # Effective Radius Calculation
    r = R - (Fz/kv)
    a = sqrt(R**2 - r**2 + 1e-6)
    Pmax = 3*Fz/(4*a*w + 1e-6)

    # Slip Ratio Definition
    eps = 1e-3
    alpha = 0.5*(alpha + ca.sqrt(alpha**2 + eps)) - 0.5*((alpha - alpha_max) + ca.sqrt((alpha - alpha_max)**2 + eps))
    kappa = (r*omega - (u))/(u)
    sigma_x = kappa/(1 + kappa)
    sigma_y = tan(alpha)/(1 + kappa)
    sigma = ca.sqrt(sigma_x**2 + sigma_y**2 + 1e-6)

    # Adhesion Length Estimation
    ad_len = 2*a - ((Cp*sigma*a**2)/(mu0*w*Pmax))
    #ad_len = fmax(ad_len, 0)
    #ad_len = fmin(ad_len,2*a)

    ad_len = 0.5*(ad_len + ca.sqrt(ad_len**2 + 1e-6))
    ad_len = 2*a - 0.5*((2*a - ad_len) + ca.sqrt((2*a - ad_len)**2 + 1e-6))
    
    # Force Calculation (Changed ^ to **)
    sliding_term = mu*w*Pmax*((4*ad_len/3) - ((ad_len**2)/a) + ((ad_len**3)/(3*a**2)))

    sigma   = ca.sqrt(sigma_x**2 + sigma_y**2 + eps**2)
    dir_x   = sigma_x / sigma      # stable because sigma >= sqrt(eps)
    dir_y   = sigma_y / sigma      # same

    Fx = (Cp/2)*sigma_x*(ad_len**2) + dir_x*sliding_term
    Fy = (Cp/2)*sigma_y*(ad_len**2) + dir_y*sliding_term
    Mz = (Cp/2)*sigma_y*((a/2)*ad_len**2 - ad_len**3/3) + dir_y*mu*w*Pmax*(ad_len**3/a - ad_len**4/(4*a**2) - ad_len**2)


    #Fx = (Cp/2)*sigma_x*(ad_len**2) + (sigma_x/sigma)*sliding_term
    #Fy = (Cp/2)*sigma_y*(ad_len**2) + (sigma_y/sigma)*sliding_term
    #Mz = (Cp/2)*sigma_y*(((a/2)*ad_len**2) - ((ad_len**3)/3)) + (sigma_y/sigma)*mu*w*Pmax*(((ad_len**3)/a) - ((ad_len**4)/(4*a**2)) - ad_len**2)
    
    # Combined Brush Tire Model Outputs
    return Fx, Fy, Mz, r

#====================================================================================================
#Registry
#====================================================================================================
tire_models = {"brush":combined_brush_tire_model}

def get_tire_model(name: str):
    if name not in tire_models:
        raise ValueError(f"Unknown tire model: {name}")
    return tire_models[name]

## End
#===================================================================================================================================================================================================================