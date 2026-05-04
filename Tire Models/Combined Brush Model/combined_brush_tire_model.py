#====================================================================================================================================================================================
# This project utilizes CasADi for symbolic differentiation and optimization.
# Andersson, J. A. E., et al. (2019). "CasADi: a software framework for 
# nonlinear optimization and optimal control." Mathematical Programming Computation.
#====================================================================================================================================================================================

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
import json
import math
from casadi import *

os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')

def combined_brush_tire_model(u,omega,alpha,Fz,vehicle):
    mu = vehicle.mu
    mu0 = vehicle.mu0
    R = vehicle.R
    w = vehicle.w
    Cp = vehicle.Cp
    kv = vehicle.kv
    # 1. Effective Radius Calculation
    r = R - (Fz/kv)
    a = sqrt(R**2 - r**2)
    Pmax = 3*Fz/(4*a*w)

    # 2. Slip Ratio Definition
    kappa = (r*omega - u)/u
    sigma_x = kappa/(1 + kappa)
    sigma_y = tan(alpha)/(1 + kappa)
    sigma = sqrt(sigma_x**2 + sigma_y**2)

    # 3. Adhesion Length Estimation
    ad_len = 2*a - ((Cp*sigma*a**2)/(mu0*w*Pmax))
    ad_len = fmax(ad_len, 0)
    ad_len = fmin(ad_len,2*a)

    # 4. Force Calculation (Changed ^ to **)
    sliding_term = mu*w*Pmax*((4*ad_len/3) - ((ad_len**2)/a) + ((ad_len**3)/(3*a**2)))
    
    Fx = (Cp/2)*sigma_x*(ad_len**2) + (sigma_x/sigma)*sliding_term
    Fy = (Cp/2)*sigma_y*(ad_len**2) + (sigma_y/sigma)*sliding_term
    
    Mz = (Cp/2)*sigma_y*(((a/2)*ad_len**2) - ((ad_len**3)/3)) + (sigma_y/sigma)*mu*w*Pmax*(((ad_len**3)/a) - ((ad_len**4)/(4*a**2)) - ad_len**2)
    
# Combined Brush Tire Model Outputs
    return Fx, Fy, Mz, r

## End
#===================================================================================================================================================================================================================