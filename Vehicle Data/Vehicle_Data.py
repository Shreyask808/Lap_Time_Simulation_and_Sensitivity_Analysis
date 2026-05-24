#====================================================================================================================================================================================
# This project utilizes CasADi for symbolic differentiation and optimization.
# Andersson, J. A. E., et al. (2019). "CasADi: a software framework for 
# nonlinear optimization and optimal control." Mathematical Programming Computation.
#====================================================================================================================================================================================
import os
import sys
import numpy as np
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from types import SimpleNamespace
import json
import math
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')

# ======================================================================================================================================================================================================================
# Vehicle Modeling Data
# ======================================================================================================================================================================================================================
# Vehicle Dimension Data
vehicle = SimpleNamespace()
vehicle.l = 3.6                                                             # Wheelbase in m
vehicle.d = 1.65                                                            # Distance between rear axle and cg in m
vehicle.hcg = 0.28                                                          # Height of the cg from ground in m
vehicle.hcp = 0.20                                                          # Height of the cp from ground in m
vehicle.a = 0.5                                                             # Position of Center of Pressure (cp) behind cg in m
vehicle.Wf = 1.7                                                            # Front Trackwidth in m
vehicle.Wr = 1.7                                                            # Rear Trackwidth in m
vehicle.Weq = (vehicle.Wf + vehicle.Wr)/2                                   # Equivalent Trackwidth of the vehicle in m
vehicle.Droll = 0.7                                                         # Lateral Load Distribution at the front axle [.]
vehicle.delta_max = 0.4                                                     # Max Steering Angle in rad
vehicle.delta_min = -0.4                                                    # Min Steering Angle in rad
vehicle.A = 1.2                                                             # Frontal Area of the vehicle in m^2
vehicle.Cd = 1.00                                                           # Drag Coefficient of the Vehicle [.] 
vehicle.Cl = -3.5                                                           # Lift Coefficient of the Vehicle [.]

#=====================================================================================================================================================================================================================
# Vehicle Mass and Interia Properties
vehicle.m = 800                                                             # Vehicle Mass in kg
vehicle.Iz = 900                                                            # Vehicle Moment of Interia about cg in kg.m^2 
# Vehicle Axle Interia
vehicle.Ifl = 1.5                                                           # Front Left Wheel Moment of Interia in kg.m^2
vehicle.Ifr = 1.5                                                           # Front Right Wheel Moment of Interia in kg.m^2
vehicle.Irl = 1.5                                                           # Rear Left Wheel Moment of Interia in kg.m^2
vehicle.Irr = 1.5                                                           # Rear Right Wheel Moment of Interia in kg.m^2

#=====================================================================================================================================================================================================================
# Vehicle Tire Data
vehicle.R = 0.360                                                           # Tire radius in m
vehicle.w = 0.305                                                           # Tire width in m
vehicle.Cp = 4e6                                                            # Tire stiffness per unit length in N/m^2
vehicle.mu0 = 2.5                                                           # Static friction coefficient [.]
vehicle.mu = 2.0                                                            # Sliding friction coefficient [.]
vehicle.kv = 2e5                                                            # Vertical stiffness of the tire in N.m
vehicle.alpha_max = math.radians(30.0)                                      # Maximum tire slip angle in rad
vehicle.Crr = 0                                                             # Tire rolling resistance coefficient [.]

#====================================================================================================================================================================================================================
# Vehicle Powertrain Limits
vehicle.peakdrivingpower = 880e3                                            # Vehicle peak power in W
vehicle.peakbrakingpower = -1500e3                                          # Vehicle peak power in W
vehicle.peakdrivingtorque = 1000                                            # Vehicle peak driving torque per tire in N.m
vehicle.peakbrakingtorque = -3500                                           # Vehicle peak braking torque per tire in N.m

#=====================================================================================================================================================================================================================
# Miscellaneous Data
vehicle.rho = 1.2                                                           # Density of Air in kg/m^3
vehicle.g = 9.81                                                            # Gravitational Constant in m/s^2
vehicle.umax = ((2*vehicle.peakdrivingpower)/(vehicle.rho*vehicle.Cd*vehicle.A))**(1/3) # Maximum Longitudinal Velocity in m/s
vehicle.vmax = vehicle.umax*np.tan(vehicle.alpha_max)                       # Maximum Lateral Velocity in m/s

#=====================================================================================================================================================================================================================
# Save Track Data
root = tk.Tk()
root.withdraw()
files = [('json Files','*.json'), ('All Files', '*.*')]
full_save_path = filedialog.asksaveasfilename(
    defaultextension=".json",
    filetypes=files,
    title="Save Vehicle Data"
)    

if full_save_path:
    data_to_save = {
        k: v.tolist() if isinstance(v, np.ndarray) else v 
        for k, v in vars(vehicle).items()
    }
    try:
        with open(full_save_path, "w") as f:
            json.dump(data_to_save, f, indent=4)
        
        print(f"Vehicle Data saved successfully to: {full_save_path}")
    except Exception as e:
        print(f"Error saving file: {e}")
else:
    print("Save cancelled.")

## End
#===================================================================================================================================================================================================================