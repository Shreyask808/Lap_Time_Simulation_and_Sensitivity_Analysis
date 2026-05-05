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
vehicle.d = 1.6                                                             # Distance between rear axle and cg in m
vehicle.h = 0.3                                                             # Height of the cg from ground in m
vehicle.a = 0.5                                                             # Position of Center of Pressure (cp) behind cg in m
vehicle.Wf = 2                                                              # Front Trackwidth in m
vehicle.Wr = 2                                                              # Rear Trackwidth in m
vehicle.Weq = (vehicle.Wf + vehicle.Wr)/2                                   # Equivalent Trackwidth of the vehicle in m
vehicle.Droll = 0.5                                                         # Lateral Load Distribution at the front axle [.]
vehicle.delta_max = 0.5                                                     # Max Steering Angle in rad
vehicle.delta_min = -0.5                                                    # Min Steering Angle in rad
vehicle.A = 1.5                                                             # Frontal Area of the vehicle in m^2
vehicle.Cd = 0.22                                                           # Drag Coefficient of the Vehicle [.] 
vehicle.Cl = -0.5                                                           # Lift Coefficient of the Vehicle [.]

#=====================================================================================================================================================================================================================
# Vehicle Mass and Interia Properties
vehicle.m = 660                                                             # Vehicle Mass in kg
vehicle.Iz = 450                                                            # Vehicle Moment of Interia about cg in kg.m^2 
# Vehicle Axle Interia
vehicle.Ifl = 5                                                             # Front Left Wheel Moment of Interia in kg.m^2
vehicle.Ifr = 5                                                             # Front Right Wheel Moment of Interia in kg.m^2
vehicle.Irl = 5                                                             # Rear Left Wheel Moment of Interia in kg.m^2
vehicle.Irr = 5                                                             # Rear Right Wheel Moment of Interia in kg.m^2

#=====================================================================================================================================================================================================================
# Vehicle Tire Data
vehicle.R = 0.334                                                           # Tire radius in m
vehicle.w = 0.225                                                           # Tire width in m
vehicle.Cp = 4e6                                                            # Tire stiffness per unit length in N/m^2
vehicle.mu0 = 0.9                                                           # Static friction coefficient [.]
vehicle.mu = 0.6                                                            # Sliding friction coefficient [.]
vehicle.kv = 1e6                                                            # Vertical stiffness of the tire in N.m
vehicle.alpha_max = d=math.radians(10)                                      # Maximum tire slip angle in rad
vehicle.Crr = 0                                                             # Tire rolling resistance coefficient [.]

#====================================================================================================================================================================================================================
# Vehicle Powertrain Limits
vehicle.peakpower = 450e3                                                   # Vehicle peak power in W
vehicle.peakdrivingtorque = 1000                                            # Vehicle peak driving torque per tire in N.m
vehicle.peakbrakingtorque = - 2000                                          # Vehicle peak braking torque per tire in N.m

#=====================================================================================================================================================================================================================
# Miscellaneous Data
vehicle.rho = 1.2                                                           # Density of Air in kg/m^3
vehicle.g = 9.81                                                            # Gravitational Constant in m/s^2
vehicle.umax = 80                                                           # Maximum Longitudinal Velocity in m/s

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