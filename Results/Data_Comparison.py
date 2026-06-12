import os
import sys
import matplotlib
import numpy as np
import pandas as pd
import tkinter as tk
from tkinter import filedialog
import matplotlib.pyplot as plt
from types import SimpleNamespace
import json
import pickle

matplotlib.use('Qt5Agg')

os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')

# ======================================================================
# Load Telemetry Data
# ======================================================================

power_change = np.array([-15,-10,-5,-4,-3,-2,-1,0,1,2,3,4,5,10,15])
drag_eff = np.array([2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5])
df_eff = np.array([2.5,2.6,2.7,2.8,2.9,3,3.1,3.2,3.3,3.4,3.5])
mass_change = np.array([-10,-8,-6,-4,-2,0,2,4,6,8,10])

monza_power = np.array([78.465,77.745,77.085,76.96,76.838,76.718,76.6,76.484,76.369,76.258,76.149,76.042,75.936,75.428,74.962])
barcelona_power = np.array([71.348,70.909,70.524,70.451,70.379,70.31,70.242,70.176,70.112,70.049,69.987,69.927,69.868,69.591,69.336])
monza_power = monza_power/np.max(monza_power)
barcelona_power = barcelona_power/np.max(barcelona_power)

monza_drag = np.array([76.899,76.481,76.092,75.731,75.395,75.0795,74.787,74.514,74.257,74.0155,73.787])
barcelona_drag = np.array([71.805,71.584,71.384,71.201,71.033,70.879,70.737,70.605,70.483,70.369,70.262])
monza_drag = monza_drag/np.max(monza_drag)
barcelona_drag = barcelona_drag/np.max(barcelona_drag)

monza_df = np.array([76.483,76.194,75.909,75.63,75.354,75.083,74.813,74.548,74.285,74.024,73.769])
barcelona_df = np.array([72.36,71.902,71.453,71.017,70.59,70.176,69.776,69.392,69.019,68.664,68.324])
monza_df = monza_df/np.max(monza_df)
barcelona_df = barcelona_df/np.max(barcelona_df)

monza_mass = np.array([76.22,76.279,76.33,76.381,76.432,76.483,76.534,76.585,76.636,76.686,76.737])
barcelona_mass = np.array([71.02,71.083,71.146,71.208,71.27,71.332,71.394,71.457,71.519,71.582,71.644])
monza_mass = monza_mass/np.max(monza_mass)
barcelona_mass = barcelona_mass/np.max(barcelona_mass)


root = tk.Tk()
root.withdraw()
root.attributes('-topmost', True)

file_path = filedialog.askopenfilename(
    title="Select Telemetry Data",
    filetypes=[("JSON files", "*.json")]
)

if not file_path:
    raise FileNotFoundError("Telemetry file not selected")

with open(file_path, 'r') as f:
    telemetrydata = json.load(f)

print("\nLoaded Telemetry Data")
print("=" * 80)

# Show all available channels
print("\nAvailable Channels:\n")

for key, value in telemetrydata.items():

    if isinstance(value, list):
        print(f"{key:20s} | list | length = {len(value)}")
    elif isinstance(value, dict):
        print(f"{key:20s} | dict")
    else:
        print(f"{key:20s} | {type(value).__name__}")

# Convert to attribute-style access
telemetry = SimpleNamespace(**telemetrydata)
arc_length = np.array(telemetrydata["tel"]["distance"])
speed = telemetrydata["tel"]["speed"]

# Optimized Racing line
# Import Optimized Data
root = tk.Tk()
root.withdraw()
root.attributes('-topmost',True)

file_path1 = filedialog.askopenfilename(title="Select the Optimized Data",filetypes=[("pickle files","*.pkl"),])
if file_path1:
    with open(file_path1, 'rb') as f:
        Optimal_Solution = pickle.load(f)
    print("Loaded Optimal Solution")
else:
    raise FileNotFoundError("Optimal Solution file not found")


fig, (ax1,ax2) = plt.subplots(2, 1, figsize=(10, 12), sharex=True)
ax2.plot(arc_length - 236,speed,color = 'orange',label = '2025 Spanish GP Pole - Piastri (Lap 13 - 71.546 sec)') #- 236
ax2.plot(Optimal_Solution.arc_s,Optimal_Solution.u_opt*(18/5),color = 'black',label = 'Optimal Speed Profile Generated (70.176 sec)')
ax2.set_xlabel('Centerline Arc Length [m]')
ax2.set_ylabel('Vehicle Speed [kmph]')
ax2.set_xlim([Optimal_Solution.arc_s[0],Optimal_Solution.arc_s[-1]])
ax2.grid(True,alpha=0.3)
ax2.set_title('Barcelona Catalunya Speed Profile Comparison')
ax2.legend()
plt.tight_layout()
plt.show()


fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,figsize=(10,12))
# Peak Power Sensitivity
ax1.plot(power_change,monza_power,color='red',label='Monza')
ax1.plot(power_change,barcelona_power,color='blue',label='Barcelona Catalunya')
ax1.set_xlabel('Change in Engine Peak Power [%]')
ax1.set_ylabel('Normalized Lap Time [.]')
ax1.grid(True,alpha=0.3)
ax1.set_title('Peak Power Sensitivity (Baseline Peak Power = 730 kW)')
ax1.legend()

# Vehicle Weight Sensitivity
ax2.plot(mass_change,monza_mass,color='red',label='Monza')
ax2.plot(mass_change,barcelona_mass,color='blue',label='Barcelona Catalunya')
ax2.set_xlabel('Change in Vehicle Mass [kg]')
ax2.set_ylabel('Normalized Lap Time [.]')
ax2.set_title('Vehicle Mass Sensitivity (Baseline Vehicle Mass = 798 kg)')
ax2.legend()
ax2.grid(True,alpha=0.3)

# Drag Sensitivity
ax3.plot(drag_eff,monza_drag,color='red',label='Monza')
ax3.plot(drag_eff,barcelona_drag,color='blue',label='Barcelona Catalunya')
ax3.set_xlabel('Aerodynamic Efficiency ($- C_l/C_d$) [.]')
ax3.set_ylabel('Normalized Lap Time [.]')
ax3.set_title('Drag Sensitivity (Constant Lift Coefficient ($C_l$))')
ax3.legend()
ax3.grid(True,alpha=0.3)

# Downforce Sensitivity
ax4.plot(df_eff,monza_df,color='red',label='Monza')
ax4.plot(df_eff,barcelona_df,color='blue',label='Barcelona Catalunya')
ax4.set_xlabel('Aerodynamic Efficiency ($- C_l/C_d$)[.]')
ax4.set_ylabel('Normalized Lap Time [.]')
ax4.set_title('Downforce Sensitivity (Constant Drag Coefficient ($C_d$))')
ax4.legend()
ax4.grid(True,alpha=0.3)

plt.tight_layout()
plt.show()