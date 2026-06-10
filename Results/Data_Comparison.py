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
ax2.plot(arc_length,speed,color = 'blue',label = '2025 Italian GP Pole - Verstappen (Lap 17 - 78.792 sec)') #- 236
ax2.plot(Optimal_Solution.arc_s,Optimal_Solution.u_opt*(18/5),color = 'black',label = 'Optimal Speed Profile Generated (76.484 sec)')
ax2.set_xlabel('Centerline Arc Length [m]')
ax2.set_ylabel('Vehicle Speed [kmph]')
ax2.set_xlim([Optimal_Solution.arc_s[0],Optimal_Solution.arc_s[-1]])
ax2.grid(True,alpha=0.3)
ax2.set_title('Monza Speed Profile Comparison')
ax2.legend()
plt.tight_layout()
plt.show()
