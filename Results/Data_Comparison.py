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
print(list(telemetrydata["tel"].keys()))

arc_length = telemetrydata["tel"]["distance"]
speed = telemetrydata["tel"]["speed"]

fig, (ax1) = plt.subplots(1, 1, figsize=(10, 12), sharex=True)
ax1.plot(arc_length,speed,color = 'black')
ax1.set_xlabel('Centerline Arc Length [m]')
ax1.set_ylabel('Vehicle Speed [kmph]')
ax1.grid(True,alpha=0.3)
plt.tight_layout()
plt.show()
