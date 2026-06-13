#====================================================================================================================================================================================
# Legendre-Gauss-Radau (LGR) Pseudospectral Optimal Control
# Reference: Limebeer & Perantoni (2015), JDSMC 137, DOI: 10.1115/1.4029466
#====================================================================================================================================================================================

import os
import sys
import matplotlib
import numpy as np
import pandas as pd
import tkinter as tk
import scipy.ndimage as ndimage
from scipy.ndimage import gaussian_filter1d
from scipy.special import roots_jacobi
from tkinter import filedialog
import matplotlib.pyplot as plt
from types import SimpleNamespace
import json
import pickle
import casadi as ca
from casadi import *
matplotlib.use('Qt5Agg')
import plotly.graph_objects as go
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Tire Models'))
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')

