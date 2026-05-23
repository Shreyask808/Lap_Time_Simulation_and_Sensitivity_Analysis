#=========================================================================================================================================================================================================================
# This project utilizes CasADi for symbolic differentiation and optimization.
# Andersson, J. A. E., et al. (2019). "CasADi: a software framework for 
# nonlinear optimization and optimal control." Mathematical Programming Computation.
#=========================================================================================================================================================================================================================

import os
import sys
import numpy as np
import pandas as pd
import tkinter as tk
import scipy.ndimage as ndimage
from scipy.ndimage import gaussian_filter1d
from tkinter import filedialog
import matplotlib.pyplot as plt
from types import SimpleNamespace
import json
import math
import casadi as ca
from casadi import *
import matplotlib
matplotlib.use('Qt5Agg')
import plotly.graph_objects as go
sys.path.append(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'Tire Models'))
from tire_models import get_tire_model
os.system('cls' if os.name == 'nt' else 'clear')
plt.close('all')
 
def Seven_DOF_Handling_Model_2D(vehicle,track_data,N,name,x0_ini):
