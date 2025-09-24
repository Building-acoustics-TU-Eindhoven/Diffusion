"""
Simulate the energy decay in a cuboid room using the Finite Volume diffusion
equation.
"""
# %%
import os
import pandas as pd
import json
import numpy as np
import matplotlib.pyplot as plt
from Diffusion_Module_ADE.Auralization.Auralization import run_auralization

script_dir = os.path.dirname(os.path.abspath(__file__))
anechoic_file_path = os.path.join(script_dir, 'anechoic_file.wav') # Full path to the file
results_fvm_path = os.path.join(script_dir, 'resultsFVM.pkl') # Full path to the file

# %%
results = run_auralization(anechoic_file_path,results_fvm_path)
