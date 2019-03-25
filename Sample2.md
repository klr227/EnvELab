```python
from aguaclara.core.units import unit_registry as u
import aguaclara.research.environmental_processes_analysis as epa
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
import collections
import os
from pathlib import Path
from random import *
import aguaclara.core.utility as ut
##Load a data file for a reactor with baffles trial 4 (4 baffles with gap).
trial_path = ['https://raw.githubusercontent.com/klr227/EnvELab/master/Reactor_Characteristics/baffle_test_1.xls','https://raw.githubusercontent.com/klr227/EnvELab/master/Reactor_Characteristics/baffle_test_2.xls','https://raw.githubusercontent.com/klr227/EnvELab/master/Reactor_Characteristics/baffle_test_3.xls','https://raw.githubusercontent.com/klr227/EnvELab/master/Reactor_Characteristics/baffle_test_4.xls']
len(trial_path)
trial_firstrow = epa.notes(trial_path[0]).last_valid_index() + 1
trial_time_data = (epa.column_of_time(trial_path[0],trial_firstrow,-1)).to(u.s)
trial_concentration_data = epa.column_of_data(trial_path[0],trial_firstrow,1,-1,'mg/L')

#I noticed that the initial concentration measured by the photometer wasn't
#zero. This suggests that there may have been a small air bubble in the
#photometer or perhaps there was some other anomaly that was causing the
#photometer to read a concentration that was higher than was actually present in
#the reactor. To correct for this I subtracted the initial concentration reading
#from all of the data. This was based on the assumption that the concentration
#measurement error persisted for the entire experiment.#
density_H2O = 997 *u.g/u.L
trial_concentration_data = trial_concentration_data - trial_concentration_data[0]
trial_mass = [137,1537,1557,1598] * u.g
trial_V = trial_mass[0]/density_H2O
Q = 380 * u.mL/u.min
trial_theta_hydraulic = (trial_V/Q).to(u.s)
trial_C_bar_guess = np.max(trial_concentration_data)/2
#use solver to get the CMFR parameters
trial_CMFR = epa.Solver_CMFR_N(trial_time_data, trial_concentration_data, trial_theta_hydraulic, trial_C_bar_guess)
trial_CMFR_model = (trial_CMFR.C_bar*epa.E_CMFR_N(trial_time_data/trial_CMFR.theta, trial_CMFR.N)).to(u.mg/u.L)

trial_AD = epa.Solver_AD_Pe(trial_time_data, trial_concentration_data, trial_theta_hydraulic, trial_C_bar_guess)
trial_AD_model = (trial_AD.C_bar*epa.E_Advective_Dispersion((trial_time_data/trial_AD.theta).to_base_units(), trial_AD.Pe)).to(u.mg/u.L)
trial_AD_model

trial_firstrow
trial_time_data
trial_concentration_data
trial_mass
trial_theta_hydraulic
trial_C_bar_guess
