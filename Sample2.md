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

N0 = 100

diazonin_half = 6.20 * u.days
diaz_exp_time = [0,15,30,60,90,180] * u.days
diaz_N_remaining = [100,18.7,1.5,0.30,0.01,0.01]
diazonin_t = np.linspace(0,180) * u.days
diazonin_Nt = N0 * (.5**((diazonin_t/diazonin_half).magnitude))
plt.plot(diazonin_t,diazonin_Nt)
plt.plot(diaz_exp_time,diaz_N_remaining,'ro')
plt.xlabel('Time (days)')
plt.ylabel('Diazonin Remaining %')
plt.legend(['Diazonin Decay Function','Diazonin Decomposition Data'])
plt.show()


benefin_half = 37.58 * u.days
benefin_exp_time = [0,30,60,90,180] * u.days
benefin_N_remaining = [100,57.5,29,14.5,7.5]
benefin_t = np.linspace(0,500) * u.days
benefin_Nt = N0 * (.5**((benefin_t/benefin_half).magnitude))
plt.plot(benefin_t,benefin_Nt)
plt.plot(benefin_exp_time,benefin_N_remaining,'ro')
plt.xlabel('Time (days)')
plt.ylabel('Benefin Remaining %')
plt.legend(['Benefin Decay Function','Benefin Decomposition Data'])
plt.show()

terbutol_half = 186.92 * u.days
terbutol_exp_time = [0,30,60,90,180] * u.days
terbutol_N_remaining = [100,94.3,80.2,73.7,51.3]
terbutol_t = np.linspace(0,3000) * u.days
terbutol_Nt = N0 * (.5**((terbutol_t/terbutol_half).magnitude))
plt.plot(terbutol_t,terbutol_Nt)
plt.plot(terbutol_exp_time,terbutol_N_remaining,'ro')
plt.xlabel('Time (days)')
plt.ylabel('Terbutol Remaining %')
plt.legend(['Terbutol Decay Function','Terbutol Decomposition Data'])
plt.show()
