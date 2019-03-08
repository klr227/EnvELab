```python
from aguaclara.core.units import unit_registry as u
u.define('equivalent = mole = eq')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import aguaclara.research.environmental_processes_analysis as epa
from scipy import optimize
import math

data_file_path_1 = "https://raw.githubusercontent.com/klr227/EnvELab/master/Gas_Transfer/350_Flow_Rate.tsv"
df_1 = pd.read_csv(data_file_path_1,delimiter='\t')
time_350 = epa.column_of_time(data_file_path_1,0).to(u.s)
DO_350 = df_1.iloc[1:,2].values * u.mg/u.L

data_file_path_4 = "https://raw.githubusercontent.com/klr227/EnvELab/master/Gas_Transfer/475_Flow_Rate.tsv"
df_4 = pd.read_csv(data_file_path_4,delimiter='\t')
time_475 = epa.column_of_time(data_file_path_4,0).to(u.s)
DO_475 = df_4.iloc[1:,2].values * u.mg/u.L

data_file_path_2 = "https://raw.githubusercontent.com/klr227/EnvELab/master/Gas_Transfer/525_Flow_Rate.tsv"
df_2 = pd.read_csv(data_file_path_2,delimiter='\t')
time_525 = epa.column_of_time(data_file_path_2,0).to(u.s)
DO_525 = df_2.iloc[1:,2].values * u.mg/u.L

data_file_path_3 = "https://raw.githubusercontent.com/klr227/EnvELab/master/Gas_Transfer/575_Flow_Rate.tsv"
df_3 = pd.read_csv(data_file_path_3,delimiter='\t')
time_575 = epa.column_of_time(data_file_path_3,0).to(u.s)
DO_575 = df_3.iloc[1:,2].values * u.mg/u.L

data_file_path_5 = "https://raw.githubusercontent.com/klr227/EnvELab/master/Gas_Transfer/650_Flow_Rate.tsv"
df_5 = pd.read_csv(data_file_path_5,delimiter='\t')
time_650 = epa.column_of_time(data_file_path_5,0).to(u.s)
DO_650 = df_5.iloc[1:,2].values * u.mg/u.L


columns = df_1.columns
columns
fig, ax = plt.subplots()
plt.plot(time_350,DO_350,'y')
plt.plot(time_475,DO_475,'b')
plt.plot(time_525,DO_525,'r')
plt.plot(time_575,DO_575,'g')
plt.plot(time_650,DO_650,'k')
plt.xlabel('Time(seconds)')
plt.ylabel('Dissolved Oyxgen (mg/L)')
ax.legend(['350 microM/s','475 microM/s', '525 microM/s', '575 microM/s','650 microM/s'])
#plt.savefig('/Users/kenrivero/Documents/EnvELab/Gas_Transfer/DO_vs_time')
plt.show()

temp = 22 * u.degC
pressure = 1 * u.atm

C_star = epa.O2_sat(pressure,temp)
C_star
