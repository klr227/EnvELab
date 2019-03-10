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

# This code is included because there is a bug in the version of this code that is in epa.

def aeration_data(DO_column, dirpath):
    """This function extracts the data from folder containing tab delimited
    files of aeration data. The file must be the original tab delimited file.
    All text strings below the header must be removed from these files.
    The file names must be the air flow rates with units of micromoles/s.
    An example file name would be "300.xls" where 300 is the flowr ate in
    micromoles/s. The function opens a file dialog for the user to select
    the directory containing the data.

    Parameters
    ----------
    DO_column : int
        index of the column that contains the dissolved oxygen concentration
        data.

    dirpath : string
        path to the directory containing aeration data you want to analyze

    Returns
    -------
    filepaths : string list
        all file paths in the directory sorted by flow rate

    airflows : numpy array
        sorted array of air flow rates with units of micromole/s attached

    DO_data : numpy array list
        sorted list of numpy arrays. Thus each of the numpy data arrays can
        have different lengths to accommodate short and long experiments

    time_data : numpy array list
        sorted list of numpy arrays containing the times with units of seconds

    Examples
    --------

    """
    #return the list of files in the directory
    filenames = os.listdir(dirpath)
    #extract the flowrates from the filenames and apply units
    airflows = ((np.array([i.split('.', 1)[0] for i in filenames])).astype(np.float32))
    #sort airflows and filenames so that they are in ascending order of flow rates
    idx = np.argsort(airflows)
    airflows = (np.array(airflows)[idx])*u.umole/u.s
    filenames = np.array(filenames)[idx]

    filepaths = [os.path.join(dirpath, i) for i in filenames]
    #DO_data is a list of numpy arrays. Thus each of the numpy data arrays can have different lengths to accommodate short and long experiments
    # cycle through all of the files and extract the column of data with oxygen concentrations and the times
    DO_data=[epa.column_of_data(i,0,DO_column,-1,'mg/L') for i in filepaths]
    time_data=[(epa.column_of_time(i,0,-1)).to(u.s) for i in filepaths]
    aeration_collection = collections.namedtuple('aeration_results','filepaths airflows DO_data time_data')
    aeration_results = aeration_collection(filepaths, airflows, DO_data, time_data)
    return aeration_results

# The column of data containing the dissolved oxygen concentrations
DO_column = 2
dirpath = "/Users/kenrivero/Documents/EnvELab/Gas_Transfer/Aeration"
filepaths, airflows, DO_data, time_data = aeration_data(DO_column,dirpath)


# Plot the raw data
#for i in range(airflows.size):
#  plt.plot(time_data[i], DO_data[i],'-')
#plt.xlabel(r'$time (s)$')
#plt.ylabel(r'Oxygen concentration $\left ( \frac{mg}{L} \right )$')
#plt.legend(airflows.magnitude)
#plt.show()

#delete data that is less than 2 or greater than 6 mg/L
DO_min = 2 * u.mg/u.L
DO_max = 6 * u.mg/u.L
for i in range(airflows.size):
  idx_start = (np.abs(DO_data[i]-DO_min)).argmin()
  idx_end = (np.abs(DO_data[i]-DO_max)).argmin()
  time_data[i] = time_data[i][idx_start:idx_end] - time_data[i][idx_start]
  DO_data[i] = DO_data[i][idx_start:idx_end]

########################Part 2 ############################

index = []
for i in range(0,airflows.size,4):
  plt.plot(time_data[i], DO_data[i],'-')
  index.append(i)
plt.xlabel(r'$time (s)$')
plt.ylabel(r'Oxygen concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(airflows[index].magnitude)
plt.savefig('/Users/kenrivero/Documents/EnvELab/Gas_Transfer/DO_vs_time_5graphs')
plt.show()
#number3: Calculate C_star with the assumptions of temperature at 22 degrees celsius and pressure of 1 atm.
temp = 22 * u.degC
pressure = 1 * u.atm

C_star = epa.O2_sat(pressure,temp)

#number 4: estimate K_vl using linear regression for each data set
#ln((C*-C)/(C*- C_0)) = -k_vl*(t-t0)
#put all data in an array
#find the C_0 values
#This will define the variables that we will use in the for loop
C_0 =[]
t_0 = []

#This for loop retrieves all the initial dissolved oxygen concentrations for each airflow
for data in DO_data:
  C_0.append(data[0])
#This  forloop retrieves all the intial times for each airflow.
#NOTE: All initial times are 0.
for data in time_data:
  t_0.append(data[0])

#This will initialize the variable for kvl which will be used in the for loop
k_vl = [] #initialize k_vl values

#This for loop calculates the kvl value using linear regression
for i in range(airflows.size):
  C_values = DO_data[:][i]
  t_values = time_data[:][i]
  y = np.log((C_star-C_values)/(C_star-C_0_values[i]))
  x = t_values - t_0[0]
  slope, intercept, r_value,p_value,std_err =stats.linregress(x,y)
  k_vl.append(slope)


#number 5: plot model and actual data -- will need to create a model
#This random function is used to index any of the flow rates randomnly
rand = randint(0,airflows.size)
#Time array is made without units because our k_vl array doesn't have units.
t = np.linspace(0,100)
#Calculate for the concentration of dissolved oxygen using equation (103) from the lab manual and solving for C.
#Used the C_0 values from a random airflow and K_vl values from the same airflow
C = C_star - ((C_star-C_0[rand])* np.exp(-k_vl[rand]*t))

plt.plot(t,C,'g')
plt.plot(time_data[:][rand],DO_data[:][rand])
plt.legend(['Model Concentration','Experimental Concentration %i'])
plt.xlabel('Time (seconds)')
plt.ylabel('Concentration of Oxygen (mg/L)')
plt.savefig('/Users/kenrivero/Documents/EnvELab/Gas_Transfer/ModelandExperimentData_Cgraph')
plt.show()
#number 6: plot k_vl as a function of air flow rate

plt.scatter(airflows, k_vl)
plt.xlabel('Air Flow Rate in micromoles/second')
plt.ylabel('k_vl values')
plt.savefig('/Users/kenrivero/Documents/EnvELab/Gas_Transfer/k_vl_vs_Airflow_Rates')
plt.show()
#number 7: plot OTE as a function of air flow rate with oxygen deficit (C*-C) set at 6 mg/L
f_O2 = 0.21
V = 750 *u.mL
O2_deficit = 6 * u.mg/u.L
MWO2 = 15.999*2 * u.g/u.mole
OTE = []

#This for loop will calculate the OTE
for k in range(airflows.size):
  value = V*k_vl[k]*(O2_deficit)/(f_O2*airflows[k]*MWO2)
  OTE.append(value)

plt.scatter(airflows,OTE)
plt.xlabel('Air Flow Rate in micromoles/second')
plt.ylabel('Oxygen Transfer Efficiency')
plt.savefig('/Users/kenrivero/Documents/EnvELab/Gas_Transfer/OTE_vs_Airflow_Rate')
plt.show()

#number 8: Comment on the oxygen transfer efficiency and the trend or trends that you observe.
#Oxygen transfer efficiency decreases as air flow increases
#number 9: Propose a change to the experimental apparatus that would increase the efficiency.

#number 10: Verify that your report and graphs meet the requirements.
