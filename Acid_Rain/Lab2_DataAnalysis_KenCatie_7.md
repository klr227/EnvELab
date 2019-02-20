Lab 2: Acid Rain
Data Analysis
Catherine Johnson (caj92)
Ken Rivero Rivera (klr227)

Question 1:

If enough NaHCO3 is added to the lake to maintain the original 50 ueq/L ANC for more than three hydraulic residence times, with the stir bar turned off, we would expect that not all of the NaHCO3 would dissolve in the lake, and some would sink the bottom of the lake. Because some of the NaHCO3 is at the bottom of the lake, and the influent acid rain is not being well mixed, parts of the lake are not contributing to the ANC, and it is overall a lower ANC, and the lake will become acidic at a faster rate.

Question 2: What are some of the complicating factors you might find in attempting to remediate a lake using CaCO3?
The solubility of NaHCO3 is significantly higher than that of CaCO3. Their solubility in water is 96 g/L and 0.013 g/L, respectively, in 20 degrees Celsius water. This would make it different to fully dissolve CaCO3 if used to remediate a lake. We used CaCO3 in the second part of our lab and noticed that even with a stir bar, there were still particles of CaCO3 in the water.
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

#Given Values
k1 = 10**-6.3
k2 = 10**-10.3
kH = 10**-1.5 *u.mol/u.L/u.atmosphere
PCo2 = 10 **-3.5 *u.atmosphere
kW = 10**-14

#Section 1: This section contains the data analysis for the lake with NaHCO3

#  Part 1: Plot measured pH of the lake versus dimensionless hydraulic residence time (t/θ).

#Calculate residence time
mass_tankwater = 3766 * u.g
density_water = 0.997 * u.kg/u.L
vol_tankwater = (mass_tankwater/density_water).to(u.L)
flow_rate = (101*u.mL/(25*u.s)).to(u.L/u.s)
theta = vol_tankwater/flow_rate

# import pH data from spreadsheet located on GitHub
data_file_path_1 = "https://raw.githubusercontent.com/klr227/EnvELab/master/Acid_Rain/acid_rain_lab1.tsv"
#df = pd.read_csv(data_file_path,delimiter='\t')

# Use the epa.notes function to find what you've commented in the data file.
lakepHnotes = epa.notes(data_file_path_1)
# set the start index of the data file to one past the note indicating the start.
start_1 = 240

# Extract the pH data starting from where pump starts comment is located. pH data is in column 1
lakepH = epa.column_of_data(data_file_path_1,start_1,1)
# Extract the corresponding time data and convert to seconds
time = epa.column_of_time(data_file_path_1,start_1).to(u.s)

# Calculate dimensionless hydraulic residence time
res = time/theta

#Now plot the graph
fig, ax = plt.subplots()
ax.plot(res,lakepH,'r')
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('pH')
plt.title('NaHCO3 pH vs. Dimensionless Hyd. Residence Time')
ax.legend(['pH'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/pHgraph')
plt.show()


# Part 2: CMFR graph
pHin = 3
C_initial_open = 50*(10**-6) * u.eq/u.L
C_influent = (-10)**-pHin* u.eq/u.L
Output_conc_CMFR = epa.CMFR(res, C_initial_open, C_influent)

fig, ax = plt.subplots()
ax.plot(res,Output_conc_CMFR,'r')
ax.legend(['Conservative ANC'])
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('ANC')
plt.title('NaHCO3 Exp. ANC in lake output(CMFR) vs hyd. residence time')
ax.legend(['Conservative ANC'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/Output_conc_CMFR')
plt.show()


# Part 3: ANC in a closed system
#Calculate total Carbonate concentration in system in mol/L
massNaHCO3 = 0.625 *u.g
MmNaHCO3 = 84.007 * u.g/u.mol
Total_Carbonates_conc = massNaHCO3/MmNaHCO3/vol_tankwater
ANC_out_closed = epa.ANC_closed(lakepH,Total_Carbonates_conc)


#Part 4: ANC in open system
ANC_out_open = epa.ANC_open(lakepH)
fig, ax = plt.subplots()
plt.plot(res,ANC_out_closed,'r')
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('ANC')
plt.plot(res,ANC_out_open,'b')
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('ANC')
plt.title('NaHCO3 Exp. ANC of closed/open System vs Hyd. Res. Time')
ax.legend(['ANC closed system','ANC open system'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/ANC_openclosed')
plt.show()


#Section 2: This section contains information for the CaCO3 solution
# import pH data from spreadsheet
data_file_path_2 = "https://raw.githubusercontent.com/klr227/EnvELab/master/Acid_Rain/acid_rain_lab1_pt2.tsv"
#df = pd.read_csv(data_file_path,delimiter='\t')

# Use the epa.notes function to find what you've commented in the data file.
lakepHnotes_2 = epa.notes(data_file_path_2)
# set the start index of the data file to one past the note indicating the start.
start_2 = 40

# Extract the pH data starting from where pump starts comment is located. pH data is in column 1
lakepH_2 = epa.column_of_data(data_file_path_2,start_2,1)
# Extract the corresponding time data and convert to seconds
time_2 = epa.column_of_time(data_file_path_2,start_2).to(u.s)

# Calculate dimensionless hydraulic residence time
res_2 = time_2/theta

#Now plot the graph
fig, ax = plt.subplots()
ax.plot(res_2,lakepH_2,'r')
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('pH')
plt.title('CaCO3 pH vs. Dimensionless Hyd. Residence Time')
ax.legend(['pH'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/pHgraph_2')
plt.show()


# Part 2: CMFR graph
pHin = 3
C_initial_open_2 = 50*(10**-6) * u.eq/u.L
C_influent = (-10)**-pHin* u.eq/u.L
Output_conc_CMFR_2 = epa.CMFR(res_2, C_initial_open_2, C_influent)

fig, ax = plt.subplots()
ax.plot(res_2,Output_conc_CMFR_2,'r')
ax.legend(['Conservative ANC'])
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('ANC')
plt.title('CaCO3 Exp. ANC in lake output (CMFR) vs Hyd residence time')
ax.legend(['Conservative ANC'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/Output_conc_CMFR_2')
plt.show()


# Part 3: ANC in a closed system
#Calculate total Carbonate concentration in system in mol/L
massCaCO3 = 0.371 *u.g
MmCaCO3 = 100.0869 * u.g/u.mol
Total_Carbonates_conc_2 = massCaCO3/MmCaCO3/vol_tankwater
ANC_out_closed_2 = epa.ANC_closed(lakepH_2,Total_Carbonates_conc_2)


#Part 4: ANC in open system
ANC_out_open_2 = epa.ANC_open(lakepH_2)
fig, ax = plt.subplots()
plt.plot(res_2,ANC_out_closed_2,'r')
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('ANC')
plt.plot(res_2,ANC_out_open_2,'b')
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('ANC')
plt.title('CaCO3 Exp. ANC of closed/open System vs Hyd Res Time')
ax.legend(['ANC closed system','ANC open system'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/ANC_openclosed_2')
plt.show()
```
![pHgraph](https://github.com/klr227/EnvELab/blob/master/Acid_Rain/pHgraph.png)
Figure 1. Measured pH of the lake versus dimensionless hydraulic residence time for NaHCO3
![pHgraph_2](https://github.com/klr227/EnvELab/blob/master/Acid_Rain/pHgraph_2.png)
Figure 2. Measured pH of the lake versus dimensionless hydraulic residence time for CaCO3
![ANC_openclosed](https://github.com/klr227/EnvELab/blob/master/Acid_Rain/ANC_openclosed.png)
Figure 3. Expected ANC in the lake effluent versus the hydraulic residence time for open and closed systems for NaHCO3
![ANC_openclosed_2](https://github.com/klr227/EnvELab/blob/master/Acid_Rain/ANC_openclosed_2.png)
Figure 4. Expected ANC in the lake effluent versus the hydraulic residence time for open and closed systems for CaCO3
![Output_conc_CMFR](https://github.com/klr227/EnvELab/blob/master/Acid_Rain/Output_conc_CMFR.png)
Figure 5. Expected ANC in the lake effluent versus the hydraulic residence time for CMFR for CaCO3
![Output_conc_CMFR_2](https://github.com/klr227/EnvELab/blob/master/Acid_Rain/Output_conc_CMFR_2.png)
Figure 6. Expected ANC in the lake effluent versus the hydraulic residence time for CMFR for CaCO3
