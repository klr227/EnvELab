Below is the code for Lab's 2 and 3 data analysis.
```python
#Import python packages.
from aguaclara.core.units import unit_registry as u
u.define('equivalent = mole = eq')
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy import stats
import aguaclara.research.environmental_processes_analysis as epa
from scipy import optimize
import math

#$$$$$$$$$$$$$$$$$$$$$Data Analysis for Lab 2 $$$$$$$$$$$$$$$$$$$$

#Given Values
k1 = 10**-6.3
k2 = 10**-10.3
kH = 10**-1.5 *u.mol/u.L/u.atmosphere
PCo2 = 10 **-3.5 *u.atmosphere
kW = 10**-14

#Section 1: This section contains the data analysis for the lake with NaHCO3
#Part 1: Plot measured pH of the lake versus dimensionless hydraulic residence time (t/θ).

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
ax.legend(['pH'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/Acid_Rain/pHgraph')
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
ax.legend(['Conservative ANC'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/Acid_Rain/Output_conc_CMFR')
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
plt.plot(res,ANC_out_open,'b')
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('ANC')
ax.legend(['ANC closed system','ANC open system','Measured ANC'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/Acid_Rain/ANC_openclosed')
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
ax.legend(['pH'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/Acid_Rain/pHgraph_2')
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
ax.legend(['Conservative ANC'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/Acid_Rain/Output_conc_CMFR_2')
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
plt.plot(res_2,ANC_out_open_2,'b')
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('ANC')
ax.legend(['ANC closed system','ANC open system'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/Acid_Rain/ANC_openclosed_2')
plt.show()

#$$$$$$$$$$$$$$$$$$$$$ Data Analysis for Lab 3 $$$$$$$$$$$$$$$$$$$$
#Import datafile from GitHub
data_file_path =('https://raw.githubusercontent.com/klr227/EnvELab/master/Acid_Rain/lab3_0minsample.tsv')
df = pd.read_csv(data_file_path,delimiter='\t')

# The column headers can be access by using the list command
#list(df)
#columns = df.columns
#print(columns)

#Part 1 ***************************************************************
#Assign variables to the data in the dataframe (Titrant volume (mL), pH, sample volume (mL), titrant normality, equivalent volume (mL), ANC (eq/L), R^2). F1 is assigned later in Part 2.
titrant_vol = df.iloc[:,0].values * u.mL
pH = df.iloc[:,1].values
Sample_vol = df.iloc[0,3] * u.mL
Titrant_norm = df.iloc[0,4]
equivalent_vol = df.iloc[0,5]
ANC = df.iloc[0,6] * u.eq/u.L
R2 = df.iloc[0,7]

#Plot pH as a function of titrant volume. Label the equivalent volume of titrant. Label the 2 regions of the graph where pH changes slowly with the dominant reaction that is occurring.
fig, ax = plt.subplots()
ax.plot(titrant_vol,pH,'r')
plt.xlabel('Titrant volume (mL)')
plt.ylabel('pH')
ax.legend(['pH'])
ax.grid(True)
ax.annotate('First Reaction', xy=(.35, 7), xytext=(1, 7.5),
            arrowprops=dict(facecolor='black', shrink=1))
ax.annotate('Second Reaction', xy=(1.75, 4.5), xytext=(.25, 4),
            arrowprops=dict(facecolor='black', shrink=0.05))
plt.savefig('/Users/kenrivero/Documents/EnvELab/Acid_Rain/pHvsTitvol')
plt.show()


#Part 2 ****************************************************************
#Assign variables to the data in the dataframe (F1 and titrant volume). F1 and titrant_vol_F1: Includes points that are not needed for the linear regression
F1 = df.iloc[0:7,2].values
titrant_vol_F1 = df.iloc[0:7,0].values * u.mL
#Gran_point and titrant_vol_GP: Includes the points that ARE needed for the linear regression.
Gran_point = df.iloc[8:,2].values
titrant_vol_GP = df.iloc[8:,0].values * u.mL

# Use the stats package to do the linear regression.
slope, intercept, r_value, p_value, std_err = stats.linregress(titrant_vol_GP,Gran_point)
# Attach correct units to slope which would be the units of Gran Point (dimensionless) divided by units of titrant_vol (mL)
slope = slope *  1/u.mL
#Create a figure and plot the data and the line from the linear regression.
fig, ax = plt.subplots()
# plot the data not needed for linear regression in red circles and the data NEEDED for linear regression in blue circles.
ax.plot(titrant_vol_F1, F1, 'ro',)
ax.plot(titrant_vol_GP, Gran_point, 'bo',)

#plot the linear regression as a black line
ax.plot(titrant_vol_GP, (slope * titrant_vol_GP) + intercept, 'k-', )
#Define the x-axis
plt.xlabel('Titrant volume (mL)')
# Add axis labels using the column labels from the dataframe
ax.legend(['Non-Linear Regression Points', 'Linear regression Points', 'Linear Regression Line'])
ax.grid(True)

# Here I save the file to my local hard drive.
plt.savefig('/Users/kenrivero/Documents/EnvELab/Acid_Rain/Gran_Plot')
plt.show()

#Calculate equivalent volume using slope from stats package and ANC from ProCoDA. This calculation is done based on the lab manual that states that if F1 is plotted as a function of Vt the result is a straight line with slope = N_t/Vs. Using equation 41, the equivalent volume can be calculated by ANC divided by slope. Our units will be off because ANC is in units of eq/L and slope is in units of [1/mL].

calc_EV = ANC / slope
calc_EV.to(u.eq)*1000

#Part 3 **************************************************************

ANC_0 = 0.001826
ANC_5 = 0.000922 * u.eq/u.L
ANC_10 = 0.000276 * u.eq/u.L
ANC_15 = -0.00008 * u.eq/u.L


time_array = [0,5,10,15] * u.min
time_array = (time_array).to(u.s)
hydraulic_res_time = time_array/theta

fig, ax = plt.subplots()
plt.plot(res,ANC_out_closed,'r')
plt.plot(res,ANC_out_open,'b')
plt.scatter(0,ANC_0,color = 'green')
plt.scatter(hydraulic_res_time[1],ANC_5,color = 'green')
plt.scatter(hydraulic_res_time[2],ANC_10,color = 'green')
plt.scatter(hydraulic_res_time[3],ANC_15,color = 'green')
plt.xlabel('Hydraulic Residence Time')
plt.ylabel('ANC')
ax.legend(['ANC closed system','ANC open system','Measured ANC'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/Acid_Rain/ANC_openclosedmeasuredpoints')
plt.show()
