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

#Import datafile from GitHub
data_file_path =('https://raw.githubusercontent.com/klr227/EnvELab/master/Acid_Rain/lab3_0minsample.tsv')
df = pd.read_csv(data_file_path,delimiter='\t')

# The column headers can be access by using the list command
list(df)
columns = df.columns
print(columns)

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
plt.title('pH vs. Titrant volume (mL)')
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

#Title of Gran Plot
plt.title('Gran Plot')

# Add axis labels using the column labels from the dataframe
ax.legend(['Non-Linear Regression Points', 'Linear regression Points', 'Linear Regression Line'])
ax.grid(True)

# Here I save the file to my local hard drive.
plt.savefig('/Users/kenrivero/Documents/EnvELab/Acid_Rain/Gran_Plot')
plt.show()

#Calculate equivalent volume using slope from stats package and ANC from ProCoDA. This calculation is done based on the lab manual that states that if F1 is plotted as a function of Vt the result is a straight line with slope = N_t/Vs. Using equation 41, the equivalent volume can be calculated by ANC divided by slope. Our units will be off because ANC is in units of eq/L and slope is in units of [1/mL].

calc_EV = ANC / slope
calc_EV.to(u.eq)*1000
