#Reactor Characteristics  
##Group 7:  
Ken Rivero Rivera - 14 hours  
Catherine Johnson - 6 hours  
***
##Introduction:  
In water and wastewater treatment plants, chlorine is used to disinfect water in contact tanks. In order to maximize the effectiveness of the disinfecting process, the time the water spends in the tank must also be maximized. This allows for a better chance that bacteria and other pathogens in the water will be neutralized. This lab experiment is a tracer study, which are used to understand the characteristics of different reactor apparatuses.  
***
##Objectives:  
The purpose of this lab is to maximize the time that water spends in the reactor and to understand how the addition of barriers and baffles effects flow and dispersion in a reactor.  
##Procedures:  
A 1 foot long reactor was set up to receive water from a reservoir via peristaltic pump with a flow rate of 380 mg/L (100rpm). The peristaltic pump was also used to retrieve water in the tank near the outlet and send through the photometer. The first experiment was modeling a CMFR. This was an open tank with a stir bar. The remaining experiments were a series of baffle configurations that can be seen below.  
####Configuration 1:  
One baffle with small holes (5mm diameter) was placed in the middle of the reactor. These holes were arranged such that there were 6 holes horizontally per row, and 4 holes vertically per column. Only the bottom two rows were submerged in the water in this experiment, meaning the flow was through 12 holes.   

![](https://github.com/klr227/EnvELab/blob/master/Reactor_Characteristics/Images/IMG_1047.JPG
)  
####Configuration 2:  
Two baffles were placed such that the reactor consisted of three equal sized sub-sections, each 4 inches in length. The baffles each had a gap on one side of the tank, and the baffles were arranged so that the gaps were on opposite sides.  

![](https://github.com/klr227/EnvELab/blob/master/Reactor_Characteristics/Images/IMG_5439.JPG
)  
####Configuration 3:  
Three baffles were placed such that the reactor consisted of four equal sized sub-sections, each 3 inches in length. The baffles each had a gap on one side of the tank, and the baffles were arranged so that the gaps were on opposite sides from the adjacent baffle gaps.  

![](https://github.com/klr227/EnvELab/blob/master/Reactor_Characteristics/Images/IMG_5448.JPG)  

####Configuration 4:  
Four baffles were placed such that the reactor consisted of five equal sized sub-sections, each 6 centimeters in length. The baffles each had a gap on one side of the tank, and the baffles were arranged so that the gaps were on opposite sides from the adjacent baffle gaps.  

![](https://github.com/klr227/EnvELab/blob/master/Reactor_Characteristics/Images/IMG_5452.JPG)  
***
##Results and Discussion:  
####1. Use multivariable nonlinear regression to obtain the best fit between the experimental data and the two models by minimizing the sum of the squared errors. Use epa.Solver_AD_Pe and epa.Solver_CMFR_N. These functions will minimize the error by varying the values of average residence time, (mass of tracer/reactor volume), and either the number of CMFR in series or the Peclet number.

####2. Generate a plot showing the experimental data as points and the model results as thin lines for each of your experiments. Explain which model fits best and discuss those results based on your expectations.  

![](https://github.com/klr227/EnvELab/blob/master/Reactor_Characteristics/One_baffle_graph_model.png)  
**Graph of CMFR and AD models compared with data from experiment with one baffle with 5 mm diameter holes**  

![](https://github.com/klr227/EnvELab/blob/master/Reactor_Characteristics/Trial_2_graph_model.png
)  
**Graph of CMFR and AD models compared with data from experiment with two baffles**
![](https://github.com/klr227/EnvELab/blob/master/Reactor_Characteristics/Trial_3_graph_model.png
)  
**Graph of CMFR and AD models compared with data from experiment with three baffles**  
![](https://github.com/klr227/EnvELab/blob/master/Reactor_Characteristics/Trial_4_Graph_model.png
)   
**Graph of CMFR and AD models compared with data from experiment with four baffles**   

It is difficult to say which model is better because they both fit the data points better in separate regions.  


####3. Compare the trends in the estimated values of N and Pe across your set of experiments. How did your chosen reactor modifications effect dispersion?  

CMFR: Pe = n/a , N= 1.00  
One Baffle w/ 5mm holes: Pe = 1.741 , N = 1.801  
Two Baffles: Pe = 2.288, N = 2.17  
Three Baffles: Pe = 2.464, N = 2.210  
Four Baffles: Pe = 4.762; N = 3.37  


####4. Report the values of t⋆ at F = 0.1 for each of your experiments. Do they meet your expectations?  
One Baffle w/ 5mm holes: t⋆ = 0.4396  
Two Baffles: t⋆ = 0.4637  
Three Baffles: t⋆ = 0.4857  
Four Baffles: t⋆ = 0.5463  

Yes they meet our expectations.  
####5. Evaluate whether there is any evidence of “dead volumes” or “short circuiting” in your reactor.  
The difference of the hydraulic residence time calculated with the flow rate and volume and the hydraulic residence time from the AD model is evidence of "dead volumes" in our reactors. The hydraulic residence time of the AD model is lower than the other one which can be attributed to dead volumes located in the reactor.  
####6. Make a recommendation for the design of a full scale chlorine contact tank. As part of your recommendation discuss the parameter you chose to vary as part of your experimentation and what the optimal value was determined to be  

We would recommend a chlorine contact tank with several baffles between the influent and the effluent. We found that there is greater dispersion with an increase in baffles. The higher number of baffles create a longer flow path thereby increasing the residence time for which the chlorine is in contact with contaminants.
***
##Conclusion:  
We tested the residence of a reactor with four different configurations. We determined that the higher number of baffles increases residence time in which chlorine is in contact with contaminant. Also by decreasing the width of the flow path and increasing it's length, we decrease the amount of dead volume in a reactor.
##Suggestions:  
1. The first suggestion we have is straightforward - include in the instructions that the bigger tube in the peristaltic pump is for the influent not the effluent. A lot of people mixed those two around. Explicit direction in this step would save some time.  
2. The photometer is used to determine the concentration of red dye at the outlet, however it is not necessarily measuring what exactly leaves the reactor. In the CMFR simulation, it would because that is well mixed. However, when using the baffles, the flow is not well mixed. Therefore, the concentration of what is leaving in real time is not necessarily what is next to the outlet. Perhaps there is a way to attach the photometer to the outlet in the future.  
***
##Appendix:
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


#Concentration of tracer and volume of tracer added
C_tr = 30 * u.mg/u.L
V_tr = 300 * u.uL
#Density of water
density_H2O = 997 *u.g/u.L

#CMFR
CMFR_path = "https://raw.githubusercontent.com/klr227/EnvELab/master/Reactor_Characteristics/red_dye_30_trial1.tsv"
CMFR_firstrow = epa.notes(CMFR_path).last_valid_index() + 1
CMFR_firstrow = CMFR_firstrow + 10
CMFR_time_data = (epa.column_of_time(CMFR_path,CMFR_firstrow,-1)).to(u.s)
CMFR_concentration_data = epa.column_of_data(CMFR_path,CMFR_firstrow,1,-1,'mg/L')
CMFR_concentration_data

#Used measured mass of water in the tank in CMFR experiment
CMFR_mass = 1646 *u.g
CMFR_V = CMFR_mass/density_H2O
CMFR_Q = 380 * u.mL/u.min

#here we set estimates that we will use as starting values for the curve fitting
CMFR_theta_hydraulic = (CMFR_V/CMFR_Q).to(u.s)
CMFR_C_bar_guess = np.max(CMFR_concentration_data)

#The Solver_CMFR_N will return the initial tracer concentration,
#residence time, and number of reactors in series.
#This experiment was for a single reactor and so we expect N to be 1!
CMFR_CMFR = epa.Solver_CMFR_N(CMFR_time_data, CMFR_concentration_data, CMFR_theta_hydraulic, CMFR_C_bar_guess)
#use dot notation to get the 3 elements of the tuple that are in CMFR.

print('The model estimated mass of tracer injected was',CMFR_CMFR.C_bar*CMFR_V)
print('The model estimate of the number of reactors in series was', CMFR_CMFR.N)
print('The tracer residence time was',CMFR_CMFR.theta)
print('The ratio of tracer to hydraulic residence time was',(CMFR_CMFR.theta/CMFR_theta_hydraulic).magnitude)

#create a model curve given the curve fit parameters.

CMFR_CMFR_model = CMFR_CMFR.C_bar * epa.E_CMFR_N(CMFR_time_data/CMFR_CMFR.theta,CMFR_CMFR.N)
plt.plot(CMFR_time_data.to(u.min), CMFR_concentration_data,'ro')
plt.plot(CMFR_time_data.to(u.min), CMFR_CMFR_model,'b')

plt.xlabel(r'$time (min)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model'])
plt.savefig('Reactor_Characteristics/CMFR.png', bbox_inches = 'tight')
plt.show()



#Load a data file for a reactor with baffles.
one_baffle_path = 'https://raw.githubusercontent.com/klr227/EnvELab/master/Reactor_Characteristics/baffle_test_1.xls'
one_baffle_firstrow = epa.notes(one_baffle_path).last_valid_index() + 1
one_baffle_time_data = (epa.column_of_time(one_baffle_path,one_baffle_firstrow,-1)).to(u.s)
one_baffle_concentration_data = epa.column_of_data(one_baffle_path,one_baffle_firstrow,1,-1,'mg/L')

#I noticed that the initial concentration measured by the photometer wasn't
#zero. This suggests that there may have been a small air bubble in the
#photometer or perhaps there was some other anomaly that was causing the
#photometer to read a concentration that was higher than was actually present in
#the reactor. To correct for this I subtracted the initial concentration reading
#from all of the data. This was based on the assumption that the concentration
#measurement error persisted for the entire experiment.#

one_baffle_concentration_data = one_baffle_concentration_data - one_baffle_concentration_data[0]
one_baffle_mass = 137 *u.g
one_baffle_V = one_baffle_mass/density_H2O
one_baffle_Q = 380 * u.mL/u.min
one_baffle_theta_hydraulic = (one_baffle_V/one_baffle_Q).to(u.s)
one_baffle_C_bar_guess = np.max(one_baffle_concentration_data)/2
#use solver to get the CMFR parameters
one_baffle_CMFR = epa.Solver_CMFR_N(one_baffle_time_data, one_baffle_concentration_data, one_baffle_theta_hydraulic, one_baffle_C_bar_guess)
one_baffle_CMFR.C_bar
one_baffle_CMFR.N
one_baffle_CMFR.theta.to(u.s)

#Create the CMFR model curve based on the scipy.optimize curve_fit
#parameters. We do this with dimensions so that we can plot both models and
#the data on the same graph. If we did this in dimensionless space it wouldn't
#be possible to plot everything on the same plot because the values used to
#create dimensionless time and dimensionless concentration are different for
#the two models.
one_baffle_CMFR_model = (one_baffle_CMFR.C_bar*epa.E_CMFR_N(one_baffle_time_data/one_baffle_CMFR.theta, one_baffle_CMFR.N)).to(u.mg/u.L)

#use solver to get the advection dispersion parameters
one_baffle_AD = epa.Solver_AD_Pe(one_baffle_time_data, one_baffle_concentration_data, one_baffle_theta_hydraulic, one_baffle_C_bar_guess)
one_baffle_AD.C_bar
one_baffle_AD.Pe
one_baffle_AD.theta

#Create the advection dispersion model curve based on the solver parameters
one_baffle_AD_model = (one_baffle_AD.C_bar*epa.E_Advective_Dispersion((one_baffle_time_data/one_baffle_AD.theta).to(u.dimensionless), one_baffle_AD.Pe)).to(u.mg/u.L)

#Plot the data and the two model curves.
plt.plot(one_baffle_time_data.to(u.s), one_baffle_concentration_data.to(u.mg/u.L),'ro')
plt.plot(one_baffle_time_data.to(u.s), one_baffle_CMFR_model,'b')
plt.plot(one_baffle_time_data.to(u.s), one_baffle_AD_model,'g')
plt.xlabel(r'$time (seconds)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model', 'AD Model'])
plt.savefig('Reactor_Characteristics/One_baffle_graph_model.png', bbox_inches = 'tight')
plt.show()

E1 = ((one_baffle_concentration_data*one_baffle_V)/(C_tr*V_tr)).to(u.dimensionless)
t1_star = one_baffle_time_data/one_baffle_AD.theta
F1 = []
for i in range(t1_star.size):
  integration = np.trapz(E1[0:i],t1_star[0:i])
  F1.append(integration)

F1max = F1[t1_star.size-1]
F1_10 = .1*F1max
index_1 = 0
#while value < F1_10:

for i in range(t1_star.size):
  if F1[i] < F1_10:
    index_1 = i + 1
  else:
    break

t1_star_10 = t1_star[index_1]


fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('t⋆')
ax1.set_ylabel('E', color=color)
ax1.plot(t1_star, E1, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('F', color=color)  # we already handled the x-label with ax1
ax2.plot(t1_star, F1, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.annotate('t⋆ at F = 0.1' , xy=(t1_star_10, F1[index_1]), xytext=(1, 7.5),
            arrowprops=dict(facecolor='black', shrink=1))
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('Reactor_Characteristics/E_F_Trial1_Graph.png', bbox_inches = 'tight')
plt.show()

#Load a data file for a reactor with baffles trial 2 (1 baffle with gap).

trial2_path = 'https://raw.githubusercontent.com/klr227/EnvELab/master/Reactor_Characteristics/baffle_test_2.xls'
trial2_firstrow = epa.notes(trial2_path).last_valid_index() + 1
trial2_time_data = (epa.column_of_time(trial2_path,trial2_firstrow,-1)).to(u.s)
trial2_concentration_data = epa.column_of_data(trial2_path,trial2_firstrow,1,-1,'mg/L')

#I noticed that the initial concentration measured by the photometer wasn't
#zero. This suggests that there may have been a small air bubble in the
#photometer or perhaps there was some other anomaly that was causing the
#photometer to read a concentration that was higher than was actually present in
#the reactor. To correct for this I subtracted the initial concentration reading
#from all of the data. This was based on the assumption that the concentration
#measurement error persisted for the entire experiment.#

trial2_concentration_data = trial2_concentration_data - trial2_concentration_data[0]
trial2_mass = 1537 *u.g
trial2_V = trial2_mass/density_H2O
trial2_Q = 380 * u.mL/u.min
trial2_theta_hydraulic = (trial2_V/trial2_Q).to(u.s)
trial2_C_bar_guess = np.max(trial2_concentration_data)/2
#use solver to get the CMFR parameters
trial2_CMFR = epa.Solver_CMFR_N(trial2_time_data, trial2_concentration_data, trial2_theta_hydraulic, trial2_C_bar_guess)
trial2_CMFR.C_bar
trial2_CMFR.N
trial2_CMFR.theta.to(u.s)

#Create the CMFR model curve based on the scipy.optimize curve_fit
#parameters. We do this with dimensions so that we can plot both models and
#the data on the same graph. If we did this in dimensionless space it wouldn't
#be possible to plot everything on the same plot because the values used to
#create dimensionless time and dimensionless concentration are different for
#the two models.
trial2_CMFR_model = (trial2_CMFR.C_bar*epa.E_CMFR_N(trial2_time_data/trial2_CMFR.theta, trial2_CMFR.N)).to(u.mg/u.L)

#use solver to get the advection dispersion parameters
trial2_AD = epa.Solver_AD_Pe(trial2_time_data, trial2_concentration_data, trial2_theta_hydraulic, trial2_C_bar_guess)
trial2_AD.C_bar
trial2_AD.Pe
trial2_AD.theta

#Create the advection dispersion model curve based on the solver parameters
trial2_AD_model = (trial2_AD.C_bar*epa.E_Advective_Dispersion((trial2_time_data/trial2_AD.theta).to_base_units(), trial2_AD.Pe)).to(u.mg/u.L)

#Plot the data and the two model curves.
plt.plot(trial2_time_data.to(u.s), trial2_concentration_data.to(u.mg/u.L),'ro')
plt.plot(trial2_time_data.to(u.s), trial2_CMFR_model,'b')
plt.plot(trial2_time_data.to(u.s), trial2_AD_model,'g')
plt.xlabel(r'$time (seconds)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model', 'AD Model'])
plt.savefig('Reactor_Characteristics/Trial_2_graph_model.png', bbox_inches = 'tight')
plt.show()

E2 = ((trial2_concentration_data*trial2_V)/(C_tr*V_tr)).to(u.dimensionless)
t2_star = trial2_time_data/trial2_AD.theta

F2 = []
for i in range(t2_star.size):
  integration = np.trapz(E2[0:i],t2_star[0:i])
  F2.append(integration)

F2max = F2[t2_star.size-1]
F2_10 = .1*F2max
index_2 = 0
#while value < F1_10:

for i in range(t2_star.size):
  if F2[i] < F2_10:
    index_2 = i + 1
  else:
    break

t2_star_10 = t2_star[index_2]

fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('t⋆')
ax1.set_ylabel('E', color=color)
ax1.plot(t2_star, E2, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('F', color=color)  # we already handled the x-label with ax1
ax2.plot(t2_star, F2, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.annotate('t⋆ at F = 0.1' , xy=(t2_star_10, F2[index_2]), xytext=(1, 7.5),
            arrowprops=dict(facecolor='black', shrink=1))
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('Reactor_Characteristics/E_F_Trial2_Graph.png', bbox_inches = 'tight')
plt.show()
#trial 3 - 3 baffles
trial3_path = 'https://raw.githubusercontent.com/klr227/EnvELab/master/Reactor_Characteristics/baffle_test_3.xls'
trial3_firstrow = epa.notes(trial3_path).last_valid_index() + 1
trial3_time_data = (epa.column_of_time(trial3_path,trial3_firstrow,-1)).to(u.s)
trial3_concentration_data = epa.column_of_data(trial3_path,trial3_firstrow,1,-1,'mg/L')

#I noticed that the initial concentration measured by the photometer wasn't
#zero. This suggests that th`ere may have been a small air bubble in the
#photometer or perhaps there was some other anomaly that was causing the
#photometer to read a concentration that was higher than was actually present in
#the reactor. To correct for this I subtracted the initial concentration reading
#from all of the data. This was based on the assumption that the concentration
#measurement error persisted for the entire experiment.#

trial3_concentration_data = trial3_concentration_data - trial3_concentration_data[0]
trial3_mass = 1557 *u.g
trial3_V = trial3_mass/density_H2O
trial3_Q = 380 * u.mL/u.min
trial3_theta_hydraulic = (trial3_V/trial3_Q).to(u.s)
trial3_C_bar_guess = np.max(trial3_concentration_data)/4
trial3_C_bar_guess = np.average(trial3_concentration_data)* u.mg/u.L
#use solver to get the CMFR parameters
trial3_CMFR = epa.Solver_CMFR_N(trial3_time_data, trial3_concentration_data, trial3_theta_hydraulic, trial3_C_bar_guess)
trial3_CMFR.C_bar
trial3_CMFR.N
trial3_CMFR.theta.to(u.s)

#Create the CMFR model curve based on the scipy.optimize curve_fit
#parameters. We do this with dimensions so that we can plot both models and
#the data on the same graph. If we did this in dimensionless space it wouldn't
#be possible to plot everything on the same plot because the values used to
#create dimensionless time and dimensionless concentration are different for
#the two models.
trial3_CMFR_model = (trial3_CMFR.C_bar*epa.E_CMFR_N(trial3_time_data/trial3_CMFR.theta, trial3_CMFR.N)).to(u.mg/u.L)

#use solver to get the advection dispersion parameters
trial3_AD = epa.Solver_AD_Pe(trial3_time_data, trial3_concentration_data, 200*u.s, trial3_C_bar_guess)
trial3_AD.C_bar
trial3_AD.Pe
trial3_AD.theta


#Create the advection dispersion model curve based on the solver parameters
trial3_AD_model = (trial3_AD.C_bar*epa.E_Advective_Dispersion((trial3_time_data/trial3_AD.theta).to_base_units(), trial3_AD.Pe)).to(u.mg/u.L)

#Plot the data and the two model curves.
plt.plot(trial3_time_data.to(u.s), trial3_concentration_data.to(u.mg/u.L),'ro')
plt.plot(trial3_time_data.to(u.s), trial3_CMFR_model,'b')
plt.plot(trial3_time_data.to(u.s), trial3_AD_model,'g')
plt.xlabel(r'$time (seconds)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model', 'AD Model'])
plt.savefig('Reactor_Characteristics/Trial_3_graph_model.png', bbox_inches = 'tight')
plt.show()

E3 = (trial3_concentration_data*trial3_V/(C_tr*V_tr)).to(u.dimensionless)
t3_star = trial3_time_data/trial3_AD.theta
F3 = []
for i in range(t3_star.size):
  integration = np.trapz(E3[0:i],t3_star[0:i])
  F3.append(integration)

F3max = F3[t3_star.size-1]
F3_10 = .1*F3max
index_3 = 0
#while value < F1_10:

for i in range(t3_star.size):
  if F3[i] < F3_10:
    index_3 = i + 1
  else:
    break

t3_star_10 = t3_star[index_3]

fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('t⋆')
ax1.set_ylabel('E', color=color)
ax1.plot(t3_star, E3, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('F', color=color)  # we already handled the x-label with ax1
ax2.plot(t3_star, F3, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.annotate('t⋆ at F = 0.1' , xy=(t3_star_10, F3[index_3]), xytext=(1, 7.5),
            arrowprops=dict(facecolor='black', shrink=1))
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('Reactor_Characteristics/E_F_Trial3_Graph.png', bbox_inches = 'tight')
plt.show()
##Load a data file for a reactor with baffles trial 4 (4 baffles with gap).
trial4_path = 'https://raw.githubusercontent.com/klr227/EnvELab/master/Reactor_Characteristics/baffle_test_4.xls'
trial4_firstrow = epa.notes(trial4_path).last_valid_index() + 1
trial4_time_data = (epa.column_of_time(trial4_path,trial4_firstrow,-1)).to(u.s)
trial4_concentration_data = epa.column_of_data(trial4_path,trial4_firstrow,1,-1,'mg/L')

#I noticed that the initial concentration measured by the photometer wasn't
#zero. This suggests that there may have been a small air bubble in the
#photometer or perhaps there was some other anomaly that was causing the
#photometer to read a concentration that was higher than was actually present in
#the reactor. To correct for this I subtracted the initial concentration reading
#from all of the data. This was based on the assumption that the concentration
#measurement error persisted for the entire experiment.#

trial4_concentration_data = trial4_concentration_data - trial4_concentration_data[0]
trial4_mass = 1598 *u.g
trial4_V = trial4_mass/density_H2O
trial4_Q = 380 * u.mL/u.min
trial4_theta_hydraulic = (trial4_V/trial4_Q).to(u.s)
trial4_C_bar_guess = np.max(trial4_concentration_data)/5
#use solver to get the CMFR parameters
trial4_CMFR = epa.Solver_CMFR_N(trial4_time_data, trial4_concentration_data, trial4_theta_hydraulic, trial4_C_bar_guess)
trial4_CMFR.C_bar
trial4_CMFR.N
trial4_CMFR.theta.to(u.s)

#Create the CMFR model curve based on the scipy.optimize curve_fit
#parameters. We do this with dimensions so that we can plot both models and
#the data on the same graph. If we did this in dimensionless space it wouldn't
#be possible to plot everything on the same plot because the values used to
#create dimensionless time and dimensionless concentration are different for
#the two models.
trial4_CMFR_model = (trial4_CMFR.C_bar*epa.E_CMFR_N(trial4_time_data/trial4_CMFR.theta, trial4_CMFR.N)).to(u.mg/u.L)

#use solver to get the advection dispersion parameters
trial4_AD = epa.Solver_AD_Pe(trial4_time_data, trial4_concentration_data, trial4_theta_hydraulic, trial4_C_bar_guess)
trial4_AD.C_bar
trial4_AD.Pe
trial4_AD.theta

#Create the advection dispersion model curve based on the solver parameters
trial4_AD_model = (trial4_AD.C_bar*epa.E_Advective_Dispersion((trial4_time_data/trial4_AD.theta).to_base_units(), trial4_AD.Pe)).to(u.mg/u.L)

#Plot the data and the two model curves.
plt.plot(trial4_time_data.to(u.s), trial4_concentration_data.to(u.mg/u.L),'ro')
plt.plot(trial4_time_data.to(u.s), trial4_CMFR_model,'b')
plt.plot(trial4_time_data.to(u.s), trial4_AD_model,'g')
plt.xlabel(r'$time (seconds)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model', 'AD Model'])
plt.savefig('Reactor_Characteristics/Trial_4_Graph_model.png', bbox_inches = 'tight')
plt.show()

E4 = (trial4_concentration_data*trial4_V/(C_tr*V_tr)).to(u.dimensionless)
t4_star = trial4_time_data/trial4_AD.theta
F4 = []
for i in range(t4_star.size):
  timedata = t4_star[0:i]
  integration = np.trapz(E4[0:i],timedata)
  F4.append(integration)

F4max = F4[t4_star.size-1]
F4_10 = .1*F4max
index_4= 0
#while value < F1_10:

for i in range(t4_star.size):
  if F4[i] < F4_10:
    index_4 = i + 1
  else:
    break

t4_star_10 = t4_star[index_4]

fig, ax1 = plt.subplots()
color = 'tab:red'
ax1.set_xlabel('t⋆')
ax1.set_ylabel('E', color=color)
ax1.plot(t4_star, E4, color=color)
ax1.tick_params(axis='y', labelcolor=color)

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

color = 'tab:blue'
ax2.set_ylabel('F', color=color)  # we already handled the x-label with ax1
ax2.plot(t4_star, F4, color=color)
ax2.tick_params(axis='y', labelcolor=color)
ax2.annotate('t⋆ at F = 0.1' , xy=(t4_star_10, F4[index_4]), xytext=(1, 7.5),
            arrowprops=dict(facecolor='black', shrink=1))
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('Reactor_Characteristics/E_F_Trial4_Graph.png', bbox_inches = 'tight')
plt.show()

print('THE FOLLOWING INFORMATION IS FOR OUR EXPERIMENT WITH 1 BAFFLE WITH SEVERAL HOLES EQUALLY SPACED')
print('The model estimated mass of tracer injected was',one_baffle_AD.C_bar*one_baffle_V )
print('The model estimate of the Peclet number was', one_baffle_AD.Pe)
print('The tracer residence time was',one_baffle_AD.theta)
print('The ratio of tracer to hydraulic residence time was',(one_baffle_AD.theta/one_baffle_theta_hydraulic).magnitude)
print('The values of t⋆ at F = 0.1 for this experiment is',t1_star_10)
print('')
print('THE FOLLOWING INFORMATION IS FOR OUR EXPERIMENT WITH 1 BAFFLES WITH A GAP ON ONE END')
print('The model estimated mass of tracer injected was',trial2_AD.C_bar*trial2_V)
print('The model estimate of the Peclet number was', trial2_AD.Pe)
print('The tracer residence time was',trial2_AD.theta)
print('The ratio of tracer to hydraulic residence time was',(trial2_AD.theta/trial2_theta_hydraulic).magnitude)
print('The values of t⋆ at F = 0.1 for this experiment is',t2_star_10)
print('')
print('THE FOLLOWING INFORMATION IS FOR OUR EXPERIMENT WITH 3 BAFFLES')
print('The model estimated mass of tracer injected was',trial3_AD.C_bar*trial3_V)
print('The model estimate of the Peclet number was', trial3_AD.Pe)
print('The tracer residence time was',trial3_AD.theta)
print('The ratio of tracer to hydraulic residence time was',(trial3_AD.theta/trial3_theta_hydraulic).magnitude)
print('The values of t⋆ at F = 0.1 for this experiment is',t3_star_10)
print('')
print('THE FOLLOWING INFORMATION IS FOR OUR EXPERIMENT WITH 4 BAFFLES')
print('The model estimated mass of tracer injected was',trial4_AD.C_bar*trial4_V)
print('The model estimate of the Peclet number was', trial4_AD.Pe)
print('The tracer residence time was',trial4_AD.theta)
print('The ratio of tracer to hydraulic residence time was',(trial4_AD.theta/trial4_theta_hydraulic).magnitude)
print('The values of t⋆ at F = 0.1 for this experiment is',t4_star_10)
