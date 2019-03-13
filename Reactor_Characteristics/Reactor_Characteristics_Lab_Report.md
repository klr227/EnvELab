#Reactor Characteristics
##Group 7:
Ken Rivero Rivera - 1 hours
Catherine Johnson - 1 hours

##Introduction:
In water and wastewater treatment plants, chlorine is used to disinfect water in contact tanks. In order to maximize the effectiveness of the disinfecting process, the time the water spends in the tank must also be maximized. This allows for a better chance that bacteria and other pathogens in the water will be neutralized. This lab experiment is a tracer study, which are used to understand the characteristics of different reactor apparatuses.

##Objectives:
The purpose of this lab is to maximize the time that water spends in the reactor and to understand how the addition of barriers and baffles effects flow and dispersion in a reactor.

##Procedures:
A 1 foot long reactor was set up to receive water from a reservoir via peristaltic pump with a flow rate of 380 mg/L (100rpm). The peristaltic pump was also used to retrieve water in the tank near the outlet and send through the photometer. The first experiment was modelling a CMFR. This was an open tank with a stir bar. The remaining experiments were a series of baffle configurations that can be seen below.
###Configuration 1:
###Configuration 2:
Two baffles were placed such that the reactor consisted of three equal sized sub-sections, each 4 inches in length. The baffles each had a gap on one side of the tank, and the baffles were arranged so that the
###Configuration 3:
###Configuration 4:

##Results and Discussion:
###1. Use multivariable nonlinear regression to obtain the best fit between the experimental data and the two models by minimizing the sum of the squared errors. Use epa.Solver_AD_Pe and epa.Solver_CMFR_N. These functions will minimize the error by varying the values of average residence time, (mass of tracer/reactor volume), and either the number of CMFR in series or the Peclet number.

###2. Generate a plot showing the experimental data as points and the model results as thin lines for each of your experiments. Explain which model fits best and discuss those results based on your expectations.

###3. Compare the trends in the estimated values of N and Pe across your set of experiments. How did your chosen reactor modifications effect dispersion?

###4. Report the values of t⋆ at F = 0.1 for each of your experiments. Do they meet your expectations?

###5. Evaluate whether there is any evidence of “dead volumes” or “short circuiting” in your reactor.

###6. Make a recommendation for the design of a full scale chlorine contact tank. As part of your recommendation discuss the parameter you chose to vary as part of your experimentation and what the optimal value was determined to be

##Conclusion:

##Suggestions:
1. The first suggestion we have is straightforward - include in the instructions that the bigger tube in the peristaltic pump is for the influent not the effluent. A lot of people mixed those two around. Explicit direction in this step would save some time.
2. The photometer is used to determine the concentration of red dye at the outlet, however it is not necessarily measuring what exactly leaves the reactor. In the CMFR simulation, it would because that is well mixed. However, when using the baffles, the flow is not well mixed. Therefore, the concentration of what is leaving in real time is not necessarily what is next to the outlet. Perhaps there is a way to attach the photometer to the outlet in the future.

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

CMFR_path = "https://raw.githubusercontent.com/klr227/EnvELab/master/Reactor_Characteristics/red_dye_30_trial1.tsv"

epa.notes(CMFR_path)

CMFR_firstrow = epa.notes(CMFR_path).last_valid_index() + 1
CMFR_firstrow

CMFR_firstrow = CMFR_firstrow + 10

CMFR_time_data = (epa.column_of_time(CMFR_path,CMFR_firstrow,-1)).to(u.s)

CMFR_concentration_data = epa.column_of_data(CMFR_path,CMFR_firstrow,1,-1,'mg/L')

#I don't actually know the values that follow. I'm guessing.
#You should use real measured values!#
CMFR_V = 1.35*u.L
CMFR_Q = 380 * u.mL/u.min

#here we set estimates that we will use as starting values for the curve fitting
CMFR_theta_hydraulic = (CMFR_V/CMFR_Q).to(u.s)
CMFR_C_bar_guess = np.max(CMFR_concentration_data)

#The Solver_CMFR_N will return the initial tracer concentration,
#residence time, and number of reactors in series.
#This experiment was for a single reactor and so we expect N to be 1!
CMFR_CMFR = epa.Solver_CMFR_N(CMFR_time_data, CMFR_concentration_data, CMFR_theta_hydraulic, CMFR_C_bar_guess)
CMFR_CMFR
#use dot notation to get the 3 elements of the tuple that are in CMFR.

print('The model estimated mass of tracer injected was',ut.round_sf(CMFR_CMFR.C_bar*CMFR_V ,2) )
print('The model estimate of the number of reactors in series was', CMFR_CMFR.N)
print('The tracer residence time was',ut.round_sf(CMFR_CMFR.theta ,2))
print('The ratio of tracer to hydraulic residence time was',(CMFR_CMFR.theta/CMFR_theta_hydraulic).magnitude)

#create a model curve given the curve fit parameters.

CMFR_CMFR_model = CMFR_CMFR.C_bar * epa.E_CMFR_N(CMFR_time_data/CMFR_CMFR.theta,CMFR_CMFR.N)
plt.plot(CMFR_time_data.to(u.min), CMFR_concentration_data.to(u.mg/u.L),'r')
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
#photometer or perhaps there was some other anomoly that was causing the
#photometer to read a concentration that was higher than was actually present in
#the reactor. To correct for this I subtracted the initial concentration reading
#from all of the data. This was based on the assumption that the concentration
#measurement error persisted for the entire experiment.#

one_baffle_concentration_data = one_baffle_concentration_data - one_baffle_concentration_data[0]
one_baffle_V = 1.637*u.L
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

print('The model estimated mass of tracer injected was',ut.round_sf(one_baffle_AD.C_bar*one_baffle_V ,2) )
print('The model estimate of the Peclet number was', one_baffle_AD.Pe)
print('The tracer residence time was',ut.round_sf(one_baffle_AD.theta ,2))
print('The ratio of tracer to hydraulic residence time was',(one_baffle_AD.theta/one_baffle_theta_hydraulic).magnitude)

#Create the advection dispersion model curve based on the solver parameters
one_baffle_AD_model = (one_baffle_AD.C_bar*epa.E_Advective_Dispersion((one_baffle_time_data/one_baffle_AD.theta).to_base_units(), one_baffle_AD.Pe)).to(u.mg/u.L)

#Plot the data and the two model curves.
plt.plot(one_baffle_time_data.to(u.s), one_baffle_concentration_data.to(u.mg/u.L),'r')
plt.plot(one_baffle_time_data.to(u.s), one_baffle_CMFR_model,'b')
plt.plot(one_baffle_time_data.to(u.s), one_baffle_AD_model,'g')
plt.xlabel(r'$time (min)$')
plt.ylabel(r'Concentration $\left ( \frac{mg}{L} \right )$')
plt.legend(['Measured dye','CMFR Model', 'AD Model'])
plt.savefig('Reactor_Characteristics/Dispersion.png', bbox_inches = 'tight')
plt.show()
