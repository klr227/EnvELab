Reactor Characteristics Prelaboratory
Ken Rivero-Rivera (klr227)
Catherine Johnson (caj92)

https://github.com/klr227/EnvELab/blob/master/Reactor_Characteristics/Pe_vs_t_star.png

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

import aide_design
from aide_design import physchem as pc
import aguaclara.core.constants as constants

#Question 1 ***************************************************************
#Retrieve vena contracta orifice ratio from AguaClara
VC_ORIFICE_RATIO = constants.VC_ORIFICE_RATIO
#The head_orifice function needs the parameters in the converted units. For some reason the function also doesn't accept the values with units so I took the magnitude of each after the converstion
flow_rate = (380 * u.mL/u.min).to(u.m**3/u.s).magnitude
baffle_diameter = (1 * u.mm).to(u.m).magnitude
orifice_HL = pc.head_orifice(baffle_diameter,VC_ORIFICE_RATIO,flow_rate)
print('The head loss through an orifice is',orifice_HL)

#Reactor baffles with 6 holes 1mm in diameter do not make a good design for this reactor  because there is a head loss of 8.353 meters which is large for our lab. (I think this is correct)

#Question 2 of prelaboratory exercise
#Created a linspace array from 0.01 to 3. As mentioned in the lab manual, the adjective dispersion equation is undefined when t_star = 0.
t_star = np.linspace(0.01,3,num = 50)
#Peclet numbers given in the prelab
Pe_1 = 1
Pe_10 = 10
Pe_100 = 100
#Calculate the exit age distribution
E_t_star_1 = epa.E_Advective_Dispersion(t_star,Pe_1)
E_t_star_10 = epa.E_Advective_Dispersion(t_star,Pe_10)
E_t_star_100 = epa.E_Advective_Dispersion(t_star,Pe_100)

#Plot exit age distribution vs t_star for the given Peclet numbers in the prelab.
fig, ax = plt.subplots()
plt.plot(t_star,E_t_star_1,'y')
plt.plot(t_star,E_t_star_10,'r')
plt.plot(t_star,E_t_star_100,'g')
plt.xlabel('t star')
plt.ylabel('E(t_star)')
ax.legend(['Pe = 1','Pe = 10', 'Pe = 100'])
plt.savefig('/Users/kenrivero/Documents/EnvELab/Reactor_Characteristics/Pe_vs_t_star')
plt.show()
