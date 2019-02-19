Gas Transfer Pre-Laboratory
Ken Rivero-Rivera
Catherine Johnson
Group 7

Question 1: Calculate the mass of sodium sulfite needed to reduce all the dissolved oxygen in 750 mL of pure water in equilibrium with the atmosphere and at 22∘C.

Using Python, we calculated the mass of sodium sulfite needed to reduce all the dissolved oxygen is 52.54 mg.

Question 2: Describe your expectations for dissolved oxygen concentration as a function of time during an aeration experiment. Assume you have added enough sodium sulfite to consume all of the oxygen at the start of the experiment. What would the shape of the curve look like?

The dissolved oxygen concentration should increase with time during an aeration experiment. The shape of the curve should be concave up and increasing. Eventually, it will reach the saturation point and plateau.

Question 3: Why is k̂ v,l not zero when the gas flow rate is zero? How can oxygen transfer into the reactor even when no air is pumped into the diffuser?

k(v,l) is not zero when the gas flow rate is zero because there will still be gas transfer with the oxygen in the atmosphere. This is explained by Henry's law.

Question 4: Describe your expectations for k̂ v,l as a function of gas flow rate. Do you expect a straight line? Why?

k(v,l) is independent of the glass flow rate. k(v,l) is dependent on interface surface area, volume, oxygen diffusion coefficient, and the thickness of the laminar boundary layer.

Question 5: A dissolved oxygen probe was placed in a small vial in such a way that the vial was sealed. The water in the vial was sterile. Over a period of several hours the dissolved oxygen concentration gradually decreased to zero. Why?

Based on equation (111) in the lab manual, the probe will consume the oxygen existing in solution and convert it to H2O.

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

#Calculate the saturated oxygen concentration with the given temperature. Assume pressure is 1 atmosphere (equilibrium with the atmosphere)
temp = 22 * u.degC
pressure = 1 * u.atm
DO_conc = epa.O2_sat(pressure,temp)

#Calculate the mass of oxygen present in solution of 750 mL
volume = 750 * u.mL
DO_amount = DO_conc * volume

#Molar mass of O2 and sodium sulfite
MM_O2 = 32 * u.g/u.mol
MM_Na2SO3 = 126 * u.g/u.mol
#Molar ratio sodium sulfite : O2
molar_ratio = 2

#Calculate the mass of sodium sulfite needed to reduce all the dissolved oxygen
Na2SO3_amount = (DO_amount*molar_ratio*MM_Na2SO3)/MM_O2
#Convert to milligrams
Na2SO3_amount.to(u.mg)
