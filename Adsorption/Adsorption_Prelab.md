#Adsorption Prelab
##Ken Rivero-Rivera
##Catherine Johnson
##Group 7

####A carbon column is packed with 15 cm of activated carbon and then used to remove 50 mg/L of red dye #40. The approach velocity is 1 mm/s, the porosity is 0.4, and the bulk density of the activated carbon is 0.5 g/cm3. How long will it take for the mass transfer zone to travel to the bottom of the carbon column?

It takes 33.35 hours for the mass transfer zone to travel to the bottom of the carbon column.
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

L_column = 15 * u.cm
RD_conc = 50 * u.mg/u.L
appr_vel = 1 * u.mm/u.s
porosity = 0.4
bulk_density = 0.5 * u.g/u.cm**3
q0 = 0.08 * u.g/u.g

t_mtz = ((L_column*porosity/appr_vel) + (L_column*q0*bulk_density/(appr_vel*RD_conc))).to(u.hr)
t_mtz
