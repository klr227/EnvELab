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

data_file_path
df = pd.read_csv(data_file_path,delimiter='\t')

#if you want to see what is in the dataframe you can print it!
print(df)
# The column headers can be access by using the list command
list(df)
columns = df.columns
print(columns)
