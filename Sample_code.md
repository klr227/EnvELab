```python
from aguaclara.core.units import unit_registry as u
import aguaclara.research.environmental_processes_analysis as epa
import aguaclara.research.procoda_parser as pp
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def aeration_data(DO_column, dirpath):
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
  DO_data=[pp.column_of_data(i,0,DO_column,-1,'mg/L') for i in filepaths]
  time_data=[(pp.column_of_time(i,0,-1)).to(u.s) for i in filepaths]
  aeration_collection = collections.namedtuple('aeration_results','filepaths airflows DO_data time_data')
  aeration_results = aeration_collection(filepaths, airflows, DO_data, time_data)
  return aeration_results

  # The column of data containing the dissolved oxygen concentrations
  DO_column = 2
  dirpath = "/Users/kenrivero/Documents/EnvELab/Gas_Transfer/Aeration"
  filepaths, airflows, DO_data, time_data = aeration_data(DO_column,dirpath)
