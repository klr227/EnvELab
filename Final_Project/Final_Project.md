#Coagulant Adsorption Final Project  
###Ken Rivero - Rivera & Catherine Johnson  

##Introduction  
Water has a lot of important features. Because of these unique features, water is vital to the survival of life on Earth, it is found all over the world, and it is a universal solvent. It’s ability to be such a powerful solvent has both helped and hurt. It can carry essential substances to cells around the body of living creatures. But it can also be harmful, as the things it carries can be toxic pollutants, especially in an increasingly polluted world. Thus, there needs to be a way to remove pollutants from water - and one of those ways is through adsorption.

**Information about coagulant here**

##Objectives  
The purpose of this experiment was to evaluate the performance of coagulant as an adsorbent. The coagulant was tested at different concentrations and was assessed for two different substances - Red Dye #40 and Red Dye #3.  
##Procedures  
A schematic of the experiment can be seen in the figure below.

![](https://github.com/klr227/EnvELab/tree/master/Final_Project/Pictures)   
**insert figure here**  
The column was first filled with sand to ensure the correct volume of sand was obtained. Then, the sand and coagulant and 100mL of water was added to a beaker for mixing. While mixing, pH was monitored to ensure a pH above 6 was achieved. This was a necessary step because coagulant is largely ineffective in acidic conditions.
Coa
##Results and Discussion  
##Conclusions and Suggestions  

This lab yielded some interesting results, including the adsorption of Red Dye #40 with 1500 μL, which performed worse than the 745 and 945 μL runs. This could be attributed to a fluke, or a result of human error in measurement or calibrating the instruments. However, the team did note all of the pHs of sand and coagulant mixture before addition to the column, and 1500

```python
from aguaclara.core.units import unit_registry as u
import aguaclara.research.environmental_processes_analysis as epa
import aguaclara.core.physchem as pc
import aguaclara.core.utility as ut
import numpy as np
import matplotlib.pyplot as plt
import collections
import os
from pathlib import Path
import pandas as pd



def adsorption_data(C_column, dirpath):
    """This function extracts the data from folder containing tab delimited
    files of adsorption data. The file must be the original tab delimited file.

    Parameters
    ----------
    C_column : int
        index of the column that contains the dissolved oxygen concentration
        data.
    dirpath : string
        path to the directory containing aeration data you want to analyze
    Returns
    -------
    filepaths : string list
        all file paths in the directory sorted by flow rate
    time_data : numpy array list
        sorted list of numpy arrays containing the times with units of seconds
    Examples
    --------
    """
    #return the list of files in the directory
    metadata = pd.read_csv(dirpath + '/metadata.csv', delimiter=',')
    filenames = metadata['file name']
    #extract the flowrates from the filenames and apply units
    #sort airflows and filenames so that they are in ascending order of flow rates


    filepaths = [dirpath + '/' + i for i in filenames]
    #C_data is a list of numpy arrays. Thus each of the numpy data arrays can have different lengths to accommodate short and long experiments
    # cycle through all of the files and extract the column of data with oxygen concentrations and the times
    C_data=[epa.column_of_data(i,epa.notes(i).last_valid_index() + 1,C_column,-1,'mg/L') for i in filepaths]
    time_data=[(epa.column_of_time(i,epa.notes(i).last_valid_index() + 1,-1)).to(u.s) for i in filepaths]

    adsorption_collection = collections.namedtuple('adsorption_results','metadata filenames C_data time_data')
    adsorption_results = adsorption_collection(metadata, filenames, C_data, time_data)
    return adsorption_results


C_column = 1
dirpath1 = "https://raw.githubusercontent.com/klr227/EnvELab/master/Final_Project/Final_3_Data"

#flow rate
flow_rate = 100*u.mL/80*u.s


metadata1, filenames1, C_data1, time_data1 = adsorption_data(C_column,dirpath1)

C_0_red3 = 70 * u.mg/u.L
for i in range(np.size(filenames1)):
  C_data1[i]=C_data1[i]-C_data1[i][0]
mylegend1 = []
for i in range(np.size(filenames1)):
  plt.plot(time_data1[i], C_data1[i]/C_0_red3,'-');
  mylegend1.append(str(metadata1['Coagulant added (uL)'][i]) + ' µL')

plt.xlabel(r'Time (seconds)');
plt.ylabel(r'Red dye #3 $\left ( \frac{C}{C0} \right )$');
plt.legend(mylegend1);
plt.savefig('Final_Project/Pictures/Final_Dye_3.png');
plt.show()

dirpath2 = "https://raw.githubusercontent.com/klr227/EnvELab/master/Final_Project/Final_40_Data"
metadata2, filenames2, C_data2, time_data2 = adsorption_data(C_column,dirpath2)

C_0_red40 = 100 * u.mg/u.L
for i in range(np.size(filenames2)):
  C_data2[i]=C_data2[i]-C_data2[i][0]
mylegend2 = []
for i in range(np.size(filenames2)):
  plt.plot(time_data2[i], C_data2[i]/C_0_red40,'-');
  mylegend2.append(str(metadata2['Coagulant added (uL)'][i]) + ' µL')

plt.xlabel(r'Time (seconds)');
plt.ylabel(r'Red dye #40$\left ( \frac{C}{C0} \right )$');
plt.legend(mylegend2);
plt.savefig('Final_Project/Pictures/Red_Dye_40.png');
plt.show()
