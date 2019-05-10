#Coagulant Adsorption Final Project  
###Ken Rivero - Rivera & Catherine Johnson  

##Introduction  
Water has a lot of important features. Because of these unique features, water is vital to the survival of life on Earth, it is found all over the world, and it is a universal solvent. It’s ability to be such a powerful solvent has both helped and hurt. It can carry essential substances to cells around the body of living creatures. But it can also be harmful, as the things it carries can be toxic pollutants, especially in an increasingly polluted world. Thus, there needs to be a way to remove pollutants from water - and one of those ways is through adsorption.

**Information about coagulant here**

##Objectives  
The purpose of this experiment was to evaluate the performance of coagulant as an adsorbent. The coagulant was tested at different concentrations and was assessed for two different substances - Red Dye #40 and Red Dye #3.  
##Procedures  
A schematic of the experiment can be seen in the figure below.

![](https://github.com/klr227/EnvELab/blob/master/Final_Project/Pictures/AddingSandColumn.JPG)   
**insert figure here**  
The column was first filled with sand to ensure the correct volume of sand was obtained. Then, the sand and coagulant and 100mL of water was added to a beaker for mixing. While mixing, pH was monitored to ensure a pH above 6 was achieved. This was a necessary step because coagulant is largely ineffective in acidic conditions.
The aluminum in the coagulant is acidic, so a basic solution of 0.5M NaOH was added until a pH of at least 6 was reached.

After a few failed attempts at different flow rates, 1.25 mL/s was determined to work for both red dye #40 and red dye #3 enough to yield useful results. Both dyes were run through a sand-only column to obtain a baseline for the experiment. After, the team began testing coagulant. We started with 445 μL, which is equivalent to 22 moles of Aluminum per meter cubed. This concentration was suggested as ideal in “Enhanced Particle Capture through Aluminum Hydroxide Addition to Pores in Sand Media”. We then increased the volumes used to obtain data for comparisons. As you can see, we used both dyes for all volumes except the 1500 microliters of coagulant. This is because we quickly realized that red dye #3 was removed far more effectively by sand and coagulant than red dye #40 was. Where the red dye #40 experiments took minutes or less, the red dye #3 experiments took a few hours. Thus, we figured that the last experiment would not be worth the time.

##Results and Discussion  
##Conclusions and Suggestions  

This lab yielded some interesting results, including the adsorption of Red Dye #40 with 1500 μL, which performed worse than the 745 and 945 μL runs. This could be attributed to a fluke, or a result of human error in measurement or calibrating the instruments. However, the team did note all of the pHs of sand and coagulant mixture before addition to the column, and the 1500 μL run had the lowest pH at 6.1. The team did not add any more NaOH after this point because nearly 4.5 mL had already been added, the most of any other run, and flocs had begun to form in the mixture. It would be an interesting experiment in the future to see if the performance of coagulant as an adsorbent has correlation with the pH it starts out at.  

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
flow_rate = 100*u.mL/(80*u.s)
flow_rate

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
