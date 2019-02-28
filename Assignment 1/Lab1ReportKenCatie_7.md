1) What is the value of the extinction coefficient, e?
  1) A = ebc where A is absorbance, e is the extinction coefficient, b is the linear path of light distance (which is approximately 1cm in this lab), and c is concentration.
  Therefore, e = A/c for this lab and is equal to the slope of the linear regression line, which is 0.0421 L/mg.
1) Did you use interpolation or extrapolation to get the concentration of the unknown?
  1) we used interpolation to get the concentration of the unknown. The concentration had a voltage of -0.312 V, which corresponds to a
  absorbance of 0.6822 and a concentration of 15.131 mg/L based on our linear regression line.
1) What measurement controls the accuracy of the density measurement for the NaCl solution?
  1) The measurement that most controls the accuracy of the density measurement for the NaCl solution is the measurement of the NaCl mass on the scale.
  Any deviation from the mass required to make 100mL of 1M solution would result in an inaccurate density. Other measurements can effect this as well, but this
  is the step with the room for the most error.
1) What density did you expect (see prelab 2)?
  1) The expected density, as calculated in the prelab, was 1037.8 g/L. However, the measured density was 1037.56 g/L.
1) Approximately what should the accuracy be for the density measurement?
  1) The percent error in the measured density was 0.02%, meaning the accuracy of the density measurement was 99.98%.
1) Donít forget to write a brief paragraph on conclusions and on suggestions using Markdown.
  1) This lab exercise was helpful for learning laboratory techniques, and how to properly and accurately use lab equipment to take measurements.
  The objectives of this lab were to gain proficiency in lab practices, such as calibrating scales, digital pipetting, preparing standards and dilutions,
  and finding the concentrations of unknowns based on measurable properties (in this lab that property was absorbance). This lab succeeded in its objectives,
  but some of the equipment was hard to use. The photometer/pump configuration allowed both air and liquid through the tubing, and the air pockets in the
  tubing created some problems with measuring voltage and absorbance. While we were able to "knock" the air out of the photometer by slapping the apparatus,
  there is most likely a better way to pump air than a manual pump that would result in less error. Perhaps a mechanical pump would be better in the future,
  if possible.

``` python
  from aguaclara.core.units import unit_registry as u
  import numpy as np
  import matplotlib.pyplot as plt
  import pandas as pd
  from scipy import stats
  #The data file path is to my local hard drive
  data_file_path ="/Users/kenrivero/Documents/EnvELab/SampleData.tsv"

  #Now we create a pandas dataframe with the data in the file
  df = pd.read_csv(data_file_path,delimiter='\t')
  df = df.drop([7,8])
  #if you want to see what is in the dataframe you can print it!
  print(df)

  # The column headers can be access by using the list command
  list(df)
  columns = df.columns
  print(columns)

  #I'll use the iloc method to get the x and y values.
  x = df.iloc[:,0].values * u.mg/u.L
  y = df.iloc[:,1].values * u.V

  #Below are the dark and blank voltages used to calibrate the photometer
  Vdark = -1.30825 * u.V
  Vblank = 3.4841 * u.V

  #Calculate the absorbance
  Absorb = -np.log10((y-Vdark)/(Vblank-Vdark))
  # Use the stats package to do the linear regression.
  slope, intercept, r_value, p_value, std_err = stats.linregress(x,Absorb)


  # Attach correct units to slope. In this lab, the Extinction coefficient ,e, is equal to the slope calculated.
  slope = slope /x.units
  slope

  # Here are the calculations to determine the concentration of the unknown solution. First the absorbance is calculated and then the concentration is determined by using the slope and intercept.
  UnknownV = -0.312 * u.V
  UnknownA = -np.log10((UnknownV-Vdark)/(Vblank-Vdark))
  UnknownC = (UnknownA - intercept) / slope
  UnknownA

  #Create a figure and plot the data and the line from the linear regression.
  fig, ax = plt.subplots()
  # plot the data as red circles
  ax.plot(x, Absorb, 'ro', )
  #plot the linear regression as a black line
  ax.plot(x, slope*x + intercept, 'k-', )


  # Add axis labels using the column labels from the dataframe
  ax.set(xlabel=list(df)[0])
  ax.set_ylabel(ylabel= 'Absorbance')
  ax.legend(['Measured', 'Linear regression'])
  ax.grid(True)

  # Here I save the file to my local hard drive.
  plt.savefig('/Users/kenrivero/Documents/EnvELab/Lab1linreggraph')
  plt.show()
```
