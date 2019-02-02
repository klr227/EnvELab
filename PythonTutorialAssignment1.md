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

  #if you want to see what is in the dataframe you can print it!
  print(df)
  # The column headers can be access by using the list command
  list(df)
  columns = df.columns
  print(columns)

  #I'll use the iloc method to get the x and y values.
  x = df.iloc[:,0].values * u.mg/u.L
  x
  y = df.iloc[:,1].values * u.V
  y

  # Use the stats package to do the linear regression.
  slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

  #Give intercept the units of the y values
  intercept = intercept * y.units
  intercept
  # Attach correct units to slope
  slope = slope * y.units/x.units
  slope
  #Create a figure and plot the data and the line from the linear regression.
  fig, ax = plt.subplots()
  # plot the data as red circles
  ax.plot(x, y, 'ro', )

  #plot the linear regression as a black line
  ax.plot(x, slope * x + intercept, 'k-', )

  # Add axis labels using the column labels from the dataframe
  ax.set(xlabel=list(df)[0])
  ax.set(ylabel=list(df)[1])
  ax.legend(['Measured', 'Linear regression'])
  ax.grid(True)

  # Here I save the file to my local hard drive.
  plt.savefig('/Users/kenrivero/Documents/EnvELab/linreggraph')
  plt.show()
```

  ![linreggraph]('https://github.com/klr227/EnvELab/blob/master/linreggraph.png')
Figure 1. This is a linear regression line for the data collected in the first lab.

  $$y= slope * x + intercept = (-0.01994 LV/mg) * x + 1.56527 V$$
