# E1: Data processing and analysis for in-situ NOx data at Loobos flux tower and Wekerom sites air monitoring site, this script solves 
# the first section of RQ 1 of the thesis, which is to compare the NOx data measurements at both sites.
# Author: Zhiyu Wu, Date: 05/02/2026
# ------ Step 1: Import necessary libraries------
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# ------ Step 2: Load the datasets and select relevant columns, create output folder ------
import os
# Read the csv file containing dataset, separated by spaces, skipping first three lines
data = pd.read_csv('Data/VeluweInSitu/NO2&meteodata.csv' , sep='\s+', skiprows=3)
# Merge Time and Timestamp columns into a single datetime column, convert to datatime format, delete timestamp column
data['Time'] = pd.to_datetime(data['Time'] + ' ' + data['Timestamp'], format='%Y-%m-%d %H:%M:%S')
data.drop(columns=['Timestamp'], inplace=True)
# Select columns containing time and NOx data, which are 'NO2_LGR_30m' and 'NO2_Wek_30m'
data = data[['Time', 'NO2_LGR_30m', 'NO2_Wek_30m']]
# Rename the columns for better readability
data.rename(columns={'NO2_LGR_30m': 'NO2_Loobos', 'NO2_Wek_30m': 'NO2_Wekerom'}, inplace=True)
# Create output folder if it does not exist
output_folder = 'Output/RQ1'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)


# ------ Step 3: Data cleaning ------
# Filter out all the NA data represented by -9999.999 and values less than 0, as NO2 concentration cannot be negative
data = data[(data['NO2_Loobos'] != -9999.999) & (data['NO2_Wekerom'] != -9999.999)]
data = data[(data['NO2_Loobos'] >= 0) & (data['NO2_Wekerom'] >= 0)]

# ------ Step 4: Data analysis and visualization ------
# Calculate basic statistics for both sites
loobos_stats = data['NO2_Loobos'].describe()
wekerom_stats = data['NO2_Wekerom'].describe()
print("Loobos NO2 Statistics:\n", loobos_stats)
print("Wekerom NO2 Statistics:\n", wekerom_stats)

# Plot a box plot of NO2 concentrations for both sites, save the plot as a png file
plt.figure(figsize=(10, 6))
plt.boxplot([data['NO2_Loobos'], data['NO2_Wekerom']], labels=['Loobos', 'Wekerom'])
plt.ylabel('$NO_2$ Concentration ($\mu g/m^3$)')
plt.title('Box Plot of $NO_2$ Concentrations at Loobos and Wekerom Sites')
plt.grid(True)
plt.savefig('Output/RQ1/NO2_BoxPlot.png')
plt.show()

# Calculate Pearson and Spearman correlation coefficient, plot a scatter plot to visualize the correlation
pearson_correlation = data['NO2_Loobos'].corr(data['NO2_Wekerom'])
print("Pearson Correlation Coefficient between Loobos and Wekerom NO2 concentrations: ", pearson_correlation)
spearman_correlation = data['NO2_Loobos'].corr(data['NO2_Wekerom'], method='spearman')
print("Spearman Correlation Coefficient between Loobos and Wekerom NO2 concentrations: ", spearman_correlation)
plt.figure(figsize=(10, 6))
plt.scatter(data['NO2_Loobos'], data['NO2_Wekerom'], alpha=0.5)
plt.xlabel('$NO_2$ Concentration at Loobos ($\mu g/m^3$)')
plt.ylabel('$NO_2$ Concentration at Wekerom ($\mu g/m^3$)')
plt.title('Scatter Plot of $NO_2$ Concentrations at Loobos vs Wekerom')
plt.grid(True)
plt.savefig('Output/RQ1/NO2_ScatterPlot.png')
plt.show()

# Calculate the mean difference between the two sites, and plot a histogram of the differences
data['NO2_Difference'] = data['NO2_Loobos'] - data['NO2_Wekerom']
mean_difference = data['NO2_Difference'].mean()
print("Mean Difference in $NO_2$ Concentration between Loobos and Wekerom: ", mean_difference)
plt.figure(figsize=(10, 6))
plt.hist(data['NO2_Difference'], bins=30, edgecolor='black')
plt.xlabel('Difference in $NO_2$ Concentration (Loobos - Wekerom) ($\mu g/m^3$)')
plt.ylabel('Frequency')
plt.title('Histogram of $NO_2$ Concentration Differences between Loobos and Wekerom')
plt.grid(True)
plt.savefig('Output/RQ1/NO2_Difference_Histogram.png')
plt.show()
