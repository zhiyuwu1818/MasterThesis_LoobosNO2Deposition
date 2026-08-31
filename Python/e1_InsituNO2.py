# E1: Data processing and analysis for in-situ NO2 data at Loobos flux tower and Wekerom sites air monitoring site, this script solves 
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
# Select columns containing time and NO2 data, which are 'NO2_LGR_30m' and 'NO2_Wek_30m'
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

# ------ Step 5: Diurnal and daily variation heatmaps (temperature and NO2) ------
# select the temperature data from the original dataset, which is 'TA_1_1_1', rename it to 'Temperature', and merge it with the NO2 data
temperature_data = pd.read_csv('Data/VeluweInSitu/NO2&meteodata.csv' , sep='\s+', skiprows=3)
temperature_data['Time'] = pd.to_datetime(temperature_data['Time'] + ' ' + temperature_data['Timestamp'], format='%Y-%m-%d %H:%M:%S')
temperature_data.drop(columns=['Timestamp'], inplace=True)
temperature_data = temperature_data[['Time', 'TA_1_1_1']]
temperature_data.rename(columns={'TA_1_1_1': 'Temperature'}, inplace=True)
data = pd.merge(data, temperature_data, on='Time', how='inner')

# Create a copy with date and hour columns
heatmap_data = data.copy()
heatmap_data = heatmap_data[(heatmap_data['Temperature'] != -9999.999) & (heatmap_data['Temperature'] >= -50)]
heatmap_data['Date'] = heatmap_data['Time'].dt.date
heatmap_data['Hour'] = heatmap_data['Time'].dt.hour + heatmap_data['Time'].dt.minute / 60

# Pivot tables: rows = hour, columns = date
temp_pivot = heatmap_data.pivot_table(
    index='Hour', columns='Date', values='Temperature', aggfunc='mean'
)
no2_pivot = heatmap_data.pivot_table(
    index='Hour', columns='Date', values='NO2_Loobos', aggfunc='mean'
)

# --- Align both pivots to the same date columns ---
all_dates = sorted(set(temp_pivot.columns) | set(no2_pivot.columns))
temp_pivot = temp_pivot.reindex(columns=all_dates)
no2_pivot  = no2_pivot.reindex(columns=all_dates)

# Shared x-axis labels
date_labels = [str(d) for d in all_dates]
n_dates = len(date_labels)

# --- Plot ---
fig, axes = plt.subplots(2, 1, figsize=(20, 10), sharex=True)
fig.subplots_adjust(hspace=0.3)

# --- Temperature heatmap ---
im1 = axes[0].imshow(
    temp_pivot.values,
    aspect='auto',
    origin='lower',
    cmap='RdYlBu_r',
    interpolation='nearest',
    extent=[0, n_dates, temp_pivot.index.min(), temp_pivot.index.max()]
)
cbar1 = plt.colorbar(im1, ax=axes[0], pad=0.01)
cbar1.set_label('Temperature (°C)', fontsize=11)
axes[0].set_ylabel('Hour (LT)', fontsize=12)
axes[0].set_title('Diurnal and Daily Variation of Temperature', fontsize=13, fontweight='bold')
axes[0].set_yticks(range(0, 25, 3))

# --- NO2 heatmap ---
no2_values = no2_pivot.values
vmin = np.nanpercentile(no2_values, 5)
vmax = np.nanpercentile(no2_values, 95)

im2 = axes[1].imshow(
    no2_values,
    aspect='auto',
    origin='lower',
    cmap='YlOrRd',
    interpolation='nearest',
    vmin=vmin,
    vmax=vmax,
    extent=[0, n_dates, no2_pivot.index.min(), no2_pivot.index.max()]
)
cbar2 = plt.colorbar(im2, ax=axes[1], pad=0.01)
cbar2.set_label('NO$_2$ (µg/m³)', fontsize=11)
axes[1].set_ylabel('Hour (LT)', fontsize=12)
axes[1].set_title('Diurnal and Daily Variation of NO$_2$ Concentration', fontsize=13, fontweight='bold')
axes[1].set_yticks(range(0, 25, 3))

# --- Shared x-axis ticks: label every ~30 days ---
tick_spacing = 30
tick_positions = list(range(0, n_dates, tick_spacing))
tick_labels = [date_labels[i] for i in tick_positions]
axes[1].set_xticks(tick_positions)
axes[1].set_xticklabels(tick_labels, rotation=45, ha='right', fontsize=9)
axes[1].set_xlabel('Date', fontsize=12)

plt.suptitle('Diurnal and Daily Variation of Temperature and NO$_2$ at Loobos',
             fontsize=14, fontweight='bold', y=1.01)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'Diurnal_Daily_Heatmap_Temp_NO2.png'),
            bbox_inches='tight', dpi=150)
plt.show()
