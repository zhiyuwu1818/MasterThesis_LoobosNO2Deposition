# F1: Inferring surface NO2 concentration from TROPOMI and SCM data, this script solves the first part of RQ 2 
# of the thesis, which is to infer the surface NO2 concentration using TROPOMI VCD and SCM vertical profiles, and 
# compare the inferred surface NO2 concentration with in-situ measurements.

# ------ Step 1: Import necessary libraries ------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import math
import os

from analysis import process_SCM_NO2_VCD
from analysis import infer_surface_NO2
from scipy import stats


# ------ Step 2: Load the datasets and pre-process, convert units, create output folder ------
# ---- TROPOMI RS NO2 VCD data for the year 2025 ----
tropomi_no2 = pd.read_csv('Data/TROPOMI/TropNO2_LoobosWekerom_3kmBuf_20CF_2025.csv')
# Choose data with the column 'location' as 'Loobos'
Loobos_no2 = tropomi_no2[tropomi_no2['location'] == 'Loobos']
# Calculate the daily average 
Loobos_no2['datetime'] = pd.to_datetime(Loobos_no2['datetime'])
TROPOMI_NO2VCD = Loobos_no2.groupby(Loobos_no2['datetime'].dt.date)['Tropospheric_NO2_value'].mean().reset_index()
TROPOMI_NO2VCD['datetime'] = pd.to_datetime(TROPOMI_NO2VCD['datetime'])
# Convert the unit from mol/m^2 to 10^15 molecules/cm^2
TROPOMI_NO2VCD['Tropospheric_NO2_value'] = TROPOMI_NO2VCD['Tropospheric_NO2_value'] * (6.022e23 / 1e4) * 1e-15
# print(TROPOMI_NO2VCD.head())
# ---- in-situ NO2 data at Loobos flux tower, from May 2025 to September 2025 ----
insitu_no2 = pd.read_csv('Data/VeluweInSitu/NO2&meteodata.csv', sep='\s+', skiprows=3)
# Merge Time and Timestamp columns into a single datetime column, convert to datetime format
insitu_no2['datetime'] = pd.to_datetime(insitu_no2['Time'].astype(str) + ' ' + insitu_no2['Timestamp'].astype(str))
# Extract hour and convert to UTC+0 (assuming original data is in local time UTC+2)
insitu_no2['Hour'] = pd.to_datetime(insitu_no2['datetime']).dt.hour
# Filter out all the NA data represented by -9999.999 in the original dataset
insitu_no2 = insitu_no2[(insitu_no2['NO2_LGR_30m'] != -9999.999)]
# Convert hour column to UTC+0
insitu_no2['Hour'] = (insitu_no2['Hour'] - 2) % 24
# Filter data for hours between 10 and 13 o'clock
daytime_no2_data = insitu_no2[(insitu_no2['Hour'] >= 10) & (insitu_no2['Hour'] <= 13)]
# Group by date and calculate average NO2 for each day
daytime_no2_data['Date'] = pd.to_datetime(daytime_no2_data['Time']).dt.date
insitu_daily_avg_no2 = (
    daytime_no2_data
    .groupby('Date')[['NO2_LGR_30m']]
    .mean()
    .reset_index()
)
insitu_daily_avg_no2.rename(columns={'NO2_LGR_30m': 'NO2_30m_Loobos'}, inplace=True)
insitu_daily_avg_no2['NO2_30m_Loobos_ugm3'] = insitu_daily_avg_no2['NO2_30m_Loobos']

# ---- Temperature data from Loobos flux dataset ----
flux_temp = pd.read_csv('Data/VeluweInSitu/Loobos_Flux_30min_23-25.csv', sep=',')
flux_temp['datetime'] = pd.to_datetime(flux_temp['Timestamp'], format='%Y-%m-%d %H:%M:%S')
flux_temp['Hour'] = pd.to_datetime(flux_temp['datetime']).dt.hour
# Filter out NA values for temperature
flux_temp = flux_temp[flux_temp['TA_1_1_1'] != -9999.999]
# Convert hour to UTC+0
flux_temp['Hour'] = (flux_temp['Hour'] - 2) % 24
# Filter for TROPOMI overpass window (10–13 UTC)
daytime_temp_data = flux_temp[(flux_temp['Hour'] >= 10) & (flux_temp['Hour'] <= 13)]
daytime_temp_data['Date'] = pd.to_datetime(daytime_temp_data['datetime']).dt.date
insitu_daily_avg_temp = (
    daytime_temp_data
    .groupby('Date')[['TA_1_1_1']]
    .mean()
    .reset_index()
)
insitu_daily_avg_temp.rename(columns={'TA_1_1_1': 'T_30m_Loobos'}, inplace=True)

# ---- Merge NO2 and temperature on Date ----
insitu_daily_avg = pd.merge(insitu_daily_avg_no2, insitu_daily_avg_temp, on='Date', how='outer')

# Select only Loobos site columns for comparison with TROPOMI
insitu_no2_Loobos_daily = insitu_daily_avg[['Date', 'NO2_30m_Loobos_ugm3']]

# ---- SCM vertical profile and physical property data for the period from 01/01/2025 to 31/08/2025, use a function to already calculate out the VCD derived from the vertical profiles ----
scm_no2_vcd, scm_prof_NO2 = process_SCM_NO2_VCD('Data/Model/NO2profiles11-14_2025.out', tropopause_height_m=12000)
# rename the column for NO2 VCD to be consistent with TROPOMI dataset
scm_no2_vcd.rename(columns={'NO2_VCD_1e15_molecules_per_cm2': 'Tropo_NO2_VCD'}, inplace=True)

# calculate the daily mean SCM VCD for the period from 01/01/2025 to 31/08/2025
scm_no2vcd = scm_no2_vcd.groupby(scm_no2_vcd['datetime'].dt.date)['Tropo_NO2_VCD'].mean().reset_index()
scm_no2vcd['datetime'] = pd.to_datetime(scm_no2vcd['datetime'])
# exclude NO2 VCD values greater than 100*10e15 molecules/cm² as outliers
scm_no2vcd = scm_no2vcd[scm_no2vcd['Tropo_NO2_VCD'] <= 100]

# Create output folder if it doesn't exist
output_folder = 'Output/RQ2'
os.makedirs(output_folder, exist_ok=True)

# ------ Step 3: Compare tropospheric NO2VCDs inferred from SCM and derived from TROPOMI ------
# Make time series comparison of the daily mean NO2 VCD from SCM and TROPOMI, save the restults figure
# Filter both datasets for the period of 01/2025 to 08/2025
mask_scm = (scm_no2vcd['datetime'] >= '2025-01-01') & (scm_no2vcd['datetime'] <= '2025-08-31')
mask_tropomi = (TROPOMI_NO2VCD['datetime'] >= '2025-01-01') & (TROPOMI_NO2VCD['datetime'] <= '2025-08-31')
scm_plot = scm_no2vcd.loc[mask_scm]
tropomi_plot = TROPOMI_NO2VCD.loc[mask_tropomi]

# Plot both datasets in one figure
plt.figure(figsize=(14, 8))
plt.plot(scm_plot['datetime'],        scm_plot['Tropo_NO2_VCD'],
         marker='o', markersize=2, label='SCM Profile')
plt.plot(tropomi_plot['datetime'],    tropomi_plot['Tropospheric_NO2_value'],
         marker='o', markersize=2, label='TROPOMI')
plt.xlabel('Date', fontsize=14)
plt.ylabel('NO$_2$ VCD (10$^{15}$ molecules/cm$^2$)', fontsize=14)
plt.title('Daily NO$_2$ VCDs at Loobos Flux Tower: SCM inferred vs TROPOMI (01/2025 - 08/2025)', fontsize=16)
plt.legend()
plt.xticks(fontsize=14, rotation=30)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'NO2_VCD_Comparison_SCM_TROPOMI.png'))

# Make another figure to show the scatter plot of SCM VCD vs TROPOMI VCD, calculate the r-squared, adding the 1:1 line, save the results figure
# Merge the two datasets on the date column to ensure we are comparing the same days
merged_vcd = pd.merge(scm_no2vcd, TROPOMI_NO2VCD, on='datetime', how='inner')
# Calculate R-squared
correlation_matrix = np.corrcoef(merged_vcd['Tropo_NO2_VCD'], merged_vcd['Tropospheric_NO2_value'])
correlation_xy = correlation_matrix[0,1]
r_squared = correlation_xy**2
print(f"R-squared between SCM VCD and TROPOMI VCD: {r_squared:.3f}")

# Calculate slope via linear regression
slope, intercept, _, _, _ = stats.linregress(merged_vcd['Tropo_NO2_VCD'], merged_vcd['Tropospheric_NO2_value'])
print(f"Slope: {slope:.3f}")

# Plot the scatter plot with 1:1 line and r-squared value in the legend
plt.figure(figsize=(8, 8))
plt.scatter(merged_vcd['Tropo_NO2_VCD'], merged_vcd['Tropospheric_NO2_value'], label=f'R² = {r_squared:.3f}\nSlope = {slope:.3f}')
plt.plot([0, max(merged_vcd['Tropo_NO2_VCD'].max(), merged_vcd['Tropospheric_NO2_value'].max())],
         [0, max(merged_vcd['Tropo_NO2_VCD'].max(), merged_vcd['Tropospheric_NO2_value'].max())],
         color='red', linestyle='--', label='1:1 Line')
plt.xlabel('SCM Inferred NO$_2$ VCD (10$^{15}$ molecules/cm$^2$)', fontsize=14)
plt.ylabel('TROPOMI NO$_2$ VCD (10$^{15}$ molecules/cm$^2$)', fontsize=14)
plt.title('SCM Inferred NO$_2$ VCD vs TROPOMI NO$_2$ VCD', fontsize=16)
plt.legend(fontsize=12)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'NO2_VCD_Scatter_SCM_TROPOMI.png'))

# ------ Step 4: Infer surface NO2 concentration using TROPOMI VCD and SCM vertical profiles, and compare with in-situ measurements ------
inferred_SL_NO2 = infer_surface_NO2(scm_no2vcd, scm_prof_NO2, TROPOMI_NO2VCD)

# Convert inferred NO2 from ppb to µg/m³ using temperature-dependent formula:
# µg/m³ = ppb × (12.187 × 46.01) / (273.15 + T_C)
# Use daily mean temperature from in-situ data for the conversion
inferred_SL_NO2 = pd.merge(inferred_SL_NO2, 
                            insitu_daily_avg[['Date', 'T_30m_Loobos']], 
                            on='Date', how='left')
inferred_SL_NO2['NO2_SL_ugm3'] = (
    inferred_SL_NO2['NO2_SL'] * (12.187 * 46.01)
) / (273.15 + inferred_SL_NO2['T_30m_Loobos'])

# save the inferred dataset to csv for record
inferred_SL_NO2.to_csv(os.path.join(output_folder, 'InferredSurfaceNO2.csv'), index=False)

# Plot the inferred surface NO2 concentration against in-situ measurements for the period 05/2025 to 08/2025, save the results figure
# Filter both datasets for the period 05/2025 to 08/2025
mask_inferred = (inferred_SL_NO2['Date'] >= pd.Timestamp('2025-05-01').date()) & \
                (inferred_SL_NO2['Date'] <= pd.Timestamp('2025-08-31').date())
mask_insitu   = (insitu_no2_Loobos_daily['Date'] >= pd.Timestamp('2025-05-01').date()) & \
                (insitu_no2_Loobos_daily['Date'] <= pd.Timestamp('2025-08-31').date())

# Plot the inferred surface NO2 concentration against in-situ measurements for the period 05/2025 to 08/2025
plt.figure(figsize=(14, 8))
plt.plot(inferred_SL_NO2[mask_inferred]['Date'], inferred_SL_NO2[mask_inferred]['NO2_SL_ugm3'],
         marker='x', markersize=8, linestyle='--', label='Inferred NO$_2$ Surface Layer')
plt.plot(insitu_no2_Loobos_daily[mask_insitu]['Date'], insitu_no2_Loobos_daily[mask_insitu]['NO2_30m_Loobos_ugm3'],
         marker='o', markersize=8, linestyle='-', label='In-situ NO$_2$ at 30m (Loobos)')
plt.xlabel('Date', fontsize=14)
plt.ylabel('NO$_2$ Concentration (µg m$^{-3}$)', fontsize=14)
plt.title('Inferred Surface NO$_2$ vs In-situ Measurements at Loobos (05/2025 - 08/2025)', fontsize=16)
plt.legend(fontsize=12)
plt.xticks(fontsize=14, rotation=30)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'Inferred_Surface_NO2_vs_InSitu_Loobos.png'))

# Plot the correlation between the inferred surface NO2 concentration and in-situ measurements, calculate the r-squared, save the results figure
# Merge the two datasets on the date column to ensure we are comparing the same days
merged_surface = pd.merge(inferred_SL_NO2[['Date', 'NO2_SL_ugm3']], 
                           insitu_no2_Loobos_daily, on='Date', how='inner')

# Calculate R-squared
correlation_matrix_surface = np.corrcoef(merged_surface['NO2_SL_ugm3'], merged_surface['NO2_30m_Loobos_ugm3'])
correlation_xy_surface = correlation_matrix_surface[0,1]
r_squared_surface = correlation_xy_surface**2
print(f"R-squared between inferred surface NO2 and in-situ NO2: {r_squared_surface:.3f}")

# Calculate slope via linear regression
slope_surface, intercept_surface, _, _, _ = stats.linregress(merged_surface['NO2_SL_ugm3'], merged_surface['NO2_30m_Loobos_ugm3'])
print(f"Slope: {slope_surface:.3f}")

# Plot the scatter plot with 1:1 line and r-squared value in the legend
plt.figure(figsize=(8, 8))
plt.scatter(merged_surface['NO2_SL_ugm3'], merged_surface['NO2_30m_Loobos_ugm3'], label=f'R² = {r_squared_surface:.3f}\nSlope = {slope_surface:.3f}')
plt.plot([0, max(merged_surface['NO2_SL_ugm3'].max(), merged_surface['NO2_30m_Loobos_ugm3'].max())],
         [0, max(merged_surface['NO2_SL_ugm3'].max(), merged_surface['NO2_30m_Loobos_ugm3'].max())],
         color='red', linestyle='--', label='1:1 Line')
plt.xlabel('Inferred NO$_2$ Surface Layer (µg m$^{-3}$)', fontsize=14)
plt.ylabel('In-situ NO$_2$ at 30m (µg m$^{-3}$)', fontsize=14)
plt.title('Inferred NO$_2$ Surface Layer vs In-situ NO$_2$ at Loobos', fontsize=16)
plt.legend(fontsize=12)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'Inferred_Surface_NO2_vs_InSitu_Loobos_Scatter.png'))