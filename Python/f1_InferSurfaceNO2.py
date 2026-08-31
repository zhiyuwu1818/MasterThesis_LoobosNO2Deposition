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

# import the function to process SCM NO2 vertical profiles and calculate VCD
from analysis import process_SCM_NO2_VCD
from analysis import infer_surface_NO2

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
# ---- in-situ NO2 data at Loobos flux tower, from April 2025 to September 2025 ----
insitu_no2 = pd.read_csv('Data/VeluweInSitu/NO2&meteodata.csv', sep='\s+', skiprows=3)
# Select the daily period of 10 o'clock to 13 o'clock, which is the typical overpass time of TROPOMI, and calculate the daily average NO2 concentration for both sites
# Merge Time and Timestamp columns into a single datetime column, convert to datatime format,
insitu_no2['datetime'] = pd.to_datetime(insitu_no2['Time'].astype(str) + ' ' + insitu_no2['Timestamp'].astype(str))
# For dataset, for each day, calculate the average NO2 pollution for the hour 10 o'clock to 13 o'clock for two locations
insitu_no2['Hour'] = pd.to_datetime(insitu_no2['datetime']).dt.hour
# Filter out all the NA data represented by -9999.999 in the original dataset
insitu_no2 = insitu_no2[(insitu_no2['NO2_LGR_30m'] != -9999.999)]
# Convert hour column to UTC +0 if necessary (assuming original data is in local time UTC+2)
insitu_no2['Hour'] = (insitu_no2['Hour'] - 2) % 24
# Filter data for hours between 10 and 13 o'clock
daytime_no2_data = insitu_no2[(insitu_no2['Hour'] >= 10) & (insitu_no2['Hour'] <= 13)]
# Group by date and calculate average NO2 pollution for each day
daytime_no2_data['Date'] = pd.to_datetime(daytime_no2_data['Time']).dt.date
insitu_daily_avg = (
    daytime_no2_data
    .groupby('Date')[['NO2_LGR_30m', 'NO2_Wek_30m']]
    .mean()
    .reset_index()
)

# Rename columns for clarity
insitu_daily_avg.rename(columns={'NO2_LGR_30m': 'NO2_30m_Loobos', 'NO2_Wek_30m': 'NO2_30m_Wekerom'}, inplace=True)
# Select only Loobos site for comparison with TROPOMI
insitu_no2_Loobos_daily = insitu_daily_avg[['Date', 'NO2_30m_Loobos']]

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
         marker='o', markersize=5, label='SCM Profile')
plt.plot(tropomi_plot['datetime'],    tropomi_plot['Tropospheric_NO2_value'],
         marker='o', markersize=5, label='TROPOMI')
plt.xlabel('Date', fontsize=14)
plt.ylabel('NO$_2$ VCD (10$^{15}$ molecules/cm$^2$)', fontsize=14)
plt.title('Daily NO$_2$ VCDs at Loobos Flux Tower: SCM inferred vs TROPOMI (01/2025 - 08/2025)', fontsize=16)
plt.legend()
plt.xticks(fontsize=14, rotation=30)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'NO2_VCD_Comparison_SCM_TROPOMI.png'))

# ------ Step 4: Infer surface NO2 concentration using TROPOMI VCD and SCM vertical profiles, and compare with in-situ measurements ------
inferred_SL_NO2 = infer_surface_NO2(scm_no2vcd, scm_prof_NO2, TROPOMI_NO2VCD)

# Plot the inferred surface NO2 concentration against in-situ measurements for the period 05/2025 to 08/2025, save the results figure
# Filter both datasets for the period 05/2025 to 08/2025
mask_inferred = (inferred_SL_NO2['Date'] >= pd.Timestamp('2025-05-01').date()) & \
                (inferred_SL_NO2['Date'] <= pd.Timestamp('2025-08-31').date())
mask_insitu   = (insitu_no2_Loobos_daily['Date'] >= pd.Timestamp('2025-05-01').date()) & \
                (insitu_no2_Loobos_daily['Date'] <= pd.Timestamp('2025-08-31').date())

# Plot the inferred surface NO2 concentration against in-situ measurements for the period 05/2025 to 08/2025
plt.figure(figsize=(14, 8))
plt.plot(inferred_SL_NO2[mask_inferred]['Date'], inferred_SL_NO2[mask_inferred]['NO2_SL'],
         marker='x', markersize=8, linestyle='--', label='Inferred NO$_2$ Surface Layer')
plt.plot(insitu_no2_Loobos_daily[mask_insitu]['Date'], insitu_no2_Loobos_daily[mask_insitu]['NO2_30m_Loobos'],
         marker='o', markersize=8, linestyle='-', label='In-situ NO$_2$ at 30m (Loobos)')
plt.xlabel('Date', fontsize=14)
plt.ylabel('NO$_2$ Concentration (ppb)', fontsize=14)
plt.title('Inferred Surface NO$_2$ vs In-situ Measurements at Loobos (05/2025 - 08/2025)', fontsize=16)
plt.legend(fontsize=12)
plt.xticks(fontsize=14, rotation=30)
plt.yticks(fontsize=14)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'Inferred_Surface_NO2_vs_InSitu_Loobos.png'))