# A2: Data acquisition for flux data at Loobos flux tower, this script is a part of data acquisition process for this thesis, 
# which we use API to download flux data related to later analysis, like footprint calculation and meteorological data analysis. 
# By using this script, we download the raw data for the period of 01/01/2023 to 31/12/2025, and merge them together, calculate
# half hourly average, to align with different measurement frequency of different variables, and save the final data for later 
# analysis.
# Author: Zhiyu Wu, Date: 28/04/2026

# ------ Step 1: Import necessary libraries ------
import datetime as dt
import os
from data import fetch_data
import glob
import pandas as pd

# ------ Step 2: Define the time range for data acquisition and station number------
start_date = dt.datetime(2023, 1, 1)
end_date = dt.datetime(2025, 12, 31)

site = 2 # 2 = Loobos

# ------ Step 3: Select interested variables, add API key information------
# Select interested variables, add API information
variables = ["SW_IN_1_1_1", "SW_OUT_1_1_1", "u*", "L", "LE", "H", "WS_2_1_1", "WD_2_1_1", "P_1_1_1", "RH_1_1_1", "TA_1_1_1"]
# Variables include: incoming shortwave radiation, outgoing shortwave radiation, friction velocity, Obukhov length, latent heat flux, sensible heat flux, wind speed, wind direction, air pressure, relative humidity and air temperature.
# Load API key from environment variable
API_KEY = os.getenv("MY_API_KEY")
# Check if API key is available
if API_KEY is None:
    raise ValueError("Missing API key. Set the MY_API_KEY environment variable.")

# ------ Step 4: Fetch the data and save as csv file ------
# Make sure save directory exists
save_dir = 'Data/VeluweInSitu'
os.makedirs(save_dir, exist_ok=True)
# ---- Monthly loop ----
current = start_date

while current < end_date:

    # Get first day of next month
    next_month = (current.replace(day=28) + dt.timedelta(days=4)).replace(day=1)

    # Prevent overshooting end_date
    request_end = min(next_month, end_date)

    # Create monthly filename
    save_filename = os.path.join(
        save_dir,
        f"Loobos_Flux_{current.strftime('%Y_%m')}.csv"
    )

    print(f"Downloading {current.date()} to {request_end.date()}")

    fetch_data(
        current,
        request_end,
        site,
        variables,
        API_KEY,
        True,
        save_filename
    )

    current = next_month
# ------ Step 5: Merge monthly files, calculate half hourly average and save the final data ------
# Get all monthly files
all_files = glob.glob(os.path.join(save_dir, "Loobos_Flux_*.csv"))
df_list = []
for file in all_files:
    df = pd.read_csv(file, skiprows=[1], low_memory=False)  # skip units row in each file
    df_list.append(df)

# Concatenate all dataframes
merged_df = pd.concat(df_list, ignore_index=True)

# Save merged dataframe to a new csv file
merged_df.to_csv(os.path.join(save_dir, "Loobos_Flux_2023_2025.csv"), index=False)

# Convert Timestamp to datetime and set as index
merged_df['Timestamp'] = pd.to_datetime(merged_df['Timestamp'], errors='coerce')
merged_df.set_index('Timestamp', inplace=True)
# Convert remaining columns to numeric
merged_df = merged_df.apply(pd.to_numeric, errors='coerce')
# Resample to half hourly average
half_hourly_df = merged_df.resample('30min').mean()
# Save the half hourly averaged data to a new csv file
half_hourly_df.to_csv(os.path.join(save_dir, "Loobos_Flux_30min_23-25.csv"), index=True)

# Separate the timestamp into Time (date) and timestamp (time) columns for later use in footprint calculation.
half_hourly_df['Time'] = half_hourly_df.index.date
half_hourly_df['Timestamp'] = half_hourly_df.index.time
half_hourly_df.to_csv(os.path.join(save_dir, "Loobos_Flux_30min_23-25.csv"), index=False)

# ------ Step 6: Clean up original monthly files and intermediate files ------
# Delete the original monthly files and the merged file to save space
for file in all_files:
    os.remove(file)
os.remove(os.path.join(save_dir, "Loobos_Flux_2023_2025.csv"))