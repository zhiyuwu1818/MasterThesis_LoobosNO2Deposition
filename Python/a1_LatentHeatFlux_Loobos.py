# A1: Data aquisition and preprocessing for latent heat flux data at Loobos flux tower, in this script, we first download the latent heat flux data from MAQ website, then 
# we will merge it with old flux tower's latent heat flux dataset, and finally delete all the original monthly files to save space. The final merged dataset will be used 
# for the following analysis and modeling.

# Author: Zhiyu Wu, Date: 24/02/2026
# ------ Step 1: Import necessary libraries ------
import datetime as dt
import os
from data import fetch_data
import pandas as pd
import numpy as np
import glob

# ------ Step 2: Define the time range for data acquisition and station number------
start_date = dt.datetime(2021, 9, 1)
end_date   = dt.datetime(2025, 12, 31)

site = 2   # 1=Veenkampen, 2=Loobos, 3=Amsterdam, for this research we choose Loobos flux tower

# ------ Step 3: Select interested variables, add API key information------
variables = [ 'LE'] # latent heat flux variable, unit is W/m^2
API_KEY = os.getenv("MY_API_KEY") # make sure to set the API key in your environment variables, for example, in terminal, you can get the API key from MAQ website.

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
        f"Loobos_LE_{current.strftime('%Y_%m')}.csv"
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

# ------ Step 5: Read old flux dataset downloaded beforehand------
df_old = pd.read_csv(os.path.join(save_dir, 'FluxTowerOld/NL-Loo_EC.csv'))
# Select only the relevant columns: 'TIMESTAMP_END' and 'LE_1_1_1' (Latent Heat Flux)
df_old = df_old[['TIMESTAMP_END', 'LE_1_1_1']]
# Convert 'TIMESTAMP_END' to datetime format
df_old['TIMESTAMP_END'] = pd.to_datetime(df_old['TIMESTAMP_END'], format='%Y%m%d%H%M')
# Rename columns for clarity and consistency
df_old.rename(columns={'TIMESTAMP_END': 'Timestamp', 'LE_1_1_1': 'LE'}, inplace=True)
# Minus 1 hour to align with UTC time (original data in central european time and does not consider daylight saving time)
df_old['Timestamp'] = df_old['Timestamp'] - pd.Timedelta(hours=1)
# Select only the data from Jan 2020
df_old = df_old[df_old['Timestamp'] >= '2020-01-01']
# Discard rows with NA values (-9999) in 'LE'
df_old = df_old[df_old['LE'] != -9999]
# Check the processed data
print(df_old.head())

# ------ Step 6: Merge the new downloaded dataset with old dataset ------
# For the period 202109-202512, Merge all monthly files for Loobos_LE into one file, skipping the second row with units

# Find all files containing "Loobos_LE"
file_pattern = os.path.join(save_dir, '*Loobos_LE*.csv')
files = glob.glob(file_pattern)

# Check if files were found
if not files:
    raise FileNotFoundError("No files containing 'Loobos_LE' were found.")

print(f"Found {len(files)} files.")

# Read and combine all files
df_list = []

for file in sorted(files):
    print(f"Reading {file}")
    
    df = pd.read_csv(file, skiprows=[1])  # Skip the second row with units
    
    # Convert timestamp column to datetime
    df['Timestamp'] = pd.to_datetime(df['Timestamp'], utc=True)
    
    df_list.append(df)

# Merge all into one DataFrame
merged_df = pd.concat(df_list, ignore_index=True)

# Remove possible duplicates
merged_df = merged_df.drop_duplicates(subset='Timestamp')

# Sort by timestamp
merged_df = merged_df.sort_values('Timestamp')

# Reset index
merged_df = merged_df.reset_index(drop=True)

# Merge old dataset and new merged dataset based on timestamp
# First, ensure both datasets have the same timestamp format and are in the same timezone (UTC)
df_old['Timestamp'] = pd.to_datetime(df_old['Timestamp'], utc=True)
# Merge the datasets on 'Timestamp', using an inner join and if there are timestamps both datasets have values, we keep the values from the new merged dataset (Loobos_LE)
merged_df = pd.merge(df_old, merged_df, on='Timestamp', how='outer', suffixes=('_old', '_new'))
# If both datasets have values for a timestamp, keep the value from the new merged dataset (Loobos_LE)
merged_df['LE'] = merged_df['LE_new'].fillna(merged_df['LE_old'])
# Drop the old and new columns
merged_df = merged_df.drop(columns=['LE_old', 'LE_new'])
# Save final merged file
output_file = os.path.join(save_dir, 'Loobos_LatentHeatFlux_2020_2025.csv')
merged_df.to_csv(output_file, index=False)

print("Done! Merged file saved as:")
print(output_file)

# Check the first few rows of the merged DataFrame
print(merged_df.head())

# --- Step 7: Delete the original monthly files to save space ---
# Delete the source files after merging
for file in files:
    os.remove(file)
    print(f"Deleted {file}")