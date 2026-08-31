# B1: Calculate flux footprint data at Loobos flux tower, in this script, we will use the meteological data at Loobos flux tower
# to calculate the flux footprint data for September 2025, which will be used for later analysis in this thesis. The function to 
# calculate flux footprint is from  Kljun et al. (2015), and we created a function to calculate monthly average footprint.
# Author: Zhiyu Wu, Date: 02/03/2026

# import numpy and pandas
import numpy as np
import pandas as pd

# Import function to calculate monthly footprint climatology
from analysis.calc_period_FFP import cal_period_FFP
# Import function to export footprint shapefile
from gis.export_footprint import export_footprint_shapefile 

# Read in-situ meteorology data for footprint calculation
insitu_data = pd.read_csv('Data/VeluweInSitu/Loobos_Flux_30min_23-25.csv')
# Ensure Timestamp column is string for compatibility with cal_period_FFP
insitu_data['Time'] = pd.to_datetime(insitu_data['Timestamp']).dt.date
insitu_data['Timestamp'] = pd.to_datetime(insitu_data['Timestamp']).dt.time
insitu_data['Timestamp'] = insitu_data['Timestamp'].astype(str)

# Define measurement height (m) and default boundary layer height (m) for Loobos site, the value of boundary layer height 
# will be adjusted based on the actual conditions during the measurement period later
zmt = 38 - 22 * 2 / 3
ht = 1500 # default boundary layer height = 1500 m during the day

# Define the location of the flux tower for exporting the shapefile at correct location later
lat = 52.167778
lon = 5.743889

# Yearly variations
# Loop over each year
for year in range(2023, 2026):
    start_date = f'{year}-01-01'
    end_date   = f'{year}-12-31'

    # Loop over different rs contour values (90% down to 50%)
    for rs in [0.9, 0.8, 0.7, 0.6, 0.5]:
        # Calculate annual footprint climatology
        annual_footprint = cal_period_FFP(start_date, end_date, zmt, ht, rs, insitu_data)

        # Export the footprint shapefile for later analysis in GIS software
        shapefile_name = f'Loobos_FFP_{year}_rs{int(rs*100)}.shp'
        export_footprint_shapefile(annual_footprint, str(year), rs, lat, lon, shapefile_name)

# Make seasonal footprint for the year from 12/2024 to 11/2025 
# For seasonal variations, we can define the seasons as follows:
# Spring: March 1 - May 31
# Summer: June 1 - August 31
# Autumn: September 1 - November 30
# Winter: December 1 - February 28/29
seasons = {
    'Spring_2025': ('2025-03-01', '2025-05-31'),
    'Summer_2025': ('2025-06-01', '2025-08-31'),
    'Autumn_2025': ('2025-09-01', '2025-11-30'),
    'Winter_2024': ('2024-12-01', '2025-02-28')
}
# Loop over each season
for season, (start_date, end_date) in seasons.items():
    # Loop over different rs contour values (90% down to 50%)
    for rs in [0.9, 0.8, 0.7, 0.6, 0.5]:
        # Calculate seasonal footprint climatology
        seasonal_footprint = cal_period_FFP(start_date, end_date, zmt, ht, rs, insitu_data)

        # Export the footprint shapefile for later analysis in GIS software
        output_folder = f'../Data/Footprint'
        shapefile_name = f'Loobos_FFP_{season}_rs{int(rs*100)}.shp'
        export_footprint_shapefile(seasonal_footprint, season, rs, lat, lon, output_folder)
