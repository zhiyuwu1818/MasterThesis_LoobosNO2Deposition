# D1: Data analysis of vegetation indices and latent heat flux, in this script, we will analyze the relationship between five vegetation indices (NDVI, 
# kNDVI, CRswir, CIre and REP) and latent heat flux at the Loobos flux tower site for the period of 2020-2025, and save the plots.

# Author: Zhiyu Wu, Date: 06/03/2026

# ------ Step 1: Import necessary libraries, assign (and create if necessary) output folder ------
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

output_folder = 'Output/RQ1'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# ------ Step 2: Load the vegetation indices and latent heat flux datasets ------
# Vegetation indices csv, which contains columns of 'Date', 'NDVI', 'kNDVI', 'CRswir', 'CIre' and 'REP'
VI_data = pd.read_csv('Data/Sentinel2/VegetationIndices_TimeSeries.csv', parse_dates=['Date'])
# Latent heat flux csv, which contains columns of 'Timestamp' and 'LE'
LE_data = pd.read_csv('Data/VeluweInSitu/Loobos_LatentHeatFlux_2020_2025.csv', parse_dates=['Timestamp'])

# ------ Step 3: Data processing ------
# For vegetation indices, create centered 5-point rolling mean to smooth the data, min_periods=1 means that if there are less than 5 data points 
# at the beginning and end of the time series, it will still calculate the mean with the available data points
VI_data['NDVI_smooth'] = VI_data['NDVI'].rolling(window=5, center=True, min_periods=1).mean()
VI_data['kNDVI_smooth'] = VI_data['kNDVI'].rolling(window=5, center=True, min_periods=1).mean()
VI_data['CRSWIR_smooth'] = VI_data['CRSWIR'].rolling(window=5, center=True, min_periods=1).mean()
VI_data['CIre_smooth'] = VI_data['CIre'].rolling(window=5, center=True, min_periods=1).mean()
VI_data['REP_smooth'] = VI_data['REP'].rolling(window=5, center=True, min_periods=1).mean()
# For latent heat flux, calcualte monthly average
LE_data['Month'] = LE_data['Timestamp'].dt.to_period('M')
LE_monthly = LE_data.groupby('Month')['LE'].mean().reset_index()
LE_monthly['Month'] = LE_monthly['Month'].dt.to_timestamp() # Convert back to timestamp for plotting

# ------ Step 4: Data visualization ------
# Make a function to create timeseries plot for vegetation index, smoothed vegetation index and latent heat flux, save the plot as a png file
def create_timeseries_plot(vi_data, vi_column, vi_label, le_data, le_column, le_label, filename):
    # --- Create figure ---
    fig, ax1 = plt.subplots(figsize=(12, 6))

    # =========================
    # LEFT AXIS (VEGETATION INDEX)
    # =========================
    ax1.plot(
        vi_data['Date'],
        vi_data[vi_column],
        label=vi_label,
        color='lightblue',
        alpha=0.6
    )
    ax1.plot(
        vi_data['Date'],
        vi_data[f'{vi_column}_smooth'],
        label=f'Smoothed {vi_label} (5-point rolling mean)',
        color='green',
        linewidth=3
    )
    ax1.set_xlabel('Year', fontsize=16)
    ax1.set_ylabel(vi_label, fontsize=16, color='green')
    ax1.tick_params(axis='both', labelsize=14)

    # =========================
    # RIGHT AXIS (Latent Heat)
    # =========================
    ax2 = ax1.twinx()
    ax2.plot(
        le_data['Month'],
        le_data[le_column],
        label=le_label,
        color='red',
        linewidth=2
    )
    ax2.set_ylabel('Latent Heat Flux (W/m$^2$)', fontsize=12, color='red')
    ax2.tick_params(axis='y', labelcolor='red', labelsize=12)

    # =========================
    # Title and Legend
    # =========================
    fig.suptitle(f'{vi_label} and Latent Heat Flux Time Series for NL-Loo Footprint', fontsize=12)

    # Combine legends from both axes
    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc='upper right', fontsize=8)

    fig.autofmt_xdate()
    plt.tight_layout()
    plt.savefig(filename)

# Create timeseries plots for each vegetation index and latent heat flux
create_timeseries_plot(VI_data, 'NDVI', 'NDVI', LE_monthly, 'LE', 'Latent Heat Flux', 'Output/RQ1/NDVI_LE_Timeseries.png')
create_timeseries_plot(VI_data, 'kNDVI', 'kNDVI', LE_monthly, 'LE', 'Latent Heat Flux', 'Output/RQ1/kNDVI_LE_Timeseries.png')
create_timeseries_plot(VI_data, 'CRSWIR', 'CRSWIR', LE_monthly, 'LE', 'Latent Heat Flux', 'Output/RQ1/CRSWIR_LE_Timeseries.png')
create_timeseries_plot(VI_data, 'CIre', 'CIre', LE_monthly, 'LE', 'Latent Heat Flux', 'Output/RQ1/CIre_LE_Timeseries.png')
create_timeseries_plot(VI_data, 'REP', 'REP', LE_monthly, 'LE', 'Latent Heat Flux', 'Output/RQ1/REP_LE_Timeseries.png')

# create a linear regression (y = mx + c) plot between CRSWIR and latent heat flux, calculate out the R-squared value, equation and show it in the plot, save the plot as a png file
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

# Calculate monthly average of smoothed CRSWIR
VI_data['Month'] = VI_data['Date'].dt.to_period('M')
CRSWIR_monthly = VI_data.groupby('Month')['CRSWIR_smooth'].mean().reset_index()
CRSWIR_monthly['Month'] = CRSWIR_monthly['Month'].dt.to_timestamp()

# Ensure LE_monthly also has a proper timestamp column for merging
LE_monthly['Month'] = LE_monthly['Month'].dt.to_timestamp() if hasattr(LE_monthly['Month'].dtype, 'freq') else LE_monthly['Month']

# merge on Month so both arrays have the same length ──
merged = pd.merge(CRSWIR_monthly, LE_monthly[['Month', 'LE']], on='Month', how='inner')

# Interpolate any remaining NaNs in CRSWIR after the merge
merged['CRSWIR_smooth'] = merged['CRSWIR_smooth'].interpolate(method='linear')

# Drop rows where LE is still NaN (if any)
merged = merged.dropna(subset=['CRSWIR_smooth', 'LE'])

X = merged['CRSWIR_smooth'].values.reshape(-1, 1)
y = merged['LE'].values

# Fit the linear regression model
model = LinearRegression()
model.fit(X, y)

# Predict and calculate R-squared
y_pred = model.predict(X)
r2 = r2_score(y, y_pred)

# Plot
plt.figure(figsize=(8, 6))
plt.scatter(merged['CRSWIR_smooth'], merged['LE'], color='blue', alpha=0.6, label='Data Points')
plt.plot(merged['CRSWIR_smooth'], y_pred, color='red',
         label=f'Linear Fit: y = {model.coef_[0]:.2f}x + {model.intercept_:.2f}\nR² = {r2:.2f}')
plt.xlabel('Smoothed CRSWIR', fontsize=12)
plt.ylabel('Latent Heat Flux (W/m$^2$)', fontsize=12)
plt.title('Linear Regression between Smoothed CRSWIR and Latent Heat Flux', fontsize=14)
plt.legend(fontsize=10)
plt.tight_layout()
plt.savefig('Output/RQ1/CRSWIR_LE_LinearRegression.png')

