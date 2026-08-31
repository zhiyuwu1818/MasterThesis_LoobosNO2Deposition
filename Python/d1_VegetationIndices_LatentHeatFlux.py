# D1: Data analysis of vegetation indices and latent heat flux. n this script, we will analyze the relationship between five vegetation indices (NDVI, 
# NDRE, CIred-edge, MTCI, MSAVI, NIRv, EVI (need to add CRswir) and latent heat flux at the Loobos flux tower site for the period of 2020-2025, and save the plots.

# Author: Zhiyu Wu, Date: 06/03/2026

# ------ Step 1: Import necessary libraries, assign (and create if necessary) output folder ------
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import os
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

output_folder = 'Output/RQ1'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# ------ Step 2: Load the vegetation indices and latent heat flux datasets ------
# Vegetation indices csv, which contains columns of 'Date', 'NDVI', 'NDRE', 'CIred-edge', 'MTCI', 'MSAVI', 'NIRv', 'EVI'
VIs_data = pd.read_csv('Data/Sentinel2/VegetationIndices_TimeSeries_New.csv', parse_dates=['Date'])
# Read another vegetation indices csv, select only the column of "CRswir" and "Date"
CRswir_data = pd.read_csv('Data/Sentinel2/VegetationIndices_TimeSeries.csv', parse_dates=['Date'], usecols=['Date', 'CRSWIR'])
# Merge the CRswir data with the VIs_data on the 'Date' column to create a combined DataFrame
VIs_data = pd.merge(VIs_data, CRswir_data, on='Date', how='left')
# Print out the first few rows of the vegetation indices data to check
print(VIs_data.head())

# Latent heat flux csv, which contains columns of 'Timestamp' and 'LE'
LE_data = pd.read_csv('Data/VeluweInSitu/Loobos_LatentHeatFlux_2020_2025.csv', parse_dates=['Timestamp'])

# ------ Step 3: Data processing ------
# For vegetation indices, create centered 5-point rolling mean to smooth the data, min_periods=1 means that if there are less than 5 data points 
# at the beginning and end of the time series, it will still calculate the mean with the available data points
# Create a function to calculate centered 5-point rolling mean to smooth the data, 
def smooth_vegetation_indices(df, index_col, date_col='Date', window=5):
    """
    Smooth the vegetation index using a centered rolling mean.
    Parameters:
    - df: DataFrame containing the vegetation index and date columns.
    - index_col: Name of the column containing the vegetation index to be smoothed.
    - date_col: Name of the column containing the date information (default is 'Date').
    - window: Size of the rolling window (default is 5).
    Returns:
    - DataFrame with an additional column for the smoothed vegetation index.
    The function first ensures that the date column is in datetime format, then drops any rows with missing values in the vegetation index column. 
    It calculates the centered rolling mean for the specified vegetation index and adds it as a new column to the DataFrame. The smoothed column 
    is named by appending '_smooth' to the original index column name.

    """
    # Ensure date column is datetime
    df = df.copy()
    df[date_col] = pd.to_datetime(df[date_col])

    # Drop missing values
    df = df.dropna(subset=[index_col])

    # Smooth the index
    smooth_col = f"{index_col}_smooth"
    df[smooth_col] = df[index_col].rolling(
        window=window, center=True, min_periods=1
    ).mean()
    return df
# Apply the smoothing function to each vegetation index
for index in ['NDVI', 'NDRE', 'CIredge', 'MTCI', 'MSAVI', 'NIRv', 'EVI', 'CRSWIR']:
    VIs_data = smooth_vegetation_indices(VIs_data, index)   
# print out the first few rows of the smoothed vegetation indices data to check
print(VIs_data.head())

# For latent heat flux, select day time period, then calcualte monthly average
LE_data = LE_data[(LE_data['Timestamp'].dt.hour >= 9) & (LE_data['Timestamp'].dt.hour <= 17)]
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
for vi in ['NDVI', 'NDRE', 'CIredge', 'MTCI', 'MSAVI', 'NIRv', 'EVI', 'CRSWIR']:
    create_timeseries_plot(
        VIs_data, vi, vi,
        LE_monthly, 'LE', 'Latent Heat Flux',
        os.path.join(output_folder, f'{vi}_LE_Timeseries.png')
    )

def vi_le_regression(VIs_data, LE_monthly, vi_col, output_folder=None, ax=None):
    """
    Calculate monthly average of a smoothed vegetation index, merge with monthly
    latent heat flux, perform linear regression, and plot the results.

    Parameters
    ----------
    VIs_data : pandas DataFrame
        DataFrame containing vegetation index data with a 'Date' column and
        smoothed VI columns (e.g. 'CRSWIR_smooth', 'NDVI_smooth').
    LE_monthly : pandas DataFrame
        DataFrame containing monthly latent heat flux with columns 'Month' and 'LE'.
    vi_col : str
        Name of the raw VI column (e.g. 'CRSWIR', 'NDVI').
        The smoothed column is automatically inferred as vi_col + '_smooth'.
    output_folder : str, optional
        Folder to save the output figure when plotted standalone. If None, shown interactively.
    ax : matplotlib Axes, optional
        If provided, plot into this axes (for combined figure). If None, creates its own figure.

    Returns
    -------
    merged : pandas DataFrame
        Merged DataFrame with monthly VI and LE values used for regression.
    slope : float
    intercept : float
    r2 : float
    """
    vi_smooth_col = f"{vi_col}_smooth"
    if vi_smooth_col not in VIs_data.columns:
        raise ValueError(f"Smoothed column '{vi_smooth_col}' not found. "
                         f"Available: {VIs_data.columns.tolist()}")

    # Calculate monthly average of smoothed VI
    VIs_data['Month'] = VIs_data['Date'].dt.to_period('M')
    vi_monthly        = VIs_data.groupby('Month')[vi_smooth_col].mean().reset_index()
    vi_monthly['Month'] = vi_monthly['Month'].dt.to_timestamp()

    # Ensure LE_monthly has a proper timestamp column
    le = LE_monthly.copy()
    le['Month'] = le['Month'].dt.to_timestamp() if hasattr(le['Month'].dtype, 'freq') else le['Month']

    # Merge on Month
    merged = pd.merge(vi_monthly, le[['Month', 'LE']], on='Month', how='inner')
    merged[vi_smooth_col] = merged[vi_smooth_col].interpolate(method='linear')
    merged = merged.dropna(subset=[vi_smooth_col, 'LE'])

    # Linear regression
    X     = merged[vi_smooth_col].values.reshape(-1, 1)
    y     = merged['LE'].values
    model = LinearRegression()
    model.fit(X, y)
    y_pred    = model.predict(X)
    r2        = r2_score(y, y_pred)
    slope     = model.coef_[0]
    intercept = model.intercept_
    print(f"LE = {slope:.2f} * {vi_col} + {intercept:.2f}, R² = {r2:.2f}")

    # Plot into provided ax or create a standalone figure
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(8, 6))

    ax.scatter(merged[vi_smooth_col], merged['LE'],
               color='blue', alpha=0.6, label='Data Points', zorder=3)
    ax.plot(merged[vi_smooth_col], y_pred, color='red',
            label=f'LE = {slope:.2f} * {vi_col} + {intercept:.2f}\nR² = {r2:.2f}')
    ax.set_xlabel(f'Smoothed {vi_col}', fontsize=14)
    ax.set_ylabel('Latent Heat Flux (W/m$^2$)', fontsize=14)
    ax.set_title(f'LE vs {vi_col}', fontsize=15)
    ax.legend(fontsize=12)
    ax.tick_params(axis='both', labelsize=13)

    if standalone:
        plt.tight_layout()
        if output_folder:
            plt.savefig(os.path.join(output_folder, f'{vi_col}_LE_LinearRegression.png'))
        else:
            plt.show()

    return merged, slope, intercept, r2

vi_list = ['CRSWIR', 'NDVI', 'NIRv', 'MSAVI', 'EVI', 'NDRE', 'CIredge', 'MTCI']

n_cols = 2
n_rows = math.ceil(len(vi_list) / n_cols)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(22, 20))
axes = axes.flatten()

le_regression_results = []

for i, vi in enumerate(vi_list):
    merged, slope, intercept, r2 = vi_le_regression(
        VIs_data, LE_monthly, vi_col=vi, ax=axes[i]
    )
    le_regression_results.append({'VI': vi, 'slope': slope, 'intercept': intercept, 'R2': r2})

# Hide unused subplots
for j in range(len(vi_list), len(axes)):
    axes[j].set_visible(False)

plt.suptitle('Latent Heat Flux vs Vegetation Indices (Monthly)', fontsize=18, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'LE_vs_all_VI_regression.png'), bbox_inches='tight')
plt.show()

# Summary table
le_regression_results_df = pd.DataFrame(le_regression_results)
print(le_regression_results_df)
# Save the summary table as a csv file
le_regression_results_df.to_csv(os.path.join(output_folder, 'LE_VI_Regression_Results.csv'), index=False)