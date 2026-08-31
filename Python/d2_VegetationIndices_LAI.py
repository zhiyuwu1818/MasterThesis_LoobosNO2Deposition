# D2: Data analysis of vegetation indices and leaf area index (LAI). In this script, we will analyze the relationship between seven vegetation indices 
# (NDVI, NDRE, CIred-edge, MTCI, MSAVI, NIRv, EVI (need to add CRswir)) and LAI at the Loobos flux tower site for the period of 2020-2025, and save the plots.

# Author: Zhiyu Wu, Date: 11/05/2026

# ------ Step 1: Import necessary libraries, assign (and create if necessary) output folder ------
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
import math
from scipy.stats import linregress

output_folder = 'Output/RQ1'
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

# ------ Step 2: Load the vegetation indices and LAI datasets ------
# Vegetation indices csv, which contains columns of 'Date', 'NDVI', 'NDRE', 'CIred-edge', 'MTCI', 'MSAVI' and 'NIRv‘'
VIs_data = pd.read_csv('Data/Sentinel2/VegetationIndices_TimeSeries_New.csv', parse_dates=['Date'])
# Read another vegetation indices csv, select only the column of "CRswir" and "Date"
CRswir_data = pd.read_csv('Data/Sentinel2/VegetationIndices_TimeSeries.csv', parse_dates=['Date'], usecols=['Date', 'CRSWIR'])
# Merge the CRswir data with the VIs_data on the 'Date' column to create a combined DataFrame
VIs_data = pd.merge(VIs_data, CRswir_data, on='Date', how='left')
# Read the in-situ LAI csv, select rows where LAI_STATISTIC is "Mean"
lai_insitu = pd.read_csv('Data/VeluweInSitu/Loobos_LAI_23-25.csv')
lai_mean = lai_insitu[lai_insitu["LAI_STATISTIC"] == "Mean"].copy()

# ------ Step 3: Data processing ------
# LAI in-situ data
# Delete entry when the period between LAI_DATE_START and LAI_DATE_END is more than 30 days when there is LAI_DATE_START and LAI_DATE_END
lai_mean["DATE_START"] = pd.to_datetime(lai_mean["LAI_DATE_START"], format="%Y%m%d", errors="coerce")
lai_mean["DATE_END"] = pd.to_datetime(lai_mean["LAI_DATE_END"], format="%Y%m%d", errors="coerce")
lai_mean["PERIOD_DAYS"] = (lai_mean["DATE_END"] - lai_mean["DATE_START"]).dt.days
lai_mean = lai_mean[~((lai_mean["PERIOD_DAYS"] > 30) & lai_mean["DATE_START"].notna() & lai_mean["DATE_END"].notna())]

# Convert date column to datetime
lai_mean["Date"] = pd.to_datetime(lai_mean["LAI_DATEu"], format="%Y%m%d")


# Convert LAI to numeric (force invalid values to NaN)
lai_mean["LAI"] = pd.to_numeric(lai_mean["LAI"], errors="coerce")

# Remove rows where LAI could not be converted
lai_mean = lai_mean.dropna(subset=["LAI"])

# sort by date
lai_mean = lai_mean.sort_values("Date")

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

# ------ Step 4: Linear regression and data visualization ------

# Create a function to process the VIs for respective LAI datapoint, and create linear regression plot.
def lai_vi_regression(lai_mean, smoothed_vi_df, vi_col, window_days=15, output_folder=None, ax=None):
    """
    For each LAI measurement, calculate the average of a smoothed vegetation index
    within a time window, then perform linear regression and plot the results.

    Parameters
    ----------
    lai_mean : pandas DataFrame
        DataFrame containing LAI measurements with columns 'Date' and 'LAI'.
    smoothed_vi_df : pandas DataFrame
        DataFrame containing both raw and smoothed VI columns with a 'Date' column.
    vi_col : str
        Name of the raw vegetation index column (e.g. 'NIRv', 'NDVI').
        The smoothed column is automatically inferred as vi_col + '_smooth'.
    window_days : int, optional
        Number of days before and after each LAI date to include in the window. Default is 15.
    output_folder : str, optional
        Folder to save the output figure when plotted individually. If None, shown interactively.
    ax : matplotlib Axes, optional
        If provided, plot into this axes (for combined figure). If None, creates its own figure.

    Returns
    -------
    lai_vi : pandas DataFrame
        DataFrame with LAI and the window-averaged vegetation index values.
    slope : float
    intercept : float
    r_squared : float
    """

    # Infer smoothed column name automatically
    vi_smooth_col = f"{vi_col}_smooth"
    if vi_smooth_col not in smoothed_vi_df.columns:
        raise ValueError(f"Smoothed column '{vi_smooth_col}' not found. "
                         f"Available: {smoothed_vi_df.columns.tolist()}")

    # Work on a copy to avoid modifying the original
    lai_vi         = lai_mean.copy()
    lai_vi[vi_col] = np.nan

    # For each LAI measurement, calculate the average VI within the time window
    for i, row in lai_vi.iterrows():
        lai_date     = row["Date"]
        window_start = lai_date - pd.Timedelta(days=window_days)
        window_end   = lai_date + pd.Timedelta(days=window_days)
        vi_window    = smoothed_vi_df[
            (smoothed_vi_df["Date"] >= window_start) &
            (smoothed_vi_df["Date"] <= window_end)
        ]
        if not vi_window.empty:
            lai_vi.at[i, vi_col] = vi_window[vi_smooth_col].mean()

    # Remove rows where VI could not be calculated
    lai_vi = lai_vi.dropna(subset=[vi_col])

    # Linear regression
    slope, intercept, r_value, p_value, std_err = linregress(lai_vi[vi_col], lai_vi['LAI'])
    r_squared = r_value ** 2
    print(f"LAI = {slope:.2f} * {vi_col} + {intercept:.2f}, R² = {r_squared:.2f}")

    # Plot into provided ax or create a new figure
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(8, 6))

    x_fit = np.linspace(lai_vi[vi_col].min(), lai_vi[vi_col].max(), 100)
    y_fit = slope * x_fit + intercept

    ax.scatter(lai_vi[vi_col], lai_vi['LAI'], label='Data Points', zorder=3)
    ax.plot(x_fit, y_fit, linestyle='--',
            label=f'LAI = {slope:.2f} * {vi_col} + {intercept:.2f}')
    ax.text(0.05, 0.95, f'R² = {r_squared:.2f}',
            transform=ax.transAxes, fontsize=14, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.5))
    ax.set_xlabel(f'{vi_col} ({window_days}-day window avg)', fontsize=14)
    ax.set_ylabel('LAI', fontsize=14)
    ax.set_title(f'LAI vs {vi_col}', fontsize=15)
    ax.legend(fontsize=12)
    ax.tick_params(axis='both', labelsize=13)

    if standalone:
        plt.tight_layout()
        if output_folder:
            plt.savefig(os.path.join(output_folder, f'LAI_vs_{vi_col}_regression.png'))
        else:
            plt.show()

    return lai_vi, slope, intercept, r_squared

vi_list = ['NDVI', 'NIRv', 'MSAVI', 'EVI', 'NDRE', 'CIredge', 'MTCI', 'CRSWIR']

# Create a grid of subplots — 4 rows x 2 columns for 8 (7) VIs
n_cols = 2
n_rows = math.ceil(len(vi_list) / n_cols)
fig, axes = plt.subplots(n_rows, n_cols, figsize=(20, 25))
axes = axes.flatten()

regression_results = []
lai_vi_results = {}  # <-- add this to store matched dataframes for later use in time series plotting

for i, vi in enumerate(vi_list):
    lai_vi, slope, intercept, r2 = lai_vi_regression(
        lai_mean, VIs_data, vi_col=vi,
        window_days=15, ax=axes[i]
    )
    regression_results.append({'VI': vi, 'slope': slope, 'intercept': intercept, 'R2': r2})
    lai_vi_results[vi] = lai_vi  # <-- store the matched dataframe

# Hide any unused subplots
for j in range(len(vi_list), len(axes)):
    axes[j].set_visible(False)

plt.suptitle('LAI vs Vegetation Indices (15-day window)', fontsize=18, fontweight='bold', y=1.02)
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'LAI_vs_all_VI_regression.png'), bbox_inches='tight')
plt.show()

# Summary table
regression_results_df = pd.DataFrame(regression_results)
print(regression_results_df)
# save the regression results to a csv file
regression_results_df.to_csv(os.path.join(output_folder, 'LAI_VI_regression_results.csv'), index=False)

# ------ Step 5: Time series of predicted LAI vs observed LAI (full smoothed VI) ------
top3_vis = regression_results_df.sort_values('R2', ascending=False)['VI'].head(3).tolist()
plot_vis = top3_vis + ['NDVI'] if 'NDVI' not in top3_vis else top3_vis

fig, axes = plt.subplots(2, 2, figsize=(18, 12), sharey=True)
axes = axes.flatten()

for ax, vi in zip(axes, plot_vis):
    slope     = regression_results_df.loc[regression_results_df['VI'] == vi, 'slope'].values[0]
    intercept = regression_results_df.loc[regression_results_df['VI'] == vi, 'intercept'].values[0]

    # Use the full smoothed VI time series to predict LAI continuously
    vi_smooth_col = f'{vi}_smooth'
    full_series = VIs_data[['Date', vi_smooth_col]].dropna().sort_values('Date').copy()
    full_series['LAI_predicted'] = slope * full_series[vi_smooth_col] + intercept

    # Observed LAI points (from the matched dataframe)
    lai_vi = lai_vi_results[vi].copy().sort_values('Date')

    # Label NDVI as reference panel
    title_suffix = ' (reference)' if vi == 'NDVI' and vi not in top3_vis else ''

    # Plot continuous predicted LAI
    ax.plot(full_series['Date'], full_series['LAI_predicted'],
            color='steelblue', linewidth=1.5, label=f'Predicted LAI ({vi})')

    # Overlay observed LAI scatter points
    ax.scatter(lai_vi['Date'], lai_vi['LAI'],
               color='darkorange', zorder=5, s=50, label='Observed LAI')

    ax.set_title(f'{vi}{title_suffix}', fontsize=14, fontweight='bold')
    ax.set_xlabel('Date', fontsize=13)
    ax.set_ylabel('LAI', fontsize=13)
    ax.tick_params(axis='both', labelsize=11)
    ax.tick_params(axis='x', rotation=30)
    ax.legend(fontsize=12)
    ax.grid(alpha=0.4)

plt.suptitle('Time Series of Predicted and Observed LAI (2020–2025)',
             fontsize=16, fontweight='bold')
plt.tight_layout()
plt.savefig(os.path.join(output_folder, 'TimeSeries_Predicted_Observed_LAI_FullSeries.png'),
            bbox_inches='tight')
plt.show()