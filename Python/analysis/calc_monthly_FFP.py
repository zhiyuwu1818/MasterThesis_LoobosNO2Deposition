def cal_monthly_FFP(year_month, zmt, ht, rs_new, insitu_data):
    """
    Calculate the monthly flux footprint climatology for the period of 9:00-17:00 of each day based on Kljun et al. (2015) with inputs of zm (measurement height above displacement height)
    and h (boundary layer height) and in-situ measurements for a given month and output the results.
    Input
    -----
    year_month: str
        Year and month in 'YYYY-MM' format, e.g., '2023-01  ' for January 2023
    zmt: float
        Measurement height above displacement height (m)
    ht: float
        Boundary layer height (m)
    rs_new: list of float     
    Percentage of source area for which to provide contours, must be between 10% and 90%. 
                       Can be either a single value (e.g., "80") or a list of values (e.g., "[10, 20, 30]")
                       Expressed either in percentages ("80") or as fractions of 1 ("0.8"). 
                       Default is [10:10:80]. Set to "None" for no output of percentages
    insitu_data: pandas DataFrame
        DataFrame containing in-situ measurements with columns: 'Time', 'Timestamp', 'WS_2_1_1', 'u*', 'L', 'WD_2_1_1'
    Output
    ------
    FFP: dict
        Dictionary containing the footprint climatology results
    References
    ----------
    Kljun, N., Calanca, P., Rotach, M. W., & Schmid, H. P. (2015). A simple parameterisation for flux footprint predictions
    https://doi.org/10.5194/bg-12-3691-2015
    Notes
    -----   
    Created: 4 Feb 2026 Zhiyu Wu together with Copilot AI
    copyright reserved
    """
    from . import calc_footprint_FFP_climatology as myfootprint_climatology
    # Filter data for the specified month and time range (9:00-17:00)  
    data_day = insitu_data[
        (insitu_data['Time'].str.startswith(year_month)) &
        (insitu_data['Timestamp'] >= '09:00:00') &
        (insitu_data['Timestamp'] <= '17:00:00')
    ]
    # Remove rows with missing or invalid data
    data_day = data_day[
        (data_day['WS_2_1_1'] > -9999) &
        (data_day['u*'] > -9999) &
        (data_day['L'] > -9999) &
        (data_day['WD_2_1_1'] > -9999)
    ]
    # Number of data points
    N = len(data_day)
    print(f"Number of data points for footprint calculation in {year_month}: {N}")
    # Convert constant value inputs to lists of length N
    zm_new = [zmt] * N
    h_new = [ht] * N
    # Convert other changing in-situ measurement columns to lists
    umean_new = data_day['WS_2_1_1'].tolist()
    ustar_new = data_day['u*'].tolist()
    sigmav_new = [1.3 * u for u in ustar_new]
    ol_new = data_day['L'].tolist()
    wind_dir_new = data_day['WD_2_1_1'].tolist()
    # Calculate the footprint climatologyß
    print(f"Calculating footprint for {year_month}...")
    FFP = myfootprint_climatology.FFP_climatology (zm = zm_new, z0 = None, umean = umean_new, h=h_new, ol=ol_new, sigmav=sigmav_new,
                                                   ustar=ustar_new, wind_dir=wind_dir_new, smooth_data=1, fig=True, rs=rs_new)
    print(f"Footprint calculation for {year_month} completed.")
    return FFP
