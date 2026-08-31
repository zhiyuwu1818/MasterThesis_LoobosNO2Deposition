def cal_period_FFP(start_date, end_date, zmt, ht, rs_new, insitu_data):
    """
    Calculate the flux footprint climatology for a set period, for the time range 9:00-17:00 of each day,
    based on Kljun et al. (2015) with inputs of zm (measurement height above displacement height)
    and h (boundary layer height) and in-situ measurements for a given period and output the results.
    Input
    -----
    start_date: str
        Start date in 'YYYY-MM-DD' format, e.g., '2023-01-01' for 1 January 2023
    end_date: str
        End date in 'YYYY-MM-DD' format, e.g., '2023-03-31' for 31 March 2023 (inclusive)
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
    Modified: 29 Apr 2026 — generalised from monthly to arbitrary date range by Claude
    copyright reserved
    """
    from . import calc_footprint_FFP_climatology as myfootprint_climatology
    import pandas as pd

    # Parse start and end dates
    start_dt = pd.Timestamp(start_date)
    end_dt   = pd.Timestamp(end_date)

    # Convert the 'Time' column to datetime for range filtering (coerce errors to NaT)
    time_as_dt = pd.to_datetime(insitu_data['Time'], errors='coerce')

    # Filter data for the specified date range and time window (9:00–17:00)
    data_day = insitu_data[
        (time_as_dt >= start_dt) &
        (time_as_dt <= end_dt) &
        (insitu_data['Timestamp'] >= '09:00:00') &
        (insitu_data['Timestamp'] <= '17:00:00')
    ]

    # Remove rows with missing or invalid data
    data_day = data_day[
        (data_day['WS_2_1_1'] > -9999) &
        (data_day['u*']        > -9999) &
        (data_day['L']         > -9999) &
        (data_day['WD_2_1_1'] > -9999)
    ]

    # Number of data points
    N = len(data_day)
    print(f"Number of data points for footprint calculation "
          f"({start_date} to {end_date}): {N}")

    if N == 0:
        raise ValueError(
            f"No valid data found between {start_date} and {end_date}. "
            "Check your date range and input data."
        )

    # Convert constant value inputs to lists of length N
    zm_new = [zmt] * N
    h_new  = [ht]  * N

    # Convert in-situ measurement columns to lists
    umean_new    = data_day['WS_2_1_1'].tolist()
    ustar_new    = data_day['u*'].tolist()
    sigmav_new   = [1.3 * u for u in ustar_new]
    ol_new       = data_day['L'].tolist()
    wind_dir_new = data_day['WD_2_1_1'].tolist()

    # Calculate the footprint climatology
    print(f"Calculating footprint for {start_date} to {end_date}...")
    FFP = myfootprint_climatology.FFP_climatology(
        zm        = zm_new,
        z0        = None,
        umean     = umean_new,
        h         = h_new,
        ol        = ol_new,
        sigmav    = sigmav_new,
        ustar     = ustar_new,
        wind_dir  = wind_dir_new,
        smooth_data = 1,
        fig       = True,
        rs        = rs_new
    )
    print(f"Footprint calculation for {start_date} to {end_date} completed.")
    return FFP