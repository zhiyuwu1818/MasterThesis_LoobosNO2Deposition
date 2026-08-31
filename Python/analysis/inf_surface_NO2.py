def infer_surface_NO2(NO2_VCD, SCM_prof_NO2, TROPOMI_NO2VCD):
    """
    Infer daily tropospheric NO2 surface layer concentrations by scaling
    SCM surface layer values with TROPOMI VCD observations.

    Parameters
    ----------
    NO2_VCD : pandas DataFrame
        SCM-derived NO2 VCD with columns 'datetime' and 'Tropo_NO2_VCD'.
    SCM_prof_NO2 : pandas DataFrame
        SCM NO2 vertical profiles in wide format with 'datetime' and height level columns
        ordered from lowest to highest (e.g. '10m', '20m', ...).
    TROPOMI_NO2VCD : pandas DataFrame
        TROPOMI overpass data with columns 'datetime' and 'Tropospheric_NO2_value'.

    Returns
    -------
    inferred_SL_NO2 : pandas DataFrame
        DataFrame with columns:
            'Date'       : date of valid TROPOMI overpass
            'NO2_SL'     : inferred NO2 surface layer concentration
            'SCM_NO2VCD' : daily mean SCM NO2 VCD used as reference
    Created: 7 May 2026 Zhiyu Wu together with Copilot AI and Claude
    copyright reserved
    """
    import pandas as pd

   # Compute daily mean VCD from SCM data
    NO2_VCD['Date']  = pd.to_datetime(NO2_VCD['datetime']).dt.date
    daily_mean_VCD   = NO2_VCD.set_index('Date')['Tropo_NO2_VCD']

    # Compute daily mean vertical profile, explicitly exclude non-height columns
    SCM_prof_NO2['Date'] = pd.to_datetime(SCM_prof_NO2['datetime']).dt.date
    height_cols          = [c for c in SCM_prof_NO2.columns if c.startswith('height_')]
    daily_mean_profile   = SCM_prof_NO2.groupby('Date')[height_cols].mean()
    surface_col          = daily_mean_profile.columns[2]
    print(f"Using '{surface_col}' as the surface layer for NO2_SL_model.")

    # Restrict to the date range covered by SCM VCD data
    scm_start = daily_mean_VCD.index.min()
    scm_end   = daily_mean_VCD.index.max()

    # Get dates with valid TROPOMI overpasses within the SCM date range
    TROPOMI_NO2VCD['Date'] = pd.to_datetime(TROPOMI_NO2VCD['datetime']).dt.date
    valid_tropomi_dates    = TROPOMI_NO2VCD[
        (TROPOMI_NO2VCD['Date'] >= scm_start) &
        (TROPOMI_NO2VCD['Date'] <= scm_end)
    ]['Date'].unique()
    print(f"Found {len(valid_tropomi_dates)} days with valid TROPOMI overpasses within SCM range ({scm_start} to {scm_end}).")

    # Loop only over valid TROPOMI dates within the SCM range
    results = []

    for date in valid_tropomi_dates:
        date_str = str(date)

        # Skip if no SCM profile available for this date
        if date not in daily_mean_profile.index:
            print(f"No SCM profile found for {date_str}, skipping.")
            continue

        scm_no2vcd   = float(daily_mean_VCD[date])
        NO2_SL_model = float(daily_mean_profile.loc[date, surface_col])
        tropomi_day  = TROPOMI_NO2VCD[TROPOMI_NO2VCD['Date'] == date]

        # Infer NO2 surface concentration using daily mean SCM VCD as reference
        NO2 = (tropomi_day['Tropospheric_NO2_value'].values[0] / scm_no2vcd) * NO2_SL_model

        results.append({
            'Date'       : date,
            'NO2_SL'     : NO2,
            'SCM_NO2VCD' : scm_no2vcd
        })

    inferred_SL_NO2 = pd.DataFrame(results).sort_values('Date').reset_index(drop=True)
    print(f"Completed: {len(inferred_SL_NO2)} days processed.")
    print(inferred_SL_NO2)

    return inferred_SL_NO2