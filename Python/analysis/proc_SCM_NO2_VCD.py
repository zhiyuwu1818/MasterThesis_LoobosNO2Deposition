def process_SCM_NO2_VCD(scm_file_path, tropopause_height_m=12000):
    """
    Read a SCM output file containing physical properties and NO2 vertical profiles,
    compute layer thicknesses from pressure, and calculate the NO2 vertical column
    density (VCD) for each timestep up to the tropopause.

    Parameters
    ----------
    scm_file_path : str
        Path to the SCM .out file (e.g., 'Data/Model/NO2profiles11-14_2025.out')
    tropopause_height_m : float, optional
        Height in meters below which layers are included in the VCD integration.
        Default is 12000 m.

    Returns
    -------
    NO2_VCD : pandas DataFrame
        DataFrame with columns:
            'Datetime'                      : timestep
            'NO2_VCD_1e15_molecules_per_cm2': tropospheric NO2 VCD per timestep
    SCM_prof_NO2 : pandas DataFrame
        Wide-format DataFrame with 'Datetime' and height level columns
        (e.g. 'height_10m', 'height_20m', ...) containing NO2 values in ppb.

    Notes
    -----
    The VCD is calculated as:
        VCD = sum( no2_ppb * 1e-9 * Rhoa * 1e3 / Zmair * Zavo * dz * 1e-4 * 1e-15 )
    where dz is derived from half-level pressure reconstruction using the
    hydrostatic equation.
    Created: 6 May 2026 Zhiyu Wu together with Copilot AI and Claude
    copyright reserved
    """
    import pandas as pd
    import numpy as np

    G      = 9.81      # gravitational acceleration [m/s²]
    Zavo   = 6.022e23  # Avogadro's number [molecules/mol]
    Zmair  = 28.97     # Molar mass of air [g/mol]

    # ------------------------------------------------------------------ #
    # 1. Read file: physical properties (first 3 rows) + NO2 profiles
    # ------------------------------------------------------------------ #
    with open(scm_file_path, "r") as f:
        physical_properties = [f.readline().strip() for _ in range(3)]
        raw_lines = [physical_properties[2]] + [line.strip() for line in f]

    # ------------------------------------------------------------------ #
    # 2. Parse physical properties into a DataFrame
    # ------------------------------------------------------------------ #
    rows  = {}
    units = {}

    for prop in physical_properties:
        tokens = prop.split()
        name   = tokens[0]

        if tokens[1].startswith('['):
            if not tokens[1].endswith(']'):
                unit        = tokens[1] + ' ' + tokens[2]
                value_start = 3
            else:
                unit        = tokens[1]
                value_start = 2
        else:
            unit        = ''
            value_start = 1

        units[name]  = unit
        rows[name]   = [float(v) for v in tokens[value_start:]]

    df_phys = pd.DataFrame(rows).rename(columns={'Datetime': 'Height_m'})

    # ------------------------------------------------------------------ #
    # 3. Compute layer thicknesses from half-level pressure reconstruction
    # ------------------------------------------------------------------ #
    press  = df_phys['Press'].values
    rhoa   = df_phys['Rhoa'].values
    height = df_phys['Height_m'].values
    n      = len(press)

    p_half    = np.zeros(n + 1)
    p_half[0] = press[0] + rhoa[0] * G * height[0]        # surface anchor
    for i in range(1, n):
        p_half[i] = (press[i - 1] + press[i]) / 2         # midpoint average
    p_half[n] = press[-1] - rhoa[-1] * G * (height[-1] - height[-2]) / 2  # top

    dz = (p_half[:-1] - p_half[1:]) / (rhoa * G)          # hydrostatic dz

    df_phys['p_half_lower']       = p_half[:-1]
    df_phys['p_half_upper']       = p_half[1:]
    df_phys['dz_m_from_pressures'] = dz
    df_phys['accumulated_dz_m']   = np.cumsum(dz)

    # ------------------------------------------------------------------ #
    # 4. Parse NO2 vertical profiles into long-format DataFrame
    # ------------------------------------------------------------------ #
    split_lines  = [line.split() for line in raw_lines]
    merged_lines = [[r[0] + " " + r[1]] + r[2:] for r in split_lines]

    header       = ["Datetime"] + merged_lines[0]
    df_no2_wide  = pd.DataFrame(merged_lines[1:], columns=header)
    df_no2_wide.columns = ["Datetime"] + [f"height_{c}m" for c in df_no2_wide.columns[1:]]
    df_no2_wide["Datetime"] = pd.to_datetime(df_no2_wide["Datetime"], format="%d-%m-%y %H:%M")
    df_no2_wide = df_no2_wide.rename(columns={'Datetime': 'datetime'})

    height_cols = [c for c in df_no2_wide.columns if c.startswith('height_')]
    df_no2_wide[height_cols] = df_no2_wide[height_cols].apply(pd.to_numeric, errors='coerce')

    df_no2_long = (
        df_no2_wide
        .melt(id_vars='datetime', var_name='height', value_name='no2_ppb')
        .assign(height_m=lambda x: x['height'].str.extract(r"(\d+)").astype(float))
        .drop(columns=['height'])
        .sort_values(['datetime', 'height_m'])
    )

    # ------------------------------------------------------------------ #
    # 5. Merge NO2 profiles with physical properties on height
    # ------------------------------------------------------------------ #
    merged = pd.merge_asof(
        df_no2_long.sort_values('height_m'),
        df_phys[['Height_m', 'Rhoa', 'dz_m_from_pressures']].sort_values('Height_m'),
        left_on='height_m',
        right_on='Height_m',
        direction='nearest',
        tolerance=1.0
    )[['datetime', 'height_m', 'no2_ppb', 'Rhoa', 'dz_m_from_pressures']]

    # ------------------------------------------------------------------ #
    # 6. Filter to troposphere and calculate VCD per timestep
    # ------------------------------------------------------------------ #
    merged_tropo = merged[merged['height_m'] <= tropopause_height_m]

    results = []
    for dt, df_dt in merged_tropo.groupby('datetime'):
        VCD = np.sum(
            df_dt['no2_ppb']  * 1e-9          # ppb → mol/mol
            * df_dt['Rhoa']   * 1e3            # kg/m³ → g/m³
            / Zmair                            # g/mol → mol/m³
            * Zavo                             # mol → molecules/m³
            * df_dt['dz_m_from_pressures']     # → molecules/m²
            * 1e-4                             # /m² → /cm²
            * 1e-15                            # → 10¹⁵ molecules/cm²
        )
        results.append({'datetime': dt, 'NO2_VCD_1e15_molecules_per_cm2': VCD})

    NO2_VCD = pd.DataFrame(results).sort_values('datetime').reset_index(drop=True)
    NO2_VCD = NO2_VCD.rename(columns={'datetime': 'datetime'})
    print(f"Processed {len(NO2_VCD)} timesteps from {NO2_VCD['datetime'].min()} to {NO2_VCD['datetime'].max()}")
    print(NO2_VCD.head())

    return NO2_VCD, df_no2_wide