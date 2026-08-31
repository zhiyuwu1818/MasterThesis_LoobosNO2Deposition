import numpy as np
import matplotlib.pyplot as plt


def pollution_rose(
    df,
    pollutant_col,
    dir_col,
    ws_col,
    *,
    pollutant_name='Pollutant',
    pollutant_unit='µg m⁻³',
    bins=np.array([1, 5, 10, 15, np.inf]),
    dir_binwidth=30,
    calm_wind_cut=0.2,
    max_pct=25,
    title=None,
    figsize=(6, 8)
):
    """
    Plot a pollution rose (directional concentration distribution).

    Parameters
    ----------
    df : pandas.DataFrame
        Input dataframe containing pollutant, wind direction, and wind speed.
    pollutant_col : str
        Column name of pollutant concentration.
    dir_col : str
        Column name of wind direction (degrees).
    ws_col : str
        Column name of wind speed (m s⁻¹).
    pollutant_name : str, optional
        Name of the pollutant (e.g. 'NO₂').
    pollutant_unit : str, optional
        Unit of the pollutant (e.g. 'ppb', 'µg m⁻³').
    bins : array-like, optional
        Upper bounds of pollutant bins.
    dir_binwidth : int, optional
        Wind direction bin width in degrees.
    calm_wind_cut : float, optional
        Wind speed threshold below which data are excluded.
    max_pct : float, optional
        Maximum radial percentage for plot scaling.
    title : str, optional
        Plot title.
    figsize : tuple, optional
        Figure size.

    Returns
    -------
    matplotlib.axes._subplots.PolarAxesSubplot
        Polar axes of the pollution rose.
    """

    # ----------------------------
    # Remove calm winds
    # ----------------------------
    df_valid = df[df[ws_col] > calm_wind_cut].copy()

    # Extract values
    wd = df_valid[dir_col].values
    poll = df_valid[pollutant_col].values

    # ----------------------------
    # Direction binning
    # ----------------------------
    width = dir_binwidth
    half = width / 2

    theta = np.radians(
        (np.floor((wd + half) / width) * width) % 360
    )

    # ----------------------------
    # Pollutant bins
    # ----------------------------
    upper_bounds = bins
    lower_bounds = np.r_[0, bins[:-1]]

    colors = plt.cm.jet(np.linspace(0, 1, len(bins)))

    # ----------------------------
    # Plot
    # ----------------------------
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, polar=True)

    ax.set_theta_offset(np.radians(90))
    ax.set_theta_direction(-1)

    unique_theta = np.unique(theta)
    n_total = len(theta)

    for ang in unique_theta:
        mask_angle = theta == ang
        p_vals = poll[mask_angle]

        bottom = 0
        for lb, ub, color in zip(lower_bounds, upper_bounds, colors):
            count = np.sum((p_vals >= lb) & (p_vals < ub))
            pct = count / n_total * 100

            ax.bar(
                ang,
                pct,
                width=np.radians(width),
                bottom=bottom,
                color=color,
                edgecolor='k'
            )
            bottom += pct

    # ----------------------------
    # Radial axis
    # ----------------------------
    ticks = np.linspace(0, max_pct, 5)[1:]
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{int(t)}%" for t in ticks])
    ax.set_rmax(max_pct)

    # ----------------------------
    # Title
    # ----------------------------
    if title is None:
        title = f'{pollutant_name} Pollution Rose'

    ax.set_title(title, fontsize=12, y=1.08)

    # ----------------------------
    # Legend (unit-aware)
    # ----------------------------
    labels = [
        f"{lb}–{ub} {pollutant_unit}"
        if ub != np.inf
        else f">{lb} {pollutant_unit}"
        for lb, ub in zip(lower_bounds, upper_bounds)
    ]

    patches = [
        plt.Rectangle((0, 0), 1, 1, color=c)
        for c in colors
    ]

    ax.legend(
        patches,
        labels,
        loc="lower center",
        bbox_to_anchor=(0.5, -0.3),
        ncol=2,
        frameon=False
    )

    plt.tight_layout()
    return ax
