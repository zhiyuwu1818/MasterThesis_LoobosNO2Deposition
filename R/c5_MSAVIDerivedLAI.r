# ----------------------------------------------------------------------------
# C4: Seasonal and full-series LAI derivation from Sentinel-2 MSAVI
#
# Description:
#   Loads seasonal (winter, spring, summer, autumn) and full Sentinel-2 MSAVI
#   raster stacks (2020-2025) for the Loobos footprint area. Filters out
#   layers/dates with a spatial mean MSAVI outside the range [0.25, 0.4] to
#   exclude invalid or noisy observations. Computes the temporal mean MSAVI
#   per pixel and derives Leaf Area Index (LAI) using the empirical
#   regression LAI = 12.484 * MSAVI - 1.753. Exports LAI rasters (per season
#   and full time series), generates MSAVI time series plots colour-coded
#   by season, and saves the extracted time series as CSV files.

# ----- Step 1: preparation -----

# install and load packages
# install.packages("terra")
# install.packages("ggplot2")
# install.packages("tibble")
library(terra)
library(ggplot2)
library(tibble)

# Helper function to filter layers by their spatial mean
filter_layers_by_mean <- function(raster_stack, min_val = 0.25, max_val = 0.4) {
  # Calculate the global mean for each layer/date
  layer_means <- global(raster_stack, mean, na.rm = TRUE)$mean
  # Determine which layers fall within the acceptable range
  valid_layers <- layer_means >= min_val & layer_means <= max_val
  # Return only the layers that match the condition
  return(raster_stack[[valid_layers]])
}

# ----- Step 2: load dataset, select layers -----
msavi_layers_winter = rast("Data/Sentinel2/MSAVI_winter_2020_2025.tif")
msavi_layers_spring = rast("Data/Sentinel2/MSAVI_spring_2020_2025.tif")
msavi_layers_summer = rast("Data/Sentinel2/MSAVI_summer_2020_2025.tif")
msavi_layers_autumn = rast("Data/Sentinel2/MSAVI_autumn_2020_2025.tif")

# **NEW FILTER METHOD**: Drop entire days if their spatial average is outside [0.25, 0.4]
msavi_layers_winter <- filter_layers_by_mean(msavi_layers_winter)
msavi_layers_spring <- filter_layers_by_mean(msavi_layers_spring)
msavi_layers_summer <- filter_layers_by_mean(msavi_layers_summer)
msavi_layers_autumn <- filter_layers_by_mean(msavi_layers_autumn)


# ----- Step 3: calculate average VI values for each pixel across the filtered time series -----
mean_msavi_winter <- app(msavi_layers_winter, mean, na.rm = TRUE)
mean_msavi_spring <- app(msavi_layers_spring, mean, na.rm = TRUE)
mean_msavi_summer <- app(msavi_layers_summer, mean, na.rm = TRUE)
mean_msavi_autumn <- app(msavi_layers_autumn, mean, na.rm = TRUE)


# ----- Step 4: derive LAI values for each pixel -----
# MSAVI: calculate by the equation: LAI = 12.484 * MSAVI - 1.753
lai_msavi_winter <- 12.484 * mean_msavi_winter - 1.753
lai_msavi_spring <- 12.484 * mean_msavi_spring - 1.753
lai_msavi_summer <- 12.484 * mean_msavi_summer - 1.753
lai_msavi_autumn <- 12.484 * mean_msavi_autumn - 1.753

# ----- Step 5: export the results -----
writeRaster(lai_msavi_winter, "Data/Sentinel2/LAI_MSAVI_mean_2020_2025_RDNew_winter.tif", overwrite = TRUE)
writeRaster(lai_msavi_spring, "Data/Sentinel2/LAI_MSAVI_mean_2020_2025_RDNew_spring.tif", overwrite = TRUE)
writeRaster(lai_msavi_summer, "Data/Sentinel2/LAI_MSAVI_mean_2020_2025_RDNew_summer.tif", overwrite = TRUE)
writeRaster(lai_msavi_autumn, "Data/Sentinel2/LAI_MSAVI_mean_2020_2025_RDNew_autumn.tif", overwrite = TRUE)


# ----- Step 6: calculate mean MSAVI per layer (each layer = one date) -----
extract_ts <- function(r, season_label = NULL) {
  ts <- global(r, mean, na.rm = TRUE)
  df <- data.frame(
    Date  = as.Date(rownames(ts), format = "%Y%m%d"),
    MSAVI = ts$mean
  )
  if (!is.null(season_label)) {
    df$Season <- season_label
  }
  return(df)
}

ts_winter <- extract_ts(msavi_layers_winter, "Winter")
ts_spring <- extract_ts(msavi_layers_spring, "Spring")
ts_summer <- extract_ts(msavi_layers_summer, "Summer")
ts_autumn <- extract_ts(msavi_layers_autumn, "Autumn")

# ----- Step 7: merge all seasons and sort by date -----
ts_all <- rbind(ts_winter, ts_spring, ts_summer, ts_autumn)
ts_all <- ts_all[order(ts_all$Date), ]


# ----- Step 8: plot MSAVI time series, coloured by season -----
season_colors <- c("Winter" = "steelblue", "Spring" = "mediumseagreen",
                   "Summer" = "goldenrod",  "Autumn" = "sienna")

ggplot(ts_all, aes(x = Date, y = MSAVI, color = Season)) +
  geom_line() +
  geom_point(size = 1.5) +
  scale_color_manual(values = season_colors) +
  labs(title = "Mean MSAVI Time Series by Season (2020–2025)",
       x = "Date", y = "Mean MSAVI", color = "Season") +
  theme_minimal()


# ----- Step 9: export merged time series to CSV -----
write.csv(ts_all, "Data/Sentinel2/MSAVI_TimeSeries_seasonal_2020_2025.csv", row.names = FALSE)


# ----- Step 10: load MSAVI stacks for the same footprint, calculate the LAI results -----
VI_layers = rast("Data/Sentinel2/VIsforLAI_10m_2020_2025_RDNew.tif")
msavi_layers = VI_layers[[grep("MSAVI", names(VI_layers))]]

# **NEW FILTER METHOD**: Filter the full stack by image average
msavi_layers <- filter_layers_by_mean(msavi_layers)

# derive LAI for each pixel
lai_msavi_layers <- 12.484 * msavi_layers - 1.753

# save LAI layers
writeRaster(lai_msavi_layers, "Data/Sentinel2/LAI_MSAVI_10m_2020_2025_RDNew.tif", overwrite = TRUE)

# Process final time-series data frame
plot_msavi <- extract_ts(msavi_layers)

# Plot time-series of MSAVI (all points will naturally be between 0.25 and 0.4 now)
ggplot(plot_msavi, aes(x = Date, y = MSAVI)) +
  geom_line(color = "blue") +
  geom_point(color = "red") +
  labs(title = "Mean MSAVI Time Series (Filtered)", x = "Date", y = "Mean MSAVI") +
  theme_minimal()

# Save to csv
write.csv(plot_msavi, "Data/Sentinel2/Mean_MSAVI_TimeSeries_2020_2025_RDNew.csv", row.names = FALSE)
