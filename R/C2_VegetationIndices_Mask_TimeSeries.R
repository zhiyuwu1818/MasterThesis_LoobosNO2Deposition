# C2: Footprint-Based Multi-Index Vegetation Processing and Time Series Extraction
#
# This script performs footprint-level vegetation monitoring using
# multiple Sentinel-2 spectral vegetation indices.
#
# The workflow includes:
# 1. Computing multi-temporal mean vegetation NDVI
# 2. Deriving a vegetation mask based on a percentile threshold
#    (10th percentile of NDVI) to isolate persistent vegetation
# 3. Filtering physically unrealistic and low-quality index values for each index with values
# 4. Extracting spatially averaged time series for each index within the defined research footprint, applying the previously derived vegetation mask.
# 5. Exporting cleaned multi-index time series for subsequent
#    atmospheric and flux comparison analyses
#
# Indices included:
# - NDVI   : Normalized Difference Vegetation Index
# - kNDVI  : Kernel NDVI (nonlinear vegetation response)
# - CR-SWIR: Continuum-Removed SWIR Index (moisture / stress sensitivity)
# - CIre   : Red-Edge Chlorophyll Index
# - REP    : Red-Edge Position (4-Point Linear Interpolation)
#
# Author: Zhiyu Wu
# Date: 04/03/2026

# ----- Step 1: preparation--------

# install and load packages
install.packages("terra")
install.packages("ggplot2")
install.packages("tibble")
library(terra)
library(ggplot2)
library(tibble)

# ----- Step 2: Select different layers for vegetation bands--------
# Load the dataset, r10 is the one with 10 m resolution, r20 is the one with 20 m resolution
r10 <- rast("Data/Sentinel2/MultiIndex_10m_Daily_2020_2025.tif")
r20 <- rast("Data/Sentinel2/MultiIndex_20m_Daily_2020_2025.tif")
# Select bands containing vegetation indices information
ndvi_layers <- r10[[grep("_NDVI$", names(r10))]]
kndvi_layers <- r10[[grep("_kNDVI$", names(r10))]]
crswir_layers <- r20[[grep("_CRSWIR$", names(r20))]]
cire_layers <- r20[[grep("_CIre$", names(r20))]]
rep_layers <- r20[[grep("_REP$", names(r20))]]

# ----- Step 3: calculate mean NDVI--------
mean_raster <- app(ndvi_layers, mean, na.rm = TRUE)

writeRaster(
  mean_raster,
  "Data/Sentinel2/NDVI_Mean_RDNew.tif",
  overwrite = TRUE
)

# ----- Step 4: Create 10th percentile mask--------

# Compute the 10th Percentile Threshold
threshold <- global(mean_raster, quantile, probs = 0.10, na.rm = TRUE)[1,1]

# Create a binary raster as a mask
veg_mask_10 <- mean_raster > threshold

# Check the histogram
hist(values(mean_raster), main="Mean NDVI Distribution")
abline(v = threshold, col="red", lwd=2)

# Save mask as GeoTIFF
writeRaster(
  veg_mask_10,
  "Data/Sentinel2/NDVI_Mask_10percent_RDNew.tif",
  overwrite = TRUE
)

# ----- Step 5: Apply mask to the original NDVI stacks, calculate average NDVI--------
# Filter out values < 0.2 across the entire stack 'ndvi_layers'
# This replaces any value below 0.2 with NA
ndvi_filtered <- ndvi_layers
ndvi_filtered[ndvi_filtered < 0.2] <- NA

# Apply 10th percentile mask on top of the filtered stack
# We ensure the mask uses NA for the 'False' areas
veg_mask_final <- ifel(veg_mask_10 == 0, NA, 1)
masked_ndvi <- mask(ndvi_filtered, veg_mask_final)

# Calculate the global mean for each layer
mean_NDVI_value <- global(masked_ndvi, mean, na.rm = TRUE)

# Check NAs values
print(mean_NDVI_value)

# ----- Step 6: Plot time series of average NDVI -----
#  Convert the results into a data frame and move row names to a column
plot_ndvi <- as.data.frame(mean_NDVI_value)
plot_ndvi$date_str <- rownames(plot_ndvi)

# Convert the string date (e.g., "20200118") to a proper Date object
# Adjust the format "%Y%m%d" if your strings look different
plot_ndvi$Date <- as.Date(plot_ndvi$date_str, format="%Y%m%d")

# Delete all the datapoints less than 0.6 or more than 0.9, as it is likely influenced by clouds etc.
plot_ndvi <- plot_ndvi[plot_ndvi$mean >= 0.6 & plot_ndvi$mean <= 0.9, ]

# Create the time-series plot
ggplot(plot_ndvi, aes(x = Date, y = mean)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(color = "forestgreen") +
  labs(title = "Average NDVI Time Series",
       x = "Date",
       y = "Mean NDVI") +
  theme_minimal()

# ----- Step 7: Apply mask to the kNDVI layer, calculate average kNDVI -----
# Filter out values < 0.04 across the entire stack 'kndvi_layers' (this value corresponds to ndvi<0.2)
# This replaces any value below 0.04 with NA
kndvi_filtered <- kndvi_layers
kndvi_filtered[kndvi_filtered < 0.04] <- NA

# Apply 10th percentile mask on top of the filtered stack
# We ensure the mask uses NA for the 'False' areas
masked_kndvi <- mask(kndvi_filtered, veg_mask_final)

# Calculate the global mean for each layer
mean_kNDVI_value <- global(masked_kndvi, mean, na.rm = TRUE)

# Check NAs values
print(mean_kNDVI_value)

# ----- Step 8: Plot time series of average kNDVI -----
#  Convert the results into a data frame and move row names to a column
plot_kndvi <- as.data.frame(mean_kNDVI_value)
plot_kndvi$date_str <- rownames(plot_kndvi)

# Convert the string date (e.g., "20200118") to a proper Date object
# Adjust the format "%Y%m%d" if your strings look different
plot_kndvi$Date <- as.Date(plot_kndvi$date_str, format="%Y%m%d")

# Delete all the datapoints less than 0.35 or more than 0.67 (corresponds to NDVI thresholds), as it is likely influenced by clouds etc.
plot_kndvi <- plot_kndvi[plot_kndvi$mean >= 0.35 & plot_kndvi$mean <= 0.67, ]

# Create the time-series plot
ggplot(plot_kndvi, aes(x = Date, y = mean)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(color = "forestgreen") +
  labs(title = "Average kNDVI Time Series",
       x = "Date",
       y = "Mean kNDVI") +
  theme_minimal()

# ----- Step 9: Apply mask to the CR-SWIR layer, calculate average CR-SWIR -----

# Replace values smaller than 0 and larger than 2 with NA
crswir_filtered <- crswir_layers
crswir_filtered[crswir_filtered < 0] <- NA
crswir_filtered[crswir_filtered > 2] <- NA

# Resample mask to the same resolution as crswir_layers
veg_mask_final_resampled <- resample(veg_mask_final, crswir_layers[[1]], method = "near")

# Apply 10th percentile mask on top of the filtered stack
# We ensure the mask uses NA for the 'False' areas
masked_crswir <- mask(crswir_filtered, veg_mask_final_resampled)

# Calculate the global mean for each layer
mean_CRSWIR_value <- global(masked_crswir, mean, na.rm = TRUE)

# Check NAs values
print(mean_CRSWIR_value)

# ----- Step 10: Plot time series of average CR-SWIR -----
plot_crswir <- as.data.frame(mean_CRSWIR_value)
plot_crswir$date_str <- rownames(plot_crswir)

# Convert the string date (e.g., "20200118") to a proper Date object
# Adjust the format "%Y%m%d" if your strings look different
plot_crswir$Date <- as.Date(plot_crswir$date_str, format="%Y%m%d")

# Create the time-series plot
ggplot(plot_crswir, aes(x = Date, y = mean)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(color = "forestgreen") +
  labs(title = "Average CRswir Time Series",
       x = "Date",
       y = "Mean CRswir value") +
  theme_minimal()

# ----- Step 11: Apply mask to the CIre layer, calculate average CIre -----

# Replace values smaller than 0 and larger than 10 with NA
cire_filtered <- cire_layers
cire_filtered[cire_filtered < 0] <- NA
cire_filtered[cire_filtered > 10] <- NA

# Apply 10th percentile mask on top of the filtered stack
# We ensure the mask uses NA for the 'False' areas
masked_cire <- mask(cire_filtered, veg_mask_final_resampled)

# Calculate the global mean for each layer
mean_CIRE_value <- global(masked_cire, mean, na.rm = TRUE)

# Check NAs values
print(mean_CIRE_value)

# ----- Step 12: Plot time series of average CIre -----
plot_cire <- as.data.frame(mean_CIRE_value)
plot_cire$date_str <- rownames(plot_cire)

# Convert the string date (e.g., "20200118") to a proper Date object
# Adjust the format "%Y%m%d" if your strings look different
plot_cire$Date <- as.Date(plot_cire$date_str, format="%Y%m%d")

# Delete all the datapoints less than 0 or more than 4, as it is likely influenced by clouds etc.
plot_cire <- plot_cire[plot_cire$mean >= 0 & plot_cire$mean <= 4, ]

# Create the time-series plot
ggplot(plot_cire, aes(x = Date, y = mean)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(color = "forestgreen") +
  labs(title = "Average CIre Time Series",
       x = "Date",
       y = "Mean CIre") +
  theme_minimal()

# ----- Step 13: Apply mask to the REP layer, calculate average REP -----
# Replace values smaller than 600 and larger than 1000 with NA
rep_filtered <- rep_layers
rep_filtered[rep_filtered < 600] <- NA
rep_filtered[rep_filtered > 1000] <- NA

# Apply 10th percentile mask on top of the filtered stack
# We ensure the mask uses NA for the 'False' areas
masked_rep <- mask(rep_filtered, veg_mask_final_resampled)

# Calculate the global mean for each layer
mean_REP_value <- global(masked_rep, mean, na.rm = TRUE)

# Check NAs values
print(mean_REP_value)

# ----- Step 14: Plot time series of average REP -----
plot_rep <- as.data.frame(mean_REP_value)
plot_rep$date_str <- rownames(plot_rep)

# Convert the string date (e.g., "20200118") to a proper Date object
# Adjust the format "%Y%m%d" if your strings look different
plot_rep$Date <- as.Date(plot_rep$date_str, format="%Y%m%d")

# Create the time-series plot
ggplot(plot_rep, aes(x = Date, y = mean)) +
  geom_line(color = "darkgreen", size = 1) +
  geom_point(color = "forestgreen") +
  labs(title = "Average REP Time Series",
       x = "Date",
       y = "Mean REP") +
  theme_minimal()

# ----- Step 15: Save the result-----
# Rename the datasets to prepare for the merge
# For each datasets, delete "date_str" column
plot_ndvi <- plot_ndvi[, c("Date", "mean")]
colnames(plot_ndvi)[2] <- "NDVI"
plot_kndvi <- plot_kndvi[, c("Date", "mean")]
colnames(plot_kndvi)[2] <- "kNDVI"
plot_crswir <- plot_crswir[, c("Date", "mean")]
colnames(plot_crswir)[2] <- "CRSWIR"
plot_cire <- plot_cire[, c("Date", "mean")]
colnames(plot_cire)[2] <- "CIre"
plot_rep <- plot_rep[, c("Date", "mean")]
colnames(plot_rep)[2] <- "REP"
# Merge the datasets into one
plot_data <- merge(plot_ndvi, plot_kndvi, by = "Date", all = TRUE)
plot_data <- merge(plot_data, plot_crswir, by = "Date", all = TRUE)
plot_data <- merge(plot_data, plot_cire, by = "Date", all = TRUE)
plot_data <- merge(plot_data, plot_rep, by = "Date", all = TRUE)
# Delete rows if the "Date" column is NA
plot_data <- plot_data[!is.na(plot_data$Date), ]
# Export the result in a csv file
# Save the intermediate csv file, replace the old one if a file with the same name already exist.
write.csv(
  plot_data,
  "Data/Sentinel2/VegetationIndices_TimeSeries.csv",
  row.names = FALSE
)

