###############################################
# C3: Vegetation Index Time Series Extraction
# Loobos Flux Tower – 100 m Buffer
#
# Indices:
# NDVI, NIRv, WDVI, MSAVI, EVI, NDRE, CIredge, MTCI
#
# Author: Zhiyu Wu
# Updated: 2026
###############################################

# -------------------------------
# Step 1: Load libraries
# -------------------------------

library(terra)
library(ggplot2)

# -------------------------------
# Step 2: Load Sentinel-2 stacks
# -------------------------------

r10 <- rast("Data/Sentinel2/VIs_10m_2020_2025_RDNew.tif")
r20 <- rast("Data/Sentinel2/VIs_20m_2020_2025_RDNew.tif")

# Extract layers
ndvi_layers <- r10[[grep("_NDVI$", names(r10))]]
nirv_layers <- r10[[grep("_NIRv$", names(r10))]]
wdvi_layers <- r10[[grep("_WDVI$", names(r10))]]
msavi_layers <- r10[[grep("_MSAVI$", names(r10))]]
evi_layers <- r10[[grep("_EVI$", names(r10))]]

ndre_layers <- r20[[grep("_NDRE$", names(r20))]]
ciredge_layers <- r20[[grep("_CIredge$", names(r20))]]
mtci_layers <- r20[[grep("_MTCI$", names(r20))]]

# -------------------------------
# Step 3: Create Loobos tower buffer
# -------------------------------

loobos <- vect(
  data.frame(lon = 5.743889, lat = 52.167778),
  geom = c("lon","lat"),
  crs = "EPSG:4326"
)

loobos_rd <- project(loobos, crs(ndvi_layers))

loobos_buffer <- buffer(loobos_rd, width = 100)

# -------------------------------
# Step 4: Crop stacks to buffer
# -------------------------------

crop_mask <- function(r){
  x <- crop(r, loobos_buffer)
  x <- mask(x, loobos_buffer)
  return(x)
}

ndvi_crop <- crop_mask(ndvi_layers)
nirv_crop <- crop_mask(nirv_layers)
wdvi_crop <- crop_mask(wdvi_layers)
msavi_crop <- crop_mask(msavi_layers)
evi_crop <- crop_mask(evi_layers)

ndre_crop <- crop_mask(ndre_layers)
ciredge_crop <- crop_mask(ciredge_layers)
mtci_crop <- crop_mask(mtci_layers)

# -------------------------------
# Step 4b: Align spatial resolution
# Resample 20 m indices to 10 m grid
# -------------------------------

ndre_crop <- resample(ndre_crop, ndvi_crop, method = "bilinear")

ciredge_crop <- resample(ciredge_crop, ndvi_crop, method = "bilinear")

mtci_crop <- resample(mtci_crop, ndvi_crop, method = "bilinear")

# -------------------------------
# Step 5: Load external vegetation mask
# -------------------------------

# Load the mask file
veg_mask_raw <- rast("Data/Sentinel2/NDVI_Mask_10percent_RDNew.tif")

# Ensure the mask matches the extent and CRS of our cropped data
# (Assuming the .tif is already in RD New/EPSG:28992)
veg_mask_final <- crop(veg_mask_raw, ndvi_crop)

# Optional: If the mask uses 0 for non-vegetation, convert 0 to NA
# to ensure the mask() function in Step 6 treats them as 'exclude' areas.
veg_mask_final <- ifel(veg_mask_final == 0, NA, 1)

# -------------------------------
# Step 6: Apply vegetation mask
# -------------------------------

masked_ndvi <- mask(ndvi_crop, veg_mask_final)
masked_nirv <- mask(nirv_crop, veg_mask_final)
masked_wdvi <- mask(wdvi_crop, veg_mask_final)
masked_msavi <- mask(msavi_crop, veg_mask_final)
masked_evi <- mask(evi_crop, veg_mask_final)
masked_ndre <- mask(ndre_crop, veg_mask_final)
masked_ciredge <- mask(ciredge_crop, veg_mask_final)
masked_mtci <- mask(mtci_crop, veg_mask_final)

# -------------------------------
# Step 7: Extract mean time series
# -------------------------------

extract_ts <- function(r){

  # Calculate mean per layer within the masked area
  ts <- global(r, mean, na.rm = TRUE)

  df <- as.data.frame(ts)

  df$date_str <- rownames(df)

  df$Date <- as.Date(df$date_str, format="%Y%m%d")

  return(df)
}

# Compute mean time series first
plot_ndvi <- extract_ts(masked_ndvi)
plot_nirv <- extract_ts(masked_nirv)
plot_wdvi <- extract_ts(masked_wdvi)
plot_msavi <- extract_ts(masked_msavi)
plot_evi <- extract_ts(masked_evi)
plot_ndre <- extract_ts(masked_ndre)
plot_ciredge <- extract_ts(masked_ciredge)
plot_mtci <- extract_ts(masked_mtci)

# -------------------------------
# Step 8: Filter unrealistic values on the mean time series
# -------------------------------

# NDVI
plot_ndvi$NDVI <- plot_ndvi$mean # rename the column
plot_ndvi$NDVI[plot_ndvi$NDVI < 0.62 | plot_ndvi$NDVI > 0.9] <- NA
plot_ndvi <- plot_ndvi[, c("Date","NDVI")] # assign new column names

# NIRv
plot_nirv$NIRv <- plot_nirv$mean
plot_nirv$NIRv[plot_nirv$NIRv < 0.07 | plot_nirv$NIRv > 0.8] <- NA
plot_nirv <- plot_nirv[, c("Date","NIRv")]

# WDVI
plot_wdvi$WDVI <- plot_wdvi$mean
plot_wdvi$WDVI[plot_wdvi$WDVI < 0.05 | plot_wdvi$WDVI > 0.25] <- NA
plot_wdvi <- plot_wdvi[, c("Date","WDVI")]

# MSAVI
plot_msavi$MSAVI <- plot_msavi$mean
plot_msavi$MSAVI[plot_msavi$MSAVI < 0.25 | plot_msavi$MSAVI > 0.4] <- NA
plot_msavi <- plot_msavi[, c("Date","MSAVI")]

# EVI
plot_evi$EVI <- plot_evi$mean
plot_evi$EVI[plot_evi$EVI < 0.3 | plot_evi$EVI > 0.6] <- NA
plot_evi <- plot_evi[, c("Date","EVI")]

# NDRE
plot_ndre$NDRE <- plot_ndre$mean
plot_ndre$NDRE[plot_ndre$NDRE < 0.3 | plot_ndre$NDRE > 0.8] <- NA
plot_ndre <- plot_ndre[, c("Date","NDRE")]

# CIredge
plot_ciredge$CIredge <- plot_ciredge$mean
plot_ciredge$CIredge[plot_ciredge$CIredge < 1 | plot_ciredge$CIredge > 4] <- NA
plot_ciredge <- plot_ciredge[, c("Date","CIredge")]

# MTCI
plot_mtci$MTCI <- plot_mtci$mean
plot_mtci$MTCI[plot_mtci$MTCI < 2 | plot_mtci$MTCI > 4] <- NA
plot_mtci <- plot_mtci[, c("Date","MTCI")]

# -------------------------------
# Step 9: Merge datasets
# -------------------------------

plot_data <- merge(plot_ndvi, plot_nirv, by = "Date", all = TRUE)
plot_data <- merge(plot_data, plot_wdvi, by = "Date", all = TRUE)
plot_data <- merge(plot_data, plot_msavi, by = "Date", all = TRUE)
plot_data <- merge(plot_data, plot_evi, by = "Date", all = TRUE)
plot_data <- merge(plot_data, plot_ndre, by = "Date", all = TRUE)
plot_data <- merge(plot_data, plot_ciredge, by = "Date", all = TRUE)
plot_data <- merge(plot_data, plot_mtci, by = "Date", all = TRUE)

plot_data <- plot_data[!is.na(plot_data$Date), ]

# -------------------------------
# Step 10: Plot NDVI time series
# -------------------------------

ggplot(plot_data, aes(x = Date, y = NDVI)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_point(color = "forestgreen") +
  labs(
    title = "NDVI Time Series (100 m Buffer, Loobos)",
    x = "Date",
    y = "Mean NDVI"
  ) +
  theme_minimal()

# -------------------------------
# Step 11: Save CSV
# -------------------------------

write.csv(
  plot_data,
  "Data/Sentinel2/VegetationIndices_TimeSeries_New.csv",
  row.names = FALSE
)
