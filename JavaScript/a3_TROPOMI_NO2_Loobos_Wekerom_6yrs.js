/************************************
# A2: Data acquisition for Sentinel-5P tropospheric and 
# total NO2 column density at Loobos and Wekerom (3 km buffers)
#
# This script extracts:
#  - Tropospheric NO₂ column density
#  - Total vertical column density (VCD) NO₂
#
# Source: COPERNICUS/S5P/OFFL/L3_NO2
# Period: 2020–2025
# Cloud filter: cloud_fraction ≤ 0.2
#
# Author: Zhiyu Wu
# Date: 02/03/2026
************************************/


/************************************
 * 1. CREATE GEOMETRIES & BUFFERS
 ************************************/

var LoobosGeometry  = ee.Geometry.Point([5.743889, 52.167778]);
var WekeromGeometry = ee.Geometry.Point([5.708333, 52.111667]);

var bufferRadius = 3000;

var LoobosBuffer  = LoobosGeometry.buffer(bufferRadius);
var WekeromBuffer = WekeromGeometry.buffer(bufferRadius);


/************************************
 * 2. DATE RANGE (2020–2025)
 ************************************/

var startDate = '2020-01-01';
var endDate   = '2025-12-31';


/************************************
 * 3. LOAD SENTINEL-5P NO2 DATA
 ************************************/

var bandNameTropo  = 'tropospheric_NO2_column_number_density';
var bandNameTotal  = 'NO2_column_number_density';
var bandNameCloud  = 'cloud_fraction';

var NO2 = ee.ImageCollection('COPERNICUS/S5P/OFFL/L3_NO2')
  .filterDate(startDate, endDate)
  .map(function(image) {

    var cloudMask = image.select(bandNameCloud).lte(0.2);
    var dateTime = ee.Date(image.get('system:time_start'));

    return image
      .updateMask(cloudMask)
      .set({
        'date': dateTime.format('YYYY-MM-dd'),
        'datetime': dateTime.format('YYYY-MM-dd HH:mm'),
        'year': dateTime.get('year'),
        'time_start': image.get('system:time_start')
      });
  });


/************************************
 * 4. FUNCTION TO EXTRACT MEAN VALUES
 ************************************/

function extractStats(image, region, locationName) {

  // Mean tropospheric NO2
  var meanTropo = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: 1000,
    maxPixels: 1e9
  }).get(bandNameTropo);

  // Mean total NO2 VCD
  var meanTotal = image.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: 1000,
    maxPixels: 1e9
  }).get(bandNameTotal);

  // Mean cloud fraction (unmasked)
  var meanCloud = image
    .unmask()
    .select(bandNameCloud)
    .reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: 1000,
      maxPixels: 1e9
    }).get(bandNameCloud);

  return ee.Feature(null, {
    'date': image.get('date'),
    'datetime': image.get('datetime'),
    'year': image.get('year'),
    'time_start': image.get('time_start'),
    'location': locationName,
    'Tropospheric_NO2_value': meanTropo,
    'Total_NO2_VCD_value': meanTotal,
    'cloud_fraction': meanCloud
  });
}


/************************************
 * 5. APPLY EXTRACTION
 ************************************/

var LoobosData = NO2.map(function(image) {
  return extractStats(image, LoobosBuffer, 'Loobos');
});

var WekeromData = NO2.map(function(image) {
  return extractStats(image, WekeromBuffer, 'Wekerom');
});

var combinedData = ee.FeatureCollection(LoobosData)
  .merge(ee.FeatureCollection(WekeromData))
  .filter(ee.Filter.notNull(['Tropospheric_NO2_value']));

print(combinedData.limit(10));


/************************************
 * 6. EXPORT
 ************************************/

Export.table.toDrive({
  collection: combinedData,
  description: 'Trop_and_Total_NO2_LoobosWekerom_3kmBuf_20CF_2020_2025',
  folder: 'MGIthesis',
  fileFormat: 'CSV',
  selectors: [
    'date',
    'datetime',
    'year',
    'time_start',
    'location',
    'Tropospheric_NO2_value',
    'Total_NO2_VCD_value',
    'cloud_fraction'
  ]
});