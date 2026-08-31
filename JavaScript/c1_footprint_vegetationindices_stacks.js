/************************************
# Footprint-Based Sentinel-2 Multi-Index Time Series Extraction
# Period: 2020–2025
# Cloud filter (AOI-based): cloud_fraction ≤ 0.5
#
# Indices:
# 10 m → NDVI, NIRv, WDVI, MSAVI, EVI
# 20 m → NDRE, CIredge, MTCI
#
# Input: Flux footprint polygon (aoi)
#
# Author: Zhiyu Wu
# Date: 26/03/2026
************************************/

// ==========================================================
// 1. Study Area
// ==========================================================
Map.centerObject(aoi, 7);
Map.addLayer(aoi, {color: 'green'}, 'AOI');

var startDate = '2020-01-01';
var endDate   = '2025-12-31';


// ==========================================================
// 2. Load Sentinel-2
// ==========================================================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(startDate, endDate);


// ==========================================================
// 3. Cloud filtering
// ==========================================================
function addCloudFraction(image) {

  var scl = image.select('SCL');

  var cloudMask = scl.eq(3)
                    .or(scl.eq(8))
                    .or(scl.eq(9))
                    .or(scl.eq(10));

  var cloudFrac = cloudMask.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: aoi,
    scale: 20,
    maxPixels: 1e13
  }).get('SCL');

  return image.set('cloud_frac_aoi', cloudFrac);
}

s2 = s2.map(addCloudFraction)
       .filter(ee.Filter.lt('cloud_frac_aoi', 0.5));


// ==========================================================
// 4. Preprocessing
// ==========================================================
function preprocess(image) {

  var scl = image.select('SCL');

  var mask = scl.neq(3)
                .and(scl.neq(8))
                .and(scl.neq(9))
                .and(scl.neq(10))
                .and(scl.neq(11));

  var bands = [
    'B2',        // Blue
    'B4','B8',   // 10 m
    'B5','B6',   // red-edge
    'B8A'
  ];

  var scaled = image.select(bands).divide(10000);

  return scaled.updateMask(mask)
               .copyProperties(image, ['system:time_start']);
}


// ==========================================================
// 5. Add Indices
// ==========================================================
function addIndices(image) {

  var blue = image.select('B2');
  var red  = image.select('B4');
  var nir  = image.select('B8');
  var nirN = image.select('B8A');
  var re1  = image.select('B5');
  var re2  = image.select('B6');

  var s = 1.2;

  // --- 10 m indices ---
  var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');

  var nirv = ndvi.multiply(nir).rename('NIRv');

  var wdvi = nir.subtract(red.multiply(s)).rename('WDVI');

  var msavi = nir.multiply(2).add(1)
    .subtract(
      nir.multiply(2).add(1).pow(2)
        .subtract(nir.subtract(red).multiply(8))
        .sqrt()
    )
    .divide(2)
    .rename('MSAVI');

  var evi = nir.subtract(red)
    .multiply(2.5)
    .divide(
      nir.add(red.multiply(6))
         .subtract(blue.multiply(7.5))
         .add(1)
    )
    .rename('EVI');

  // --- 20 m indices ---
  var ndre = nirN.subtract(re1).divide(nirN.add(re1)).rename('NDRE');

  var cire = nirN.divide(re1).subtract(1).rename('CIredge');

  var mtci = re2.subtract(re1)
                .divide(re1.subtract(red))
                .rename('MTCI');

  return image.addBands([
    ndvi, nirv, wdvi, msavi, evi,
    ndre, cire, mtci
  ]);
}


// ==========================================================
// 6. Build collection
// ==========================================================
var vegCollection = s2
  .map(preprocess)
  .map(addIndices)
  .select([
    'NDVI','NIRv','WDVI','MSAVI','EVI',
    'NDRE','CIredge','MTCI'
  ]);


// ==========================================================
// 7. Daily composites
// ==========================================================
function dailyComposite(collection) {

  var days = ee.List(
    collection.aggregate_array('system:time_start')
      .map(function(t){
        return ee.Date(t).format('YYYY-MM-dd');
      })
  ).distinct();

  var dailyImages = days.map(function(dayStr){

    var day = ee.Date(dayStr);
    var nextDay = day.advance(1,'day');

    return collection.filterDate(day, nextDay).mean()
      .set('system:time_start', day.millis());
  });

  return ee.ImageCollection(dailyImages);
}

var vegDaily = dailyComposite(vegCollection);


// ==========================================================
// 8. Separate collections
// ==========================================================
var collection10m = vegDaily.select([
  'NDVI','NIRv','WDVI','MSAVI','EVI'
]);

var collection20m = vegDaily.select([
  'NDRE','CIredge','MTCI'
]);


// ==========================================================
// 9. Convert to stacks
// ==========================================================
var timestamps = vegDaily.aggregate_array('system:time_start');

// --- 10 m ---
var stack10m = collection10m.toBands();

var bandNames10 = ee.List.sequence(0, timestamps.size().subtract(1))
  .map(function(i){
    var dateStr = ee.Date(timestamps.get(i)).format('YYYYMMdd');
    return [
      dateStr.cat('_NDVI'),
      dateStr.cat('_NIRv'),
      dateStr.cat('_WDVI'),
      dateStr.cat('_MSAVI'),
      dateStr.cat('_EVI')
    ];
  }).flatten();

stack10m = stack10m.rename(bandNames10).toFloat();


// --- 20 m ---
var stack20m = collection20m.toBands();

var bandNames20 = ee.List.sequence(0, timestamps.size().subtract(1))
  .map(function(i){
    var dateStr = ee.Date(timestamps.get(i)).format('YYYYMMdd');
    return [
      dateStr.cat('_NDRE'),
      dateStr.cat('_CIredge'),
      dateStr.cat('_MTCI')
    ];
  }).flatten();

stack20m = stack20m.rename(bandNames20).toFloat();


// ==========================================================
// 10. Export
// ==========================================================
Export.image.toDrive({
  image: stack10m.clip(aoi),
  description: 'VIs_10m_Stack_RDNew',
  folder: 'MGIThesis',
  fileNamePrefix: 'VIs_10m_2020_2025_RDNew',
  region: aoi,
  crs: 'EPSG:28992',
  scale: 10,
  maxPixels: 1e13
});

Export.image.toDrive({
  image: stack20m.clip(aoi),
  description: 'VIs_20m_Stack_RDNew',
  folder: 'MGIThesis',
  fileNamePrefix: 'VIs_20m_2020_2025_RDNew',
  region: aoi,
  crs: 'EPSG:28992',
  scale: 20,
  maxPixels: 1e13
});