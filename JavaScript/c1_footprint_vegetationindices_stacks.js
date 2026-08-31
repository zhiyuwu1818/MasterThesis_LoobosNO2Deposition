/************************************
# C1 – Footprint-Based Sentinel-2 Multi-Index Time Series Extraction
# Period: 2020–2025
# Cloud filter (AOI-based): cloud_fraction ≤ 0.5
# Indices: NDVI, CR_SWIR, CIre, REP, kNDVI
#
# Input: Flux footprint polygon (aoi)
# Dual Resolution Output:
# - 10 m stack (NDVI, kNDVI)
# - 20 m stack (CRSWIR, CIre, REP)
# Author: Zhiyu Wu
# Date: 03/03/2026
************************************/

// ==========================================================
// 1. Study Area
// ==========================================================
Map.centerObject(aoi, 7);
Map.addLayer(aoi, {color: 'green'}, 'AOI');

var startDate = '2020-01-01';
var endDate   = '2025-12-31';

// ==========================================================
// 2. Load Sentinel-2 SR Harmonized
// ==========================================================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate(startDate, endDate);

// ==========================================================
// 3. AOI-Based Cloud Fraction (SCL)
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

print('Usable S2 images:', s2.size());

// ==========================================================
// 4. Cloud Masking + Normalization
// ==========================================================
function preprocess(image) {

  var scl = image.select('SCL');

  var mask = scl.neq(3)
                .and(scl.neq(8))
                .and(scl.neq(9))
                .and(scl.neq(10))
                .and(scl.neq(11));

  var bands = [
    'B4','B8',        // NDVI
    'B8A','B11','B12',// CR-SWIR
    'B5','B6','B7'    // Red-edge
  ];

  var scaled = image.select(bands).divide(10000);

  return scaled.updateMask(mask)
               .copyProperties(image, ['system:time_start']);
}

// ==========================================================
// 5. Add Vegetation Indices
// ==========================================================
function addIndices(image) {

  var red  = image.select('B4');
  var nir  = image.select('B8');
  var nirN = image.select('B8A');
  var sw1  = image.select('B11');
  var sw2  = image.select('B12');
  var re1  = image.select('B5');
  var re2  = image.select('B6');
  var re3  = image.select('B7');

  // NDVI
  var ndvi = nir.subtract(red).divide(nir.add(red)).rename('NDVI');

  // kNDVI
  var kndvi = ndvi.pow(2).tanh().rename('kNDVI');

  // CR-SWIR
  // The denominator calculates the "continuum line" value at 1614nm (B11)
  // by interpolating between 865nm (B8A) and 2202nm (B12).
  var denom = nirN.add(
      sw2.subtract(nirN)
        .divide(2202 - 865)  // Step 1: Divide by the total distance (x2 - x1)
        .multiply(1614 - 865) // Step 2: Multiply by the distance to the target (x - x1)
  );

  var crswir = sw1.divide(denom).rename('CR_SWIR');

  // CIre
  var cire = re3.divide(re1).subtract(1).rename('CIre');

  // REP (4PLI)
  var rep = red.add(re3).divide(2)
    .subtract(re1)
    .divide(re2.subtract(re1))
    .multiply(40)
    .add(700)
    .rename('REP');

  return image.addBands([ndvi, kndvi, crswir, cire, rep]);
}

// ==========================================================
// 6. Build Multi-Index Collection
// ==========================================================
var vegCollection = s2
  .map(preprocess)
  .map(addIndices)
  .select(['NDVI','kNDVI','CR_SWIR','CIre','REP']);

print('Vegetation collection:', vegCollection.size());

// ==========================================================
// 7. Daily Composites
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

    var daily = collection.filterDate(day, nextDay).mean()
      .set('system:time_start', day.millis())
      .set('date', dayStr);

    return daily;
  });

  return ee.ImageCollection(dailyImages);
}

var vegDaily = dailyComposite(vegCollection);
print('Daily images:', vegDaily.size());

// ==========================================================
// 8. Separate 10 m and 20 m Index Collections
// ==========================================================

// 10 m indices (native resolution)
var collection10m = vegDaily.select(['NDVI','kNDVI']);

// 20 m indices (native resolution)
var collection20m = vegDaily.select(['CR_SWIR','CIre','REP']);

// ----------------------------------------------------------
// 9. Convert Each Collection to Multi-Band Stack
// ----------------------------------------------------------

// ---------- 10 m STACK ----------
var stack10m = collection10m.toBands();

// Rename bands → YYYYMMDD_index
var timestamps = vegDaily.aggregate_array('system:time_start');

var bandNames10 = ee.List.sequence(0, timestamps.size().subtract(1))
  .map(function(i){
    var t = ee.Number(timestamps.get(i));
    var dateStr = ee.Date(t).format('YYYYMMdd');
    return [
      dateStr.cat('_NDVI'),
      dateStr.cat('_kNDVI')
    ];
  }).flatten();

stack10m = stack10m.rename(bandNames10).toFloat();


// ---------- 20 m STACK ----------
var stack20m = collection20m.toBands();

var bandNames20 = ee.List.sequence(0, timestamps.size().subtract(1))
  .map(function(i){
    var t = ee.Number(timestamps.get(i));
    var dateStr = ee.Date(t).format('YYYYMMdd');
    return [
      dateStr.cat('_CRSWIR'),
      dateStr.cat('_CIre'),
      dateStr.cat('_REP')
    ];
  }).flatten();

stack20m = stack20m.rename(bandNames20).toFloat();


// ==========================================================
// 10. Export 10 m Stack (RD New)
// ==========================================================
Export.image.toDrive({
  image: stack10m.clip(aoi),
  description: 'MultiIndex_10m_Stack_RDNew',
  folder: 'MGIThesis',
  fileNamePrefix: 'MultiIndex_10m_Daily_2020_2025',
  region: aoi,
  crs: 'EPSG:28992',
  scale: 10,
  maxPixels: 1e13
});


// ==========================================================
// 11. Export 20 m Stack (RD New)
// ==========================================================
Export.image.toDrive({
  image: stack20m.clip(aoi),
  description: 'MultiIndex_20m_Stack_RDNew',
  folder: 'MGIThesis',
  fileNamePrefix: 'MultiIndex_20m_Daily_2020_2025',
  region: aoi,
  crs: 'EPSG:28992',
  scale: 20,
  maxPixels: 1e13
});