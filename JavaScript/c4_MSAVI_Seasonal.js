/************************************
# C4: Footprint-Based Sentinel-2 Multi-Index Time Series Extraction
# Period: 2020–2025, by Season
# Cloud filter (AOI-based): cloud_fraction ≤ 0.5
#
# Indices: MSAVI (10 m)
#
# Input: Seasonal flux footprint polygons
#        (aoi_spring, aoi_summer, aoi_autumn, aoi_winter)
#
# Season definition:
#   Spring : March–May
#   Summer : June–August
#   Autumn : September–November
#   Winter : December–February
#
# Author: Zhiyu Wu
# Date: 27/05/2026
************************************/

// ==========================================================
// 1. Season definitions
// ==========================================================
var seasons = [
  { name: 'spring', aoi: aoi_spring, months: [3, 4, 5] },
  { name: 'summer', aoi: aoi_summer, months: [6, 7, 8] },
  { name: 'autumn', aoi: aoi_autumn, months: [9, 10, 11] },
  { name: 'winter', aoi: aoi_winter, months: [12, 1, 2] }
];

var startYear = 2020;
var endYear   = 2025;


// ==========================================================
// 2. Cloud masking
// ==========================================================
function addCloudFraction(image, aoi) {
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


// ==========================================================
// 3. Preprocessing
// ==========================================================
function preprocess(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(3)
                .and(scl.neq(8))
                .and(scl.neq(9))
                .and(scl.neq(10))
                .and(scl.neq(11));
  var scaled = image.select(['B4', 'B8']).divide(10000);
  return scaled.updateMask(mask)
               .copyProperties(image, ['system:time_start']);
}


// ==========================================================
// 4. Compute MSAVI
// ==========================================================
function addMSAVI(image) {
  var red = image.select('B4');
  var nir = image.select('B8');
  var msavi = nir.multiply(2).add(1)
    .subtract(
      nir.multiply(2).add(1).pow(2)
        .subtract(nir.subtract(red).multiply(8))
        .sqrt()
    )
    .divide(2)
    .rename('MSAVI');
  return image.addBands(msavi);
}


// ==========================================================
// 5. Daily composites
// ==========================================================
function dailyComposite(collection) {
  var days = ee.List(
    collection.aggregate_array('system:time_start')
      .map(function(t) {
        return ee.Date(t).format('YYYY-MM-dd');
      })
  ).distinct();

  var dailyImages = days.map(function(dayStr) {
    var day = ee.Date(dayStr);
    var nextDay = day.advance(1, 'day');
    return collection.filterDate(day, nextDay).mean()
      .set('system:time_start', day.millis());
  });

  return ee.ImageCollection(dailyImages);
}


// ==========================================================
// 6. Build seasonal date ranges for 2020–2025
// ==========================================================
function getSeasonalDateRanges(months, startYear, endYear) {
  var ranges = [];
  for (var y = startYear; y <= endYear; y++) {
    if (months[0] === 12) {
      // Winter: Dec(y-1)–Feb(y), but also include Dec(y) for completeness
      // Dec of previous year to Feb of current year
      ranges.push({ start: (y - 1) + '-12-01', end: y + '-03-01' });
    } else {
      var startMonth = String(months[0]).length === 1 ? '0' + months[0] : String(months[0]);
      var endMonth   = months[months.length - 1] + 1;
      var endMonthStr = String(endMonth).length === 1 ? '0' + endMonth : String(endMonth);
      ranges.push({ start: y + '-' + startMonth + '-01', end: y + '-' + endMonthStr + '-01' });
    }
  }
  return ranges;
}


// ==========================================================
// 7. Process and export each season
// ==========================================================
seasons.forEach(function(season) {

  var aoi    = season.aoi;
  var ranges = getSeasonalDateRanges(season.months, startYear, endYear);

  // Collect all seasonal images across years
  var allImages = ee.ImageCollection([]);

  ranges.forEach(function(range) {
    var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
      .filterBounds(aoi)
      .filterDate(range.start, range.end)
      .map(function(image) { return addCloudFraction(image, aoi); })
      .filter(ee.Filter.lt('cloud_frac_aoi', 0.5))
      .map(preprocess)
      .map(addMSAVI)
      .select('MSAVI');

    allImages = allImages.merge(s2);
  });

  // Daily composites
  var daily = dailyComposite(allImages);

  // Build band-named stack
  var timestamps = daily.aggregate_array('system:time_start');

  var bandNames = ee.List.sequence(0, timestamps.size().subtract(1))
    .map(function(i) {
      return ee.Date(timestamps.get(i)).format('YYYYMMdd').cat('_MSAVI');
    });

  var stack = daily.toBands().rename(bandNames).toFloat();

  // Export
  Export.image.toDrive({
    image: stack.clip(aoi),
    description: 'MSAVI_' + season.name + '_2020_2025',
    folder: 'MGIThesis',
    fileNamePrefix: 'MSAVI_' + season.name + '_2020_2025',
    region: aoi,
    crs: 'EPSG:28992',
    scale: 10,
    maxPixels: 1e13
  });

  print('Export submitted for season: ' + season.name);
});

// ==========================================================
// 8. Visualise AOIs
// ==========================================================
Map.centerObject(aoi_spring, 7);
Map.addLayer(aoi_spring, {color: 'green'},  'AOI Spring');
Map.addLayer(aoi_summer, {color: 'yellow'}, 'AOI Summer');
Map.addLayer(aoi_autumn, {color: 'orange'}, 'AOI Autumn');
Map.addLayer(aoi_winter, {color: 'blue'},   'AOI Winter');