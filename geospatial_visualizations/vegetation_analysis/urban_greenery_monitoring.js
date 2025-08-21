// step 1: define the study area: reims boundary
Map.centerObject(reims, 13)

Map.addLayer(reims.style({
  color: '000000',
  fillColor: '00000000',
  width: 4
}), {}, 'Reims boundary');

// step 2: load the harmonized sentinel-2 image collection
var s2 = ee.ImageCollection('COPERNICUS/S2_HARMONIZED')
  .filterBounds(reims)
  .filterDate('2017-01-01', '2025-08-17')
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 20));

// step 3: cloud masking using QA60
function maskClouds(image) {
  var qa = image.select('QA60');
  var mask = qa.bitwiseAnd(1 << 10).eq(0)
               .and(qa.bitwiseAnd(1 << 11).eq(0));
  return image.updateMask(mask)
              .copyProperties(image, image.propertyNames());
}
var s2Masked = s2.map(maskClouds);

// step 4: calculate NDVI, EVI, NDBI indices
function addIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var evi = image.expression(
    '2.5 * ((NIR - RED) / (NIR + 6 * RED - 7.5 * BLUE + 1))', {
      'NIR': image.select('B8'),
      'RED': image.select('B4'),
      'BLUE': image.select('B2')
    }).rename('EVI');
  var ndbi = image.normalizedDifference(['B11', 'B8']).rename('NDBI');
  return image.addBands([ndvi, evi, ndbi])
              .copyProperties(image, image.propertyNames());
}
var s2WithIndices = s2Masked.map(addIndices);

// step 5: annual composites of indices
var years = ee.List.sequence(2017, 2025);
var annualStats = ee.ImageCollection.fromImages(
  years.map(function(y) {
    var start = ee.Date.fromYMD(y, 1, 1);
    var end = start.advance(1, 'year');
    var meanImage = s2WithIndices.filterDate(start, end)
      .mean()
      .select(['NDVI', 'EVI', 'NDBI'])
      .clip(reims)
      .set({
        'year': y,
        'system:time_start': ee.Date.fromYMD(y, 6, 1).millis()
      });
    return meanImage;
  })
);

// Number of images per year
var imageCountsPerYear = years.map(function(y) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var end = start.advance(1, 'year');
  
  var count = s2WithIndices.filterDate(start, end).size();
  
  return ee.Feature(null, {
    year: y,
    image_count: count
  });
});

var imageCountsFC = ee.FeatureCollection(imageCountsPerYear);

print('Years:', years);
print('Image counts per year:', imageCountsFC.aggregate_array('image_count'));

var chart = ui.Chart.feature.byFeature(imageCountsFC, 'year', 'image_count')
  .setOptions({
    title: 'Number of Images per Year',
    hAxis: {title: 'Year'},
    vAxis: {title: 'Image Count'},
    lineWidth: 2,
    pointSize: 5,
  });

print(chart);

// step 6: visualize NDVI / NDBI for 2024
var ndviVis = {min: 0, max: 0.8, palette: ['white', 'green']};
var ndbiVis = {min: -0.5, max: 0.5, palette: ['blue', 'gray', 'brown']};

var ndvi2024 = annualStats.filter(ee.Filter.eq('year', 2024)).first();
Map.addLayer(ndvi2024.select('NDVI'), ndviVis, 'NDVI 2024');
Map.addLayer(ndvi2024.select('NDBI'), ndbiVis, 'NDBI 2024');

// step 7: NDVI / NDBI annual time series charts
print(ui.Chart.image.series({
  imageCollection: annualStats.select('NDVI'),
  region: reims,
  reducer: ee.Reducer.mean(),
  scale: 10
}).setOptions({
  title: 'Mean annual NDVI in Reims (2017–2025)',
  hAxis: {title: 'Year'},
  vAxis: {title: 'NDVI'}
}));

print(ui.Chart.image.series({
  imageCollection: annualStats.select('NDBI'),
  region: reims,
  reducer: ee.Reducer.mean(),
  scale: 10
}).setOptions({
  title: 'Mean annual NDBI in Reims (2017–2025)',
  hAxis: {title: 'Year'},
  vAxis: {title: 'NDBI'}
}));

// step 8: NDVI linear trend (regression)
var ndviTimeSeries = s2WithIndices.select('NDVI').map(function(img) {
  var year = ee.Number(ee.Date(img.date()).get('year')).toFloat();
  var yearImage = ee.Image.constant(year).rename('year').toFloat();
  return img.addBands(yearImage);
});

var linearFit = ndviTimeSeries
  .select(['year', 'NDVI'])
  .map(function(img) {
    return img.toFloat();
  })
  .reduce(ee.Reducer.linearFit());

Map.addLayer(linearFit.select('scale').clip(reims), {
  min: -0.05, max: 0.05,
  palette: ['red', 'orange', 'yellow', 'white', 'lime', 'green', 'darkgreen']
}, 'NDVI trend (clipped)');

// step 9: classification (3 classes: vegetation / built-up / other)
var greenThreshold = 0.3;    // NDVI threshold for vegetation
var builtNdbiThreshold = 0.05;  // NDBI threshold for built-up
var builtEviThreshold = 0.1;   // low EVI threshold for built-up

function classifyVegBuiltOther(img) {
  var ndvi = img.select('NDVI');
  var evi = img.select('EVI');
  var ndbi = img.select('NDBI');

  var veg = ndvi.gt(greenThreshold);
  var built = ndbi.gt(builtNdbiThreshold).and(evi.lt(builtEviThreshold)).and(veg.not());

  var classified = ee.Image(0)
    .where(veg, 1)
    .where(built, 2)
    .rename('landcover');

  return classified.set('year', img.get('year'));
}

var classifiedAnnual = annualStats.map(classifyVegBuiltOther);

// step 10: compute % for each land cover class
var landCoverStats = classifiedAnnual.map(function(img) {
  var total = img.reduceRegion({
    reducer: ee.Reducer.count(),
    geometry: reims,
    scale: 10,
    maxPixels: 1e9
  }).get('landcover');

  var veg = img.eq(1).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: reims,
    scale: 10,
    maxPixels: 1e9
  }).get('landcover');

  var built = img.eq(2).reduceRegion({
    reducer: ee.Reducer.sum(),
    geometry: reims,
    scale: 10,
    maxPixels: 1e9
  }).get('landcover');

  var totalNum = ee.Number(total);
  return ee.Feature(null, {
    year: img.get('year'),
    vegetation: ee.Number(veg).divide(totalNum).multiply(100),
    built: ee.Number(built).divide(totalNum).multiply(100),
    other: ee.Number(100).subtract(
      ee.Number(veg).divide(totalNum).multiply(100)
        .add(ee.Number(built).divide(totalNum).multiply(100))
    )
  });
});

// step 11: bar chart for 3-class land cover evolution
var landCoverBarChart = ui.Chart.feature.byFeature({
  features: landCoverStats,
  xProperty: 'year',
  yProperties: ['vegetation', 'built', 'other']
}).setChartType('ColumnChart')
  .setOptions({
    title: 'Land cover evolution: vegetation, built-up, and others (%)',
    hAxis: {title: 'Year'},
    vAxis: {title: '% of area'},
    colors: ['green', 'black', 'gray'],
    bar: {groupWidth: '75%'},
    isStacked: true
  });

print(landCoverBarChart);

// step 12: display land cover classification for 2024
var classified2024 = classifiedAnnual.filter(ee.Filter.eq('year', 2024)).first();
var classified2024Masked = classified2024.updateMask(classified2024.neq(0));

Map.addLayer(classified2024Masked, {
  min: 1,
  max: 2,
  palette: ['green', 'black']
}, 'Land cover 2024 (vegetation / built-up)');

// step 13: detect vegetation change year to year
var changes = ee.ImageCollection(years.slice(1).map(function(year) {
  year = ee.Number(year);
  var prevYear = year.subtract(1);

  var imgPrev = classifiedAnnual.filter(ee.Filter.eq('year', prevYear)).first();
  var imgCurr = classifiedAnnual.filter(ee.Filter.eq('year', year)).first();

  var vegGained = imgCurr.eq(1).and(imgPrev.neq(1)).rename('veg_gain');
  var vegLost = imgCurr.neq(1).and(imgPrev.eq(1)).rename('veg_loss');

  var gainVis = vegGained.updateMask(vegGained)
                         .set('year', year)
                         .set('type', 'gain');

  var lossVis = vegLost.updateMask(vegLost)
                       .set('year', year)
                       .set('type', 'loss');

  return gainVis.addBands(lossVis)
                .set('year', year)
                .set('system:time_start', ee.Date.fromYMD(year, 6, 1).millis());
}));

// step 14: display vegetation gains/losses for a selected year
var changeYear = 2024;
var changeImg = changes.filter(ee.Filter.eq('year', changeYear)).first();

Map.addLayer(changeImg.select('veg_gain'), {
  palette: ['lime'],
}, 'Vegetation gain in ' + changeYear);

Map.addLayer(changeImg.select('veg_loss'), {
  palette: ['red'],
}, 'Vegetation loss in ' + changeYear);

