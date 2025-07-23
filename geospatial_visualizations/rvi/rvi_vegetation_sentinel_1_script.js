//  RVI (Radar Vegetation Index) Monitoring Using Sentinel-1 SAR Imagery (Google Earth Engine Tutorial)

// 1. Defining the Area Of Interest
var StudyAreaCoords = [ 
  [1.0451427170895977,47.5985680423624],
  [1.5134349534177227,47.5985680423624],
  [1.5134349534177227,47.85721657974307],
  [1.0451427170895977,47.85721657974307],
  [1.0451427170895977,47.5985680423624]
];

// 2. Convert AOI into a polygon geometry
var StudyArea = ee.Geometry.Polygon(StudyAreaCoords);

// 3. Center the map on the study area
Map.centerObject(StudyArea, 10);
Map.addLayer(StudyArea, {}, 'Study Area (France)')

// 4. Load Sentinel-1 SAR Image Collection
var s1Collection = ee.ImageCollection("COPERNICUS/S1_GRD")
  .filterDate('2020-01-01', '2024-12-31')
  .filterBounds(StudyArea)
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV')) // Ensure VV polarization is available
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH')) // Ensure VH polarization is available
  .filter(ee.Filter.eq('instrumentMode', 'IW'))                            // Select only Interferometric Wide (IW) mode images
  .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))               // Filter ascending orbit images
  .select(['VV', 'VH']);                                                   // Select VV and VH bands for analysis

// 4. Print details of the filtered Sentinel-1 Dataset
print('Filtered Sentinel-1 Collection:', s1Collection);

// 5. Creating RVI (Radar Vegetation Index) for each image in the collection
var rviCollection = s1Collection.map(function(image) {
  var sigma = ee.Image(10).pow(image.divide(10)); // Convert dB values to linear scale
  
  // 5.1. Compute RVI using the formula: RVI = (4 * VH) / (VH + VV)
  var rvi = sigma.expression('(4 * vh) / (vh + vv)', {
    'vv': sigma.select('VV'),
    'vh': sigma.select('VH')
  }).rename('RVI');

  // 5.2. Apply a 30-meter focal median filter for noise reduction
  var rviSmoothed = rvi.focalMedian(30, 'square', 'meters');

  // 5.3. Copy time properties and return the processed image
  return rviSmoothed.copyProperties(image, ['system:time_start', 'system:time_end']);
});

// 6. Display the mean RVI image
Map.addLayer(rviCollection.mean().clip(StudyArea), {min: 0, max: 1, palette: ['blue', 'green', 'yellow']}, 'Mean RVI');

// 7. Define a smaller region for time-series analysis
var SubRegCoords = [
    [1.2637569807407312,47.71299528412205],
    [1.4312984846469812,47.71299528412205],
    [1.4312984846469812,47.8113096174144],
    [1.2637569807407312,47.8113096174144],
    [1.2637569807407312,47.71299528412205]
];

// 8. Convert to a geometry polygon
var SubReg = ee.Geometry.Polygon(SubRegCoords);
Map.addLayer(SubReg, {}, 'Sub-Region (France)')

// 9. Generate RVI time series for the subregion
print(
  ui.Chart.image.series({
    imageCollection: rviCollection,
    region: SubReg,
    reducer: ee.Reducer.mean(),
    scale: 10,
    xProperty: 'system:time_start'
  }).setOptions({
    title: 'Radar Vegetation Index (RVI) Time Series',
    vAxis: {title: 'RVI'},
    hAxis: {title: 'Year'},
    series: {0: {color: 'green'}},
    pointSize: 3
  })
);

// 10. Generate Sentinel-1 backscatter time series for VV & VH
print(
  ui.Chart.image.series({
    imageCollection: s1Collection,
    region: SubReg,
    reducer: ee.Reducer.mean(),
    scale: 10,
    xProperty: 'system:time_start'
  }).setOptions({
    title: 'Sentinel-1 VV & VH Backscatter Time Series',
    vAxis: {title: 'Backscatter (dB)'},
    hAxis: {title: 'Year'},
    series: {0: {color: 'blue'}, 1: {color: 'orange'}},
    pointSize: 3
  })
);

// Extract RVI images for 2023-2024
var rvi_2023_2024 = rviCollection.filterDate('2023-01-01', '2024-12-31').mean();

// Add RVI layer for 2023-2024
Map.addLayer(rvi_2023_2024.clip(StudyArea), {min: 0, max: 1, palette: ['blue', 'green', 'yellow']}, 'RVI 2023-2024');

// Export the processed RVI image
Export.image.toDrive({
  image: rvi_2023_2024.clip(StudyArea),
  description: 'RVI_2023_2024',
  region: StudyArea,
  scale: 30,
  crs: 'EPSG:4326',
  folder: 'RVI_Analysis',
  maxPixels: 1e13
});
