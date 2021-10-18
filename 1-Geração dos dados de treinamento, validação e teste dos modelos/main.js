// IMPORTS

var ARDLib = require('users/saraiva/packages:ARDLib.js');

// END IMPORTS

// SETTINGS

/* To train using LEM Dataset */
var startDate = '2017-06-01';
var endDate = '2018-06-01';
var datasetPath = 'users/saraiva/puc_artificial_intelligence/datasets/LEM_2017_2018';
var datasetMask;
var dates = ['06_2017','07_2017','08_2017','09_2017','10_2017','11_2017','12_2017', '01_2018',
                '02_2018','03_2018','04_2018','05_2018', '06_2018'];
var points_qtd = 20000;
var seed = 1;

/* To train using LEM+ Dataset */
var startDate = '2019-10-01';
var endDate = '2020-10-01';
var datasetPath = 'users/saraiva/puc_artificial_intelligence/datasets/LEM_2019_2020';
var datasetMask;
var dates = ['10_2019','11_2019','12_2019','01_2020','02_2020','03_2020','04_2020', '05_2020', '06_2020',
            '07_2020','08_2020','09_2020'];
var points_qtd = 20000;
var seed = 1;

/* To predict and evaluation */
// var startDate = '2020-10-01';
// var endDate = '2021-10-01';
// var datasetPath = 'users/saraiva/puc_artificial_intelligence/rois/extremo_oeste_bahiano_shape';
// var datasetMask = 'users/saraiva/puc_artificial_intelligence/extremo_oeste_baiano_raster';
// var dates = [];
// var points_qtd = 0;
// var seed = 1;

var outputPath = 'users/saraiva/puc_artificial_intelligence/mosaics';
var trainingProportion = 0.8;
var scale = 30;
var cloudCover = 80;
var lemDataset = ee.FeatureCollection(datasetPath);

var roi = lemDataset.geometry();

Map.addLayer(roi, {}, 'ROI');

var classes = ee.Dictionary({
  'not identified': -1,
  'soybean': 1,
  'maize': 2,
  'corn': 2, // LEM 2019/2020
  'cotton': 3,
  'coffee': 4,
  'beans': 5,
  'wheat': -1,
  'sorghum': 6,
  'millet': 7,
  'eucalyptus': 8,
  'pasture': 9,
  'hay': 10,
  'grass': 11,
  'crotalari': -1,
  'crotalaria': -1, // LEM 2019/2020
  'maize+crotalari': -1,
  'cerrado': 12,
  'conversion area': 13,
  'uncultivated soil': 13,
  'ncc': -1,
  'brachiaria': -1, // LEM 2019/2020,
});

// END SETTINGS

// TIME SERIES PREPROCESSING
var addNDVI = function(image){
  var ndvi = image.normalizedDifference(['NIR', 'RED']).rename('NDVI');
  return ee.Image(image.addBands(ndvi).toFloat().copyProperties(image, image.propertyNames()));
};

var l8Collection = ee.ImageCollection(ARDLib.getCollection('LANDSAT_8', roi, startDate, endDate, cloudCover));
var s2Collection = ee.ImageCollection(ARDLib.getCollection('SENTINEL_2', roi, startDate, endDate, cloudCover));

var collection = l8Collection.merge(s2Collection)
  .select(['NIR', 'SWIR1', 'RED'])
  .map(addNDVI)
  .sort('system:time_start');

Export.table.toDrive({
  collection: ee.FeatureCollection([ee.Feature(null, {
    'images_count': collection.size(),
    'images_size_bytes':  collection.aggregate_sum('system:asset_size')
  })])
});

var smoothed = ARDLib.smoorthing(collection, 8);

var products = ARDLib.get16DayProducts(collection, startDate, endDate, 'NDVI');

var stackedImage = products.toBands();

if(datasetMask !== undefined){
  stackedImage = stackedImage.updateMask(ee.Image(datasetMask));
}

var filename = startDate + '_' + endDate;
var mosaicPath = outputPath + '/' + filename + '_mosaic';


Export.image.toDrive({
  image: stackedImage.multiply(10000).int16(),
  description: filename + '_drive',
  folder: 'lem_dataset_output',
  fileNamePrefix: filename + '_mosaic',
  scale: 30,
  region: roi.bounds(),
  crs: 'EPSG:3857',
  maxPixels: 1E13
});

Export.image.toAsset({
  image: stackedImage.multiply(10000).int16(),
  description: filename + '_asset',
  assetId: mosaicPath,
  scale: 30,
  region: roi.bounds(),
  crs: 'EPSG:3857',
  maxPixels: 1E13
});

// END TIME SERIES PREPROCESSING



// LABELS PREPROCESSING

lemDataset = lemDataset
  .select(dates)
  .map(function(feature){
    var propertyNames = feature.propertyNames().remove('system:index');

    var properties = feature.toDictionary(propertyNames);

    var newProperties = propertyNames.iterate(function(propertyName, allProperties){
      propertyName = ee.String(propertyName);
      allProperties = ee.Dictionary(allProperties);


      var propertyValue = allProperties.getString(propertyName).toLowerCase();

      var newPropertyValue = classes.getNumber(propertyValue);

      return allProperties.set(propertyName, newPropertyValue);
    }, properties);

    newProperties = ee.Dictionary(newProperties);

    var validSample = newProperties.values().indexOf(-1).eq(-1);

    newProperties = newProperties
      .set('validSample', validSample);

    return feature.setMulti(newProperties);
  }).filterMetadata('validSample', 'equals', 1);

lemDataset = lemDataset.randomColumn({seed: 1}).sort('random');

var trainingSamples = lemDataset.size().multiply(trainingProportion).int();

var trainingDataset = ee.FeatureCollection(lemDataset.limit(trainingSamples, 'random', true));
var validationDataset = ee.FeatureCollection(lemDataset.limit(lemDataset.size().subtract(trainingSamples), 'random', false));

var exportedMosaic = ee.Image(mosaicPath);

buildSamples('training', trainingDataset, parseInt(points_qtd * trainingProportion));
buildSamples('validation', validationDataset, parseInt(points_qtd * (1 - trainingProportion)));

function buildSamples(title, dataset, maxPoints){
  var points = ee.FeatureCollection.randomPoints({
    region: dataset.geometry().buffer(2 * -scale),
    points: maxPoints,
    seed: seed,
    maxError: ee.ErrorMargin(scale)
  }).map(function(feature){
    var lemFeature = dataset.filterBounds(feature.geometry()).first();
    var properties = lemFeature.toDictionary(lemFeature.propertyNames());
    properties = properties.set('id', feature.id());
    return feature.setMulti(properties);
  });

  // END LABELS PREPROCESSING

  // MERGE LABELS + TIME SERIES

  var samples = ee.FeatureCollection(points.map(function(feature) {

    var imageProperties = exportedMosaic
      .reduceRegion({
        reducer: ee.Reducer.median(),
        geometry: feature.geometry(),
        scale: 30
      });

    imageProperties = ee.FeatureCollection(imageProperties.map(function(key, value){
      var parts = ee.String(key).split('_');
      return ee.Feature(null, {
        date: ee.Date.parse('yyyyMMdd', ee.String(parts.get(1))),
        band: parts.get(-1),
        value: ee.Number(value)
      });
    }).values());

    var samples = imageProperties
      .aggregate_array('date')
      .distinct()
      .map(function(date){
        date = ee.Date(date);

        var features = imageProperties
          .filterMetadata('date', 'equals', date)
          .filter(ee.Filter.notNull(['band', 'value']));

        var bandNames = features.aggregate_array('band');
        var bandValues = features.aggregate_array('value');

        var properties = ee.Dictionary
          .fromLists(bandNames, bandValues)
          .set('DATE', date);

        return ee.Feature(null, feature
          .select(['id', date.format("MM_yyyy")], ['ID', 'CLASS'])
          .toDictionary(['ID', 'CLASS'])
          .combine(properties));

      });

    return ee.Feature(null, {'samples': samples});

  }).aggregate_array('samples').flatten());

  var sampleFilename = [filename, title,'part', seed].join('_');

  Export.table.toDrive({
    collection: samples,
    description: sampleFilename,
    fileNamePrefix: sampleFilename,
    folder: 'lem_dataset_output',
    fileFormat: 'csv'
  });

  // END MERGE LABELS + TIME SERIES


  Map.onClick(function(coords){
    var point = ee.Geometry.Point([coords.lon, coords.lat]);

    var rawL8Collection = ee.ImageCollection(ARDLib.satelliteCollections['LANDSAT_8'])
      .filterDate(startDate, endDate)
      .filterBounds(point)
      .map(ARDLib.padronizeProperties)
      .map(ARDLib.padronizeBandNames)
      .map(ARDLib.maskImage)
      .map(addNDVI);

    var rawS2Collection = ee.ImageCollection(ARDLib.satelliteCollections['SENTINEL_2'])
      .filterDate(startDate, endDate)
      .filterBounds(point)
      .map(ARDLib.padronizeProperties)
      .map(ARDLib.padronizeBandNames)
      .map(ARDLib.maskImage)
      .map(addNDVI);

    var l8Collection = ee.ImageCollection(ARDLib.getCollection('LANDSAT_8', point, startDate, endDate, cloudCover));
    var s2Collection = ee.ImageCollection(ARDLib.getCollection('SENTINEL_2', point, startDate, endDate, cloudCover));

    var collection = l8Collection.merge(s2Collection);

    var bandPassAjustmentCol = collection
      .select(['NIR', 'SWIR1', 'RED'])
      .map(addNDVI)
      .sort('system:time_start');

    var smoothed = ARDLib.smoorthing(bandPassAjustmentCol, 8);

    var formatProperties = function(newBandName){
      return function(list){
        list = ee.List(list);
        var ndvi = list.get(-1);
        var time = ee.Date(list.get(-2));

        var properties = {'TIME': time};
        properties[newBandName] = ndvi;
        return ee.Feature(null, properties);

      };
    };

    var rawL8TimeSeries = rawL8Collection.select(['NDVI'])
      .getRegion({geometry: point, scale: scale})
      .slice(1)
      .map(formatProperties('Landsat 8'));

    var rawS2TimeSeries = rawS2Collection.select(['NDVI'])
      .getRegion({geometry: point, scale: scale})
      .slice(1)
      .map(formatProperties('Sentinel 2A/2B'));

    var bandPassTimeSeries = bandPassAjustmentCol.select(['NDVI'])
      .getRegion({geometry: point, scale: scale})
      .slice(1)
      .map(formatProperties('Série Ajustada'));

    var normTimeSeries = smoothed.select(['NDVI'])
      .getRegion({geometry: point, scale: scale})
      .slice(1)
      .map(formatProperties('Série Ajustada + Média Móvel (16 dias)'));

    var merged = rawL8TimeSeries
      .cat(rawS2TimeSeries)
      .cat(bandPassTimeSeries)
      .cat(normTimeSeries);

    var chart = ui.Chart.feature.byFeature({
      features: merged,
      xProperty: 'TIME',
      yProperties: ['Landsat 8', 'Sentinel 2A/2B',
                    'Série Ajustada', 'Série Ajustada + Média Móvel (16 dias)']
    });

    var options = {
      series: {
        0: {targetAxisIndex:0, lineWidth: 0, pointSize: 5 },
        1: {targetAxisIndex:0, lineWidth: 0, pointSize: 5 },
        2: {targetAxisIndex:0, lineWidth: 3, pointSize: 0},
        3: {targetAxisIndex:0, lineWidth: 3, pointSize: 0},
      },
      interpolateNulls: true,
      hAxis: {
        title: 'Data da Imagem',
        textStyle: {fontSize: 12},
        format: 'dd/MM/yyyy',
      },
      vAxis: {
        title: 'NDVI',
        textStyle: {fontSize: 12},
      }
    };

    chart.setOptions(options);

    print(chart);

  });
}
