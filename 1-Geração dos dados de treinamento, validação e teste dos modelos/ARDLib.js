// https://github.com/ndminhhus/geeguide

var landsatSpacecrafts = ['LANDSAT_5', 'LANDSAT_7', 'LANDSAT_8'];
var copernicusSpacecrafts = ['Sentinel-2A', 'Sentinel-2B'];

var satelliteCollections = {
  'LANDSAT_5': 'LANDSAT/LT05/C02/T1_L2',
  'LANDSAT_7': 'LANDSAT/LE07/C02/T1_L2',
  'LANDSAT_8': 'LANDSAT/LC08/C02/T1_L2',
  'SENTINEL_2': 'COPERNICUS/S2'
};

var satelliteBandReescale = ee.Dictionary({
  'LANDSAT_5': {
    slope: 0.0000275, 
    offset: -0.2
  },
  'LANDSAT_7': {
    slope: 0.0000275, 
    offset: -0.2
  },
  'LANDSAT_8': {
    slope: 0.0000275, 
    offset: -0.2
  },
  'Sentinel-2A': {
    slope: 0.0001, 
    offset: 0
  },
  'Sentinel-2B': {
    slope: 0.0001, 
    offset: 0
  },
});

var bandNames = ee.Dictionary({
  'LANDSAT_5': {
    'REF': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
    'QA': 'QA_PIXEL'
  },
  'LANDSAT_7': {
    'REF': ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'],
    'QA': 'QA_PIXEL'
  },
  'LANDSAT_8': {
    'REF': ['SR_B2', 'SR_B3','SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'],
    'QA': 'QA_PIXEL'
  },
  'Sentinel-2A': {
    'REF': ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'],
    'QA': 'QA60'
  },
  'Sentinel-2B': {
    'REF': ['B2', 'B3', 'B4', 'B8', 'B11', 'B12'],
    'QA': 'QA60'
  }
});

var newBandNames = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2', 'BQA'];

var satelliteZenith = ee.Dictionary({
  'LANDSAT_5': 7.5,
  'LANDSAT_7': 7.5,
  'LANDSAT_8': 7.5,
  'Sentinel-2A': 10.3,
  'Sentinel-2B': 10.3,
});

var satelliteBandPassAdjustment = ee.Dictionary({
  'Sentinel-2A': { 
    'BLUE': {slope: 1.0443, offset: -0.0644},
    'GREEN': {slope: 1.1553, offset: -0.0388},
    'RED': {slope: 1.0278, offset: -0.0200},
    'NIR': {slope: 1.1393, offset: -0.0054},
    'SWIR1': {slope: 1.0257, offset: -0.0067},
    'SWIR2': {slope: 1.0514, offset: -0.0078},
    'BQA': {slope: 1, offset: 0 }
  },
  'Sentinel-2B': {
    'BLUE': {slope: 1.0443, offset: -0.0644},
    'GREEN': {slope: 1.1553, offset: -0.0388},
    'RED': {slope: 1.0278, offset: -0.0200},
    'NIR': {slope: 1.1393, offset: -0.0054},
    'SWIR1': {slope: 1.0257, offset: -0.0067},
    'SWIR2': {slope: 1.0514, offset: -0.0078},
    'BQA': {slope: 1, offset: 0 }
  }
});


function padronizeProperties(image){
  image = ee.Image(image);
  
  var landsatImage = image
    .set('SPACECRAFT', image.get('SPACECRAFT_ID'))
    .set('CLOUD_COVER', image.get('CLOUD_COVER'))
    .set('SOLAR_ZENITH_ANGLE', image.get('SUN_ELEVATION'))
    .set('SOLAR_AZIMUTH_ANGLE', image.get('SUN_AZIMUTH'));
  
  var copernicusImage = image
    .set('SPACECRAFT', image.get('SPACECRAFT_NAME'))
    .set('CLOUD_COVER', image.get('CLOUDY_PIXEL_PERCENTAGE'))
    .set('SOLAR_ZENITH_ANGLE', image.get('MEAN_SOLAR_ZENITH_ANGLE'))
    .set('SOLAR_AZIMUTH_ANGLE', image.get('MEAN_SOLAR_AZIMUTH_ANGLE'));
  
  
  return ee.Image(ee.Algorithms.If(
      ee.List(landsatSpacecrafts).containsAll([image.get('SPACECRAFT_ID')]),
      landsatImage,
      ee.Algorithms.If(
        ee.List(copernicusSpacecrafts).containsAll([image.get('SPACECRAFT_NAME')]),
        copernicusImage,
        image
      )
  ));
}

function padronizeBandNames(image){
  
  var spacecraft = image.get('SPACECRAFT');

  var oldBandNames = ee.Dictionary(bandNames.get(spacecraft));
  var bandScales = ee.Dictionary(satelliteBandReescale.get(spacecraft));
  
  var reescaledREFBands = image.select(oldBandNames.get('REF'))
    .multiply(bandScales.getNumber('slope'))
    .add(bandScales.getNumber('offset'));
  var qaBand = image.select([oldBandNames.get('QA')]);
  
  var renamedBands = reescaledREFBands.addBands(qaBand).rename(newBandNames);

  return ee.Image(renamedBands.copyProperties(image, image.propertyNames()));
}

function bandPassAjustment(image){ 
  var spacecraft = image.get('SPACECRAFT');
  
  var bandNames = image.bandNames();
  
  function convertToImage(e, l){
    e = ee.Number.parse(e);
    l = ee.List(l);
    return l.add(ee.Image.constant(e)); 
  }
  
  function convert(image){
    return ee.ImageCollection
      .fromImages(image.iterate(convertToImage, ee.List([])))
      .toBands().rename(bandNames);
  }
  
  var bandPass = ee.Dictionary(satelliteBandPassAdjustment.get(spacecraft, {}));
  
  var slope = bandPass.map(function(key, value){
    return ee.Dictionary(value).get('slope');
  }).values(newBandNames);
  
  var offset = bandPass.map(function(key, value){
    return ee.Dictionary(value).get('offset');
  }).values(newBandNames);
	
  var adjustedImage = ee.Image(ee.Algorithms.If(ee.Algorithms.IsEqual(bandPass, {}), 
    image, 
    image.multiply(convert(slope)).add(convert(offset))));
    
  return image.addBands(adjustedImage, null, true);
}


function brdfCorrect(rawImage) {
  // Script to correct Landsat data for BRDF effects using a c-factor approach as published in:
  // D.P. Roy, H.K. Zhang, J. Ju, J.L. Gomez-Dans, P.E. Lewis, C.B. Schaaf, Q. Sun, J. Li, H. Huang, V. Kovalskyy, 
  // A general method to normalize Landsat reflectance data to nadir BRDF adjusted reflectance, Remote Sensing of Environment, Volume 176, April 2016, Pages 255-271
  // Interpreted and coded here by Daniel Wiell and Erik Lindquist of the United Nations Food and Agriculture Organization

  var roi = rawImage.geometry(ee.ErrorMargin(300));

  var constants = {
    pi: Math.PI
  };

  var coefficientsByBand = {
    'BLUE': {fiso: 0.0774, fgeo: 0.0079, fvol: 0.0372},
    'GREEN': {fiso: 0.1306, fgeo: 0.0178, fvol: 0.0580},
    'RED': {fiso: 0.1690, fgeo: 0.0227, fvol: 0.0574},
    'NIR': {fiso: 0.3093, fgeo: 0.0330, fvol: 0.1535},
    'SWIR1': {fiso: 0.3430, fgeo: 0.0453, fvol: 0.1154},
    'SWIR2': {fiso: 0.2658, fgeo: 0.0387, fvol: 0.0639}
  };
  
  var corners = findCorners();

  // viewAngles
  var maxDistanceToSceneEdge = 210000;
  var maxSatelliteZenith = satelliteZenith.getNumber(rawImage.get('SPACECRAFT'));
  var upperCenter = pointBetween(corners.upperLeft, corners.upperRight);
  var lowerCenter = pointBetween(corners.lowerLeft, corners.lowerRight);
  var slope = slopeBetween(lowerCenter, upperCenter);
  var slopePerp = ee.Number(-1).divide(slope);
  
  // viewAz
  var viewAz = ee.Number(constants.pi / 2).subtract((slopePerp).atan());

  var leftLine = toLine(corners.upperLeft, corners.lowerLeft);
  var rightLine = toLine(corners.upperRight, corners.lowerRight);
  var leftDistance = ee.FeatureCollection(leftLine).distance(maxDistanceToSceneEdge);
  var rightDistance = ee.FeatureCollection(rightLine).distance(maxDistanceToSceneEdge);
  // Map.addLayer(rightDistance)
  
  var viewZenith = rightDistance.multiply(maxSatelliteZenith.multiply(2))
    .divide(rightDistance.add(leftDistance))
    .subtract(maxSatelliteZenith);

  // Map.addLayer(viewZenith, {}, 'viewZenith');
    
  //  viewZen
  var viewZen = viewZenith.multiply(constants.pi).divide(180);
  // Map.addLayer(viewZen, {}, 'viewZen');
  
  // Solar position
  // Ported from http://pythonfmask.org/en/latest/_modules/fmask/landsatangles.html
  var date = rawImage.date();
  var secondsInHour = 3600;
  
  var longDeg = ee.Image.pixelLonLat().select('longitude');
  
  var latRad = ee.Image.pixelLonLat().select('latitude').multiply(constants.pi).divide(180);
    
  var hourGMT = ee.Number(date.getRelative('second', 'day')).divide(secondsInHour);
  
  // Julian Date Proportion
  var jdp = ee.Number.parse(date.getFraction('year'));
  // print(jdp);
    
  // Julian Date Proportion in Radians
  var jdpr = jdp.multiply(2 * constants.pi);
  // print(jdpr);
  
  var meanSolarTime = toImage("hourGMT + longDeg / 15", 
    {hourGMT: hourGMT, longDeg: longDeg});
  // print(meanSolarTime);
  var localSolarDiff = ee.Image().expression('(0.000075 + 0.001868 * cos( jdpr ) - 0.032077 * sin( jdpr ) - 0.014615 * cos(2 * jdpr ) - 0.040849 * sin(2 * jdpr )) * 12 * 60 / pi', 
    {jdpr: jdpr, pi: constants.pi});
  // print(localSolarDiff);
  
  var  trueSolarTime = toImage('meanSolarTime + localSolarDiff / 60 - 12', 
    {meanSolarTime: meanSolarTime, localSolarDiff: localSolarDiff});
  // print(trueSolarTime);
  
  var angleHour = toImage('trueSolarTime * 15 * pi / 180', 
    {trueSolarTime: trueSolarTime, pi: constants.pi});
  // print(angleHour);
  
  var delta = toImage('0.006918 - 0.399912 * cos( jdpr ) + 0.070257 * sin( jdpr ) - 0.006758 * cos(2 * jdpr ) + 0.000907 * sin(2 * jdpr ) - 0.002697 * cos(3 * jdpr ) + 0.001480 * sin(3 * jdpr )', 
    {jdpr: jdpr});
  // print(delta);
    
  var sunZen = toImage('acos(sin(latRad) * sin(delta) + cos(latRad) * cos(delta) * cos(angleHour))', 
    {latRad: latRad, delta: delta, angleHour: angleHour});
  // Map.addLayer(sunZen, {}, 'sunZen');

  
  var sinSunAzSW = toImage('cos(delta) * sin(angleHour) / sin(sunZen)', 
    {delta: delta, angleHour: angleHour, sunZen: sunZen})
    .clamp(-1, 1);
  // print(sinSunAzSW);
  
  var cosSunAzSW = toImage('(-cos(latRad) * sin(delta) + sin(latRad) * cos(delta) * cos(angleHour)) / sin(sunZen)', 
    {latRad: latRad, delta: delta, angleHour: angleHour, sunZen: sunZen});
  // print(cosSunAzSW);
    
  var sunAzSW = toImage(' asin(sinSunAzSW) ', 
    {sinSunAzSW: sinSunAzSW});
  // print(sunAzSW);
  
  
  sunAzSW =  sunAzSW.where(cosSunAzSW.lte(0), toImage('pi - sunAzSW', {pi: constants.pi, sunAzSW: sunAzSW}));
  // print(sunAzSW);
    
  sunAzSW = sunAzSW.where(cosSunAzSW.gt(0).and(sinSunAzSW.lte(0)), toImage("2 * pi + sunAzSW", {pi: constants.pi, sunAzSW: sunAzSW}));
  // print(sunAzSW);
  
  var sunAz = sunAzSW.add(constants.pi);
  // Map.addLayer(sunAz)
  // print(sunAz);
  
  sunAz = sunAz.where(sunAz.gt(2 * constants.pi), sunAz.subtract(2 * constants.pi));
  // print(sunAz);

  // sunZenOut
      // https://nex.nasa.gov/nex/static/media/publication/HLS.v1.0.UserGuide.pdf
  var centerLat = ee.Number.parse(roi.bounds(ee.ErrorMargin(30)).centroid(ee.ErrorMargin(30)).coordinates().get(0)).multiply(constants.pi).divide(180);
  // print(centerLat);

  var sunZenOut = ee.Image().expression(
    '(31.0076 - 0.1272 * centerLat + 0.01187 * pow(centerLat, 2) + 2.40E-05 * pow(centerLat, 3) - 9.48E-07 * pow(centerLat, 4) - 1.95E-09 * pow(centerLat, 5) + 6.15E-11 * pow(centerLat, 6)) * pi/180',
    {centerLat: centerLat, pi: constants.pi});
  // print(sunZenOut);

  var relativeSunViewAz = sunAz.subtract(viewAz);
  // print(relativeSunViewAz);
  
  var kvol = rossThick(sunZen, viewZen, relativeSunViewAz);
  // print(kvol)
  // Map.addLayer(kvol, {}, "kvol", false);
  
  var kvol0 = rossThick(sunZenOut, 0, 0);
  // print(kvol0);
  // Map.addLayer(kvol0, {}, "kvol0", false);
  
  var cFactors = getCFactors(kvol, kvol0);
  // print(cFactors);
  // Map.addLayer(cFactors)
  
  var bqa = rawImage.select("BQA").int();
  
  var ajustedImage = rawImage.select('BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2').multiply(cFactors);
  
  ajustedImage = rawImage.addBands(ajustedImage, null, true);
  
  return ajustedImage.addBands(bqa, null, true);

  function rossThick(sunZen, viewZen, relativeSunViewAz) {
    var args = {
      sunZen: sunZen, 
      viewZen: viewZen, 
      relativeSunViewAz: relativeSunViewAz,
      pi: Math.PI
    };
    var cosPhaseAngle = cosPhaseAngleFunc(sunZen, viewZen, relativeSunViewAz);
    args.cosPhaseAngle = cosPhaseAngle;
    
    var phaseAngle = toImage(' acos(cosPhaseAngle) ', {cosPhaseAngle: cosPhaseAngle});
    args.phaseAngle = phaseAngle;
    
    var t = toImage('((pi/2 - phaseAngle) * cosPhaseAngle + sin(phaseAngle)) / (cos(sunZen) + cos(viewZen)) - pi/4', args);
    return t;
  }

  function cosPhaseAngleFunc(sunZen, viewZen, relativeSunViewAz) {
    var args = {
      sunZen: sunZen,
      viewZen: viewZen,
      relativeSunViewAz: relativeSunViewAz
    };
    
    return toImage('cos(sunZen) * cos(viewZen) + sin(sunZen) * sin(viewZen) * cos(relativeSunViewAz)', args)
      .clamp(-1, 1);
  }

  function getCFactors(kvol, kvol0) {
    var cFactors = [];
    for (var bandName in coefficientsByBand){
      var cFactor = getCFactor(bandName, coefficientsByBand[bandName], kvol, kvol0);
      cFactors.push(cFactor);
    }
    return ee.Image.cat(cFactors);
  }

  function getCFactor(bandName, coefficients, kvol, kvol0) {
    var kgeo = -1;
    var kgeo0 = -1;
    
    var brdf = brdfFunc(kvol, kgeo, coefficients);
    var brdf0 = brdfFunc(kvol0, kgeo0, coefficients);
    // Map.addLayer(brdf, {}, "brdf", false);
    // Map.addLayer(brdf0, {}, "brdf0", false);
    
    var cFactor = brdf0.divide(brdf).rename(bandName);
    // Map.addLayer(cFactor, {},"", false);
    return cFactor;
  }

  function brdfFunc(kvolBand, kgeoBand, coefficients) {
    var args = merge(coefficients, {
      // kvol: 'i.' + kvolBand,
      kvol: kvolBand.multiply(3),     // check this multiplication factor.  Is there an 'optimal' value?  Without a factor here, there is not enough correction.
      kgeo: kgeoBand
    });
    
    return toImage('fiso + fvol * kvol + fgeo * kvol', args);
  }

  function findCorners() {
    var footprint = ee.Geometry(rawImage.get('system:footprint'));
    var bounds = ee.List(footprint.bounds().coordinates().get(0));
    var coords = footprint.coordinates();
    
    var xs = coords.map(function (item) {
      return x(item);
    });
    
    var ys = coords.map(function (item) {
      return y(item);
    });

    function findCorner(targetValue, values) {
      var diff = values.map(function (value) {
        return ee.Number(value).subtract(targetValue).abs();
      });
      
      var minValue = diff.reduce(ee.Reducer.min());
      var idx = diff.indexOf(minValue);
      return coords.get(idx);
    }

    var lowerLeft = findCorner(x(bounds.get(0)), xs);
    var lowerRight = findCorner(y(bounds.get(1)), ys);
    var upperRight = findCorner(x(bounds.get(2)), xs);
    var upperLeft = findCorner(y(bounds.get(3)), ys);
    
    return {
      upperLeft: upperLeft,
      upperRight: upperRight,
      lowerRight: lowerRight,
      lowerLeft: lowerLeft
    };
    
  }

  function x(point) {
    return ee.Number(ee.List(point).get(0));
  }

  function y(point) {
    return ee.Number(ee.List(point).get(1));
  }

  function pointBetween(pointA, pointB) {
    return ee.Geometry.LineString([pointA, pointB]).centroid(ee.ErrorMargin(30)).coordinates();
  }

  function slopeBetween(pointA, pointB) {
    return ((y(pointA)).subtract(y(pointB))).divide((x(pointA)).subtract(x(pointB)));
  }

  function toLine(pointA, pointB) {
    return ee.Geometry.LineString([pointA, pointB]);
  }

// ************** COMMON HELPERS **************

  function toImage(band, args) {
    if ((typeof band) === 'string') {
        return ee.Image().expression(band, args).float();
    }else{
      return ee.Image(band).float(); 
    }
  }
  
  function merge(o1, o2) {
    function addAll(target, toAdd) {
      for (var key in toAdd) target[key] = toAdd[key];
    }

    var result = {};
    addAll(result, o1);
    addAll(result, o2);
    return result;
  }
}

function terrainCorrection(image) {
  // Credits: https://mygeoblog.com/2018/10/17/terrain-correction-in-gee
  
  var scale = 300;
 
  // get terrain layers
  var dem = ee.Image('USGS/SRTMGL1_003');
  var degree2radian = 0.01745;
 
  var imageCondition = illuminationCondition(image);
  var correctedImage = illuminationCorrection(imageCondition);

  return correctedImage.select(newBandNames);
  
   
  ////////////////////////////////////////////////////////////////////////////////
  // Function to calculate illumination condition (IC). Function by Patrick Burns and Matt Macander
  function illuminationCondition(img){
 
    // Extract image metadata about solar position
    var SZ_rad = ee.Image.constant(ee.Number(img.get('SOLAR_ZENITH_ANGLE'))).multiply(3.14159265359).divide(180);
    var SA_rad = ee.Image.constant(ee.Number(img.get('SOLAR_AZIMUTH_ANGLE')).multiply(3.14159265359).divide(180));
  
    // Creat terrain layers
    var slp = ee.Terrain.slope(dem);
    var slp_rad = ee.Terrain.slope(dem).multiply(3.14159265359).divide(180);
    var asp_rad = ee.Terrain.aspect(dem).multiply(3.14159265359).divide(180);
   
    // Calculate the Illumination Condition (IC)
    // slope part of the illumination condition
    var cosZ = SZ_rad.cos();
    var cosS = slp_rad.cos();
    var slope_illumination = cosS.expression("cosZ * cosS",
                                            {'cosZ': cosZ,
                                             'cosS': cosS.select('slope')});
    // aspect part of the illumination condition
    var sinZ = SZ_rad.sin();
    var sinS = slp_rad.sin();
    var cosAziDiff = (SA_rad.subtract(asp_rad)).cos();
    var aspect_illumination = sinZ.expression("sinZ * sinS * cosAziDiff",
                                             {'sinZ': sinZ,
                                              'sinS': sinS,
                                              'cosAziDiff': cosAziDiff});
    // full illumination condition (IC)
    var ic = slope_illumination.add(aspect_illumination);
   
    // Add IC to original image
    var img_plus_ic = ee.Image(img.addBands(ic.rename('IC')).addBands(cosZ.rename('cosZ')).addBands(cosS.rename('cosS')).addBands(slp.rename('slope')));
    return img_plus_ic.clip(img.geometry());
  }
 
  // Function to apply the Sun-Canopy-Sensor + C (SCSc) correction method to each
  // image. Function by Patrick Burns and Matt Macander

  function illuminationCorrection(img){
      var props = img.toDictionary();
      var st = img.get('system:time_start');
   
      var img_plus_ic = img;
      var mask1 = img_plus_ic.select('NIR').gt(-0.1);
      var mask2 = img_plus_ic.select('slope').gte(5)
                              .and(img_plus_ic.select('IC').gte(0))
                              .and(img_plus_ic.select('NIR').gt(-0.1));
      var img_plus_ic_mask2 = ee.Image(img_plus_ic.updateMask(mask2));
   
      // Specify Bands to topographically correct
      var bandList = newBandNames.filter(function(item) {
          return item !== 'BQA';
      });
      var compositeBands = img.bandNames();
      var nonCorrectBands = img.select(compositeBands.removeAll(bandList));
   
      var geom = ee.Geometry(img.get('system:footprint')).bounds().buffer(10000);
   
      function apply_SCSccorr(band){
        var method = 'SCSc';
        var out = img_plus_ic_mask2
          .select('IC', band)
          .reduceRegion({
            reducer: ee.Reducer.linearFit(), // Compute coefficients: a(slope), b(offset), c(b/a)
            geometry: ee.Geometry(img.geometry().buffer(-5000)),
            scale: scale,
            maxPixels: 1E13
          }); 
          
        if (out === null || out === undefined ){
          return img_plus_ic_mask2.select(band);
        }else{
          var out_a = ee.Number(out.get('scale'));
          var out_b = ee.Number(out.get('offset'));
          
          var out_c = ee.Number(ee.Algorithms.If(ee.Algorithms.IsEqual(out_a, null), 
                                1, ee.Algorithms.If(ee.Algorithms.IsEqual(out_b, null), 
                                1, out_b.divide(out_a))));
          // Apply the SCSc correction
          var SCSc_output = img_plus_ic_mask2.expression(
            "((image * (cosB * cosZ + cvalue)) / (ic + cvalue))", {
            'image': img_plus_ic_mask2.select(band),
            'ic': img_plus_ic_mask2.select('IC'),
            'cosB': img_plus_ic_mask2.select('cosS'),
            'cosZ': img_plus_ic_mask2.select('cosZ'),
            'cvalue': out_c
          });
     
          return SCSc_output;
        }
      }
   
      var img_SCSccorr = ee.Image(bandList.map(apply_SCSccorr)).addBands(img_plus_ic.select('IC'));
      var bandList_IC = ee.List([bandList, 'IC']).flatten();
      img_SCSccorr = img_SCSccorr.unmask(img_plus_ic.select(bandList_IC)).select(bandList);
   
      return ee.Image(img_SCSccorr
        .addBands(nonCorrectBands)
        .setMulti(props)
        .copyProperties(image, image.propertyNames()));
    }
 
}  

function landsatCloudMask(image){
  // Cloud Mask
  var landsatQABits = ee.List([
    [0, 0, 0], // Designated Fill
    [6, 6, 1], // Clear
    [8, 9, 1], //  Cloud Confidence
    [10, 11, 1], // Cloud Shadow Confidence
    [14, 15, 1] // Cirrus Confidence
  ]);

  function getQABits(image, start, end) {
    var pattern = ee.Number(ee.List.sequence(start, end).distinct().iterate(function(i, pattern){
      i = ee.Number(i);
      pattern = ee.Number(pattern);
  
      return pattern.add(ee.Number(2).pow(i));
    }, ee.Number(0)));
  
    return image.select(0).bitwiseAnd(pattern.int()).rightShift(start);
  }
  
  var bqa = image.select('BQA');


  var inital_state = ee.Dictionary({
    'bqa': bqa,
    'mask': bqa.eq(bqa).not()
  });

  var finalState = ee.Dictionary(landsatQABits.iterate(function(bits, state){
    bits = ee.List(bits);
    state = ee.Dictionary(state);

    var bqa = ee.Image(state.get('bqa'));
    var mask = ee.Image(state.get('mask'));

    var start = bits.getNumber(0);
    var end = bits.getNumber(1);
    var desired = bits.getNumber(2);

    var blueprint = getQABits(bqa, start, end).eq(desired);

    return ee.Dictionary({
        'bqa': bqa,
        'mask': mask.blend(blueprint.not().or(mask))
    });

  }, inital_state));

  var isClouds = ee.Image(finalState.get('mask'));
  return isClouds; 
}

function sentinelCloudMask(image){
  var sentinelCloudCoverProbability = 40;
  // Identify dark NIR pixels (potential cloud shadow pixels).
  var NIR_DRK_THRESH = 0.15;
  var CLD_PRJ_DIST = 2;
  var BUFFER = 50;

  var s2CloudsImage = ee.Image(ee.ImageCollection('COPERNICUS/S2_CLOUD_PROBABILITY')
    .filterMetadata('system:index', 'equals', image.get('system:index'))
    .first()).unmask();
    
  // 1: clouds 0: non-clouds.
  var isClouds = s2CloudsImage.gte(sentinelCloudCoverProbability);

  var darkPixels = image.select('NIR').lt(NIR_DRK_THRESH).rename('dark_pixels');
  
  // Determine the direction to project cloud shadow from clouds (assumes UTM projection).
  var shadowAzimuth = ee.Number(90).subtract(ee.Number(image.get('SOLAR_AZIMUTH_ANGLE')));

  // Project shadows from clouds for the distance specified by the CLD_PRJ_DIST input.
  var isCloudsProj = isClouds.directionalDistanceTransform(shadowAzimuth, CLD_PRJ_DIST*10)
        .reproject({'crs': image.select(0).projection(), 'scale': 100})
        .select('distance')
        .mask()
        .rename('cloud_transform');

  // Identify the intersection of dark pixels with cloud shadow projection.
  var isShadows = isCloudsProj.multiply(darkPixels).rename('shadows');
  
  // Combine cloud and shadow mask, set cloud and shadow as value 1, else 0.
  var cloudsOrCloudShadows = isClouds.or(isShadows);
  
  // Remove small cloud-shadow patches and dilate remaining pixels by BUFFER input.
  // 20 m scale is for speed, and assumes clouds don't require 10 m precision.
  cloudsOrCloudShadows = cloudsOrCloudShadows
    .focal_min(2)
    .focal_max(BUFFER*2/20)
    .reproject({'crs': image.select([0]).projection(), 'scale': 20})
    .rename('cloudmask');
    
  return cloudsOrCloudShadows;
}

function cloudMask1(image){
  
  // 1: clouds 
  // 0: Non-clouds.
  return ee.Image(ee.Image(ee.Algorithms.If(
     ee.List(['Sentinel-2A', 'Sentinel-2B']).containsAll([image.get('SPACECRAFT')]),
    sentinelCloudMask(image),
    landsatCloudMask(image)
  )).copyProperties(image, image.propertyNames()));
}

function cloudMask2(image){
  /*
    1 - Land                  : Clear-sky land observation;
    2 - Water                 : Clear-sky water observation;
    3 - Cloud                 : Cloud detected;
    4 - Cloud Shadow          : Shadow detected.
    5 - Topographic shadow    : Shadow detected. The pixel located outside cloud projections and within estimated topographic shadow (estimated using DEM and solar elevation and azimuth);
    6 - Snow/Ice              : Snow or ice detected;
    7 - Haze                  : Dense semi-transparent clouds/fog detected;
    8 - Sensor Degradation    : Errors in image acquisition caused by sensor degradation;
  */
  
  var properties = ['BLUE', 'GREEN', 'RED', 'NIR', 'SWIR1', 'SWIR2'];
  
  var samples = ee.FeatureCollection('users/saraiva/puc_artificial_intelligence/qualityBand/samples');

  var classifier = ee.Classifier
    .smileRandomForest({numberOfTrees: 10, seed: 42})
    .train({
      features: samples, 
      classProperty: 'class', 
      inputProperties: properties
    });
    
  var qaBand = image
    .select(properties)
    .classify(classifier);
  
  var elevation = ee.Image("USGS/SRTMGL1_003");
  
  var shadowMap = ee.Terrain.hillShadow({
    image: elevation,
    azimuth: image.get('SOLAR_AZIMUTH_ANGLE'),
    zenith: image.get('SOLAR_ZENITH_ANGLE'),
    neighborhoodSize: 200,
    hysteresis: true
  });
  
  qaBand = qaBand.blend(shadowMap.not().selfMask().remap([1], [5]));
    
  var badPixels = qaBand.remap([3, 4, 7, 8], [1, 1, 1, 1], 0);
  
  // 1: clouds 0: non-clouds.
  return badPixels;
}

function maskImage(image){
  var mask1 = cloudMask1(image);
  var mask2 = cloudMask2(image);
  var combinedMask = mask1;//.or(mask2);
  return image.updateMask(combinedMask.not());
}

function normalizeImage(image){
  return  maskImage(terrainCorrection(brdfCorrect(bandPassAjustment(padronizeBandNames(padronizeProperties(image))))));
}

function getCollection(satelite, roi, startDate, endDate, cloudCover){
  
  var col = ee.ImageCollection(satelliteCollections[satelite])
      .filterDate(startDate, endDate)
      .filterBounds(roi)
      .filterMetadata('system:asset_size', 'greater_than', 1024 * 1024 * 36) // > 36MB
      .filter(ee.Filter.or(ee.Filter.lt('CLOUD_COVER', cloudCover), 
                            ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', cloudCover)))
      .map(normalizeImage);
      
  return col;
}


function smoorthing(collection, interval){
  var timeField = 'system:time_start';
  var date = ee.Date(collection.get(timeField));
  
  var join = ee.Join.saveAll({
    matchesKey: 'images'
  });
  
  var diffFilter = ee.Filter.maxDifference({
    difference: 1000 * 60 * 60 * 24 * interval,
    leftField: timeField, 
    rightField: timeField
  });
  
  var threeNeighborJoin = join.apply({
    primary: collection, 
    secondary: collection, 
    condition: diffFilter
  });
  
  var smoothed = ee.ImageCollection(threeNeighborJoin.map(function(image) {
    image = ee.Image(image);
    
    var col = ee.ImageCollection.fromImages(image.get('images'));
    
    var index = ee.Algorithms.String(ee.String('BAND_').cat(ee.String(ee.Image(image).date().format('yyyy-MM-dd'))));
    
    var smoothedImage = ee.Image(col.mean().copyProperties(image, image.propertyNames())).set('date', index)
      
    return smoothedImage;
    
  }).distinct('date'));

  return smoothed;
}

function getIntervals(startDateInput, endDateInput){
  var dates = ee.List([startDateInput.get('year'), endDateInput.get('year')])
    .distinct()
    .iterate(function(year, list){
      year = ee.Number.parse(year);
      list = ee.List(list);
      
      for(var i=1; i<=23; i++){
        var startDOY = ((i * 16) - 15);
        var endDOY = (((i+1) * 16) - 15);
        
        var startDate = ee.Date.parse('Y D', ee.String(year).cat(' ' + startDOY));
        var endDate = ee.Date.parse('Y D', ee.String(year).cat(' ' + endDOY));
        
        if(i >= 5){
          startDate = ee.Date.parse('Y D', ee.String(year).cat(' ' +(startDOY-1)));
          if(i >= 23){
            endDate = ee.Date.parse('Y D', ee.String(year.add(1)).cat((' ' + 1)));
          }else{
            endDate = ee.Date.parse('Y D', ee.String(year).cat(' ' + (endDOY-1)));
          }
        }
        
        var feature = ee.Feature(null, {
          'startDate': startDate,
          'endDate': endDate,
          'interval': i
        });
        list = list.add(feature);
      }
      return list;
    }, ee.List([]));
  
  dates = ee.List(dates)
    .filter(ee.Filter.greaterThanOrEquals('startDate', startDateInput))
    .filter(ee.Filter.lessThanOrEquals('endDate', ee.Date(endDateInput).advance(1, 'day')));
    
  return dates;
}

function get16DayProducts(collection, startDateCol, endDateCol, qualityBand){
  startDateCol = ee.Date(startDateCol);
  endDateCol = ee.Date(endDateCol);
  
  var intervals = ee.List(getIntervals(startDateCol, endDateCol));
  
  var result = ee.ImageCollection.fromImages(intervals.iterate(function(feature, products){
    feature = ee.Feature(feature);
    products = ee.List(products);
    
    
    var startD = ee.Date(feature.get('startDate'));
    var endD = ee.Date(feature.get('endDate'));
    
    var filteredNormCollection = collection.filterDate(startD, endD);

    var filteredNormMosaic = filteredNormCollection.qualityMosaic(qualityBand)
      .set({
        'system:index': ee.String('BAND_').cat(ee.String(startD.format('yyyyMMdd'))),
        'system:time_start': startD.millis(),
        'SPACECRAFT_ID': 'ARD',
        'START_DATE': startD,
        'END_DATE': endD,
        'INTERVAL': feature.get('interval')
      });

    var finalProduct = ee.Image(ee.Algorithms.If(filteredNormCollection.size().gte(1), filteredNormMosaic, ee.Image(0)));
      
    return products.add(finalProduct);
    
  }, ee.List([])));
  
  var validProducts = result.filterMetadata('SPACECRAFT_ID', 'equals', 'ARD');

  return validProducts;
}


function visualizeImage(image, roi){
  var stats = image.reduceRegion({
      reducer: ee.Reducer.percentile([2, 98]),
      geometry: roi, 
      scale: 500,
      maxPixels: 1E9
    });
  
  var imageVis = ee.Image(image.bandNames().iterate(function(bandName, currentImage){
      bandName = ee.String(bandName);
      currentImage = ee.Image(currentImage);
      
      var p2 = bandName.cat(ee.String('_p2'));
      var p98 = bandName.cat(ee.String('_p98'));
      
      return currentImage
        .addBands(currentImage.select(bandName).clamp(stats.getNumber(p2), stats.getNumber(p98)), null, true);
        
    }, image)).visualize();
  return imageVis;
} 

exports.satelliteCollections = satelliteCollections;
exports.padronizeProperties = padronizeProperties;
exports.padronizeBandNames = padronizeBandNames;
exports.bandPassAjustment = bandPassAjustment;
exports.applyBRDFCorrection = brdfCorrect;
exports.terrainCorrection = terrainCorrection;
exports.maskImage = maskImage;
exports.normalizeImage = normalizeImage;
exports.smoorthing = smoorthing;
exports.get16DayProducts = get16DayProducts;
exports.getCollection = getCollection;
exports.visualizeImage = visualizeImage;
