// =================================================================
// Coastal Blue Carbon — AOI Covariate Extraction & Prior Mapping
// Canada Coastal Ecosystems (Tidal Marsh + Seagrass)
// =================================================================
// Version: 1.0  (adapted from Forest Carbon Assessment v4.2)
//
// PURPOSE:
//   Step-by-step Google Earth Engine script to:
//   1. Define your AOI (draw on map or provide an asset path)
//   2. Load and display SoilGrids SOC as a spatial prior
//   3. Build the canonical 27-band coastal covariate stack
//   4. Export covariate raster + Google Satellite Embedding to Drive
//   5. Generate stratified soil sampling points
//   6. Export reports and sampling point files
//
// CANONICAL COVARIATE STACK — 27 bands (identical in Python notebooks):
//   Group 1 — Topography & Tidal Position (7):
//     elevation_m, slope, elevationRelMHW, twi, tpi,
//     tidal_flat_prob, coastal_dist_m
//   Group 2 — Sentinel-1 SAR (3):
//     VV_mean, VH_mean, VVVH_ratio
//   Group 3 — Sentinel-2 Optical (12):
//     B, G, R, NIR, SWIR1, SWIR2,
//     NDVI_median, LSWI_median, mNDWI_median,
//     EVI_median, SAVI_median, tidal_wetness
//   Group 4 — Climate (2):
//     MAT_C, MAP_mm
//   Group 5 — SoilGrids SOC Prior (3):
//     sg_soc_0_30cm, sg_soc_30_100cm, sg_soc_0_100cm
//
// NOTE on SoilGrids:
//   SoilGrids is trained primarily on upland mineral soils and has
//   known limitations in tidal sediment environments. It is used here
//   as a weak spatial prior only. Uncertainty inflation is applied
//   downstream in the Bayesian R workflow
//   (P3_0c_bayesian_prior_setup_bluecarbon.R).
//
// WORKFLOW:
//   Run steps ①–⑥ in order using the buttons in the side panel.
//   Each step populates global variables used by subsequent steps.
// =================================================================


// ─────────────────────────────────────────────────────────────────
// SECTION A — CONFIGURATION
// ─────────────────────────────────────────────────────────────────

// ── AOI ──────────────────────────────────────────────────────────
// Option 1: Set AOI_ASSET to your GEE feature collection asset path, e.g.:
//   var AOI_ASSET = 'projects/your-project/assets/your_coastal_boundary';
// Option 2: Leave as null — a polygon drawing tool will activate in Step ①.
var AOI_ASSET = null;

// ── Export settings ──────────────────────────────────────────────
var EXPORT_CRS    = 'EPSG:3347';            // Canada Albers Equal Area (all Canada)
var EXPORT_SCALE  = 25;                     // metres — covariate snapshot
var EMBED_SCALE   = 10;                     // metres — Google Satellite Embedding (native)
var EXPORT_FOLDER = 'BlueCarbon_CoastalBC'; // Google Drive destination folder
var PROJECT_YEAR  = '2020_2023';            // label used in export filenames

// ── Sentinel-2 date range & cloud threshold ───────────────────────
var S2_START          = '2020-01-01';
var S2_END            = '2023-12-31';
var S2_CLOUD_THRESHOLD = 20;    // max cloud cover % per image

// ── Sentinel-1 SAR date range ────────────────────────────────────
var SAR_START = '2020-01-01';
var SAR_END   = '2023-12-31';

// ── TerraClimate date range ──────────────────────────────────────
var TC_START = '2000-01-01';
var TC_END   = '2022-12-31';

// ── Google Satellite Embedding ───────────────────────────────────
// null = median composite across all available years (recommended)
// Set to an integer (e.g. 2022) to use a single year.
var EMBEDDING_YEAR = null;

// ── Sampling ─────────────────────────────────────────────────────
var N_SOIL_SAMPLES = 100;   // total stratified sampling points to generate

// ── AOI display buffer ───────────────────────────────────────────
var AOI_BUFFER_M = 10000;   // metres — used only for map layer display


// ─────────────────────────────────────────────────────────────────
// SECTION B — GLOBAL STATE VARIABLES
// (populated sequentially as steps are run)
// ─────────────────────────────────────────────────────────────────
var aoi          = null;  // ee.Geometry — project boundary
var aoi_display  = null;  // ee.Geometry — buffered for display
var sg_soc_prior = null;  // SoilGrids SOC prior (3 canonical bands)
var cov_stack    = null;  // canonical 27-band covariate image
var embed_img    = null;  // Google Satellite Embedding V1 (64 bands)
var soil_pts     = null;  // stratified sampling points (FeatureCollection)
var snapshot_btn = null;  // Step ④ button — inserted after Step ③ runs


// ─────────────────────────────────────────────────────────────────
// SECTION C — HELPER FUNCTIONS
// ─────────────────────────────────────────────────────────────────

// Mask Sentinel-2 clouds using QA60 band
function maskS2clouds(image) {
  var qa = image.select('QA60');
  var cloudBitMask  = 1 << 10;
  var cirrusBitMask = 1 << 11;
  var mask = qa.bitwiseAnd(cloudBitMask).eq(0)
               .and(qa.bitwiseAnd(cirrusBitMask).eq(0));
  return image.updateMask(mask).divide(10000)
    .copyProperties(image, ['system:time_start']);
}

// Rename Sentinel-2 bands and compute coastal spectral indices
function addS2IndicesAndRename(image) {
  var img = image.select(
    ['B2',  'B3', 'B4', 'B8',  'B11',   'B12'],
    ['B',   'G',  'R',  'NIR', 'SWIR1', 'SWIR2']
  );
  var ndvi  = img.normalizedDifference(['NIR', 'R']).rename('NDVI_median');
  var lswi  = img.normalizedDifference(['NIR', 'SWIR1']).rename('LSWI_median');
  var mndwi = img.normalizedDifference(['G', 'SWIR1']).rename('mNDWI_median');
  var evi   = img.expression(
    '2.5 * ((NIR - R) / (NIR + 6*R - 7.5*B + 1))',
    {NIR: img.select('NIR'), R: img.select('R'), B: img.select('B')}
  ).rename('EVI_median');
  var savi  = img.expression(
    '1.5 * (NIR - R) / (NIR + R + 0.5)',
    {NIR: img.select('NIR'), R: img.select('R')}
  ).rename('SAVI_median');
  // Tasseled Cap Wetness (Sentinel-2 SR coefficients, Nedkov 2017)
  // Wetness is the only TC component used — most sensitive to inundation.
  // Brightness and Greenness are excluded (more relevant to forest/upland).
  var tc_wet = img.expression(
    '0.1511*B + 0.1973*G + 0.3283*R + 0.3407*NIR + (-0.7117)*SWIR1 + (-0.4559)*SWIR2',
    {B: img.select('B'), G: img.select('G'), R: img.select('R'),
     NIR: img.select('NIR'), SWIR1: img.select('SWIR1'), SWIR2: img.select('SWIR2')}
  ).rename('tidal_wetness');
  return img.addBands([ndvi, lswi, mndwi, evi, savi, tc_wet]);
}

// Compute TWI — Topographic Wetness Index
// TWI = ln(contributing_area / tan(slope))
// High TWI = flat/concave, water-accumulating terrain (marsh hollows, channels)
function computeTWI(dem) {
  var slope_rad = ee.Terrain.slope(dem).multiply(Math.PI / 180);
  var tan_slope = slope_rad.tan().max(0.001);  // guard against log(0)
  // Approximate upslope contributing area using focal pixel count
  var contrib = dem.gte(-9999).unmask(0).reduceNeighborhood({
    reducer: ee.Reducer.sum(),
    kernel: ee.Kernel.circle({radius: 20, units: 'pixels'})
  }).max(1);
  return contrib.divide(tan_slope).log().rename('twi');
}

// Compute TPI — Topographic Position Index
// TPI = elevation - focal mean elevation (300 m window)
// Positive = ridge/levee; Negative = marsh hollow/channel; ~0 = flat marsh plain
function computeTPI(dem) {
  var focal_mean = dem.focalMean({radius: 300, units: 'meters'});
  return dem.subtract(focal_mean).rename('tpi');
}


// ─────────────────────────────────────────────────────────────────
// SECTION D — UI LAYOUT
// ─────────────────────────────────────────────────────────────────

var STATUS_LABEL = ui.Label({
  value: 'Ready. Run steps ①–⑥ in order.',
  style: {color: 'gray', fontSize: '12px', whiteSpace: 'pre'}
});

function setStatus(msg) { STATUS_LABEL.setValue(msg); }

var panel = ui.Panel({style: {width: '370px', padding: '10px'}});

panel.add(ui.Label({
  value: 'Coastal Blue Carbon — Covariate Extraction',
  style: {fontWeight: 'bold', fontSize: '15px', margin: '0 0 4px 0'}
}));
panel.add(ui.Label({
  value: 'Canada: Tidal Marsh + Seagrass  |  VM0033',
  style: {color: '#2e7d32', fontSize: '12px', margin: '0 0 8px 0'}
}));
panel.add(ui.Label('─────────────────────────────────', {color: '#ccc'}));

panel.add(ui.Label('Configuration', {fontWeight: 'bold', fontSize: '12px'}));
panel.add(ui.Label('  CRS:   EPSG:3347 (Canada Albers Equal Area)', {fontSize: '11px', color: '#555'}));
panel.add(ui.Label('  Scale: ' + EXPORT_SCALE + ' m (covariate) / ' + EMBED_SCALE + ' m (embedding)', {fontSize: '11px', color: '#555'}));
panel.add(ui.Label('  S2:    ' + S2_START + ' → ' + S2_END, {fontSize: '11px', color: '#555'}));
panel.add(ui.Label('  SAR:   ' + SAR_START + ' → ' + SAR_END, {fontSize: '11px', color: '#555'}));
panel.add(ui.Label('  TC:    ' + TC_START + ' → ' + TC_END, {fontSize: '11px', color: '#555'}));
panel.add(ui.Label('  N sampling points: ' + N_SOIL_SAMPLES, {fontSize: '11px', color: '#555'}));
panel.add(ui.Label('─────────────────────────────────', {color: '#ccc'}));

var btn1 = ui.Button('① Import AOI');
var btn2 = ui.Button('② Import Raster Priors (SoilGrids SOC)');
var btn3 = ui.Button('③ Build Coastal Covariate Stack (27 bands)');
// btn4 inserted dynamically after Step ③ completes
var btn5 = ui.Button('⑤ Generate Stratified Sampling Points');
var btn6 = ui.Button('⑥ Reports and Exports');

panel.add(btn1); panel.add(btn2); panel.add(btn3);
panel.add(btn5); panel.add(btn6);
panel.add(ui.Label('─────────────────────────────────', {color: '#ccc'}));
panel.add(ui.Label('Status:', {fontWeight: 'bold', fontSize: '12px'}));
panel.add(STATUS_LABEL);

ui.root.add(panel);


// ─────────────────────────────────────────────────────────────────
// STEP 1 — IMPORT AOI
// ─────────────────────────────────────────────────────────────────
btn1.onClick(function step1_importAOI() {
  setStatus('Step ①: Loading AOI…');

  if (AOI_ASSET !== null) {
    // Load from GEE asset
    var aoiFC = ee.FeatureCollection(AOI_ASSET);
    aoi = aoiFC.union().geometry();
    aoi_display = aoi.buffer(AOI_BUFFER_M);
    Map.centerObject(aoi, 11);
    Map.addLayer(aoiFC, {color: '1a9641', fillColor: '00000000', width: 2}, 'AOI Boundary');
    setStatus('Step ① ✓ — AOI loaded from asset.\n→ Run Step ②.');

  } else {
    // Interactive drawing mode
    Map.drawingTools().setShown(true);
    Map.drawingTools().setShape('polygon');
    Map.drawingTools().layers().reset();
    setStatus('Step ① — Draw your AOI on the map.\n' +
              'Click vertices, then click the first\nvertex again to close the polygon.\n' +
              'Then click ① again to confirm.');

    var layers = Map.drawingTools().layers();
    if (layers.length() > 0) {
      aoi = layers.get(0).toGeometry();
      aoi_display = aoi.buffer(AOI_BUFFER_M);
      Map.centerObject(aoi, 11);
      setStatus('Step ① ✓ — AOI set from drawn polygon.\n→ Run Step ②.');
    }
  }

  // Canada boundary for reference (display only)
  var canada = ee.FeatureCollection('USDOS/LSIB_SIMPLE/2017')
    .filter(ee.Filter.eq('country_na', 'Canada'));
  Map.addLayer(canada, {color: 'aaaaaa', fillColor: '00000000', width: 1},
               'Canada Boundary', false);
});


// ─────────────────────────────────────────────────────────────────
// STEP 2 — IMPORT RASTER PRIORS (SoilGrids SOC)
// ─────────────────────────────────────────────────────────────────
btn2.onClick(function step2_importPriors() {
  if (!aoi) { setStatus('⚠ Run Step ① first.'); return; }
  setStatus('Step ②: Loading SoilGrids SOC prior…');

  // SoilGrids v2.0 — Soil Organic Carbon Stock (OCS), kg/m² per depth interval
  // Asset: projects/soilgrids-isric/soc_mean
  // Bands available: ocs_0-5cm_mean, ocs_5-15cm_mean, ocs_15-30cm_mean,
  //                  ocs_30-60cm_mean, ocs_60-100cm_mean  (units: dg/kg·cm)
  //
  // Conversion to kg/m²:
  //   OCS (kg/m²) = OCS_band_value / 10 * thickness_m
  //   (band value is in dg/kg; multiplying by thickness gives stock per area)
  //
  // IMPORTANT: SoilGrids was trained primarily on upland mineral soils.
  // It has high uncertainty in tidal sediment environments (saltmarsh, seagrass).
  // This layer is used ONLY as a weak spatial prior.
  // Uncertainty inflation (×1.2–1.5) is applied in the R Bayesian workflow:
  //   BlueCarbon_Workflow_V1.0/P3_0c_bayesian_prior_setup_bluecarbon.R

  var sg = ee.Image('projects/soilgrids-isric/soc_mean');

  var ocs_0_5    = sg.select('ocs_0-5cm_mean').divide(10).multiply(0.05);
  var ocs_5_15   = sg.select('ocs_5-15cm_mean').divide(10).multiply(0.10);
  var ocs_15_30  = sg.select('ocs_15-30cm_mean').divide(10).multiply(0.15);
  var ocs_30_60  = sg.select('ocs_30-60cm_mean').divide(10).multiply(0.30);
  var ocs_60_100 = sg.select('ocs_60-100cm_mean').divide(10).multiply(0.40);

  // Aggregate to 3-band canonical prior (matching Python notebook output columns)
  var sg_0_30   = ocs_0_5.add(ocs_5_15).add(ocs_15_30).rename('sg_soc_0_30cm');
  var sg_30_100 = ocs_30_60.add(ocs_60_100).rename('sg_soc_30_100cm');
  var sg_0_100  = sg_0_30.add(sg_30_100).rename('sg_soc_0_100cm');

  sg_soc_prior = sg_0_30.addBands(sg_30_100).addBands(sg_0_100);

  // Display
  var socVis = {min: 0, max: 20, palette: ['f7fbff','c6dbef','6baed6','2171b5','08306b']};
  Map.addLayer(sg_0_100.clip(aoi_display || aoi), socVis, 'SoilGrids SOC 0–100 cm (kg/m²)');
  Map.addLayer(sg_0_30.clip(aoi_display || aoi),  socVis, 'SoilGrids SOC 0–30 cm (kg/m²)', false);

  var stats = sg_0_100.reduceRegion({
    reducer: ee.Reducer.mean().combine(ee.Reducer.stdDev(), '', true),
    geometry: aoi, scale: 250, maxPixels: 1e9, bestEffort: true
  });
  stats.evaluate(function(s) {
    var mean = s.sg_soc_0_100cm_mean;
    var sd   = s.sg_soc_0_100cm_stdDev;
    print('SoilGrids SOC 0–100 cm AOI statistics (kg/m²):');
    print('  Mean: ' + (mean !== null ? mean.toFixed(2) : 'N/A'));
    print('  SD:   ' + (sd   !== null ? sd.toFixed(2)   : 'N/A'));
    print('⚠ Caveat: SoilGrids has high uncertainty in tidal sediment environments.');
    print('  Uncertainty inflation applied in R Bayesian workflow.');
    setStatus('Step ② ✓ — SoilGrids SOC prior loaded.\n→ Run Step ③.');
  });
});


// ─────────────────────────────────────────────────────────────────
// STEP 3 — BUILD COASTAL COVARIATE STACK (27 canonical bands)
// ─────────────────────────────────────────────────────────────────
btn3.onClick(function step3_buildCovariates() {
  if (!sg_soc_prior) { setStatus('⚠ Run Steps ①–② first.'); return; }
  setStatus('Step ③: Building 27-band covariate stack…\n(compositing may take a moment)');

  // ─────────────────────────────────────────
  // GROUP 1 — Topography & Tidal Position (7)
  // ─────────────────────────────────────────
  var dem = ee.Image('NASA/NASADEM_HGT/001').select('elevation').rename('elevation_m');
  var slope_img = ee.Terrain.slope(dem).rename('slope');

  // Elevation relative to Mean High Water (MHW)
  // Approximation: subtract 0.5 m as a proxy for MHW in BC coastal areas.
  // Replace 0.5 with your site-specific MHW datum offset if available.
  var elevRelMHW = dem.subtract(0.5).rename('elevationRelMHW');

  var twi = computeTWI(dem);
  var tpi = computeTPI(dem);

  // Tidal flat probability — Murray et al. 2019 Global Intertidal Change
  // Asset: UQ/murray/Intertidal/v1_1/global_intertidal
  // Value 1 = intertidal flat; reprojected to float probability layer.
  // If this asset is unavailable in your GEE session, tidal_flat_prob will
  // be set to 0 everywhere — update the asset path and rerun Step ③.
  var tidal_flat_prob;
  try {
    tidal_flat_prob = ee.Image('UQ/murray/Intertidal/v1_1/global_intertidal')
      .select('classification')
      .eq(1)
      .unmask(0)
      .rename('tidal_flat_prob')
      .float();
  } catch(e) {
    tidal_flat_prob = ee.Image(0).rename('tidal_flat_prob').float();
    print('⚠ Murray intertidal dataset unavailable (asset: UQ/murray/Intertidal/v1_1/global_intertidal).');
    print('  tidal_flat_prob set to 0. Update asset path in Section C and rerun Step ③.');
  }

  // Distance to coast — derived from JRC Global Surface Water (occurrence > 50%)
  // as a proxy for the shoreline. fastDistanceTransform returns pixel distance;
  // multiply by EXPORT_SCALE to convert to metres.
  var water_mask = ee.Image('JRC/GSW1_4/GlobalSurfaceWater')
    .select('occurrence').gt(50).unmask(0);
  var coastal_dist_m = water_mask
    .fastDistanceTransform(500, 'pixels', 'squared_euclidean')
    .sqrt()
    .multiply(EXPORT_SCALE)
    .rename('coastal_dist_m')
    .float();

  var topo_stack = dem
    .addBands(slope_img)
    .addBands(elevRelMHW)
    .addBands(twi)
    .addBands(tpi)
    .addBands(tidal_flat_prob)
    .addBands(coastal_dist_m);

  // ─────────────────────────────
  // GROUP 2 — Sentinel-1 SAR (3)
  // ─────────────────────────────
  // IW mode, ascending orbit, VV + VH polarisation.
  // Noise filter: VV > -30 dB (removes obvious noise artifacts).
  var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
    .filter(ee.Filter.date(SAR_START, SAR_END))
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .filter(ee.Filter.gt('VV', -30));

  var vv    = s1.select('VV').median().rename('VV_mean');
  var vh    = s1.select('VH').median().rename('VH_mean');
  var vvvh  = vv.subtract(vh).rename('VVVH_ratio');
  var sar_stack = vv.addBands(vh).addBands(vvvh);

  // ──────────────────────────────────────
  // GROUP 3 — Sentinel-2 Optical (12)
  // ──────────────────────────────────────
  // Summer months (May–September): peak tidal marsh biomass in BC.
  // Cloud mask + median composite. Indices: NDVI, LSWI, mNDWI, EVI, SAVI,
  // and Tasseled Cap Wetness (most inundation-sensitive TC component).
  var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
    .filter(ee.Filter.date(S2_START, S2_END))
    .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', S2_CLOUD_THRESHOLD))
    .filter(ee.Filter.calendarRange(5, 9, 'month'))  // May–September
    .map(maskS2clouds)
    .map(addS2IndicesAndRename)
    .median();

  var s2_raw = s2.select(['B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2']);
  var s2_idx = s2.select([
    'NDVI_median', 'LSWI_median', 'mNDWI_median',
    'EVI_median',  'SAVI_median', 'tidal_wetness'
  ]);
  var s2_stack = s2_raw.addBands(s2_idx);

  // ────────────────────────────
  // GROUP 4 — TerraClimate (2)
  // ────────────────────────────
  var tc   = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE')
    .filter(ee.Filter.date(TC_START, TC_END));
  var mat  = tc.select('tmmx').mean().subtract(273.15).rename('MAT_C');
  var map_mm = tc.select('pr').mean().multiply(12).rename('MAP_mm');
  var climate_stack = mat.addBands(map_mm);

  // ────────────────────────────────────────────────────────
  // GROUP 5 — SoilGrids SOC Prior (3 — already in Step 2)
  // ────────────────────────────────────────────────────────
  // sg_soc_prior already contains: sg_soc_0_30cm, sg_soc_30_100cm, sg_soc_0_100cm

  // ────────────────────────────────────────────────────
  // ASSEMBLE canonical 27-band stack
  // ────────────────────────────────────────────────────
  cov_stack = topo_stack
    .addBands(sar_stack)
    .addBands(s2_stack)
    .addBands(climate_stack)
    .addBands(sg_soc_prior);

  // GROUP 6 — Google Satellite Embedding V1 (64 bands, 10 m — separate export)
  var embed_col = ee.ImageCollection('GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL');
  if (EMBEDDING_YEAR !== null) {
    embed_col = embed_col.filter(ee.Filter.calendarRange(EMBEDDING_YEAR, EMBEDDING_YEAR, 'year'));
  }
  embed_img = embed_col.median();

  // ── Band inventory ────────────────────────────────────────────
  print('═══════════════════════════════════════════════════════');
  print('Canonical 27-band covariate stack:');
  print('  Group 1 — Topography & Tidal (7):');
  print('    elevation_m, slope, elevationRelMHW, twi, tpi, tidal_flat_prob, coastal_dist_m');
  print('  Group 2 — SAR (3):');
  print('    VV_mean, VH_mean, VVVH_ratio');
  print('  Group 3 — Sentinel-2 (12):');
  print('    B, G, R, NIR, SWIR1, SWIR2');
  print('    NDVI_median, LSWI_median, mNDWI_median, EVI_median, SAVI_median, tidal_wetness');
  print('  Group 4 — Climate (2):');
  print('    MAT_C, MAP_mm');
  print('  Group 5 — SoilGrids Prior (3):');
  print('    sg_soc_0_30cm, sg_soc_30_100cm, sg_soc_0_100cm');
  print('  Embedding (separate 10 m export):');
  print('    64 bands from GOOGLE/SATELLITE_EMBEDDING/V1/ANNUAL');
  print('═══════════════════════════════════════════════════════');

  // ── Map visualisation layers ──────────────────────────────────
  var displayRegion = aoi_display || aoi;
  Map.addLayer(cov_stack.select('NDVI_median').clip(displayRegion),
    {min: -0.2, max: 0.8, palette: ['d73027','fee090','91cf60','1a9641']}, 'NDVI');
  Map.addLayer(cov_stack.select('tidal_wetness').clip(displayRegion),
    {min: -0.35, max: 0.05, palette: ['ffffcc','41b6c4','0c2c84']}, 'Tidal Wetness', false);
  Map.addLayer(cov_stack.select('mNDWI_median').clip(displayRegion),
    {min: -0.5, max: 0.5, palette: ['a50026','ffffbf','4575b4']}, 'mNDWI', false);
  Map.addLayer(cov_stack.select('LSWI_median').clip(displayRegion),
    {min: -0.3, max: 0.5, palette: ['d73027','ffffbf','4575b4']}, 'LSWI', false);
  Map.addLayer(cov_stack.select('VV_mean').clip(displayRegion),
    {min: -25, max: -5, palette: ['000000','cccccc','ffffff']}, 'SAR VV', false);
  Map.addLayer(cov_stack.select('elevation_m').clip(displayRegion),
    {min: -5, max: 10, palette: ['0571b0','92c5de','f7f7f7','d6604d','ca0020']}, 'Elevation (m)', false);
  Map.addLayer(cov_stack.select('twi').clip(displayRegion),
    {min: 0, max: 12, palette: ['f7fbff','c6dbef','6baed6','2171b5','08306b']}, 'TWI', false);
  Map.addLayer(cov_stack.select('tpi').clip(displayRegion),
    {min: -3, max: 3, palette: ['d73027','ffffbf','1a9641']}, 'TPI', false);
  Map.addLayer(cov_stack.select('coastal_dist_m').clip(displayRegion),
    {min: 0, max: 5000, palette: ['08306b','6baed6','f7fbff']}, 'Coastal Distance (m)', false);
  Map.addLayer(cov_stack.select('tidal_flat_prob').clip(displayRegion),
    {min: 0, max: 1, palette: ['ffffff','fee090','d73027']}, 'Tidal Flat Probability', false);
  Map.addLayer(cov_stack.select('MAT_C').clip(displayRegion),
    {min: 0, max: 12, palette: ['4575b4','ffffbf','d73027']}, 'MAT (°C)', false);
  Map.addLayer(cov_stack.select('sg_soc_0_100cm').clip(displayRegion),
    {min: 0, max: 20, palette: ['f7fbff','6baed6','08306b']}, 'SoilGrids SOC 0–100 cm', false);

  // ── Insert ④ Export button after Step ③ ──────────────────────
  if (!snapshot_btn) {
    snapshot_btn = ui.Button('④ Export Covariate Snapshot to Drive');
    snapshot_btn.onClick(step4_exportSnapshot);
    var children = panel.widgets();
    var btn5_idx = -1;
    for (var i = 0; i < children.length(); i++) {
      if (children.get(i) === btn5) { btn5_idx = i; break; }
    }
    if (btn5_idx >= 0) { panel.insert(btn5_idx, snapshot_btn); }
    else               { panel.add(snapshot_btn); }
  }

  setStatus('Step ③ ✓ — 27-band stack built. Embedding ready.\n→ Run Step ④ to export to Drive.');
});


// ─────────────────────────────────────────────────────────────────
// STEP 4 — EXPORT COVARIATE SNAPSHOT (inserted dynamically)
// ─────────────────────────────────────────────────────────────────
function step4_exportSnapshot() {
  if (!cov_stack) { setStatus('⚠ Run Step ③ first.'); return; }
  setStatus('Step ④: Queuing export tasks to Drive…');

  var exportRegion = aoi_display || aoi;

  // Task 1 — 27-band covariate stack at 25 m, Canada Albers Equal Area
  Export.image.toDrive({
    image       : cov_stack.clip(aoi),
    description : 'BlueCarbon_Covariate_Snapshot_25m_' + PROJECT_YEAR,
    folder      : EXPORT_FOLDER,
    fileNamePrefix: 'BlueCarbon_Covariate_Snapshot_25m_' + PROJECT_YEAR,
    region      : exportRegion,
    scale       : EXPORT_SCALE,
    crs         : EXPORT_CRS,
    maxPixels   : 1e13,
    fileFormat  : 'GeoTIFF'
  });

  // Task 2 — Google Satellite Embedding V1 (64 bands) at 10 m, Canada Albers Equal Area
  Export.image.toDrive({
    image       : embed_img.clip(aoi),
    description : 'BlueCarbon_GoogleEmbedding_V1_10m_' + PROJECT_YEAR,
    folder      : EXPORT_FOLDER,
    fileNamePrefix: 'BlueCarbon_GoogleEmbedding_V1_10m_' + PROJECT_YEAR,
    region      : exportRegion,
    scale       : EMBED_SCALE,
    crs         : EXPORT_CRS,
    maxPixels   : 1e13,
    fileFormat  : 'GeoTIFF'
  });

  print('Export tasks queued (check Tasks panel → click Run):');
  print('  1. BlueCarbon_Covariate_Snapshot_25m_' + PROJECT_YEAR);
  print('     27 bands | ' + EXPORT_SCALE + ' m | ' + EXPORT_CRS + ' | Drive: ' + EXPORT_FOLDER);
  print('  2. BlueCarbon_GoogleEmbedding_V1_10m_' + PROJECT_YEAR);
  print('     64 bands | ' + EMBED_SCALE + ' m | ' + EXPORT_CRS + ' | Drive: ' + EXPORT_FOLDER);
  print('  Once downloaded, copy GeoTIFFs to:');
  print('  BlueCarbon_Workflow_V1.0/Pre-Analysis Data Preparation/covariates/');

  setStatus('Step ④ ✓ — 2 export tasks queued.\nCheck Tasks panel → click Run on each.\n→ Run Step ⑤ for sampling points.');
}


// ─────────────────────────────────────────────────────────────────
// STEP 5 — GENERATE STRATIFIED SOIL SAMPLING POINTS
// ─────────────────────────────────────────────────────────────────
btn5.onClick(function step5_generateSampling() {
  if (!cov_stack) { setStatus('⚠ Run Steps ①–③ first.'); return; }
  setStatus('Step ⑤: Generating ' + N_SOIL_SAMPLES + ' sampling points…');

  var sg_0_100   = sg_soc_prior.select('sg_soc_0_100cm');
  var elev       = cov_stack.select('elevation_m');
  var tidal_prob = cov_stack.select('tidal_flat_prob');
  var mndwi      = cov_stack.select('mNDWI_median');

  // Define the potential tidal zone:
  //   - Elevation between -2 m and +3 m above MHW proxy, AND
  //   - Not persistently inundated (mNDWI < 0.2 = vegetated or exposed)
  //   - OR high tidal flat probability from Murray et al.
  var tidal_zone = elev.lte(3.0).and(elev.gte(-2.0))
                       .and(mndwi.lt(0.2))
                       .or(tidal_prob.gte(0.3));

  // SOC uncertainty strata — percentile-based allocation
  var soc_pct = sg_0_100.reduceRegion({
    reducer: ee.Reducer.percentile([33, 67]),
    geometry: aoi, scale: 250, maxPixels: 1e9, bestEffort: true
  });

  soc_pct.evaluate(function(pct) {
    var p33 = pct['sg_soc_0_100cm_p33'] || 5;
    var p67 = pct['sg_soc_0_100cm_p67'] || 15;

    // 3 strata within tidal zone (only tidal zone sampled for blue carbon)
    var strata_img = ee.Image(0)
      .where(tidal_zone.and(sg_0_100.lt(p33)),              1)  // low SOC uncertainty
      .where(tidal_zone.and(sg_0_100.gte(p33).and(sg_0_100.lt(p67))), 2)  // med
      .where(tidal_zone.and(sg_0_100.gte(p67)),             3)  // high SOC uncertainty
      .rename('stratum');

    var n_per_stratum = Math.ceil(N_SOIL_SAMPLES / 3);

    soil_pts = strata_img.stratifiedSample({
      numPoints: n_per_stratum,
      classBand: 'stratum',
      region: aoi,
      scale: EXPORT_SCALE,
      seed: 42,
      geometries: true
    });

    // Add human-readable stratum labels
    var labels = ee.Dictionary({
      '1': 'Tidal_LowUncertainty',
      '2': 'Tidal_MedUncertainty',
      '3': 'Tidal_HighUncertainty'
    });
    soil_pts = soil_pts.map(function(f) {
      return f.set('stratum_label',
        labels.get(ee.String(f.get('stratum').toInt())));
    });

    Map.addLayer(soil_pts, {color: 'ff0000'}, 'Sampling Points (' + N_SOIL_SAMPLES + ')');

    print('Stratified sampling points generated:');
    print('  Total: ~' + (n_per_stratum * 3) + ' points across 3 tidal strata');
    print('  Stratum 1: Tidal, Low SOC uncertainty   (SOC < p33 = ' + p33.toFixed(1) + ' kg/m²)');
    print('  Stratum 2: Tidal, Med SOC uncertainty   (p33 – p67)');
    print('  Stratum 3: Tidal, High SOC uncertainty  (SOC > p67 = ' + p67.toFixed(1) + ' kg/m²)');

    setStatus('Step ⑤ ✓ — Sampling points generated.\n→ Run Step ⑥ to export reports.');
  });
});


// ─────────────────────────────────────────────────────────────────
// STEP 6 — REPORTS AND EXPORTS
// ─────────────────────────────────────────────────────────────────
btn6.onClick(function step6_reportsAndExports() {
  if (!cov_stack) { setStatus('⚠ Run Steps ①–③ first.'); return; }
  setStatus('Step ⑥: Computing statistics and queuing exports…');

  // ── AOI-wide covariate statistics ─────────────────────────────
  var stats = cov_stack.reduceRegion({
    reducer: ee.Reducer.mean()
      .combine(ee.Reducer.stdDev(), '', true)
      .combine(ee.Reducer.minMax(), '', true),
    geometry: aoi,
    scale: EXPORT_SCALE * 4,  // coarser for speed
    maxPixels: 1e9,
    bestEffort: true
  });

  stats.evaluate(function(s) {
    var bands = [
      'elevation_m', 'slope', 'elevationRelMHW', 'twi', 'tpi',
      'tidal_flat_prob', 'coastal_dist_m',
      'VV_mean', 'VH_mean', 'VVVH_ratio',
      'NDVI_median', 'LSWI_median', 'mNDWI_median', 'EVI_median', 'SAVI_median', 'tidal_wetness',
      'MAT_C', 'MAP_mm',
      'sg_soc_0_30cm', 'sg_soc_30_100cm', 'sg_soc_0_100cm'
    ];
    print('═══════════════════════════════════════════════════════');
    print('AOI Covariate Summary:');
    bands.forEach(function(b) {
      var mean = s[b + '_mean'];
      var sd   = s[b + '_stdDev'];
      if (mean !== undefined && mean !== null) {
        print('  ' + b + ': mean = ' + mean.toFixed(3) + '  sd = ' + sd.toFixed(3));
      }
    });
    print('═══════════════════════════════════════════════════════');
  });

  // ── Export band manifest CSV ──────────────────────────────────
  // Lists all 27 band names in canonical order — use to verify
  // that Python notebook output CSVs have matching column names.
  var band_names  = [
    'elevation_m', 'slope', 'elevationRelMHW', 'twi', 'tpi',
    'tidal_flat_prob', 'coastal_dist_m',
    'VV_mean', 'VH_mean', 'VVVH_ratio',
    'B', 'G', 'R', 'NIR', 'SWIR1', 'SWIR2',
    'NDVI_median', 'LSWI_median', 'mNDWI_median', 'EVI_median', 'SAVI_median', 'tidal_wetness',
    'MAT_C', 'MAP_mm',
    'sg_soc_0_30cm', 'sg_soc_30_100cm', 'sg_soc_0_100cm'
  ];
  var band_groups = [
    'Topography', 'Topography', 'Topography', 'Topography', 'Topography',
    'Tidal', 'Tidal',
    'SAR', 'SAR', 'SAR',
    'Sentinel-2_Raw', 'Sentinel-2_Raw', 'Sentinel-2_Raw', 'Sentinel-2_Raw', 'Sentinel-2_Raw', 'Sentinel-2_Raw',
    'Sentinel-2_Index', 'Sentinel-2_Index', 'Sentinel-2_Index', 'Sentinel-2_Index', 'Sentinel-2_Index', 'Sentinel-2_Index',
    'Climate', 'Climate',
    'SoilGrids_Prior', 'SoilGrids_Prior', 'SoilGrids_Prior'
  ];

  var manifest_fc = ee.FeatureCollection(
    band_names.map(function(b, i) {
      return ee.Feature(null, {
        canonical_order: i + 1,
        band_name: b,
        group: band_groups[i]
      });
    })
  );

  Export.table.toDrive({
    collection   : manifest_fc,
    description  : 'BlueCarbon_Covariate_BandManifest_' + PROJECT_YEAR,
    folder       : EXPORT_FOLDER,
    fileNamePrefix: 'BlueCarbon_Covariate_BandManifest_' + PROJECT_YEAR,
    fileFormat   : 'CSV'
  });

  // ── Export sampling points (CSV + KML) ────────────────────────
  if (soil_pts) {
    Export.table.toDrive({
      collection   : soil_pts,
      description  : 'BlueCarbon_Sampling_Points_CSV_' + PROJECT_YEAR,
      folder       : EXPORT_FOLDER,
      fileNamePrefix: 'BlueCarbon_Sampling_Points_' + PROJECT_YEAR,
      fileFormat   : 'CSV'
    });
    Export.table.toDrive({
      collection   : soil_pts,
      description  : 'BlueCarbon_Sampling_Points_KML_' + PROJECT_YEAR,
      folder       : EXPORT_FOLDER,
      fileNamePrefix: 'BlueCarbon_Sampling_Points_' + PROJECT_YEAR,
      fileFormat   : 'KML'
    });
    print('Sampling point exports queued: CSV + KML');
  } else {
    print('⚠ Sampling points not generated — run Step ⑤ before Step ⑥ to include them.');
  }

  print('Band manifest export queued.');
  print('→ Go to Tasks panel and click Run on all queued tasks.');
  print('');
  print('Export summary:');
  print('  BlueCarbon_Covariate_Snapshot_25m_' + PROJECT_YEAR + '.tif  (27 bands, 25 m, EPSG:3347)');
  print('  BlueCarbon_GoogleEmbedding_V1_10m_' + PROJECT_YEAR + '.tif  (64 bands, 10 m, EPSG:3347)');
  print('  BlueCarbon_Sampling_Points_' + PROJECT_YEAR + '.csv/.kml    (' + N_SOIL_SAMPLES + ' points)');
  print('  BlueCarbon_Covariate_BandManifest_' + PROJECT_YEAR + '.csv  (27-band canonical list)');
  print('');
  print('After downloading, move GeoTIFFs to:');
  print('  BlueCarbon_Workflow_V1.0/Pre-Analysis Data Preparation/covariates/');

  setStatus('Step ⑥ ✓ — All exports queued.\nCheck Tasks panel → click Run.\n\nDone! Files go to Drive:\n' + EXPORT_FOLDER);
});


// ─────────────────────────────────────────────────────────────────
// INITIAL SETUP
// ─────────────────────────────────────────────────────────────────
Map.setCenter(-125.5, 50.5, 7);   // Default: central BC coast
Map.setOptions('SATELLITE');

print('Coastal Blue Carbon — AOI Covariate Extraction v1.0');
print('Canada Coastal Ecosystems: Tidal Marsh + Seagrass | VM0033');
print('─────────────────────────────────────────────────────────');
print('Set AOI_ASSET in Section A (line 1) to your boundary asset,');
print('or leave as null to draw interactively on the map.');
print('Then run steps ①–⑥ using the side panel buttons.');
