##### Creating NDVI to use as mask for Functional Trait Plots #####
library(pacman)
p_load(tidyverse, terra, sf, parallel)

data_dir <- "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/"


##### Reading in Data #####
### Airborne hyperspectral data collected in 2025
TGPPAirborne_2025_rast <- rast(paste0(data_dir, "Tallgrass2025_Airborne/AirborneClean/TIFF/1768_TGGP_Subset_REF_Mosaic_001-015_CleanBands.tif"))

### band wavelengths of airborne data
band_wavelengths  <- c(433.350006,  440.130005,  446.920013,  453.709991,  460.500000,  467.299988, 474.109985,  480.920013,  487.730011,  494.559998,  501.380005,  508.209991, 515.039978,  521.880005,  528.719971,  535.570007,  542.419983,  549.270020, 556.130005,  562.989990,  569.849976,  576.719971,  583.590027,  590.460022, 597.340027,  604.210022,  611.090027,  617.979980,  624.859985,  631.750000, 638.640015,  645.530029,  652.419983,  659.320007,  666.210022,  673.109985, 680.010010,  686.909973,  693.809998,  700.710022,  707.609985,  714.510010, 721.409973,  728.320007,  735.219971,  742.119995,  749.020020,  755.929993, 762.830017,  769.729980,  776.630005,  783.530029,  790.429993,  797.330017, 804.229980,  811.119995,  818.010010,  824.909973,  831.799988,  838.690002, 845.570007,  852.460022,  859.340027,  866.219971,  873.099976,  879.969971, 886.849976,  893.719971,  900.580017,  907.440002,  914.299988,  992.390015, 998.710022, 1005.039978, 1011.369995, 1017.700012, 1024.030029, 1030.349976, 1036.680054, 1043.010010, 1049.329956, 1055.660034, 1061.979980, 1068.310059, 1074.630005, 1080.959961, 1087.280029, 1093.599976, 1099.920044, 1106.250000, 1163.119995, 1169.439941, 1175.750000, 1182.069946, 1188.390015, 1194.699951, 1201.010010, 1207.329956, 1213.640015, 1219.949951, 1226.270020, 1232.579956, 1238.890015, 1245.199951, 1251.510010, 1257.819946, 1264.130005, 1270.430054, 1276.739990, 1283.050049, 1289.349976, 1295.660034, 1490.780029, 1497.060059, 1503.349976, 1509.630005, 1515.910034, 1522.189941, 1528.459961, 1534.739990, 1541.020020, 1547.290039, 1553.569946, 1559.839966, 1566.119995, 1572.390015, 1578.660034, 1584.930054, 1591.199951, 1597.469971, 1603.739990, 1610.000000, 1616.270020, 1622.530029, 1628.800049, 1635.060059, 1641.319946, 1647.579956, 1653.839966, 1660.099976, 1666.359985, 1672.619995, 1678.880005, 1685.130005, 1691.390015, 1697.640015, 1703.890015, 1710.140015, 1716.390015, 1722.640015, 1728.890015, 1735.140015, 1741.390015, 1747.630005, 1753.880005, 1760.119995, 1766.359985, 1772.599976, 2008.949951, 2015.150024, 2021.339966, 2027.540039, 2033.729980, 2039.920044, 2046.109985, 2052.300049, 2058.479980, 2064.669922, 2070.860107, 2077.040039, 2083.219971, 2089.399902, 2095.580078, 2101.760010, 2107.939941, 2114.110107, 2120.290039, 2126.459961, 2132.629883, 2138.800049, 2144.969971, 2151.139893, 2157.300049, 2163.469971, 2169.629883, 2175.790039, 2181.949951, 2188.110107, 2194.270020, 2200.419922, 2206.580078, 2212.729980, 2218.879883, 2225.030029, 2231.179932, 2237.330078, 2243.479980, 2249.620117, 2255.760010, 2261.899902, 2268.040039, 2274.179932, 2280.320068, 2286.449951, 2292.590088, 2298.719971, 2304.850098, 2310.979980, 2317.110107, 2323.239990, 2329.360107, 2335.479980, 2341.610107, 2347.729980, 2353.840088)

### checking which band is closest to 680 (~red wavelength) for NDVI calculation and what the wavelength is
which.min(abs(band_wavelengths - 680))
red_center <- band_wavelengths[which.min(abs(band_wavelengths - 680))]
red_center

### filtering out wavelengths to be used for red band average for NDVI calculation
red_NDVIbands <- band_wavelengths[band_wavelengths > red_center - 15 & band_wavelengths < red_center + 15]
red_NDVIbands

### checking which bands in the "clean" raster correspond with my chosen red bands
redband_index <- which(band_wavelengths %in% red_NDVIbands)
redband_index

### creating raster of just red bands
RedBands_rast <- TGPPAirborne_2025_rast[[redband_index]]

### repeating above process for NIR bands
which.min(abs(band_wavelengths - 800))
NIR_center <- band_wavelengths[which.min(abs(band_wavelengths - 800))]
NIR_center

### filtering out wavelengths to be used for NIR band average for NDVI calculation
NIR_NDVIbands <- band_wavelengths[band_wavelengths > NIR_center - 15 & band_wavelengths < NIR_center + 15]
NIR_NDVIbands

### checking which bands in the "clean" raster correspond with my chosen NIR bands
NIRband_index <- which(band_wavelengths %in% NIR_NDVIbands)
NIRband_index

### creating raster of just NIR bands
NIRBands_rast <- TGPPAirborne_2025_rast[[NIRband_index]]


### creating rasters of red and NIR reflectance averages 
RedAvg_rast <- app(RedBands_rast, fun = "mean")
NIRAvg_rast <- app(NIRBands_rast, fun = "mean")



### calculating NDVI
NDVI_rast <- (NIRAvg_rast - RedAvg_rast) / (NIRAvg_rast + RedAvg_rast)
### plotting NDVI
plot(NDVI_rast)

### saving NDVI raster
writeRaster(NDVI_rast, paste0(data_dir, "SpectralIndices/TGPP/2025/TGPP2025NDVI.tif"))
