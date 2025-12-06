##### Masking Functional Trait Layers with NDVI and Woodland area masks #####

##### Loading packages, reading in data #####
library(pacman)
p_load(tidyverse, sf, terra, tidyterra, ggspatial, parallel, mcprogress)


### local folder with all data in subdirectories
data_dir <- "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/"

### Airborne hyperspectral data collected in 2025
TGPPAirborne_2025_rast <- rast(paste0(data_dir, "Tallgrass2025_Airborne/AirborneClean/TIFF/1768_TGGP_Subset_REF_Mosaic_001-015_CleanBands.tif")) #/ 10000
# writeRaster(TGPPAirborne_2025_rast, paste0(data_dir, "Tallgrass2025_Airborne/AirborneClean/TIFF/1768_TGGP_Subset_REF_Mosaic_001-015_CleanBands.tif"), overwrite = TRUE)

### the band names are long and clunky and I don't like it so I'm just calling them "bandx" based on the order they're stacked in
band_ids <- c()
for(x in 1:nlyr(TGPPAirborne_2025_rast)){
  band_ids = c(band_ids, paste0("band", x))
}
names(TGPPAirborne_2025_rast) <- band_ids

### AOI polygon
TGPP_AOI <- read_sf(paste0(data_dir, "GIS_Files/TGPP/Northern_TGPP_2025.shp"))

### Now that extraction is done, reading in gpkg object of quadrat points with all functional traits and airborne band reflectances
TGPP_QuadPts_FunBmass_sf <- read_sf(paste0(data_dir, "GIS_Files/TGPP/ModelingData/TGPP2025BiomassFunTrait_Airborne.gpkg"))

### reading in oak woodland shapefile for masking
TGPPOakWoodland_sf <- read_sf(paste0(data_dir, "GIS_Files/TGPP/TreeMask2025/TGP_TreeMaskUTM_2025.shp"))

### reading in functional rasters
TGPPTNpct2025_rast <- rast(paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Nitrogen/MeanRaster/TGPP2025TNpct.tif"))
TGPPPpct2025_rast <- rast(paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Phosphorus/MeanRaster/TGPP2025Ppct.tif"))
TGPPKpct2025_rast <- rast(paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Potassium/MeanRaster/TGPP2025Kpct.tif"))

### reading in NDVI raster
TGPPNDVI_rast <- rast(paste0(data_dir, "SpectralIndices/TGPP/2025/TGPP2025NDVI.tif"))

plot(TGPPTNpct2025_rast)
plot(st_geometry(TGPPOakWoodland_sf), add = TRUE)

##### Masking low NDVI areas and woodland area out of functional trait rasters #####
### masking out woodland areas first
TGPPTNpct2025Mask_rast <- mask(TGPPTNpct2025_rast, TGPPOakWoodland_sf, inverse = TRUE)
TGPPPpct2025Mask_rast <- mask(TGPPPpct2025_rast, TGPPOakWoodland_sf, inverse = TRUE)
TGPPKpct2025Mask_rast <- mask(TGPPKpct2025_rast, TGPPOakWoodland_sf, inverse = TRUE)
# plot(TGPPTNpct2025Mask_rast)

### masking by NDVI values, setting all pixels with NDVI values < 0.4 to NA
TGPPTNpct2025Mask_rast <- lapp(c(TGPPTNpct2025Mask_rast, TGPPNDVI_rast), fun = function(x, y){x[y < 0.4] <- NA; return(x)})
TGPPPpct2025Mask_rast <- lapp(c(TGPPPpct2025Mask_rast, TGPPNDVI_rast), fun = function(x, y){x[y < 0.4] <- NA; return(x)})
TGPPKpct2025Mask_rast <- lapp(c(TGPPKpct2025Mask_rast, TGPPNDVI_rast), fun = function(x, y){x[y < 0.4] <- NA; return(x)})

# plot(TGPPTNpct2025Mask_rast)
# plot(TGPPPpct2025Mask_rast)
plot(TGPPKpct2025Mask_rast)
plot(st_geometry(TGPP_AOI), add = TRUE)


##### Saving Masked Rasters #####
writeRaster(TGPPTNpct2025Mask_rast, paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Nitrogen/MeanRaster/TGPP2025TNpct_Mask.tif"))
writeRaster(TGPPPpct2025Mask_rast, paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Phosphorus/MeanRaster/TGPP2025Ppct_Mask.tif"))
writeRaster(TGPPKpct2025Mask_rast, paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Potassium/MeanRaster/TGPP2025Kpct_Mask.tif"))