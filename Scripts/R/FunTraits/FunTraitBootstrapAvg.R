##### Taking Average of 50 bootstrapped TN_pct Prediction Rasters #####

##### Loading packages, reading in data #####
library(pacman)
p_load(tidyverse, sf, terra, tidyterra, parallel, ggspatial)
options(mc.cores = ceiling(detectCores() / 3))

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

##### Taking Model Average for Nitrogen #####
TNpct_rast_filelist <- list.files(path = "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/FunctionalTraits/TGPP/Outputs/Rasters/Nitrogen", pattern = "*.tif", full.names = TRUE)

### reading in all files stored in list above
TNpct_rast_list <- lapply(TNpct_rast_filelist, FUN = rast)

### storing first model prediction raster in the object that's going to store all 50 after the for loop below
TNpct_rasts <- TNpct_rast_list[[1]]

### stacking different 
for (i in 2:length(TNpct_rast_list)){
  TNpct_rasts <- c(TNpct_rasts, TNpct_rast_list[[i]])
}

### taking mean of all 50 models
TGPPTNpctMean_rast <- app(TNpct_rasts, fun = "mean")

TGPPTNpct2025_rast <- app(TGPPTNpctMean_rast, fun = function(x){x[x < 0] <- 0; return(x)})
TGPPTNpct2025_rast <- lapp(c(TGPPTNpct2025_rast, TGPPAirborne_2025_rast[[1]]), fun = function(x, y){x[y == 0] <- NA; return(x)})
min(values(TGPPTNpct2025_rast), na.rm = TRUE)
max(values(TGPPTNpct2025_rast), na.rm = TRUE)
# plot(TGPPTNpct2025_rast)
# TGPPTNpct2025Crop_rast <- crop(TGPPTNpct2025_rast, TGPP_AOI, mask = TRUE)
# min(values(TGPPTNpct2025Crop_rast), na.rm = TRUE)
# max(values(TGPPTNpct2025Crop_rast), na.rm = TRUE)
# plot(TGPPTNpct2025Crop_rast)

### creating ggplot
TGPPTNpct2025_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = TGPPTNpct2025_rast) +
  scale_fill_whitebox_c(palette = "viridi", name = "TN (%)") +#, limits = c(min(values(TGPPTNpct2025_rast), na.rm = TRUE), max(values(TGPPTNpct2025_rast), na.rm = TRUE))) +
  geom_sf(data = TGPP_QuadPts_FunBmass_sf, fill = NA, color = 'black') +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", pad_y = unit(1, units = "cm")) +#, pad_x = unit(1, units = "cm")) +
  ggtitle("Foliar Nitrogen at TGPP in 2025") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
## displaying plot
TGPPTNpct2025_plot
### saving plot
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/TGPPTNpct2025.png")

### Saving nitrogen Raster
writeRaster(TGPPTNpct2025_rast, paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Nitrogen/MeanRaster/TGPP2025TNpct.tif"))

##### Taking Model Average for Phosphorus #####
Ppct_rast_filelist <- list.files(path = "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/FunctionalTraits/TGPP/Outputs/Rasters/Phosphorus", pattern = "*.tif", full.names = TRUE)

### reading in all files stored in list above
Ppct_rast_list <- lapply(Ppct_rast_filelist, FUN = rast)

### storing first model prediction raster in the object that's going to store all 50 after the for loop below
Ppct_rasts <- Ppct_rast_list[[1]]

### stacking different 
for (i in 2:length(Ppct_rast_list)){
  Ppct_rasts <- c(Ppct_rasts, Ppct_rast_list[[i]])
}

### taking mean of all 50 models
TGPPPpctMean_rast <- app(Ppct_rasts, fun = "mean")

### removing negative values and pixels outside of AOI
TGPPPpct2025_rast <- app(TGPPPpctMean_rast, fun = function(x){x[x < 0] <- 0; return(x)})
TGPPPpct2025_rast <- lapp(c(TGPPPpct2025_rast, TGPPAirborne_2025_rast[[1]]), fun = function(x, y){x[y == 0] <- NA; return(x)})
min(values(TGPPPpct2025_rast), na.rm = TRUE)
max(values(TGPPPpct2025_rast), na.rm = TRUE)
# plot(TGPPTNpct2025_rast)
# TGPPTNpct2025Crop_rast <- crop(TGPPTNpct2025_rast, TGPP_AOI, mask = TRUE)
# min(values(TGPPTNpct2025Crop_rast), na.rm = TRUE)
# max(values(TGPPTNpct2025Crop_rast), na.rm = TRUE)
# plot(TGPPTNpct2025Crop_rast)

### creating ggplot
TGPPPpct2025_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = TGPPPpct2025_rast) +
  scale_fill_whitebox_c(palette = "viridi", name = "P (%)") +#, limits = c(min(values(TGPPPpct2025_rast), na.rm = TRUE), max(values(TGPPPpct2025_rast), na.rm = TRUE))) +
  geom_sf(data = TGPP_QuadPts_FunBmass_sf, fill = NA, color = 'black') +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", pad_y = unit(1, units = "cm")) +#, pad_x = unit(1, units = "cm")) +
  ggtitle("Foliar Phosphorus at TGPP in 2025") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
## displaying plot
TGPPPpct2025_plot
### saving plot
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/TGPPPpct2025.png")

### Saving Phosphorus Raster
writeRaster(TGPPPpct2025_rast, paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Phosphorus/MeanRaster/TGPP2025Ppct.tif"))

##### Taking Model Average for Potassium #####
Kpct_rast_filelist <- list.files(path = "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/FunctionalTraits/TGPP/Outputs/Rasters/Potassium", pattern = "*.tif", full.names = TRUE)

### reading in all files stored in list above
Kpct_rast_list <- lapply(Kpct_rast_filelist, FUN = rast)

### storing first model prediction raster in the object that's going to store all 50 after the for loop below
Kpct_rasts <- Kpct_rast_list[[1]]

### stacking different 
for (i in 2:length(Kpct_rast_list)){
  Kpct_rasts <- c(Kpct_rasts, Kpct_rast_list[[i]])
}

### taking mean of all 50 models
TGPPKpctMean_rast <- app(Kpct_rasts, fun = "mean")

TGPPKpct2025_rast <- app(TGPPKpctMean_rast, fun = function(x){x[x < 0] <- 0; return(x)})
TGPPKpct2025_rast <- lapp(c(TGPPKpct2025_rast, TGPPAirborne_2025_rast[[1]]), fun = function(x, y){x[y == 0] <- NA; return(x)})
min(values(TGPPKpct2025_rast), na.rm = TRUE)
max(values(TGPPKpct2025_rast), na.rm = TRUE)
# plot(TGPPTNpct2025_rast)
# TGPPTNpct2025Crop_rast <- crop(TGPPTNpct2025_rast, TGPP_AOI, mask = TRUE)
# min(values(TGPPTNpct2025Crop_rast), na.rm = TRUE)
# max(values(TGPPTNpct2025Crop_rast), na.rm = TRUE)
# plot(TGPPTNpct2025Crop_rast)

### creating ggplot
TGPPKpct2025_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = TGPPKpct2025_rast) +
  scale_fill_whitebox_c(palette = "viridi", name = "K (%)") + #, limits = c(min(values(TGPPKpct2025_rast), na.rm = TRUE), max(values(TGPPKpct2025_rast), na.rm = TRUE))) +
  geom_sf(data = TGPP_QuadPts_FunBmass_sf, fill = NA, color = 'black') +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", pad_y = unit(1, units = "cm")) +#, pad_x = unit(1, units = "cm")) +
  ggtitle("Foliar Potassium at TGPP in 2025") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
## displaying plot
TGPPKpct2025_plot
### saving plot
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/TGPPKpct2025.png")


### Saving Potassium Raster
writeRaster(TGPPKpct2025_rast, paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Potassium/MeanRaster/TGPP2025Kpct.tif"))



##### Taking Model Average for Protein #####
Proteinpct_rast_filelist <- list.files(path = "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/FunctionalTraits/TGPP/Outputs/Rasters/Protein", pattern = "*.tif", full.names = TRUE)

### reading in all files stored in list above
Proteinpct_rast_list <- lapply(Proteinpct_rast_filelist, FUN = rast)

### storing first model prediction raster in the object that's going to store all 50 after the for loop below
Proteinpct_rasts <- Proteinpct_rast_list[[1]]

### stacking different 
for (i in 2:length(Proteinpct_rast_list)){
  Proteinpct_rasts <- c(Proteinpct_rasts, Proteinpct_rast_list[[i]])
}

### taking mean of all 50 models
TGPPProteinpctMean_rast <- app(Proteinpct_rasts, fun = "mean")

TGPPProteinpct2025_rast <- app(TGPPProteinpctMean_rast, fun = function(x){x[x < 0] <- 0; return(x)})
TGPPProteinpct2025_rast <- lapp(c(TGPPProteinpct2025_rast, TGPPAirborne_2025_rast[[1]]), fun = function(x, y){x[y == 0] <- NA; return(x)})
min(values(TGPPProteinpct2025_rast), na.rm = TRUE)
max(values(TGPPProteinpct2025_rast), na.rm = TRUE)
# plot(TGPPTNpct2025_rast)
# TGPPTNpct2025Crop_rast <- crop(TGPPTNpct2025_rast, TGPP_AOI, mask = TRUE)
# min(values(TGPPTNpct2025Crop_rast), na.rm = TRUE)
# max(values(TGPPTNpct2025Crop_rast), na.rm = TRUE)
# plot(TGPPTNpct2025Crop_rast)

### creating ggplot
TGPPProteinpct2025_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = TGPPProteinpct2025_rast) +
  scale_fill_whitebox_c(palette = "viridi", name = "Protein (%)") + #, limits = c(min(values(TGPPKpct2025_rast), na.rm = TRUE), max(values(TGPPKpct2025_rast), na.rm = TRUE))) +
  geom_sf(data = TGPP_QuadPts_FunBmass_sf, fill = NA, color = 'black') +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", pad_y = unit(1, units = "cm")) +#, pad_x = unit(1, units = "cm")) +
  ggtitle("Protein at TGPP in 2025") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
## displaying plot
TGPPProteinpct2025_plot
### saving plot
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/TGPPProteinpct2025.png")


### Saving Protein Raster
writeRaster(TGPPProteinpct2025_rast, paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Protein/MeanRaster/TGPP2025Proteinpct.tif"))
