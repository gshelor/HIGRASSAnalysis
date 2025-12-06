##### Calculating Rao's Q with Airborne Data for TGPP 2025 data #####

##### loading packages, reading in data #####
library(pacman)
p_load(tidyverse, sf, terra, tidyterra, ggspatial, rasterdiv, parallel, doParallel)

### local folder with all data in subdirectories
data_dir <- "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/"

### reading in rasters
TGPPTNpct2025_rast <- rast(paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Nitrogen/MeanRaster/TGPP2025TNpct.tif"))
TGPPPpct2025_rast <- rast(paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Phosphorus/MeanRaster/TGPP2025Ppct.tif"))
TGPPKpct2025_rast <- rast(paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Potassium/MeanRaster/TGPP2025Kpct.tif"))
plot(TGPPTNpct2025_rast)
plot(TGPPPpct2025_rast)
plot(TGPPKpct2025_rast)

### AOI polygon
TGPP_AOI <- read_sf(paste0(data_dir, "GIS_Files/TGPP/Northern_TGPP_2025.shp"))

### quadrat points
TGPP_QuadPts_FunBmass_sf <- read_sf(paste0(data_dir, "GIS_Files/TGPP/ModelingData/TGPP2025BiomassFunTrait_Airborne.gpkg"))

### combining individual functional trait layers into 1 spatraster
TGPP2025FunTraits_rast <- c(TGPPTNpct2025_rast, TGPPPpct2025_rast, TGPPKpct2025_rast)

# TGPP2025Rao_list <- paRao(x = TGPP2025FunTraits_rast, rasterOut = TRUE, np = 8, method = "multidimension")

TGPP2025Rao_rast <- TGPP2025Rao_list[[1]][[1]]
TGPP2025Rao_rast
plot(TGPP2025Rao_rast)


##### Reading in Rao's Q from Python Function #####
RaosQNoMask_rast <- rast(paste0(data_dir, 'FunctionalDiversity/TGPP/NoMask/TGPP2025RaoQ.tif'))
RaosQMask_rast <- rast(paste0(data_dir, 'FunctionalDiversity/TGPP/TGPP2025RaoQ.tif'))
plot(RaosQNoMask_rast)
plot(RaosQMask_rast)
plot(st_geometry(TGPP_AOI), add = TRUE)
plot(st_geometry(TGPP_QuadPts_FunBmass_sf), add = TRUE)

raosq_outlier <- quantile(values(RaosQMask_rast), probs = 0.99, na.rm = TRUE)

RaosQOutlierMask_rast <- app(RaosQMask_rast, fun = function(x){x[x > raosq_outlier] <- NA; return(x)})
plot(RaosQOutlierMask_rast)


RaosQOutlierMaskAOICrop_rast <- crop(RaosQOutlierMask_rast, TGPP_AOI, mask = TRUE)
plot(RaosQOutlierMaskAOICrop_rast)

### creating ggplot of Rao's Q layer
TGPPRaosQ2025_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = RaosQOutlierMaskAOICrop_rast) +
  scale_fill_whitebox_c(palette = "viridi", name = "Rao's Q") +#, limits = c(min(values(TGPPTNpct2025_rast_Intadj), na.rm = TRUE), max(values(TGPPTNpct2025_rast_Intadj), na.rm = TRUE))) +
  geom_sf(data = TGPP_QuadPts_FunBmass_sf, fill = NA, color = 'black') +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", pad_y = unit(1, units = "cm")) +#, pad_x = unit(1, units = "cm")) +
  ggtitle("Rao's Q at TGPP in 2025") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
## displaying plot
TGPPRaosQ2025_plot
### saving plot
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalDiversity/TGPP/TGPPRaosQ2025.png")
