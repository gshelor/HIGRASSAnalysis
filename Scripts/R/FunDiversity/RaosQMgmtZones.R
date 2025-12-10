##### Examining Variation in Rao's Q by different management treatments #####

##### loading packages, reading in data #####
library(pacman)
p_load(tidyverse, sf, terra, tidyterra, ggspatial, parallel, doParallel)

### local folder with all data in subdirectories
data_dir <- "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/"

### AOI polygon
TGPP_AOI <- read_sf(paste0(data_dir, "GIS_Files/TGPP/Northern_TGPP_2025.shp"))

### quadrat points
TGPP_QuadPts_FunBmass_sf <- read_sf(paste0(data_dir, "GIS_Files/TGPP/ModelingData/TGPP2025BiomassFunTrait_Airborne.gpkg"))

### TGPP management treatment zones, updated as of Dec 8, 2025
TGPP_2025MgmtZones_sf <- read_sf(paste0(data_dir, "GIS_Files/TGPP/MgmtZones/GPKG/TGPPMgmtZones.gpkg")) |>
  mutate(last_fire = paste(burn_season, last_fire_year))
TGPP_2025TSF_sf <- read_sf(paste0(data_dir, "GIS_Files/TGPP/MgmtZones/SHP/TGPP2025TSF.shp"))
TGPP_2025BurnSeason_sf <- read_sf(paste0(data_dir, "GIS_Files/TGPP/MgmtZones/SHP/TGPP2025BurnSeason.shp"))
TGPP_2025Herbicide_sf <- read_sf(paste0(data_dir, "GIS_Files/TGPP/MgmtZones/SHP/TGPP2025Herbicide.shp"))

### Rao's Q output
RaosQMask_rast <- rast(paste0(data_dir, 'FunctionalDiversity/TGPP/TGPP2025RaoQ.tif'))
plot(RaosQMask_rast)
plot(st_geometry(TGPP_AOI), add = TRUE)
plot(st_geometry(TGPP_QuadPts_FunBmass_sf), add = TRUE)
plot(st_geometry(TGPP_2025MgmtZones_sf), add = TRUE)

raosq_outlier <- quantile(values(RaosQMask_rast), probs = 0.99, na.rm = TRUE)

RaosQOutlierMask_rast <- app(RaosQMask_rast, fun = function(x){x[x > raosq_outlier] <- NA; return(x)})
plot(RaosQOutlierMask_rast)
plot(st_geometry(TGPP_2025MgmtZones_sf), add = TRUE)

RaosQOutlierMaskAOICrop_rast <- crop(RaosQOutlierMask_rast, TGPP_AOI, mask = TRUE)
plot(RaosQOutlierMaskAOICrop_rast)
plot(st_geometry(TGPP_2025MgmtZones_sf), add = TRUE)


TGPPMgmtZonesRaosQ <- st_as_sf(terra::extract(RaosQOutlierMaskAOICrop_rast, TGPP_2025MgmtZones_sf, fun = "mean", bind = TRUE, na.rm = TRUE))

TGPPBurnSeasonRaosQ <- st_as_sf(terra::extract(RaosQOutlierMaskAOICrop_rast, TGPP_2025BurnSeason_sf, fun = "mean", bind = TRUE, na.rm = TRUE))

TGPPTSFRaosQ <- st_as_sf(terra::extract(RaosQOutlierMaskAOICrop_rast, TGPP_2025TSF_sf, fun = "mean", bind = TRUE, na.rm = TRUE))

plot(TGPPMgmtZonesRaosQ["lyr.1"])
plot(TGPPBurnSeasonRaosQ["lyr.1"])
plot(TGPP_2025TSF_sf["years_sinc"])
plot(TGPPTSFRaosQ["lyr.1"])
head(TGPPMgmtZonesRaosQ)
head(TGPPBurnSeasonRaosQ)


##### Extracting Values to look at variance of Rao's Q by time since fire #####
### empty list to store Rao's Q values and different number of years since fire
years_since_fire_list <- list()
raosq_vals <- list()
for (i in unique(TGPP_2025MgmtZones_sf$last_fire)){
  ### filtering TSF layer so I can crop the Rao's Q layer to only pixels within the polygons for a given number of years since fire
  temp_TSF_sf <- TGPP_2025MgmtZones_sf |>
    filter(last_fire == i)

  ### cropping Rao's Q to the temp layer above
  temp_RaosQ_rast <- crop(RaosQOutlierMaskAOICrop_rast, temp_TSF_sf, mask = TRUE)
  ### taking temp Raos Q values to store in list
  temp_RaosQ_vals <- values(temp_RaosQ_rast)
  ### creating vector of same length as the values with the number of years since fire as the only value
  temp_tsf <- rep(i, times = length(temp_RaosQ_vals))

  ### adding tsf and rao's q values to lists which will store vectors of values for each tsf
  years_since_fire_list[[length(years_since_fire_list) + 1]] <- temp_tsf
  raosq_vals[[length(raosq_vals) + 1]] <- temp_RaosQ_vals

}

### initializing dataframe to store tsf and associated rao's q vals
## putting first set in it automatically, will append to it in for loop below
TGPPRaosQTSFVals_df <- data.frame(tsf = years_since_fire_list[[1]], raos_q = raosq_vals[[1]])
colnames(TGPPRaosQTSFVals_df) <- c("years_since_fire", "raos_q")

for (x in 2:length(raosq_vals)){
  temp_df <- data.frame(years_since_fire = years_since_fire_list[[x]], raos_q = raosq_vals[[x]])
  colnames(temp_df) <- c("years_since_fire", "raos_q")

  TGPPRaosQTSFVals_df <- rbind(TGPPRaosQTSFVals_df, temp_df)
}

TGPPRaosQTSFVals_df <- TGPPRaosQTSFVals_df |>
  drop_na()

TGPPRaosQTSFVals_df$last_fire_factor <- as.factor(TGPPRaosQTSFVals_df$years_since_fire)

RaosQ_violinplot <- ggplot(TGPPRaosQTSFVals_df, aes(x = last_fire_factor, y = raos_q, fill = last_fire_factor)) +
  theme_bw() +
  geom_violin() +
  scale_fill_brewer(palette = "Dark2") +
  geom_boxplot(width = 0.15, fill = "white") +
  scale_x_discrete(limits=c("Summer 2018", "Summer 2020", "Summer 2022", "Spring 2023", "Summer 2023", "Spring 2024", "Summer 2024", "Spring 2025"))
  # stat_summary(fun.data = mean_sdl, mult = 1, geom = "pointrange", color = "black")

RaosQ_violinplot
