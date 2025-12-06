##### Caret-based Modeling Nitrogen Functional Traits Using TGPP Airborne Data #####
##### Loading packages #####
library(pacman)
p_load(tidyverse, here, brms, lme4, terra, sf, parallel, mcprogress, plsmod, ranger, xgboost, caret, caretEnsemble, pls, betareg, rsample, tidyterra, ggspatial)
options(mc.cores = ceiling(detectCores() / 3))

### local folder with all data in subdirectories
data_dir <- "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/"


##### Reading in Data, extracting airborne data to quadrat pts #####
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
### quadrat points
# TGPP_quadpts <- read_sf(paste0(data_dir, "GIS_Files/TGPP/QuadPoints_9PerPlot_TGPP_2025_Updated.shp"))
### reading in functional trait data
# TGPP_funtrait_df <- read_csv(paste0(data_dir, "FunctionalTraits/TGPP/TGPPFunctionalTraits2025.csv"))
# colnames(TGPP_funtrait_df) <- c("Quadrat_number", "TN_pct", "Protein_pct", "LabID2", "P_pct", "Ca_pct", "K_pct", "Mg_pct", "Na_pct", "S_pct", "Fe_ppm", "Zn_ppm", "Cu_ppm", "Mn_ppm", "B_ppm", "LabID3", "dry_weight_g")
### reading in biomass csv that I will join to functional trait data in order to join functional trait data to quadpts
# TGPP_EstBiomass_df <- read_csv(paste0(data_dir, "Biomass/Estimated/TGPP/CSVs/TGPP_2025_Avg_Height_Biomass_Litter.csv"))

##### Formatting data for modeling #####
### some non-numeric values in the functional trait spreadsheet sent to us, need to clean that up
# TGPP_funtrait_df$B_ppm = as.numeric(TGPP_funtrait_df$B_ppm)
# TGPP_funtrait_df$Na_pct = as.numeric(TGPP_funtrait_df$Na_pct)
# for (i in 1:nrow(TGPP_funtrait_df)){
#   if (is.na(TGPP_funtrait_df$Na_pct[i])){
#     TGPP_funtrait_df$Na_pct[i] = 0
#   }
#   if (is.na(TGPP_funtrait_df$B_ppm[i])){
#     TGPP_funtrait_df$B_ppm[i] = 0
#   }
# }


### binding functional trait df to biomass df, then binding that df to the quadpts
# TGPP_BiomassFunTrait_df <- full_join(TGPP_funtrait_df, TGPP_EstBiomass_df, by = "Quadrat_number")

### changing column name "plot_new" in the quad pts to match "Quadrat_name" from the biomass and joined biomass/functional trait dfs
# colnames(TGPP_quadpts)[9] <- "Quadrat_name"
### joining tabular biomass and functional trait data based on that common quadrat name column
# TGPP_QuadPts_FunBmass_sf <- full_join(TGPP_quadpts, TGPP_BiomassFunTrait_df, by = "Quadrat_name") #|>
  # drop_na()

### Extracting airborne raster values to quadrat pts
# TGPP_QuadPts_FunBmass_sf <- st_as_sf(terra::extract(TGPPAirborne_2025_rast, TGPP_QuadPts_FunBmass_sf, bind = TRUE, method = "bilinear")) #|>
  # drop_na()

### writing out gpkg with quadrat points, all biomass and functional trait data, and airborne reflectances
# write_sf(TGPP_QuadPts_FunBmass_sf, "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/GIS_Files/TGPP/ModelingData/TGPP2025BiomassFunTrait_Airborne.gpkg", append = FALSE)

### Now that extraction is done, reading in gpkg object of quadrat points with all functional traits and airborne band reflectances
TGPP_QuadPts_FunBmass_sf <- read_sf(paste0(data_dir, "GIS_Files/TGPP/ModelingData/TGPP2025BiomassFunTrait_Airborne.gpkg"))


##### Modeling Total Nitrogen based on airborne hyperspectral band reflectances #####
### splitting data into training and testing datasets
### filtering out columns not used in modeling
TGPP_Kpct_ModelingCols_df <- TGPP_QuadPts_FunBmass_sf |>
  dplyr::select(K_pct, starts_with("band")) |>
  st_drop_geometry()
### hamed asked for a csv because he refuses to believe that PLSR isn't a magical holy grail, if he asks for another one then add Quadrat_name to the selected columns and write out again
# write_csv(TGPP_Kpct_ModelingCols_df, "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/FunctionalTraits/TGPP/TGPP_Kpct_AirborneModel.csv")
set.seed(802)
TGPP_Kpct_split <- initial_split(TGPP_Kpct_ModelingCols_df, prop = 0.7, strata = "K_pct")
TGPP_Kpct_train <- training(TGPP_Kpct_split)
TGPP_Kpct_test <- testing(TGPP_Kpct_split)


##### simple linear model #####
set.seed(802)
Kpct_lm <- lm(K_pct ~ ., data = TGPP_Kpct_train)
summary(Kpct_lm)

### making predictions on test set to evaluate error
set.seed(802)
TGPP_Kpct_test$K_pct_lmpred <- predict(Kpct_lm, TGPP_Kpct_test)
TGPP_Kpct_lmRMSE = yardstick::rmse(TGPP_Kpct_test, truth = K_pct, estimate = K_pct_lmpred)
TGPP_Kpct_lmRMSE


##### PLSR Model for K_pct but using caret #####
### creating trainControl object
# Kpct_TrainCtrl <- trainControl(method = "boot", number = 50, search = "random", verboseIter = TRUE, allowParallel = TRUE, returnResamp = "final", returnData = FALSE, savePredictions = "final", predictionBounds = c(0, 100))

# ### training PLSR model with caret
# set.seed(802)
# Kpct_PLSFit <- train(form = K_pct ~ ., data = TGPP_Kpct_train, method = "pls", metric = "Rsquared", trControl = Kpct_TrainCtrl, tuneLength = 200)


# ### printing best hyperparameters
# Kpct_PLSFit$bestTune
# Kpct_PLSFinalFitResults <- Kpct_PLSFit$results |>
#   filter(ncomp == Kpct_PLSFit$bestTune$ncomp)
# Kpct_PLSFinalFitResults

# ### making predictions on test set
# set.seed(802)
# Kpct_pls_preds <- predict(Kpct_PLSFit, TGPP_Kpct_test)
# TGPP_Kpct_test$K_pct_PLSR_pred <- Kpct_pls_preds
# yardstick::rmse(TGPP_Kpct_test, K_pct, K_pct_PLSR_pred)


# Kpct_PLSFit_TrainPreds <- Kpct_PLSFit$pred
# Kpct_PLSFit_TrainPreds_Mean <- Kpct_PLSFit_TrainPreds |>
#   group_by(rowIndex) |>
#   summarise(pred_mean = mean(pred),
#             obs = mean(obs),
#             pred_mean_errorsd = sd(pred - obs))

# ### making a scatterplot of observed vs predicted with error bars representing 1 SD +- of the mean predicted value
# ggplot(data = Kpct_PLSFit_TrainPreds_Mean, mapping = aes(x = obs, y = pred_mean)) +
#   theme_bw() +
#   geom_point() +
#   geom_errorbar(aes(ymin = pred_mean - pred_mean_errorsd, ymax = pred_mean + pred_mean_errorsd)) +
#   xlab("Observed TN Value (%)") +
#   ylab("Predicted TN Value (%)") +
#   ggtitle("PLSR Predictions of TN (%) with Error Bars", subtitle = "Error Bars range between 1 standard deviation of predictions")

##### trying to bootstrap but using Leave-Group-Out CV because it gives more control over the size of the bootstrap #####

### creating trainControl object
Kpct_TrainCtrl_lgocv <- trainControl(method = "LGOCV", p = 0.7, number = 50, search = "random", verboseIter = TRUE, allowParallel = TRUE, returnResamp = "final", returnData = FALSE, savePredictions = "final", predictionBounds = c(0, 100))

### training PLSR model with caret
set.seed(802)
Kpct_PLSFit_lgocv <- train(form = K_pct ~ ., data = TGPP_Kpct_train, method = "pls", metric = "Rsquared", trControl = Kpct_TrainCtrl_lgocv, tuneLength = 200, center = TRUE, scale = FALSE)

### saving caret object as rds file
write_rds(Kpct_PLSFit_lgocv, paste0(data_dir, "FunctionalTraits/Models/RCaret/TGPPKpct2025.rds"), compress = "gz")

### extracting model results for best tune
Kpct_PLSFinalFitResults_lgocv <- Kpct_PLSFit_lgocv$results |>
  filter(ncomp == Kpct_PLSFit_lgocv$bestTune$ncomp)
Kpct_PLSFinalFitResults_lgocv

### storing model coefficients
Kpct_coefs <- coef(Kpct_PLSFit_lgocv$finalModel, intercept = TRUE)

### storing coefficients as dfs so I can multiply them by the bands in python instead since it's way faster
Kpct_coefs_df <- data.frame(coefficients = Kpct_coefs)
write_csv(Kpct_coefs_df, paste0(data_dir, "FunctionalTraits/Models/RCaret/Coefficients/TGPPKpct2025Coefs.csv"))

### making predictions on test set
set.seed(802)
Kpct_pls_lgocv_preds <- predict(Kpct_PLSFit_lgocv, TGPP_Kpct_test)
### adding predictions to test df for rmse evaluation
TGPP_Kpct_test$K_pct_PLSR_lgocv_pred <- Kpct_pls_lgocv_preds
yardstick::rmse(data = TGPP_Kpct_test, truth = K_pct, estimate = K_pct_PLSR_lgocv_pred)

### extracting predictions from each fold, 
Kpct_PLSFit_lgocv_TrainPreds <- Kpct_PLSFit_lgocv$pred
# table(Kpct_PLSFit_lgocv_TrainPreds$Resample)
Kpct_PLSFit_lgocv_TrainPreds_Mean <- Kpct_PLSFit_lgocv_TrainPreds |>
  group_by(rowIndex) |>
  summarise(pred_mean = mean(pred),
            obs = mean(obs),
            pred_mean_errorsd = sd(pred - obs))

### making a scatterplot of observed vs predicted with error bars representing 1 SD +- of the mean predicted value
KpctTrainPreds_SDPlot <- ggplot(data = Kpct_PLSFit_lgocv_TrainPreds_Mean, mapping = aes(x = obs, y = pred_mean)) +
  theme_bw() +
  geom_point() +
  geom_errorbar(aes(ymin = pred_mean - pred_mean_errorsd, ymax = pred_mean + pred_mean_errorsd)) +
  xlab("Observed TN Value (%)") +
  ylab("Predicted TN Value (%)") +
  ggtitle("PLSR Predictions of Potassium (%) with Error Bars", subtitle = "Error Bars range between 1 standard deviation of predictions")
### displaying plot
KpctTrainPreds_SDPlot
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/PLSRConfidenceIntervals/KpctSDPlot.png")

### making spatial predictions
# set.seed(802)
# TGPP_Kpct_rast <- terra::predict(TGPPAirborne_2025_rast, Kpct_PLSFit_lgocv, cores = detectCores() / 2)
### making function to make predictions by layer since predicting all at once is not working (too RAM intensive)
PLSRRasterPred_func <- function(coef_num, coefs, pred_rast){
  ### first, check correct coefficients
  ## coefficients are in order of: intercept, then band order
  ## so second coefficient corresponds with first band, and so on
  # if (length(coefs) != nlyr(pred_rast)){
  #   print("incorrect coefficient numbers, ensure correct model and raster are being used for predictions")
  #   break
  # }
  ### perform operation (since this is for PLSR models, should be similar to basic linear regression)
  ## take coefficient index, grab correct coefficient, then grab raster band at index below coefficient index, then multiply like in basic regression model
  PLSR_predlyr_rast <- coefs[coef_num] * pred_rast[[coef_num - 1]]
  
  ### return raster layer representing band reflectance multiplied by correct coefficient as output
  return(PLSR_predlyr_rast)
}


### Now calling raster prediction function
set.seed(802)
# Kpct_PLSRPredLyrs_list <- pmclapply(X = Kpct_coefs[2:length(Kpct_coefs)], FUN = PLSRRasterPred_func, mc.set.seed = FALSE, mc.cores = ceiling(detectCores() / 2), coefs = Kpct_coefs, pred_rast = TGPPAirborne_2025_rast)

### checking class of PLSRRasterPred_func outputs to make sure they are spatrasters
# class(Kpct_PLSRPredLyrs_list[[1]])

### fuck it, brute forcing the prediction since the pmclapply craps out at 7% no matter what
# TGPP_Kpct_rast <- (Kpct_coefs[2] * TGPPAirborne_2025_rast[[1]]) + (Kpct_coefs[3] * TGPPAirborne_2025_rast[[2]]) + (Kpct_coefs[4] * TGPPAirborne_2025_rast[[3]]) + (Kpct_coefs[5] * TGPPAirborne_2025_rast[[4]]) + (Kpct_coefs[6] * TGPPAirborne_2025_rast[[5]]) + (Kpct_coefs[7] * TGPPAirborne_2025_rast[[6]]) + (Kpct_coefs[8] * TGPPAirborne_2025_rast[[7]]) + (Kpct_coefs[9] * TGPPAirborne_2025_rast[[8]]) + (Kpct_coefs[10] * TGPPAirborne_2025_rast[[9]]) + (Kpct_coefs[11] * TGPPAirborne_2025_rast[[10]]) + (Kpct_coefs[12] * TGPPAirborne_2025_rast[[11]]) + (Kpct_coefs[13] * TGPPAirborne_2025_rast[[12]]) + (Kpct_coefs[14] * TGPPAirborne_2025_rast[[13]]) + (Kpct_coefs[15] * TGPPAirborne_2025_rast[[14]]) + (Kpct_coefs[16] * TGPPAirborne_2025_rast[[15]]) + (Kpct_coefs[17] * TGPPAirborne_2025_rast[[16]]) + (Kpct_coefs[18] * TGPPAirborne_2025_rast[[17]]) + (Kpct_coefs[19] * TGPPAirborne_2025_rast[[18]]) + (Kpct_coefs[20] * TGPPAirborne_2025_rast[[19]]) + (Kpct_coefs[21] * TGPPAirborne_2025_rast[[20]]) + (Kpct_coefs[22] * TGPPAirborne_2025_rast[[21]]) + (Kpct_coefs[23] * TGPPAirborne_2025_rast[[22]]) + (Kpct_coefs[24] * TGPPAirborne_2025_rast[[23]]) + (Kpct_coefs[25] * TGPPAirborne_2025_rast[[24]]) + (Kpct_coefs[26] * TGPPAirborne_2025_rast[[25]]) + (Kpct_coefs[27] * TGPPAirborne_2025_rast[[26]]) 
# ### R says the line is too long so I can't do it all at once, fuck it we ball
# TGPP_Kpct_rast <- TGPP_Kpct_rast + (Kpct_coefs[28] * TGPPAirborne_2025_rast[[27]]) + (Kpct_coefs[29] * TGPPAirborne_2025_rast[[28]]) + (Kpct_coefs[30] * TGPPAirborne_2025_rast[[29]]) + (Kpct_coefs[31] * TGPPAirborne_2025_rast[[30]]) + (Kpct_coefs[32] * TGPPAirborne_2025_rast[[31]]) + (Kpct_coefs[33] * TGPPAirborne_2025_rast[[32]]) + (Kpct_coefs[34] * TGPPAirborne_2025_rast[[33]]) + (Kpct_coefs[35] * TGPPAirborne_2025_rast[[34]]) + (Kpct_coefs[36] * TGPPAirborne_2025_rast[[35]]) + (Kpct_coefs[37] * TGPPAirborne_2025_rast[[36]]) + (Kpct_coefs[38] * TGPPAirborne_2025_rast[[37]]) + (Kpct_coefs[39] * TGPPAirborne_2025_rast[[38]]) + (Kpct_coefs[40] * TGPPAirborne_2025_rast[[39]]) + (Kpct_coefs[41] * TGPPAirborne_2025_rast[[40]]) + (Kpct_coefs[42] * TGPPAirborne_2025_rast[[41]]) + (Kpct_coefs[43] * TGPPAirborne_2025_rast[[42]]) + (Kpct_coefs[44] * TGPPAirborne_2025_rast[[43]]) + (Kpct_coefs[45] * TGPPAirborne_2025_rast[[44]]) + (Kpct_coefs[46] * TGPPAirborne_2025_rast[[45]]) + (Kpct_coefs[47] * TGPPAirborne_2025_rast[[46]]) + (Kpct_coefs[48] * TGPPAirborne_2025_rast[[47]]) + (Kpct_coefs[49] * TGPPAirborne_2025_rast[[48]]) + (Kpct_coefs[50] * TGPPAirborne_2025_rast[[49]]) + (Kpct_coefs[51] * TGPPAirborne_2025_rast[[50]]) + (Kpct_coefs[52] * TGPPAirborne_2025_rast[[51]]) + (Kpct_coefs[53] * TGPPAirborne_2025_rast[[52]]) + (Kpct_coefs[54] * TGPPAirborne_2025_rast[[53]]) + (Kpct_coefs[55] * TGPPAirborne_2025_rast[[54]]) + (Kpct_coefs[56] * TGPPAirborne_2025_rast[[55]]) + (Kpct_coefs[57] * TGPPAirborne_2025_rast[[56]]) + (Kpct_coefs[58] * TGPPAirborne_2025_rast[[57]]) + (Kpct_coefs[59] * TGPPAirborne_2025_rast[[58]]) + (Kpct_coefs[60] * TGPPAirborne_2025_rast[[59]]) + (Kpct_coefs[61] * TGPPAirborne_2025_rast[[60]]) + (Kpct_coefs[62] * TGPPAirborne_2025_rast[[61]]) + (Kpct_coefs[63] * TGPPAirborne_2025_rast[[62]]) + (Kpct_coefs[64] * TGPPAirborne_2025_rast[[63]]) + (Kpct_coefs[65] * TGPPAirborne_2025_rast[[64]]) + (Kpct_coefs[66] * TGPPAirborne_2025_rast[[65]]) + (Kpct_coefs[67] * TGPPAirborne_2025_rast[[66]]) + (Kpct_coefs[68] * TGPPAirborne_2025_rast[[67]]) + (Kpct_coefs[69] * TGPPAirborne_2025_rast[[68]]) + (Kpct_coefs[70] * TGPPAirborne_2025_rast[[69]]) + (Kpct_coefs[71] * TGPPAirborne_2025_rast[[70]]) + (Kpct_coefs[72] * TGPPAirborne_2025_rast[[71]]) + (Kpct_coefs[73] * TGPPAirborne_2025_rast[[72]]) + (Kpct_coefs[74] * TGPPAirborne_2025_rast[[73]]) + (Kpct_coefs[75] * TGPPAirborne_2025_rast[[74]]) + (Kpct_coefs[76] * TGPPAirborne_2025_rast[[75]]) + (Kpct_coefs[77] * TGPPAirborne_2025_rast[[76]]) + (Kpct_coefs[78] * TGPPAirborne_2025_rast[[77]]) + (Kpct_coefs[79] * TGPPAirborne_2025_rast[[78]]) + (Kpct_coefs[80] * TGPPAirborne_2025_rast[[79]]) + (Kpct_coefs[81] * TGPPAirborne_2025_rast[[80]]) + (Kpct_coefs[82] * TGPPAirborne_2025_rast[[81]]) + (Kpct_coefs[83] * TGPPAirborne_2025_rast[[82]]) + (Kpct_coefs[84] * TGPPAirborne_2025_rast[[83]]) + (Kpct_coefs[85] * TGPPAirborne_2025_rast[[84]]) + (Kpct_coefs[86] * TGPPAirborne_2025_rast[[85]]) + (Kpct_coefs[87] * TGPPAirborne_2025_rast[[86]]) + (Kpct_coefs[88] * TGPPAirborne_2025_rast[[87]]) + (Kpct_coefs[89] * TGPPAirborne_2025_rast[[88]]) + (Kpct_coefs[90] * TGPPAirborne_2025_rast[[89]]) + (Kpct_coefs[91] * TGPPAirborne_2025_rast[[90]]) + (Kpct_coefs[92] * TGPPAirborne_2025_rast[[91]]) + (Kpct_coefs[93] * TGPPAirborne_2025_rast[[92]]) + (Kpct_coefs[94] * TGPPAirborne_2025_rast[[93]]) + (Kpct_coefs[95] * TGPPAirborne_2025_rast[[94]]) + (Kpct_coefs[96] * TGPPAirborne_2025_rast[[95]]) + (Kpct_coefs[97] * TGPPAirborne_2025_rast[[96]]) + (Kpct_coefs[98] * TGPPAirborne_2025_rast[[97]]) + (Kpct_coefs[99] * TGPPAirborne_2025_rast[[98]]) + (Kpct_coefs[100] * TGPPAirborne_2025_rast[[99]]) + (Kpct_coefs[101] * TGPPAirborne_2025_rast[[100]]) 
# ### part 3 of adding coefficients * bands
# TGPP_Kpct_rast <- TGPP_Kpct_rast + (Kpct_coefs[102] * TGPPAirborne_2025_rast[[101]]) + (Kpct_coefs[103] * TGPPAirborne_2025_rast[[102]]) + (Kpct_coefs[104] * TGPPAirborne_2025_rast[[103]]) + (Kpct_coefs[105] * TGPPAirborne_2025_rast[[104]]) + (Kpct_coefs[106] * TGPPAirborne_2025_rast[[105]]) + (Kpct_coefs[107] * TGPPAirborne_2025_rast[[106]]) + (Kpct_coefs[108] * TGPPAirborne_2025_rast[[107]]) + (Kpct_coefs[109] * TGPPAirborne_2025_rast[[108]]) + (Kpct_coefs[110] * TGPPAirborne_2025_rast[[109]]) + (Kpct_coefs[111] * TGPPAirborne_2025_rast[[110]]) + (Kpct_coefs[112] * TGPPAirborne_2025_rast[[111]]) + (Kpct_coefs[113] * TGPPAirborne_2025_rast[[112]]) + (Kpct_coefs[114] * TGPPAirborne_2025_rast[[113]]) + (Kpct_coefs[115] * TGPPAirborne_2025_rast[[114]]) + (Kpct_coefs[116] * TGPPAirborne_2025_rast[[115]]) + (Kpct_coefs[117] * TGPPAirborne_2025_rast[[116]]) + (Kpct_coefs[118] * TGPPAirborne_2025_rast[[117]]) + (Kpct_coefs[119] * TGPPAirborne_2025_rast[[118]]) + (Kpct_coefs[120] * TGPPAirborne_2025_rast[[119]]) + (Kpct_coefs[121] * TGPPAirborne_2025_rast[[120]]) + (Kpct_coefs[122] * TGPPAirborne_2025_rast[[121]]) + (Kpct_coefs[123] * TGPPAirborne_2025_rast[[122]]) + (Kpct_coefs[124] * TGPPAirborne_2025_rast[[123]]) + (Kpct_coefs[125] * TGPPAirborne_2025_rast[[124]]) + (Kpct_coefs[126] * TGPPAirborne_2025_rast[[125]]) + (Kpct_coefs[127] * TGPPAirborne_2025_rast[[126]]) + (Kpct_coefs[128] * TGPPAirborne_2025_rast[[127]]) + (Kpct_coefs[129] * TGPPAirborne_2025_rast[[128]]) + (Kpct_coefs[130] * TGPPAirborne_2025_rast[[129]]) + (Kpct_coefs[131] * TGPPAirborne_2025_rast[[130]]) + (Kpct_coefs[132] * TGPPAirborne_2025_rast[[131]]) + (Kpct_coefs[133] * TGPPAirborne_2025_rast[[132]]) + (Kpct_coefs[134] * TGPPAirborne_2025_rast[[133]]) + (Kpct_coefs[135] * TGPPAirborne_2025_rast[[134]]) + (Kpct_coefs[136] * TGPPAirborne_2025_rast[[135]]) + (Kpct_coefs[137] * TGPPAirborne_2025_rast[[136]]) + (Kpct_coefs[138] * TGPPAirborne_2025_rast[[137]]) + (Kpct_coefs[139] * TGPPAirborne_2025_rast[[138]]) + (Kpct_coefs[140] * TGPPAirborne_2025_rast[[139]]) + (Kpct_coefs[141] * TGPPAirborne_2025_rast[[140]]) + (Kpct_coefs[142] * TGPPAirborne_2025_rast[[141]]) + (Kpct_coefs[143] * TGPPAirborne_2025_rast[[142]]) + (Kpct_coefs[144] * TGPPAirborne_2025_rast[[143]]) + (Kpct_coefs[145] * TGPPAirborne_2025_rast[[144]]) + (Kpct_coefs[146] * TGPPAirborne_2025_rast[[145]]) + (Kpct_coefs[147] * TGPPAirborne_2025_rast[[146]]) + (Kpct_coefs[148] * TGPPAirborne_2025_rast[[147]]) + (Kpct_coefs[149] * TGPPAirborne_2025_rast[[148]]) + (Kpct_coefs[150] * TGPPAirborne_2025_rast[[149]]) + (Kpct_coefs[151] * TGPPAirborne_2025_rast[[150]]) + (Kpct_coefs[152] * TGPPAirborne_2025_rast[[151]]) + (Kpct_coefs[153] * TGPPAirborne_2025_rast[[152]]) + (Kpct_coefs[154] * TGPPAirborne_2025_rast[[153]]) + (Kpct_coefs[155] * TGPPAirborne_2025_rast[[154]]) + (Kpct_coefs[156] * TGPPAirborne_2025_rast[[155]]) + (Kpct_coefs[157] * TGPPAirborne_2025_rast[[156]]) + (Kpct_coefs[158] * TGPPAirborne_2025_rast[[157]]) + (Kpct_coefs[159] * TGPPAirborne_2025_rast[[158]]) + (Kpct_coefs[160] * TGPPAirborne_2025_rast[[159]]) + (Kpct_coefs[161] * TGPPAirborne_2025_rast[[160]])
# ### part 4 of spatial predictions
# TGPP_Kpct_rast <- TGPP_Kpct_rast + (Kpct_coefs[162] * TGPPAirborne_2025_rast[[161]]) + (Kpct_coefs[163] * TGPPAirborne_2025_rast[[162]]) + (Kpct_coefs[164] * TGPPAirborne_2025_rast[[163]]) + (Kpct_coefs[165] * TGPPAirborne_2025_rast[[164]]) + (Kpct_coefs[166] * TGPPAirborne_2025_rast[[165]]) + (Kpct_coefs[167] * TGPPAirborne_2025_rast[[166]]) + (Kpct_coefs[168] * TGPPAirborne_2025_rast[[167]]) + (Kpct_coefs[169] * TGPPAirborne_2025_rast[[168]]) + (Kpct_coefs[170] * TGPPAirborne_2025_rast[[169]]) + (Kpct_coefs[171] * TGPPAirborne_2025_rast[[170]]) + (Kpct_coefs[172] * TGPPAirborne_2025_rast[[171]]) + (Kpct_coefs[173] * TGPPAirborne_2025_rast[[172]]) + (Kpct_coefs[174] * TGPPAirborne_2025_rast[[173]]) + (Kpct_coefs[175] * TGPPAirborne_2025_rast[[174]]) + (Kpct_coefs[176] * TGPPAirborne_2025_rast[[175]]) + (Kpct_coefs[177] * TGPPAirborne_2025_rast[[176]]) + (Kpct_coefs[178] * TGPPAirborne_2025_rast[[177]]) + (Kpct_coefs[179] * TGPPAirborne_2025_rast[[178]]) + (Kpct_coefs[180] * TGPPAirborne_2025_rast[[179]]) + (Kpct_coefs[181] * TGPPAirborne_2025_rast[[180]]) + (Kpct_coefs[182] * TGPPAirborne_2025_rast[[181]]) + (Kpct_coefs[183] * TGPPAirborne_2025_rast[[182]]) + (Kpct_coefs[184] * TGPPAirborne_2025_rast[[183]]) + (Kpct_coefs[185] * TGPPAirborne_2025_rast[[184]]) + (Kpct_coefs[186] * TGPPAirborne_2025_rast[[185]]) + (Kpct_coefs[187] * TGPPAirborne_2025_rast[[186]]) + (Kpct_coefs[188] * TGPPAirborne_2025_rast[[187]]) + (Kpct_coefs[189] * TGPPAirborne_2025_rast[[188]]) + (Kpct_coefs[190] * TGPPAirborne_2025_rast[[189]]) + (Kpct_coefs[191] * TGPPAirborne_2025_rast[[190]]) + (Kpct_coefs[192] * TGPPAirborne_2025_rast[[191]]) + (Kpct_coefs[193] * TGPPAirborne_2025_rast[[192]]) + (Kpct_coefs[194] * TGPPAirborne_2025_rast[[193]]) + (Kpct_coefs[195] * TGPPAirborne_2025_rast[[194]]) + (Kpct_coefs[196] * TGPPAirborne_2025_rast[[195]]) + (Kpct_coefs[197] * TGPPAirborne_2025_rast[[196]]) + (Kpct_coefs[198] * TGPPAirborne_2025_rast[[197]]) + (Kpct_coefs[199] * TGPPAirborne_2025_rast[[198]]) + (Kpct_coefs[200] * TGPPAirborne_2025_rast[[199]]) + (Kpct_coefs[201] * TGPPAirborne_2025_rast[[200]]) + (Kpct_coefs[202] * TGPPAirborne_2025_rast[[201]]) + (Kpct_coefs[203] * TGPPAirborne_2025_rast[[202]]) + (Kpct_coefs[204] * TGPPAirborne_2025_rast[[203]]) + (Kpct_coefs[205] * TGPPAirborne_2025_rast[[204]]) + (Kpct_coefs[206] * TGPPAirborne_2025_rast[[205]]) + (Kpct_coefs[207] * TGPPAirborne_2025_rast[[206]]) + (Kpct_coefs[208] * TGPPAirborne_2025_rast[[207]]) + (Kpct_coefs[209] * TGPPAirborne_2025_rast[[208]]) + (Kpct_coefs[210] * TGPPAirborne_2025_rast[[209]]) + (Kpct_coefs[211] * TGPPAirborne_2025_rast[[210]]) + (Kpct_coefs[212] * TGPPAirborne_2025_rast[[211]]) + (Kpct_coefs[213] * TGPPAirborne_2025_rast[[212]]) + (Kpct_coefs[214] * TGPPAirborne_2025_rast[[213]]) + (Kpct_coefs[215] * TGPPAirborne_2025_rast[[214]]) + (Kpct_coefs[216] * TGPPAirborne_2025_rast[[215]]) + Kpct_coefs[1]

### plotting raster as made in python using same process as above
TGPPKpct2025_rast <- rast(paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/TGPP2025Kpct.tif"))
min(values(TGPPKpct2025_rast))
max(values(TGPPKpct2025_rast))
plot(TGPPKpct2025_rast)
TGPPKpct2025Crop_rast <- crop(TGPPKpct2025_rast, TGPP_AOI, mask = TRUE)
min(values(TGPPKpct2025Crop_rast), na.rm = TRUE)
max(values(TGPPKpct2025Crop_rast), na.rm = TRUE)
plot(TGPPKpct2025Crop_rast)

### creating ggplot
TGPPKpct2025_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = TGPPKpct2025_rast) +
  scale_fill_whitebox_c(palette = "viridi", name = "TN (%)", limits = c(min(values(TGPPKpct2025_rast)), max(values(TGPPKpct2025_rast)))) +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", pad_x = unit(1, units = "cm"), pad_y = unit(1, units = "cm")) +
  ggtitle("Nitrogen Concentation at TGPP in June 2025") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
## displaying plot
TGPPKpct2025_plot
### saving plot
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/TGPPKpct2025.png")

 ##### Modeling Kpct with Random Forest #####
### fitting it with caret
## using train control object from PLSR since it's model-agnostic
set.seed(802)
# Kpct_rf_train <- train(form = K_pct ~ ., data = TGPP_Kpct_train, method = "ranger", metric = "Rsquared", trControl = Kpct_TrainCtrl, tuneLength = 200)

# Kpct_rf_train$bestTune
# Kpct_RFFinalFit_metrics <- Kpct_rf_train$results |>
#   filter(mtry == Kpct_rf_train$bestTune$mtry & splitrule == Kpct_rf_train$bestTune$splitrule & min.node.size == Kpct_rf_train$bestTune$min.node.size)
# Kpct_RFFinalFit_metrics

# Kpct_RF_FinalModel <- Kpct_rf_train$finalModel
# summary(Kpct_RF_FinalModel$variable.importance)

# ### making predictions on test set
# set.seed(802)
# Kpct_rf_caret_preds <- predict(Kpct_rf_train, TGPP_Kpct_test)
# TGPP_Kpct_test$K_pct_rf_pred_caret <- Kpct_rf_caret_preds
# yardstick::rmse(TGPP_Kpct_test, K_pct, K_pct_rf_pred_caret)


##### Modeling Kpct with xgBoost #####
### caret-based workflow
## using train control object from PLSR since it's model-agnostic
# set.seed(802)
# Kpct_xgb_train <- train(form = K_pct ~ ., data = TGPP_Kpct_train, method = "xgbTree", metric = "Rsquared", trControl = Kpct_TrainCtrl, tuneLength = 200)

# Kpct_xgb_train$bestTune
# Kpct_xgb_FinalFit <- Kpct_xgb_train$results |>
#   filter(nrounds == Kpct_xgb_train$bestTune$nrounds & max_depth == Kpct_xgb_train$bestTune$max_depth & eta == Kpct_xgb_train$bestTune$eta & gamma == Kpct_xgb_train$bestTune$gamma & colsample_bytree == Kpct_xgb_train$bestTune$colsample_bytree & min_child_weight == Kpct_xgb_train$bestTune$min_child_weight & subsample == Kpct_xgb_train$bestTune$subsample)
# Kpct_xgb_FinalFit

# ### making predictions on test set
# set.seed(802)
# Kpct_xgb_caret_preds <- predict(Kpct_xgb_train, TGPP_Kpct_test)
# TGPP_Kpct_test$K_pct_xgb_pred_caret <- Kpct_xgb_caret_preds
# print(yardstick::rmse(TGPP_Kpct_test, K_pct, K_pct_xgb_pred_caret))


##### printing R^2 and accuracy metrics for all models #####
summary(Kpct_lm)
# max(Kpct_pls_tune_r2metrics$.estimate)
# max(Kpct_rf_tune_r2metrics$.estimate)
### PLSR metrics when fit using caret
# Kpct_PLSFinalFitResults
Kpct_PLSFinalFitResults_lgocv
# Kpct_RFFinalFit
# Kpct_xgb_FinalFit
# max(Kpct_xgb_tune_r2metrics$.estimate, na.rm = TRUE)
yardstick::rmse(TGPP_Kpct_test, truth = K_pct, estimate = K_pct_lmpred)
# yardstick::rmse(TGPP_Kpct_test, truth = K_pct, estimate = K_pct_PLSR_pred)
# yardstick::rmse(TGPP_Kpct_test, truth = K_pct, estimate = K_pct_PLSR_pred)
yardstick::rmse(TGPP_Kpct_test, truth = K_pct, estimate = K_pct_PLSR_lgocv_pred)
# yardstick::rmse(TGPP_Kpct_test, truth = K_pct, estimate = K_pct_rf_pred_caret)
# yardstick::rmse(TGPP_Kpct_test, truth = K_pct, estimate = K_pct_RF_pred)
# yardstick::rmse(TGPP_Kpct_test, K_pct, K_pct_xgb_pred_caret)
# yardstick::rmse(TGPP_Kpct_test, truth = K_pct, estimate = K_pct_xgb_pred)



