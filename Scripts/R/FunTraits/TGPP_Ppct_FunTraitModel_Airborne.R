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
TGPP_Ppct_ModelingCols_df <- TGPP_QuadPts_FunBmass_sf |>
  dplyr::select(P_pct, starts_with("band")) |>
  st_drop_geometry()
### hamed asked for a csv because he refuses to believe that PLSR isn't a magical holy grail, if he asks for another one then add Quadrat_name to the selected columns and write out again
# write_csv(TGPP_Ppct_ModelingCols_df, "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/FunctionalTraits/TGPP/TGPP_Ppct_AirborneModel.csv")
set.seed(802)
TGPP_Ppct_split <- initial_split(TGPP_Ppct_ModelingCols_df, prop = 0.7, strata = "P_pct")
TGPP_Ppct_train <- training(TGPP_Ppct_split)
TGPP_Ppct_test <- testing(TGPP_Ppct_split)


##### simple linear model #####
set.seed(802)
Ppct_lm <- lm(P_pct ~ ., data = TGPP_Ppct_train)
summary(Ppct_lm)

### making predictions on test set to evaluate error
set.seed(802)
TGPP_Ppct_test$P_pct_lmpred <- predict(Ppct_lm, TGPP_Ppct_test)
TGPP_Ppct_lmRMSE = yardstick::rmse(TGPP_Ppct_test, truth = P_pct, estimate = P_pct_lmpred)
TGPP_Ppct_lmRMSE


##### PLSR Model for P_pct but using caret #####
### creating trainControl object
# Ppct_TrainCtrl <- trainControl(method = "boot", number = 50, search = "random", verboseIter = TRUE, allowParallel = TRUE, returnResamp = "final", returnData = FALSE, savePredictions = "final", predictionBounds = c(0, 100))

# ### training PLSR model with caret
# set.seed(802)
# Ppct_PLSFit <- train(form = P_pct ~ ., data = TGPP_Ppct_train, method = "pls", metric = "Rsquared", trControl = Ppct_TrainCtrl, tuneLength = 200)


# ### printing best hyperparameters
# Ppct_PLSFit$bestTune
# Ppct_PLSFinalFitResults <- Ppct_PLSFit$results |>
#   filter(ncomp == Ppct_PLSFit$bestTune$ncomp)
# Ppct_PLSFinalFitResults

# ### making predictions on test set
# set.seed(802)
# Ppct_pls_preds <- predict(Ppct_PLSFit, TGPP_Ppct_test)
# TGPP_Ppct_test$P_pct_PLSR_pred <- Ppct_pls_preds
# yardstick::rmse(TGPP_Ppct_test, P_pct, P_pct_PLSR_pred)


# Ppct_PLSFit_TrainPreds <- Ppct_PLSFit$pred
# Ppct_PLSFit_TrainPreds_Mean <- Ppct_PLSFit_TrainPreds |>
#   group_by(rowIndex) |>
#   summarise(pred_mean = mean(pred),
#             obs = mean(obs),
#             pred_mean_errorsd = sd(pred - obs))

# ### making a scatterplot of observed vs predicted with error bars representing 1 SD +- of the mean predicted value
# ggplot(data = Ppct_PLSFit_TrainPreds_Mean, mapping = aes(x = obs, y = pred_mean)) +
#   theme_bw() +
#   geom_point() +
#   geom_errorbar(aes(ymin = pred_mean - pred_mean_errorsd, ymax = pred_mean + pred_mean_errorsd)) +
#   xlab("Observed P Value (%)") +
#   ylab("Predicted P Value (%)") +
#   ggtitle("PLSR Predictions of P (%) with Error Bars", subtitle = "Error Bars range between 1 standard deviation of predictions")

##### trying to bootstrap but using Leave-Group-Out CV because it gives more control over the size of the bootstrap #####

### creating trainControl object
Ppct_TrainCtrl_lgocv <- trainControl(method = "LGOCV", p = 0.7, number = 50, search = "random", verboseIter = TRUE, allowParallel = TRUE, returnResamp = "final", returnData = FALSE, savePredictions = "final", predictionBounds = c(0, 100))

### training PLSR model with caret
set.seed(802)
Ppct_PLSFit_lgocv <- train(form = P_pct ~ ., data = TGPP_Ppct_train, method = "pls", metric = "Rsquared", trControl = Ppct_TrainCtrl_lgocv, tuneLength = 200, center = TRUE, scale = FALSE)

### saving caret object as rds file
# write_rds(Ppct_PLSFit_lgocv, paste0(data_dir, "FunctionalTraits/Models/RCaret/TGPPPpct2025.rds"), compress = "gz")

### reading saved model in so I don't need to refit the model
Ppct_PLSFit_lgocv <- read_rds(paste0(data_dir, "FunctionalTraits/Models/RCaret/TGPPPpct2025.rds"))

### extracting model results for best tune
Ppct_PLSFinalFitResults_lgocv <- Ppct_PLSFit_lgocv$results |>
  filter(ncomp == Ppct_PLSFit_lgocv$bestTune$ncomp)
Ppct_PLSFinalFitResults_lgocv

### storing model coefficients
Ppct_coefs <- coef(Ppct_PLSFit_lgocv$finalModel, intercept = TRUE)

### storing coefficients as dfs so I can multiply them by the bands in python instead since it's way faster
Ppct_coefs_df <- data.frame(coefficients = Ppct_coefs)
write_csv(Ppct_coefs_df, paste0(data_dir, "FunctionalTraits/Models/RCaret/Coefficients/TGPPPpct2025Coefs.csv"))

### making predictions on test set
set.seed(802)
Ppct_pls_lgocv_preds <- predict(Ppct_PLSFit_lgocv, TGPP_Ppct_test)
### adding predictions to test df for rmse evaluation
TGPP_Ppct_test$P_pct_PLSR_lgocv_pred <- Ppct_pls_lgocv_preds
yardstick::rmse(data = TGPP_Ppct_test, truth = P_pct, estimate = P_pct_PLSR_lgocv_pred)

### extracting predictions from each fold, 
Ppct_PLSFit_lgocv_TrainPreds <- Ppct_PLSFit_lgocv$pred
# table(Ppct_PLSFit_lgocv_TrainPreds$Resample)
Ppct_PLSFit_lgocv_TrainPreds_Mean <- Ppct_PLSFit_lgocv_TrainPreds |>
  group_by(rowIndex) |>
  summarise(pred_mean = mean(pred),
            obs = mean(obs),
            pred_mean_errorsd = sd(pred - obs))

### making a scatterplot of observed vs predicted with error bars representing 1 SD +- of the mean predicted value
PpctTrainPreds_SDPlot <- ggplot(data = Ppct_PLSFit_lgocv_TrainPreds_Mean, mapping = aes(x = obs, y = pred_mean)) +
  theme_bw() +
  geom_point() +
  geom_errorbar(aes(ymin = pred_mean - pred_mean_errorsd, ymax = pred_mean + pred_mean_errorsd)) +
  xlab("Observed P Value (%)") +
  ylab("Predicted P Value (%)") +
  ggtitle("PLSR Predictions of Phosphorous (%) with Error Bars", subtitle = "Error Bars range between 1 standard deviation of predictions")
### displaying plot
PpctTrainPreds_SDPlot
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/PLSRConfidenceIntervals/PpctSDPlot.png")


### making spatial predictions
# set.seed(802)
# TGPP_Ppct_rast <- terra::predict(TGPPAirborne_2025_rast, Ppct_PLSFit_lgocv, cores = detectCores() / 2)
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
# Ppct_PLSRPredLyrs_list <- pmclapply(X = Ppct_coefs[2:length(Ppct_coefs)], FUN = PLSRRasterPred_func, mc.set.seed = FALSE, mc.cores = ceiling(detectCores() / 2), coefs = Ppct_coefs, pred_rast = TGPPAirborne_2025_rast)

### checking class of PLSRRasterPred_func outputs to make sure they are spatrasters
# class(Ppct_PLSRPredLyrs_list[[1]])

### fuck it, brute forcing the prediction since the pmclapply craps out at 7% no matter what
# TGPP_Ppct_rast <- (Ppct_coefs[2] * TGPPAirborne_2025_rast[[1]]) + (Ppct_coefs[3] * TGPPAirborne_2025_rast[[2]]) + (Ppct_coefs[4] * TGPPAirborne_2025_rast[[3]]) + (Ppct_coefs[5] * TGPPAirborne_2025_rast[[4]]) + (Ppct_coefs[6] * TGPPAirborne_2025_rast[[5]]) + (Ppct_coefs[7] * TGPPAirborne_2025_rast[[6]]) + (Ppct_coefs[8] * TGPPAirborne_2025_rast[[7]]) + (Ppct_coefs[9] * TGPPAirborne_2025_rast[[8]]) + (Ppct_coefs[10] * TGPPAirborne_2025_rast[[9]]) + (Ppct_coefs[11] * TGPPAirborne_2025_rast[[10]]) + (Ppct_coefs[12] * TGPPAirborne_2025_rast[[11]]) + (Ppct_coefs[13] * TGPPAirborne_2025_rast[[12]]) + (Ppct_coefs[14] * TGPPAirborne_2025_rast[[13]]) + (Ppct_coefs[15] * TGPPAirborne_2025_rast[[14]]) + (Ppct_coefs[16] * TGPPAirborne_2025_rast[[15]]) + (Ppct_coefs[17] * TGPPAirborne_2025_rast[[16]]) + (Ppct_coefs[18] * TGPPAirborne_2025_rast[[17]]) + (Ppct_coefs[19] * TGPPAirborne_2025_rast[[18]]) + (Ppct_coefs[20] * TGPPAirborne_2025_rast[[19]]) + (Ppct_coefs[21] * TGPPAirborne_2025_rast[[20]]) + (Ppct_coefs[22] * TGPPAirborne_2025_rast[[21]]) + (Ppct_coefs[23] * TGPPAirborne_2025_rast[[22]]) + (Ppct_coefs[24] * TGPPAirborne_2025_rast[[23]]) + (Ppct_coefs[25] * TGPPAirborne_2025_rast[[24]]) + (Ppct_coefs[26] * TGPPAirborne_2025_rast[[25]]) + (Ppct_coefs[27] * TGPPAirborne_2025_rast[[26]]) 
# ### R says the line is too long so I can't do it all at once, fuck it we ball
# TGPP_Ppct_rast <- TGPP_Ppct_rast + (Ppct_coefs[28] * TGPPAirborne_2025_rast[[27]]) + (Ppct_coefs[29] * TGPPAirborne_2025_rast[[28]]) + (Ppct_coefs[30] * TGPPAirborne_2025_rast[[29]]) + (Ppct_coefs[31] * TGPPAirborne_2025_rast[[30]]) + (Ppct_coefs[32] * TGPPAirborne_2025_rast[[31]]) + (Ppct_coefs[33] * TGPPAirborne_2025_rast[[32]]) + (Ppct_coefs[34] * TGPPAirborne_2025_rast[[33]]) + (Ppct_coefs[35] * TGPPAirborne_2025_rast[[34]]) + (Ppct_coefs[36] * TGPPAirborne_2025_rast[[35]]) + (Ppct_coefs[37] * TGPPAirborne_2025_rast[[36]]) + (Ppct_coefs[38] * TGPPAirborne_2025_rast[[37]]) + (Ppct_coefs[39] * TGPPAirborne_2025_rast[[38]]) + (Ppct_coefs[40] * TGPPAirborne_2025_rast[[39]]) + (Ppct_coefs[41] * TGPPAirborne_2025_rast[[40]]) + (Ppct_coefs[42] * TGPPAirborne_2025_rast[[41]]) + (Ppct_coefs[43] * TGPPAirborne_2025_rast[[42]]) + (Ppct_coefs[44] * TGPPAirborne_2025_rast[[43]]) + (Ppct_coefs[45] * TGPPAirborne_2025_rast[[44]]) + (Ppct_coefs[46] * TGPPAirborne_2025_rast[[45]]) + (Ppct_coefs[47] * TGPPAirborne_2025_rast[[46]]) + (Ppct_coefs[48] * TGPPAirborne_2025_rast[[47]]) + (Ppct_coefs[49] * TGPPAirborne_2025_rast[[48]]) + (Ppct_coefs[50] * TGPPAirborne_2025_rast[[49]]) + (Ppct_coefs[51] * TGPPAirborne_2025_rast[[50]]) + (Ppct_coefs[52] * TGPPAirborne_2025_rast[[51]]) + (Ppct_coefs[53] * TGPPAirborne_2025_rast[[52]]) + (Ppct_coefs[54] * TGPPAirborne_2025_rast[[53]]) + (Ppct_coefs[55] * TGPPAirborne_2025_rast[[54]]) + (Ppct_coefs[56] * TGPPAirborne_2025_rast[[55]]) + (Ppct_coefs[57] * TGPPAirborne_2025_rast[[56]]) + (Ppct_coefs[58] * TGPPAirborne_2025_rast[[57]]) + (Ppct_coefs[59] * TGPPAirborne_2025_rast[[58]]) + (Ppct_coefs[60] * TGPPAirborne_2025_rast[[59]]) + (Ppct_coefs[61] * TGPPAirborne_2025_rast[[60]]) + (Ppct_coefs[62] * TGPPAirborne_2025_rast[[61]]) + (Ppct_coefs[63] * TGPPAirborne_2025_rast[[62]]) + (Ppct_coefs[64] * TGPPAirborne_2025_rast[[63]]) + (Ppct_coefs[65] * TGPPAirborne_2025_rast[[64]]) + (Ppct_coefs[66] * TGPPAirborne_2025_rast[[65]]) + (Ppct_coefs[67] * TGPPAirborne_2025_rast[[66]]) + (Ppct_coefs[68] * TGPPAirborne_2025_rast[[67]]) + (Ppct_coefs[69] * TGPPAirborne_2025_rast[[68]]) + (Ppct_coefs[70] * TGPPAirborne_2025_rast[[69]]) + (Ppct_coefs[71] * TGPPAirborne_2025_rast[[70]]) + (Ppct_coefs[72] * TGPPAirborne_2025_rast[[71]]) + (Ppct_coefs[73] * TGPPAirborne_2025_rast[[72]]) + (Ppct_coefs[74] * TGPPAirborne_2025_rast[[73]]) + (Ppct_coefs[75] * TGPPAirborne_2025_rast[[74]]) + (Ppct_coefs[76] * TGPPAirborne_2025_rast[[75]]) + (Ppct_coefs[77] * TGPPAirborne_2025_rast[[76]]) + (Ppct_coefs[78] * TGPPAirborne_2025_rast[[77]]) + (Ppct_coefs[79] * TGPPAirborne_2025_rast[[78]]) + (Ppct_coefs[80] * TGPPAirborne_2025_rast[[79]]) + (Ppct_coefs[81] * TGPPAirborne_2025_rast[[80]]) + (Ppct_coefs[82] * TGPPAirborne_2025_rast[[81]]) + (Ppct_coefs[83] * TGPPAirborne_2025_rast[[82]]) + (Ppct_coefs[84] * TGPPAirborne_2025_rast[[83]]) + (Ppct_coefs[85] * TGPPAirborne_2025_rast[[84]]) + (Ppct_coefs[86] * TGPPAirborne_2025_rast[[85]]) + (Ppct_coefs[87] * TGPPAirborne_2025_rast[[86]]) + (Ppct_coefs[88] * TGPPAirborne_2025_rast[[87]]) + (Ppct_coefs[89] * TGPPAirborne_2025_rast[[88]]) + (Ppct_coefs[90] * TGPPAirborne_2025_rast[[89]]) + (Ppct_coefs[91] * TGPPAirborne_2025_rast[[90]]) + (Ppct_coefs[92] * TGPPAirborne_2025_rast[[91]]) + (Ppct_coefs[93] * TGPPAirborne_2025_rast[[92]]) + (Ppct_coefs[94] * TGPPAirborne_2025_rast[[93]]) + (Ppct_coefs[95] * TGPPAirborne_2025_rast[[94]]) + (Ppct_coefs[96] * TGPPAirborne_2025_rast[[95]]) + (Ppct_coefs[97] * TGPPAirborne_2025_rast[[96]]) + (Ppct_coefs[98] * TGPPAirborne_2025_rast[[97]]) + (Ppct_coefs[99] * TGPPAirborne_2025_rast[[98]]) + (Ppct_coefs[100] * TGPPAirborne_2025_rast[[99]]) + (Ppct_coefs[101] * TGPPAirborne_2025_rast[[100]]) 
# ### part 3 of adding coefficients * bands
# TGPP_Ppct_rast <- TGPP_Ppct_rast + (Ppct_coefs[102] * TGPPAirborne_2025_rast[[101]]) + (Ppct_coefs[103] * TGPPAirborne_2025_rast[[102]]) + (Ppct_coefs[104] * TGPPAirborne_2025_rast[[103]]) + (Ppct_coefs[105] * TGPPAirborne_2025_rast[[104]]) + (Ppct_coefs[106] * TGPPAirborne_2025_rast[[105]]) + (Ppct_coefs[107] * TGPPAirborne_2025_rast[[106]]) + (Ppct_coefs[108] * TGPPAirborne_2025_rast[[107]]) + (Ppct_coefs[109] * TGPPAirborne_2025_rast[[108]]) + (Ppct_coefs[110] * TGPPAirborne_2025_rast[[109]]) + (Ppct_coefs[111] * TGPPAirborne_2025_rast[[110]]) + (Ppct_coefs[112] * TGPPAirborne_2025_rast[[111]]) + (Ppct_coefs[113] * TGPPAirborne_2025_rast[[112]]) + (Ppct_coefs[114] * TGPPAirborne_2025_rast[[113]]) + (Ppct_coefs[115] * TGPPAirborne_2025_rast[[114]]) + (Ppct_coefs[116] * TGPPAirborne_2025_rast[[115]]) + (Ppct_coefs[117] * TGPPAirborne_2025_rast[[116]]) + (Ppct_coefs[118] * TGPPAirborne_2025_rast[[117]]) + (Ppct_coefs[119] * TGPPAirborne_2025_rast[[118]]) + (Ppct_coefs[120] * TGPPAirborne_2025_rast[[119]]) + (Ppct_coefs[121] * TGPPAirborne_2025_rast[[120]]) + (Ppct_coefs[122] * TGPPAirborne_2025_rast[[121]]) + (Ppct_coefs[123] * TGPPAirborne_2025_rast[[122]]) + (Ppct_coefs[124] * TGPPAirborne_2025_rast[[123]]) + (Ppct_coefs[125] * TGPPAirborne_2025_rast[[124]]) + (Ppct_coefs[126] * TGPPAirborne_2025_rast[[125]]) + (Ppct_coefs[127] * TGPPAirborne_2025_rast[[126]]) + (Ppct_coefs[128] * TGPPAirborne_2025_rast[[127]]) + (Ppct_coefs[129] * TGPPAirborne_2025_rast[[128]]) + (Ppct_coefs[130] * TGPPAirborne_2025_rast[[129]]) + (Ppct_coefs[131] * TGPPAirborne_2025_rast[[130]]) + (Ppct_coefs[132] * TGPPAirborne_2025_rast[[131]]) + (Ppct_coefs[133] * TGPPAirborne_2025_rast[[132]]) + (Ppct_coefs[134] * TGPPAirborne_2025_rast[[133]]) + (Ppct_coefs[135] * TGPPAirborne_2025_rast[[134]]) + (Ppct_coefs[136] * TGPPAirborne_2025_rast[[135]]) + (Ppct_coefs[137] * TGPPAirborne_2025_rast[[136]]) + (Ppct_coefs[138] * TGPPAirborne_2025_rast[[137]]) + (Ppct_coefs[139] * TGPPAirborne_2025_rast[[138]]) + (Ppct_coefs[140] * TGPPAirborne_2025_rast[[139]]) + (Ppct_coefs[141] * TGPPAirborne_2025_rast[[140]]) + (Ppct_coefs[142] * TGPPAirborne_2025_rast[[141]]) + (Ppct_coefs[143] * TGPPAirborne_2025_rast[[142]]) + (Ppct_coefs[144] * TGPPAirborne_2025_rast[[143]]) + (Ppct_coefs[145] * TGPPAirborne_2025_rast[[144]]) + (Ppct_coefs[146] * TGPPAirborne_2025_rast[[145]]) + (Ppct_coefs[147] * TGPPAirborne_2025_rast[[146]]) + (Ppct_coefs[148] * TGPPAirborne_2025_rast[[147]]) + (Ppct_coefs[149] * TGPPAirborne_2025_rast[[148]]) + (Ppct_coefs[150] * TGPPAirborne_2025_rast[[149]]) + (Ppct_coefs[151] * TGPPAirborne_2025_rast[[150]]) + (Ppct_coefs[152] * TGPPAirborne_2025_rast[[151]]) + (Ppct_coefs[153] * TGPPAirborne_2025_rast[[152]]) + (Ppct_coefs[154] * TGPPAirborne_2025_rast[[153]]) + (Ppct_coefs[155] * TGPPAirborne_2025_rast[[154]]) + (Ppct_coefs[156] * TGPPAirborne_2025_rast[[155]]) + (Ppct_coefs[157] * TGPPAirborne_2025_rast[[156]]) + (Ppct_coefs[158] * TGPPAirborne_2025_rast[[157]]) + (Ppct_coefs[159] * TGPPAirborne_2025_rast[[158]]) + (Ppct_coefs[160] * TGPPAirborne_2025_rast[[159]]) + (Ppct_coefs[161] * TGPPAirborne_2025_rast[[160]])
# ### part 4 of spatial predictions
# TGPP_Ppct_rast <- TGPP_Ppct_rast + (Ppct_coefs[162] * TGPPAirborne_2025_rast[[161]]) + (Ppct_coefs[163] * TGPPAirborne_2025_rast[[162]]) + (Ppct_coefs[164] * TGPPAirborne_2025_rast[[163]]) + (Ppct_coefs[165] * TGPPAirborne_2025_rast[[164]]) + (Ppct_coefs[166] * TGPPAirborne_2025_rast[[165]]) + (Ppct_coefs[167] * TGPPAirborne_2025_rast[[166]]) + (Ppct_coefs[168] * TGPPAirborne_2025_rast[[167]]) + (Ppct_coefs[169] * TGPPAirborne_2025_rast[[168]]) + (Ppct_coefs[170] * TGPPAirborne_2025_rast[[169]]) + (Ppct_coefs[171] * TGPPAirborne_2025_rast[[170]]) + (Ppct_coefs[172] * TGPPAirborne_2025_rast[[171]]) + (Ppct_coefs[173] * TGPPAirborne_2025_rast[[172]]) + (Ppct_coefs[174] * TGPPAirborne_2025_rast[[173]]) + (Ppct_coefs[175] * TGPPAirborne_2025_rast[[174]]) + (Ppct_coefs[176] * TGPPAirborne_2025_rast[[175]]) + (Ppct_coefs[177] * TGPPAirborne_2025_rast[[176]]) + (Ppct_coefs[178] * TGPPAirborne_2025_rast[[177]]) + (Ppct_coefs[179] * TGPPAirborne_2025_rast[[178]]) + (Ppct_coefs[180] * TGPPAirborne_2025_rast[[179]]) + (Ppct_coefs[181] * TGPPAirborne_2025_rast[[180]]) + (Ppct_coefs[182] * TGPPAirborne_2025_rast[[181]]) + (Ppct_coefs[183] * TGPPAirborne_2025_rast[[182]]) + (Ppct_coefs[184] * TGPPAirborne_2025_rast[[183]]) + (Ppct_coefs[185] * TGPPAirborne_2025_rast[[184]]) + (Ppct_coefs[186] * TGPPAirborne_2025_rast[[185]]) + (Ppct_coefs[187] * TGPPAirborne_2025_rast[[186]]) + (Ppct_coefs[188] * TGPPAirborne_2025_rast[[187]]) + (Ppct_coefs[189] * TGPPAirborne_2025_rast[[188]]) + (Ppct_coefs[190] * TGPPAirborne_2025_rast[[189]]) + (Ppct_coefs[191] * TGPPAirborne_2025_rast[[190]]) + (Ppct_coefs[192] * TGPPAirborne_2025_rast[[191]]) + (Ppct_coefs[193] * TGPPAirborne_2025_rast[[192]]) + (Ppct_coefs[194] * TGPPAirborne_2025_rast[[193]]) + (Ppct_coefs[195] * TGPPAirborne_2025_rast[[194]]) + (Ppct_coefs[196] * TGPPAirborne_2025_rast[[195]]) + (Ppct_coefs[197] * TGPPAirborne_2025_rast[[196]]) + (Ppct_coefs[198] * TGPPAirborne_2025_rast[[197]]) + (Ppct_coefs[199] * TGPPAirborne_2025_rast[[198]]) + (Ppct_coefs[200] * TGPPAirborne_2025_rast[[199]]) + (Ppct_coefs[201] * TGPPAirborne_2025_rast[[200]]) + (Ppct_coefs[202] * TGPPAirborne_2025_rast[[201]]) + (Ppct_coefs[203] * TGPPAirborne_2025_rast[[202]]) + (Ppct_coefs[204] * TGPPAirborne_2025_rast[[203]]) + (Ppct_coefs[205] * TGPPAirborne_2025_rast[[204]]) + (Ppct_coefs[206] * TGPPAirborne_2025_rast[[205]]) + (Ppct_coefs[207] * TGPPAirborne_2025_rast[[206]]) + (Ppct_coefs[208] * TGPPAirborne_2025_rast[[207]]) + (Ppct_coefs[209] * TGPPAirborne_2025_rast[[208]]) + (Ppct_coefs[210] * TGPPAirborne_2025_rast[[209]]) + (Ppct_coefs[211] * TGPPAirborne_2025_rast[[210]]) + (Ppct_coefs[212] * TGPPAirborne_2025_rast[[211]]) + (Ppct_coefs[213] * TGPPAirborne_2025_rast[[212]]) + (Ppct_coefs[214] * TGPPAirborne_2025_rast[[213]]) + (Ppct_coefs[215] * TGPPAirborne_2025_rast[[214]]) + (Ppct_coefs[216] * TGPPAirborne_2025_rast[[215]]) + Ppct_coefs[1]

### plotting raster as made in python using same process as above
TGPPPpct2025_rast <- rast(paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/TGPP2025Ppct.tif"))
TGPPPpct2025_rast_adj <- app(TGPPPpct2025_rast, fun = function(x){x[x < 0] <- 0; return(x)})
TGPPPpct2025_rast_Intadj <- lapp(c(TGPPPpct2025_rast_adj, TGPPAirborne_2025_rast[[1]]), fun = function(x, y){x[y == 0] <- NA; return(x)})
min(values(TGPPPpct2025_rast))
max(values(TGPPPpct2025_rast))
plot(TGPPPpct2025_rast)
TGPPPpct2025Crop_rast <- crop(TGPPPpct2025_rast, TGPP_AOI, mask = TRUE)
min(values(TGPPPpct2025Crop_rast), na.rm = TRUE)
max(values(TGPPPpct2025Crop_rast), na.rm = TRUE)
plot(TGPPPpct2025Crop_rast)

### creating ggplot
TGPPPpct2025_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = TGPPPpct2025_rast_Intadj) +
  scale_fill_whitebox_c(palette = "viridi", name = "P (%)", limits = c(min(values(TGPPPpct2025_rast_Intadj)), max(values(TGPPPpct2025_rast_Intadj)))) +
  geom_sf(data = TGPP_QuadPts_FunBmass_sf, fill = NA, color = 'black') +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", pad_y = unit(1, units = "cm")) +
  ggtitle("Phosphorous (%) at TGPP in 2025") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
## displaying plot
TGPPPpct2025_plot
### saving plot
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/TGPPPpct2025.png")

 ##### Modeling Ppct with Random Forest #####
### fitting it with caret
## using train control object from PLSR since it's model-agnostic
set.seed(802)
# Ppct_rf_train <- train(form = P_pct ~ ., data = TGPP_Ppct_train, method = "ranger", metric = "Rsquared", trControl = Ppct_TrainCtrl, tuneLength = 200)

# Ppct_rf_train$bestTune
# Ppct_RFFinalFit_metrics <- Ppct_rf_train$results |>
#   filter(mtry == Ppct_rf_train$bestTune$mtry & splitrule == Ppct_rf_train$bestTune$splitrule & min.node.size == Ppct_rf_train$bestTune$min.node.size)
# Ppct_RFFinalFit_metrics

# Ppct_RF_FinalModel <- Ppct_rf_train$finalModel
# summary(Ppct_RF_FinalModel$variable.importance)

# ### making predictions on test set
# set.seed(802)
# Ppct_rf_caret_preds <- predict(Ppct_rf_train, TGPP_Ppct_test)
# TGPP_Ppct_test$P_pct_rf_pred_caret <- Ppct_rf_caret_preds
# yardstick::rmse(TGPP_Ppct_test, P_pct, P_pct_rf_pred_caret)


##### Modeling Ppct with xgBoost #####
### caret-based workflow
## using train control object from PLSR since it's model-agnostic
# set.seed(802)
# Ppct_xgb_train <- train(form = P_pct ~ ., data = TGPP_Ppct_train, method = "xgbTree", metric = "Rsquared", trControl = Ppct_TrainCtrl, tuneLength = 200)

# Ppct_xgb_train$bestTune
# Ppct_xgb_FinalFit <- Ppct_xgb_train$results |>
#   filter(nrounds == Ppct_xgb_train$bestTune$nrounds & max_depth == Ppct_xgb_train$bestTune$max_depth & eta == Ppct_xgb_train$bestTune$eta & gamma == Ppct_xgb_train$bestTune$gamma & colsample_bytree == Ppct_xgb_train$bestTune$colsample_bytree & min_child_weight == Ppct_xgb_train$bestTune$min_child_weight & subsample == Ppct_xgb_train$bestTune$subsample)
# Ppct_xgb_FinalFit

# ### making predictions on test set
# set.seed(802)
# Ppct_xgb_caret_preds <- predict(Ppct_xgb_train, TGPP_Ppct_test)
# TGPP_Ppct_test$P_pct_xgb_pred_caret <- Ppct_xgb_caret_preds
# print(yardstick::rmse(TGPP_Ppct_test, P_pct, P_pct_xgb_pred_caret))


##### printing R^2 and accuracy metrics for all models #####
summary(Ppct_lm)
# max(Ppct_pls_tune_r2metrics$.estimate)
# max(Ppct_rf_tune_r2metrics$.estimate)
### PLSR metrics when fit using caret
# Ppct_PLSFinalFitResults
Ppct_PLSFinalFitResults_lgocv
# Ppct_RFFinalFit
# Ppct_xgb_FinalFit
# max(Ppct_xgb_tune_r2metrics$.estimate, na.rm = TRUE)
yardstick::rmse(TGPP_Ppct_test, truth = P_pct, estimate = P_pct_lmpred)
# yardstick::rmse(TGPP_Ppct_test, truth = P_pct, estimate = P_pct_PLSR_pred)
# yardstick::rmse(TGPP_Ppct_test, truth = P_pct, estimate = P_pct_PLSR_pred)
yardstick::rmse(TGPP_Ppct_test, truth = P_pct, estimate = P_pct_PLSR_lgocv_pred)
# yardstick::rmse(TGPP_Ppct_test, truth = P_pct, estimate = P_pct_rf_pred_caret)
# yardstick::rmse(TGPP_Ppct_test, truth = P_pct, estimate = P_pct_RF_pred)
# yardstick::rmse(TGPP_Ppct_test, P_pct, P_pct_xgb_pred_caret)
# yardstick::rmse(TGPP_Ppct_test, truth = P_pct, estimate = P_pct_xgb_pred)


# hist(TGPP_BiomassFunTrait_df$P_pct, main = "P (%)")
# ### plotting predictions vs real values from the test dataset for all models
# ggplot(data = TGPP_Ppct_test, mapping = aes(x = P_pct, y = P_pct_lmpred)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   ggtitle("OLS Linear Regression Test Set Predictions vs Real P (%) Values")

# ggplot(data = TGPP_Ppct_test, mapping = aes(x = P_pct, y = P_pct_PLSR_pred)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   ggtitle("PLSR Test Set Predictions vs Real P (%) Values")


# ggplot(data = TGPP_Ppct_test, mapping = aes(x = P_pct_rf_pred_caret, y = P_pct)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   ggtitle("Random Forest Test Set Predictions vs Real P (%) Values")


# ggplot(data = TGPP_Ppct_test, mapping = aes(x = P_pct_xgb_pred_caret, y = P_pct)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   ggtitle("xgBoost Test Set Predictions vs Real P (%) Values")
# # plot(TGPP_Ppct_test$P_pct_PLSR_pred, TGPP_Ppct_test$P_pct)



