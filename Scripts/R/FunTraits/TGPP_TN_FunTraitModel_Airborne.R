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
TGPP_TNpct_ModelingCols_df <- TGPP_QuadPts_FunBmass_sf |>
  dplyr::select(TN_pct, starts_with("band")) |>
  st_drop_geometry()
### hamed asked for a csv because he refuses to believe that PLSR isn't a magical holy grail, if he asks for another one then add Quadrat_name to the selected columns and write out again
# write_csv(TGPP_TNpct_ModelingCols_df, "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/FunctionalTraits/TGPP/TGPP_TNPct_AirborneModel.csv")
set.seed(802)
TGPP_TNpct_split <- initial_split(TGPP_TNpct_ModelingCols_df, prop = 0.7, strata = "TN_pct")
TGPP_TNpct_train <- training(TGPP_TNpct_split)
TGPP_TNpct_test <- testing(TGPP_TNpct_split)

##### simple linear model #####
set.seed(802)
TNPct_lm <- lm(TN_pct ~ ., data = TGPP_TNpct_train)
summary(TNPct_lm)

### making predictions on test set to evaluate error
set.seed(802)
TGPP_TNpct_test$TN_pct_lmpred <- predict(TNPct_lm, TGPP_TNpct_test)
TGPP_TNpct_lmRMSE = yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_lmpred)
TGPP_TNpct_lmRMSE


##### PLSR Model for TN_pct but using caret #####
### creating trainControl object
# TNpct_TrainCtrl <- trainControl(method = "boot", number = 50, search = "random", verboseIter = TRUE, allowParallel = TRUE, returnResamp = "final", returnData = FALSE, savePredictions = "final", predictionBounds = c(0, 100))

# ### training PLSR model with caret
# set.seed(802)
# TNpct_PLSFit <- train(form = TN_pct ~ ., data = TGPP_TNpct_train, method = "pls", metric = "Rsquared", trControl = TNpct_TrainCtrl, tuneLength = 200)


# ### printing best hyperparameters
# TNpct_PLSFit$bestTune
# TNpct_PLSFinalFitResults <- TNpct_PLSFit$results |>
#   filter(ncomp == TNpct_PLSFit$bestTune$ncomp)
# TNpct_PLSFinalFitResults

# ### making predictions on test set
# set.seed(802)
# TNpct_pls_preds <- predict(TNpct_PLSFit, TGPP_TNpct_test)
# TGPP_TNpct_test$TN_pct_PLSR_pred <- TNpct_pls_preds
# yardstick::rmse(TGPP_TNpct_test, TN_pct, TN_pct_PLSR_pred)


# TNpct_PLSFit_TrainPreds <- TNpct_PLSFit$pred
# TNpct_PLSFit_TrainPreds_Mean <- TNpct_PLSFit_TrainPreds |>
#   group_by(rowIndex) |>
#   summarise(pred_mean = mean(pred),
#             obs = mean(obs),
#             pred_mean_errorsd = sd(pred - obs))

# ### making a scatterplot of observed vs predicted with error bars representing 1 SD +- of the mean predicted value
# ggplot(data = TNpct_PLSFit_TrainPreds_Mean, mapping = aes(x = obs, y = pred_mean)) +
#   theme_bw() +
#   geom_point() +
#   geom_errorbar(aes(ymin = pred_mean - pred_mean_errorsd, ymax = pred_mean + pred_mean_errorsd)) +
#   xlab("Observed TN Value (%)") +
#   ylab("Predicted TN Value (%)") +
#   ggtitle("PLSR Predictions of TN (%) with Error Bars", subtitle = "Error Bars range between 1 standard deviation of predictions")

##### trying to bootstrap but using Leave-Group-Out CV because it gives more control over the size of the bootstrap #####

### creating trainControl object
TNpct_TrainCtrl_lgocv <- trainControl(method = "LGOCV", p = 0.7, number = 50, search = "random", verboseIter = TRUE, allowParallel = TRUE, returnResamp = "final", returnData = FALSE, savePredictions = "final", predictionBounds = c(0, 100))

### training PLSR model with caret
set.seed(802)
TNpct_PLSFit_lgocv <- train(form = TN_pct ~ ., data = TGPP_TNpct_train, method = "pls", metric = "Rsquared", trControl = TNpct_TrainCtrl_lgocv, tuneLength = 215, center = TRUE, scale = FALSE)

### saving caret object as rds file
# write_rds(TNpct_PLSFit_lgocv, paste0(data_dir, "FunctionalTraits/Models/RCaret/TGPPTNpct2025.rds"), compress = "gz")

### reading saved model in so I don't need to refit the model
TNpct_PLSFit_lgocv <- read_rds(paste0(data_dir, "FunctionalTraits/Models/RCaret/TGPPTNpct2025.rds"))

### extracting model results for best tune
TNpct_PLSFinalFitResults_lgocv <- TNpct_PLSFit_lgocv$results |>
  filter(ncomp == TNpct_PLSFit_lgocv$bestTune$ncomp)
TNpct_PLSFinalFitResults_lgocv

### storing model coefficients
TNpct_coefs <- coef(TNpct_PLSFit_lgocv$finalModel, intercept = TRUE)

### storing coefficients as dfs so I can multiply them by the bands in python instead since it's way faster
TNpct_coefs_df <- data.frame(coefficients = TNpct_coefs)
write_csv(TNpct_coefs_df, paste0(data_dir, "FunctionalTraits/Models/RCaret/Coefficients/TGPPTNpct2025Coefs.csv"))

### making predictions on test set
set.seed(802)
TNpct_pls_lgocv_preds <- predict(TNpct_PLSFit_lgocv, TGPP_TNpct_test)
### adding predictions to test df for rmse evaluation
TGPP_TNpct_test$TN_pct_PLSR_lgocv_pred <- TNpct_pls_lgocv_preds
yardstick::rmse(data = TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_PLSR_lgocv_pred)

### extracting predictions from each fold, 
TNpct_PLSFit_lgocv_TrainPreds <- TNpct_PLSFit_lgocv$pred
# table(TNpct_PLSFit_lgocv_TrainPreds$Resample)
TNpct_PLSFit_lgocv_TrainPreds_Mean <- TNpct_PLSFit_lgocv_TrainPreds |>
  group_by(rowIndex) |>
  summarise(pred_mean = mean(pred),
            obs = mean(obs),
            pred_mean_errorsd = sd(pred - obs))

### making a scatterplot of observed vs predicted with error bars representing 1 SD +- of the mean predicted value
TNpctTrainPreds_SDPlot <- ggplot(data = TNpct_PLSFit_lgocv_TrainPreds_Mean, mapping = aes(x = obs, y = pred_mean)) +
  theme_bw() +
  geom_point() +
  geom_errorbar(aes(ymin = pred_mean - pred_mean_errorsd, ymax = pred_mean + pred_mean_errorsd)) +
  xlab("Observed TN Value (%)") +
  ylab("Predicted TN Value (%)") +
  ggtitle("PLSR Predictions of TN (%) with Error Bars", subtitle = "Error Bars range between 1 standard deviation of predictions")
### displaying plot
TNpctTrainPreds_SDPlot
### saving plot after displaying
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/PLSRConfidenceIntervals/TNpctSDPlot.png")


### making spatial predictions
# set.seed(802)
# TGPP_TNpct_rast <- terra::predict(TGPPAirborne_2025_rast, TNpct_PLSFit_lgocv, cores = detectCores() / 2)
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
# TNpct_PLSRPredLyrs_list <- pmclapply(X = TNpct_coefs[2:length(TNpct_coefs)], FUN = PLSRRasterPred_func, mc.set.seed = FALSE, mc.cores = ceiling(detectCores() / 2), coefs = TNpct_coefs, pred_rast = TGPPAirborne_2025_rast)

### plotting raster as made in python using same process as above
TGPPTNpct2025_rast <- rast(paste0(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/TGPP2025TNpct_poopypants.tif"))
TGPPTNpct2025_rast_adj <- app(TGPPTNpct2025_rast, fun = function(x){x[x < 0] <- 0; return(x)})
TGPPTNpct2025_rast_Intadj <- lapp(c(TGPPTNpct2025_rast_adj, TGPPAirborne_2025_rast[[1]]), fun = function(x, y){x[y == 0] <- NA; return(x)})
min(values(TGPPTNpct2025_rast))
max(values(TGPPTNpct2025_rast))
plot(TGPPTNpct2025_rast_Intadj)
TGPPTNpct2025Crop_rast <- crop(TGPPTNpct2025_rast_Intadj, TGPP_AOI, mask = TRUE)
min(values(TGPPTNpct2025Crop_rast), na.rm = TRUE)
max(values(TGPPTNpct2025Crop_rast), na.rm = TRUE)
plot(TGPPTNpct2025Crop_rast)

### creating ggplot
TGPPTNpct2025_plot <- ggplot() +
  theme_bw() +
  geom_spatraster(data = TGPPTNpct2025_rast_Intadj) +
  scale_fill_whitebox_c(palette = "viridi", name = "TN (%)", limits = c(min(values(TGPPTNpct2025_rast_Intadj), na.rm = TRUE), max(values(TGPPTNpct2025_rast_Intadj), na.rm = TRUE))) +
  geom_sf(data = TGPP_QuadPts_FunBmass_sf, fill = NA, color = 'black') +
  annotation_scale(location = "bl") +
  annotation_north_arrow(location = "bl", pad_y = unit(1, units = "cm")) +#, pad_x = unit(1, units = "cm")) +
  ggtitle("Total Nitrogen at TGPP in 2025") +
  theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
## displaying plot
TGPPTNpct2025_plot
### saving plot
ggsave("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/TGPPTNpct2025.png")

 ##### Modeling TNPct with Random Forest #####
### fitting it with caret
## using train control object from PLSR since it's model-agnostic
set.seed(802)
# TNpct_rf_train <- train(form = TN_pct ~ ., data = TGPP_TNpct_train, method = "ranger", metric = "Rsquared", trControl = TNpct_TrainCtrl, tuneLength = 200)

# TNpct_rf_train$bestTune
# TNpct_RFFinalFit_metrics <- TNpct_rf_train$results |>
#   filter(mtry == TNpct_rf_train$bestTune$mtry & splitrule == TNpct_rf_train$bestTune$splitrule & min.node.size == TNpct_rf_train$bestTune$min.node.size)
# TNpct_RFFinalFit_metrics

# TNpct_RF_FinalModel <- TNpct_rf_train$finalModel
# summary(TNpct_RF_FinalModel$variable.importance)

# ### making predictions on test set
# set.seed(802)
# TNpct_rf_caret_preds <- predict(TNpct_rf_train, TGPP_TNpct_test)
# TGPP_TNpct_test$TN_pct_rf_pred_caret <- TNpct_rf_caret_preds
# yardstick::rmse(TGPP_TNpct_test, TN_pct, TN_pct_rf_pred_caret)


##### Modeling TNPct with xgBoost #####
### caret-based workflow
## using train control object from PLSR since it's model-agnostic
# set.seed(802)
# TNpct_xgb_train <- train(form = TN_pct ~ ., data = TGPP_TNpct_train, method = "xgbTree", metric = "Rsquared", trControl = TNpct_TrainCtrl, tuneLength = 200)

# TNpct_xgb_train$bestTune
# TNpct_xgb_FinalFit <- TNpct_xgb_train$results |>
#   filter(nrounds == TNpct_xgb_train$bestTune$nrounds & max_depth == TNpct_xgb_train$bestTune$max_depth & eta == TNpct_xgb_train$bestTune$eta & gamma == TNpct_xgb_train$bestTune$gamma & colsample_bytree == TNpct_xgb_train$bestTune$colsample_bytree & min_child_weight == TNpct_xgb_train$bestTune$min_child_weight & subsample == TNpct_xgb_train$bestTune$subsample)
# TNpct_xgb_FinalFit

# ### making predictions on test set
# set.seed(802)
# TNpct_xgb_caret_preds <- predict(TNpct_xgb_train, TGPP_TNpct_test)
# TGPP_TNpct_test$TN_pct_xgb_pred_caret <- TNpct_xgb_caret_preds
# print(yardstick::rmse(TGPP_TNpct_test, TN_pct, TN_pct_xgb_pred_caret))


##### printing R^2 and accuracy metrics for all models #####
summary(TNPct_lm)
# max(TNpct_pls_tune_r2metrics$.estimate)
# max(TNpct_rf_tune_r2metrics$.estimate)
### PLSR metrics when fit using caret
# TNpct_PLSFinalFitResults
TNpct_PLSFinalFitResults_lgocv
# TNpct_RFFinalFit
# TNpct_xgb_FinalFit
# max(TNpct_xgb_tune_r2metrics$.estimate, na.rm = TRUE)
yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_lmpred)
# yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_PLSR_pred)
# yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_PLSR_pred)
yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_PLSR_lgocv_pred)
# yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_rf_pred_caret)
# yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_RF_pred)
# yardstick::rmse(TGPP_TNpct_test, TN_pct, TN_pct_xgb_pred_caret)
# yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_xgb_pred)


# hist(TGPP_BiomassFunTrait_df$TN_pct, main = "TN (%)")
# ### plotting predictions vs real values from the test dataset for all models
# ggplot(data = TGPP_TNpct_test, mapping = aes(x = TN_pct, y = TN_pct_lmpred)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   ggtitle("OLS Linear Regression Test Set Predictions vs Real TN (%) Values")

# ggplot(data = TGPP_TNpct_test, mapping = aes(x = TN_pct, y = TN_pct_PLSR_pred)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   ggtitle("PLSR Test Set Predictions vs Real TN (%) Values")


# ggplot(data = TGPP_TNpct_test, mapping = aes(x = TN_pct_rf_pred_caret, y = TN_pct)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   ggtitle("Random Forest Test Set Predictions vs Real TN (%) Values")


# ggplot(data = TGPP_TNpct_test, mapping = aes(x = TN_pct_xgb_pred_caret, y = TN_pct)) +
#   geom_point() +
#   geom_smooth(method = "lm") +
#   ggtitle("xgBoost Test Set Predictions vs Real TN (%) Values")
# # plot(TGPP_TNpct_test$TN_pct_PLSR_pred, TGPP_TNpct_test$TN_pct)

