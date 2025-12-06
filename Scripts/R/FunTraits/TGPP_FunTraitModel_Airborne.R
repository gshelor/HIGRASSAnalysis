##### Modeling Nitrogen Functional Traits Using TGPP Airborne Data #####

##### Loading packages #####
library(pacman)
p_load(tidyverse, here, brms, lme4, terra, sf, tidymodels, parallel, mcprogress, plsmod, mixOmics, stacks, ranger, xgboost)
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
TGPP_quadpts <- read_sf(paste0(data_dir, "GIS_Files/TGPP/QuadPoints_9PerPlot_TGPP_2025_Updated.shp"))
### reading in functional trait data
TGPP_funtrait_df <- read_csv(paste0(data_dir, "FunctionalTraits/TGPP/TGPPFunctionalTraits2025.csv"))
colnames(TGPP_funtrait_df) <- c("Quadrat_number", "TN_pct", "Protein_pct", "LabID2", "P_pct", "Ca_pct", "K_pct", "Mg_pct", "Na_pct", "S_pct", "Fe_ppm", "Zn_ppm", "Cu_ppm", "Mn_ppm", "B_ppm", "LabID3", "dry_weight_g")
### reading in biomass csv that I will join to functional trait data in order to join functional trait data to quadpts
TGPP_EstBiomass_df <- read_csv(paste0(data_dir, "Biomass/Estimated/CSVs/TGPP_2025_Avg_Height_Biomass_Litter.csv"))

##### Formatting data for modeling #####
### some non-numeric values in the functional trait spreadsheet sent to us, need to clean that up
TGPP_funtrait_df$B_ppm = as.numeric(TGPP_funtrait_df$B_ppm)
TGPP_funtrait_df$Na_pct = as.numeric(TGPP_funtrait_df$Na_pct)
for (i in 1:nrow(TGPP_funtrait_df)){
  if (is.na(TGPP_funtrait_df$Na_pct[i])){
    TGPP_funtrait_df$Na_pct[i] = 0
  }
  if (is.na(TGPP_funtrait_df$B_ppm[i])){
    TGPP_funtrait_df$B_ppm[i] = 0
  }
}


### binding functional trait df to biomass df, then binding that df to the quadpts
TGPP_BiomassFunTrait_df <- full_join(TGPP_funtrait_df, TGPP_EstBiomass_df, by = "Quadrat_number")

### changing column name "plot_new" in the quad pts to match "Quadrat_name" from the biomass and joined biomass/functional trait dfs
colnames(TGPP_quadpts)[9] <- "Quadrat_name"
### joining tabular biomass and functional trait data based on that common quadrat name column
TGPP_QuadPts_FunBmass_sf <- full_join(TGPP_quadpts, TGPP_BiomassFunTrait_df, by = "Quadrat_name") #|>
  # drop_na()

### Extracting airborne raster values to quadrat pts
TGPP_QuadPts_FunBmass_sf <- st_as_sf(terra::extract(TGPPAirborne_2025_rast, TGPP_QuadPts_FunBmass_sf, bind = TRUE)) #|>
  # drop_na()

##### Modeling Total Nitrogen based on airborne hyperspectral band reflectances #####
### splitting data into training and testing datasets
### filtering out columns not used in modeling
TGPP_TNpct_ModelingCols_df <- TGPP_QuadPts_FunBmass_sf |>
  dplyr::select(TN_pct, starts_with("band")) |>
  st_drop_geometry()
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


##### PLSR Model for TN_pct #####
### Create a recipe with the formula and preprocessing steps
TNpct_pls_recipe <- recipe(TN_pct ~ ., data = TGPP_TNpct_train) |>
  step_center(all_predictors()) |>
  step_scale(all_predictors())
### Define the PLS regression model
TNpct_pls_model <- parsnip::pls() |>
  set_engine("mixOmics") |>
  set_mode("regression") |>
  set_args(num_comp = tune(), predictor_prop = tune())


### Create a workflow
TNpct_pls_workflow <- workflow() |>
  add_model(TNpct_pls_model) |>
  add_recipe(TNpct_pls_recipe)

### Create a resampling object for cross-validation
set.seed(802)
TNpct_folds <- vfold_cv(TGPP_TNpct_train, v = 10)


### tuning hyperparameters with a random grid search
set.seed(802)
TNpct_pls_tune_result <- TNpct_pls_workflow |>
  tune_grid(resamples = TNpct_folds,
    ## number of random search combinations to try
    grid = 100,
    control = control_grid(save_pred = TRUE),
    metrics = metric_set(rsq, rmse, huber_loss)
  )



### Select the best hyperparameters based on the RMSE metric
TNpct_best_params <- select_best(TNpct_pls_tune_result, metric = "rmse")
TNpct_best_params

### looking at R2 metric for best hyperparamaters
TNpct_pls_tune_r2metrics <- as.data.frame(TNpct_pls_tune_result$.metrics) |>
  filter(.config == TNpct_best_params$.config) |>
  filter(.metric == "rsq")
TNpct_pls_tune_r2metrics$.estimate

### Finalize the workflow with the best parameters
TNpct_final_pls_workflow <- TNpct_pls_workflow |>
  finalize_workflow(TNpct_best_params)


### fitting model
set.seed(802)
TNpct_pls_fit <- fit(TNpct_final_pls_workflow, data = TGPP_TNpct_train) |>
  extract_fit_parsnip()
# TNpct_pls_fit
summary(TNpct_pls_fit)

### making predictions with PLSR model
set.seed(802)
TNpct_plsr_preds <- predict(TNpct_pls_fit, TGPP_TNpct_test)
TGPP_TNpct_test$TN_pct_PLSR_pred <- TNpct_plsr_preds$.pred
yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_PLSR_pred)


### making spatial predictions
set.seed(802)
# TN_pct_rast <- terra::predict(TGPPAirborne_2025_rast, TNpct_pls_fit)


##### Modeling TNPct with Random Forest #####
### create recipe for RF model
TNpct_rf_recipe <- recipe(TN_pct ~ ., data = TGPP_TNpct_train)

### Define the random forest model
TNpct_rf_mod <- rand_forest() |>
  set_engine("ranger", importance = "impurity", num.threads = detectCores() / 2) |>
  set_mode("regression") |>
  set_args(mtry = tune(),
          trees = tune(),
          min_n = tune())
  
### building random forest workflow with model specs and recipe
TNpct_rf_flow <- workflow() |>
  add_model(TNpct_rf_mod) |>
  add_recipe(TNpct_rf_recipe)


### tuning hyperparameters with a random grid search
TNpct_rf_tune_results <- TNpct_rf_flow |>
    tune_grid(resamples = TNpct_folds,
    ## number of random search combinations to try
    grid = 100,
    control = control_grid(save_pred = TRUE),
    metrics = metric_set(rsq, rsq_trad, rmse, huber_loss)
  )

### looking at R2 metrics
TNpct_rf_tune_r2metrics <- as.data.frame(TNpct_rf_tune_results$.metrics) |>
  filter(.metric == "rsq")
max(TNpct_rf_tune_r2metrics$.estimate)

### Select the best hyperparameters based on the RMSE metric
TNpct_rf_best_params <- select_best(TNpct_rf_tune_results, metric = "rmse")
TNpct_rf_best_params

### Finalize the workflow with the best parameters
TNpct_final_rf_flow <- TNpct_rf_flow |>
  finalize_workflow(TNpct_rf_best_params)

### fitting model
set.seed(802)
TNpct_rf_fit <- fit(TNpct_final_rf_flow, data = TGPP_TNpct_train) |>
  extract_fit_parsnip()
# TNpct_rf_fit
TNpct_rf_fit$fit$r.squared

### making predictions with Random Forest model
set.seed(802)
TNpct_rf_preds <- predict(TNpct_rf_fit, TGPP_TNpct_test)
TGPP_TNpct_test$TN_pct_RF_pred <- TNpct_rf_preds$.pred
yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_RF_pred)



##### Modeling TNPct with xgBoost #####
### create recipe for RF model
TNpct_xgb_recipe <- recipe(TN_pct ~ ., data = TGPP_TNpct_train)

### Define the random forest model
TNpct_xgb_mod <- boost_tree() |>
  set_engine("xgboost") |>
  set_mode("regression") |>
  set_args(tree_depth = tune(),
          trees = tune(),
          learn_rate = tune(),
          mtry = tune(),
          min_n = tune(),
          loss_reduction = tune(),
          sample_size = tune(),
          stop_iter = tune())
  
### building random forest workflow with model specs and recipe
TNpct_xgb_flow <- workflow() |>
  add_model(TNpct_xgb_mod) |>
  add_recipe(TNpct_xgb_recipe)


### tuning hyperparameters with a random grid search
TNpct_xgb_tune_results <- TNpct_xgb_flow |>
    tune_grid(resamples = TNpct_folds,
    ## number of random search combinations to try
    grid = 100,
    control = control_grid(save_pred = TRUE),
    metrics = metric_set(rsq, rsq_trad, rmse, huber_loss)
  )

### looking at R2 metrics
TNpct_xgb_tune_r2metrics <- as.data.frame(TNpct_xgb_tune_results$.metrics) |>
  filter(.metric == "rsq")
max(TNpct_xgb_tune_r2metrics$.estimate, na.rm = TRUE)

### Select the best hyperparameters based on the training RMSE metric
TNpct_xgb_best_params <- select_best(TNpct_xgb_tune_results, metric = "rmse")
TNpct_xgb_best_params

### Finalize the workflow with the best parameters
TNpct_final_xgb_flow <- TNpct_xgb_flow |>
  finalize_workflow(TNpct_xgb_best_params)

### fitting model
set.seed(802)
TNpct_xgb_fit <- fit(TNpct_final_xgb_flow, data = TGPP_TNpct_train) |>
  extract_fit_parsnip()
# TNpct_xgb_fit
summary(TNpct_xgb_fit$fit$best_score)

### making predictions with xgBoost model
set.seed(802)
TNpct_xgb_preds <- predict(TNpct_xgb_fit, TGPP_TNpct_test)
TGPP_TNpct_test$TN_pct_xgb_pred <- TNpct_xgb_preds$.pred

##### printing R^2 and accuracy metrics for all models #####
summary(TNPct_lm)
max(TNpct_pls_tune_r2metrics$.estimate)
max(TNpct_rf_tune_r2metrics$.estimate)
max(TNpct_xgb_tune_r2metrics$.estimate, na.rm = TRUE)
yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_lmpred)
yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_PLSR_pred)
yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_RF_pred)
yardstick::rmse(TGPP_TNpct_test, truth = TN_pct, estimate = TN_pct_xgb_pred)


hist(TGPP_BiomassFunTrait_df$TN_pct)
plot(TGPP_TNpct_test$TN_pct_RF_pred, TGPP_TNpct_test$TN_pct)
plot(TGPP_TNpct_test$TN_pct_PLSR_pred, TGPP_TNpct_test$TN_pct)
