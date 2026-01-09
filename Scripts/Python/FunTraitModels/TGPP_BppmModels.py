##### Functional Trait Model for Boron #####
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sbn
from sklearn.cross_decomposition import PLSRegression
# from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, ShuffleSplit, cross_validate
# from sklearn.linear_model import Lasso, Ridge, BayesianRidge
from sklearn.metrics import mean_squared_error
import os
import polars as pl
import geopandas as gpd
import xarray as xr
import rioxarray as rxr
# from sklearn.preprocessing import scale
import random
import time
import datetime
# import earthpy.plot as ep
from dask.distributed import LocalCluster

### measuring how long script takes to run
start = datetime.datetime.now()
# start = time.time()
### setting path to local data
data_dir = "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/"

##### reading in data #####
### reading in geopackage of functional traits, biomass variables, and extracted airborne reflectance data
TGPP_FunBmass_Airborne_gpd = gpd.read_file(os.path.join(data_dir, "GIS_Files/TGPP/ModelingData/TGPP2025BiomassFunTrait_Airborne.gpkg"))

### converting geodataframe to polars dataframe because it's easier to select the specific columns needed for modeling
## also polars very fast and very good, so yippee
TGPP_FunBmass_Airborne_df = (
    pl.from_pandas(TGPP_FunBmass_Airborne_gpd.drop(columns = ["geometry"]))
    .with_columns(pl.Series("geometry", TGPP_FunBmass_Airborne_gpd["geometry"].to_list()))
)

### reading in raster in chunks so it's hopefully less RAM intensive
TGPPAirborne_2025_rast = rxr.open_rasterio(os.path.join(data_dir, "Tallgrass2025_Airborne/AirborneClean/TIFF/1768_TGGP_Subset_REF_Mosaic_001-015_CleanBands.tif"), chunks = {'band': 215, 'x': 1500, 'y': 1500})

### splitting x (band reflectances) and y (Bppm, variable being modeled)
TGPP_Bppm_x = TGPP_FunBmass_Airborne_df.select(pl.selectors.starts_with("band"))
TGPP_Bppm_y = TGPP_FunBmass_Airborne_df.select("B_ppm")

### separating x and y into respective training and testing datasets
## uncommented to test bayesian ridge
# TGPP_Bppm_trainx, TGPP_Bppm_testx, TGPP_Bppm_trainy, TGPP_Bppm_testy = train_test_split(TGPP_Bppm_x, TGPP_Bppm_y, test_size = 0.3, random_state = 802)
# random.seed(802)
# BR = BayesianRidge()
# random.seed(802)
# BRfit = poopypantsBR.fit(TGPP_Bppm_trainx, TGPP_Bppm_trainy)

# BRfit.score(TGPP_Bppm_testx, TGPP_Bppm_testy)

### setting up Leave Group Out Cross Validation (LGOCV)
# Shuffle = ShuffleSplit(n_splits = 50, test_size = 0.3, random_state = 802)


##### Trying to replicate R/Caret's bootstrapping in Python #####
### scoring metrics to evaluate for different folds when doing cross-validation/bootstrapping
PLS_metrics = ["r2", "neg_root_mean_squared_error"]

### lists to store averages from each bootstrap
BootR2Avg_list = []
BootR2SD_list = []
BootRMSEAvg_list = []
BootRMSESD_list = []

### evaluating different numbers of components with 50-fold bootstrapping
for x in np.arange(1, 51):
    ### creating empty lists to store mean R^2 and mean test RMSE values from bootstrapping process
    boot_R2_list = []
    boot_rmse_list = []
    ### setting seed, initializing PLS model with some number of components (1 through length of for loop)
    random.seed(802)
    temp_pls = PLSRegression(n_components = x, scale = False)
    # temp_cv = cross_validate(temp_pls, X = TGPP_Bppm_trainx, y = TGPP_Bppm_trainy, cv = Shuffle, scoring = PLS_metrics, return_train_score = True, return_estimator = True, return_indices = True, n_jobs = 10)
    # temp_cv.keys()
    for i in np.arange(0, 50):
        ### separating x and y into respective training and testing datasets
        ## using the iteration of the for as the seed (adding it to 802 because that's my go-to seed, Vermont forever)
        TGPP_Bppm_trainx, TGPP_Bppm_testx, TGPP_Bppm_trainy, TGPP_Bppm_testy = train_test_split(TGPP_Bppm_x, TGPP_Bppm_y, test_size = 0.3, random_state = 802 + i)
        ### fit model
        random.seed(802)
        temp_pls_fit = temp_pls.fit(TGPP_Bppm_trainx, TGPP_Bppm_trainy)
        ### predict on bootstrapped test set
        random.seed(802)
        temp_preds = temp_pls.predict(TGPP_Bppm_testx)
        ### extracting relevant metrics (R2, RMSE)
        ## will be used to identify optimal number of components
        temp_rmse = np.sqrt(mean_squared_error(TGPP_Bppm_testy, temp_preds))
        boot_rmse_list.append(temp_rmse)
        boot_R2_list.append(temp_pls_fit.score(TGPP_Bppm_testx, TGPP_Bppm_testy))
    ### storing average R^2 and test RMSE across bootstraps for each unique number of components
    # OptimalComps_list.append(boot_R2_list.index(max(boot_R2_list)) + 1)
    BootR2Avg_list.append(np.mean(boot_R2_list))
    BootR2SD_list.append(np.std(boot_R2_list))
    # OptimalCompsRMSEInd_list.append(boot_rmse_list.index(min(boot_rmse_list)) + 1)
    BootRMSEAvg_list.append(np.mean(boot_rmse_list))
    BootRMSESD_list.append(np.std(boot_rmse_list))

### checking number of components with highest avg R^2
## adding 1 because python indexes from 0 and I'm not using 0 components
len(BootR2Avg_list)
len(boot_rmse_list)
BootR2Avg_list.index(max(BootR2Avg_list)) + 1
optimal_comps = BootR2Avg_list.index(max(BootR2Avg_list)) + 1
max(BootR2Avg_list)
BootR2SD_list[BootR2Avg_list.index(max(BootR2Avg_list))]
BootRMSEAvg_list.index(min(BootRMSEAvg_list)) + 1
min(BootRMSEAvg_list)
BootRMSESD_list[BootRMSEAvg_list.index(min(BootRMSEAvg_list))]


### plot of avg R2 values for each number of components
plt.plot(np.arange(1, 51, 1), BootR2Avg_list)
plt.xlabel("Number of Components")
# plt.xlim(left = 0, right = 50)
plt.xticks(np.arange(0, 51, 5))
plt.ylabel("R-Squared")
plt.title("Mean R-Squared Values for Each Number of Components \n Evaluated with PLSR of Boron Concentration \n (with 50 bootstraps)")
plt.close()

### plot of avg RMSE values for each number of components
plt.plot(np.arange(1, 51, 1), BootRMSEAvg_list)
plt.xlabel("Number of Components")
# plt.xlim(left = 0, right = 50)
plt.xticks(np.arange(0, 51, 5))
plt.ylabel("RMSE")
plt.title("Mean Test RMSE for Each Number of Components \n Evaluated with PLSR of Boron Concentration \n (with 50 bootstraps)")
plt.close()


##### Fitting 50 bootstraps of final model #####
### splitting x (band reflectances) and y (Bppm, variable being modeled)
## keeping quadrat number so I can use it to bind predictions from withheld sets to their corresponding points
TGPP_Bppm_x = TGPP_FunBmass_Airborne_df.select("Quadrat_name", pl.selectors.starts_with("band"))
TGPP_Bppm_y = TGPP_FunBmass_Airborne_df.select("B_ppm")
### initializing PLSR model
### fitting PLSR with the optimal number of components based on R2 and RMSE from for loop above
OptimalPLS = PLSRegression(n_components = optimal_comps, scale = False)
### making lists which will store arrays/lists of quadrats withheld, R2 and RMSE values, model coefficients, and model intercepts
TestQuadrats = []
TestPreds = []
boot_R2_list = []
boot_rmse_list = []
# ModelCoefs_list = []
# ModelIntercept_list = []

for i in np.arange(0, 50):
        ### separating x and y into respective training and testing datasets
        ## using the iteration of the for as the seed (adding it to 802 because that's my go-to seed, Vermont forever)
        TGPP_Bppm_trainx, TGPP_Bppm_testx, TGPP_Bppm_trainy, TGPP_Bppm_testy = train_test_split(TGPP_Bppm_x, TGPP_Bppm_y, test_size = 0.3, random_state = 802 + i)
        ### removing Quadrat_name from train_x and storing it in a list
        TestQuadrats.append(TGPP_Bppm_testx['Quadrat_name'])
        ### only keeping band reflectances for PLSR model
        TGPP_Bppm_trainx = TGPP_Bppm_trainx.select(pl.selectors.starts_with("band"))
        TGPP_Bppm_testx = TGPP_Bppm_testx.select(pl.selectors.starts_with("band"))
        ### fit model
        random.seed(802)
        OptimalPLS_fit = OptimalPLS.fit(TGPP_Bppm_trainx, TGPP_Bppm_trainy)
        ### predict on bootstrapped test set
        random.seed(802)
        temp_preds = OptimalPLS_fit.predict(TGPP_Bppm_testx)
        TestPreds.append(temp_preds)
        ### extracting relevant metrics (R2, RMSE)
        ## will be used to identify optimal number of components
        temp_rmse = np.sqrt(mean_squared_error(TGPP_Bppm_testy, temp_preds))
        boot_rmse_list.append(temp_rmse)
        boot_R2_list.append(OptimalPLS_fit.score(TGPP_Bppm_testx, TGPP_Bppm_testy))
        ### storing coefficients
        OptimalPLSCoefs = OptimalPLS_fit.coef_[0]
        OptimalPLSIntercept = OptimalPLS_fit.intercept_
        # OptimalPLSIntercept
        ##### Making spatial predictions #####
        ### setting up Dask client for parallelization
        # DaskCluster = LocalCluster(n_workers = 4, threads_per_worker = 1)
        # DaskClient = DaskCluster.get_client()
        # # print(DaskClient.dashboard_link)
        # ### multiplying bands by coefficient as a dask operation
        # TGPP_Bppm_array = (TGPPAirborne_2025_rast * OptimalPLSCoefs[:, None, None]
        #     ).sum(dim = "band") + OptimalPLSIntercept
        
        # ### Path to save the result
        # output_file_path = os.path.join(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Boron/", "TGPP2025Bppm_Pred" + str(i) + ".tif")
        # ### Start time of function
        # # time.strftime("%H:%M:%S", time.localtime())

        # ### Calling .compute() or .rio.to_raster() triggers Dask to compute the result chunk-by-chunk.
        # ### This performs the entire PLSR coefficients*bands math without loading the full intermediate array.
        # TGPP_Bppm_array.rio.to_raster(output_file_path, dtype=np.float32)

        # ### calling this in the hope that the time will print in the console when the above line finishes running
        # # time.strftime("%H:%M:%S", time.localtime())

        # ### closing Dask Client and cluster
        # DaskCluster.close()
        # DaskClient.close()
        ##### End of For Loop #####


### printing summary stats from bootstrapping
np.mean(boot_R2_list)
max(boot_R2_list)
min(boot_R2_list)
np.std(boot_R2_list)
np.mean(boot_rmse_list)
np.std(boot_rmse_list)

### making polars dataframe of quadrat names and corresponding prediction so I can plot mean prediction and standard error against actual y_values
### final dataframe combines all sets of withheld quadrat pts and their corresponding prediction for whichever bootstrap they were withheld from
for i in np.arange(0, len(TestQuadrats)):
    if i == 0:
        ### Storing first iteration as the main dataframe with all withheld quadrats and just binding all the later withheld sets to this one as one df
        TestPreds_df = pl.DataFrame({"Quadrat_name": TestQuadrats[i], "TestPreds": pl.from_numpy(TestPreds[i])})
    else:
        TestPreds2_df = pl.DataFrame({"Quadrat_name": TestQuadrats[i], "TestPreds": pl.from_numpy(TestPreds[i])})

        ### binding new iteration to TestPreds_df
        TestPreds_df = pl.concat([TestPreds_df, TestPreds2_df])


### taking full set of withheld bootstraps and getting averages and standard deviations
TestPredAvg_df = (
    TestPreds_df.group_by("Quadrat_name")
    .agg(
        pl.col("TestPreds").mean().alias("TestPredAvg"),
        pl.col("TestPreds").std().alias("TestPredSD")
    )
)

### joining df of validation set data made from bootstrapping with true values from original datset based on quadrat name
TestPredAvg_df = TestPredAvg_df.join(TGPP_FunBmass_Airborne_df.select("Quadrat_name", "B_ppm"), on = "Quadrat_name")

### plotting mean predictions with standard deviation error bars
plt.style.use("seaborn-v0_8-colorblind")
plt.figure(figsize=(10, 8))
# Create the scatter plot using matplotlib's errorbar for precise control
plt.errorbar(
    x = TestPredAvg_df["B_ppm"],
    y = TestPredAvg_df["TestPredAvg"],
    # yerr represents the +/- error from the mean (here, 1 standard deviation)
    yerr=TestPredAvg_df["TestPredSD"], 
    fmt='o',  # Format as circles
    color='darkblue',
    ecolor='black',
    capsize=4,
    alpha=0.6,
    label=r'Mean Prediction $\pm 1$ SD' # Use LaTeX for SD symbol
)
plt.title('Foliar Boron (ppm) \n PLSR Predicted vs. Observed Values \n (Bootstrap Averaged, k = 50)', fontsize=16)
plt.xlabel('Observed B (ppm)', fontsize = 12)
plt.ylabel('Mean Predicted B (ppm)', fontsize = 12)
# plt.xlim(min_val - padding, max_val + padding)
# plt.ylim(min_val - padding, max_val + padding)
plt.legend()
# plt.tight_layout()
# plt.show()
plt.savefig("/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/Outputs/Figures/FunctionalTraits/TGPP/2025/PLSRConfidenceIntervals/BppmSDPlot.png")
plt.close()

### measuring how long script takes to run
# end = time.time()
end = datetime.datetime.now()
print(end - start)
end - start
### reading in exported raster for plotting
# TGPP_Bppm_rast = rxr.open_rasterio(os.path.join(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/TGPP2025Bppm.tif"))

### trying to plot it
# ep.plot_bands(TGPP_Bppm_rast, cmap = 'viridis')
