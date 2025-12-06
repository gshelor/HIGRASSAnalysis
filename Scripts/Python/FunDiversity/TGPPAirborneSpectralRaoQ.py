import numpy as np
import pandas as pd
import math
from itertools import combinations
import rasterio as rio
from rasterio import plot
import matplotlib.pyplot as plt
import os
from numba import njit, prange, float64, int64, boolean, types

### local folder with all data in subdirectories
data_dir = "/home/gshelor/Documents/Schoolwork/OKSt/HIGRASS/AnalysisData/"



@njit(float64[:](float64[:], int64))
def compute_isNaN(inputArr, valueOfNan=1):
    # Return an array with valueOFNan instead of the NaN of the input array
    # NOTE: In Numba, this logic is safer and faster if rewritten for NumPy arrays.
    # The original implementation looked at the value of NaN, which is complex.
    # We redefine it to check for actual NaNs and return 1 where not NaN.
    
    # inputArr is already flattened (1D)
    outArr = np.zeros(inputArr.shape[0], dtype=np.float64)
    for i in range(inputArr.shape[0]):
        # Check if the value is not NaN. Numba requires np.isnan()
        if not np.isnan(inputArr[i]):
            outArr[i] = 1.0
    return outArr

@njit(parallel=True)
def _numba_rao_multidim_core(trastersm_3d, raoqe, w, window, na_tolerance, distance_m_flag):
    # trastersm_3d: (Bands, Padded_R, Padded_C)
    # raoqe: (R, C) - the output array to be modified in-place
    
    rows, cols = raoqe.shape
    num_bands = trastersm_3d.shape[0]
    window_pixels = window * window
    
    # Loop over original image dimensions (output matrix)
    for rw in prange(rows):
        for cl in prange(cols):
            # Coordinates in the PADDED matrix
            prw = rw + w
            pcl = cl + w
            
            # Extract window data: (Bands, Window, Window)
            window_data_3d = trastersm_3d[:, prw - w : prw + w + 1, pcl - w : pcl + w + 1]
            
            # Reshape to (Pixels_in_Window, Bands) for easier pairing (e.g., (81, 3))
            # .T is safe in Numba when chaining with reshape
            lv_data = window_data_3d.reshape(num_bands, window_pixels).T
            
            # 1. Check NaN condition (skip if too many NaN pixels)
            
            valid_pixels = 0
            for k in range(window_pixels):
                is_nan = False
                for b in range(num_bands):
                    if np.isnan(lv_data[k, b]):
                        is_nan = True
                        break
                if not is_nan:
                    valid_pixels += 1
            
            max_tolerated_nan = window_pixels * na_tolerance
            min_valid_pixels = window_pixels - max_tolerated_nan
            
            if valid_pixels < min_valid_pixels:
                # Value is NaN (default initialization), so we skip the calculation
                continue
            
            # 2. Core RAO Calculation
            
            window_sum = 0.0
            
            # Iterate over all unique pairs of pixel indices (k1, k2)
            for k1 in range(window_pixels):
                for k2 in range(k1 + 1, window_pixels):
                    
                    # Check for NaN in either pixel vector being compared (k1 or k2)
                    # If either pixel vector contains NaN in *any* band, the distance is undefined (or 0)
                    is_nan_pair = False
                    for b in range(num_bands):
                        if np.isnan(lv_data[k1, b]) or np.isnan(lv_data[k2, b]):
                            is_nan_pair = True
                            break
                    
                    if is_nan_pair:
                        continue
                        
                    # Calculate distance
                    dist = 0.0
                    
                    if distance_m_flag == 0: # Euclidean
                        squared_diff_sum = 0.0
                        for b in range(num_bands):
                            diff = lv_data[k1, b] - lv_data[k2, b]
                            squared_diff_sum += diff * diff
                        dist = np.sqrt(squared_diff_sum)
                        
                    elif distance_m_flag == 1: # Manhattan
                        abs_diff_sum = 0.0
                        for b in range(num_bands):
                            abs_diff_sum += np.abs(lv_data[k1, b] - lv_data[k2, b])
                        dist = abs_diff_sum
                        
                    window_sum += dist
            
            # Final result for the pixel: Sum of all distances * 2 / N^2
            # N = window_pixels
            raoqe[rw, cl] = window_sum * 2.0 / (window_pixels ** 2.0)
            
    # The array raoqe is modified in place
    return

# --- ORIGINAL HOST FUNCTION (Modified) ---

# Note: Keeping the original structure for non-Numba compatible parts (rasterio, pandas)
# The original distance functions (euclidean_distance, mana_distance) are removed 
# as their logic is now integrated into the JIT-compiled function.

def spectralrao(data_input, output_path, distance_m="euclidean", p=np.nan, window=9, mode="classic", na_tolerance=0.0, simplify=1, debugging=False):

    def tiff_to_numpy(tiff_input):
        # Convert Tiff input tu Numpy
        matrix1 = tiff_input.read()
        if matrix1.ndim == 3:
             # Assuming (Bands, Rows, Cols) format for Rasterio read() output
             matrix1 = matrix1[0, :, :] # Take the first band for classic mode, or handle multi-band
        
        # Reshape to (Rows, Cols) if it's (Bands, Rows, Cols) and we only need one band
        if matrix1.ndim == 3 and matrix1.shape[0] == 1:
            matrix1 = matrix1.reshape(matrix1.shape[1], matrix1.shape[2])
        elif matrix1.ndim == 3 and matrix1.shape[0] > 1:
            # For multidimension mode, we handle the list of 2D arrays later
            matrix1 = matrix1.transpose(1, 2, 0) # R, C, B for standard NumPy analysis
        
        minNum = -999
        # Assuming the nan replacement logic should only run if the value is explicitly -999
        matrix1[matrix1 == minNum] = np.nan
        return matrix1

    def export_tiff(naip_meta, output_rao, output_path):
        # Write the computation output on Tiff file
        naip_meta['count'] = 1
        naip_meta['dtype'] = "float64"
        with rio.open(output_path, 'w', **naip_meta) as dst:
            dst.write(output_rao, 1)

    # Simplified distance string to integer flag mapping for Numba
    distance_map = {"euclidean": 0, "mana": 1}
    distance_m_flag = distance_map.get(distance_m, 0) # Default to Euclidean (0)

    integer_dtypes = ("int8", "int16", "int32", "int64", "uint8", "uint16", "uint32", "uint64")
    float_dtypes = ("float_", "float16", "float32", "float64")

    # --- Setup Window Size ---
    if window % 2 == 1:
        w = int((window - 1) / 2)
    else:
        raise Exception("The size of moving window must be an odd number. Exiting...")

    # --- CLASSIC MODE (Only slight change for JIT compute_isNaN) ---
    if mode == "classic":  # one dimension input
        info = data_input.profile
        data_input = tiff_to_numpy(data_input)
        rasterm = np.copy(data_input)

        isfloat = False  
        if (rasterm.dtype.name not in integer_dtypes):
            print("Converting input data in an integer matrix...")
            isfloat = True
            mfactor = np.power(1000000, simplify)
            rasterm = rasterm * mfactor
            rasterm = rasterm.astype('int64')
        else:
            rasterm = rasterm.astype('int64')
        
        raoqe = np.zeros(shape=rasterm.shape)
        raoqe[:] = np.nan

        # Check if there are NaN value which becomes int
        nan_to_int64 = np.array([np.nan]).astype('int64')[0]

        intNaNPresence = nan_to_int64 in rasterm 
        s = pd.Series(rasterm.flatten('F'), dtype="category")
        values = [s.cat.codes[i] + 1 for i in range(len(s.cat.codes))]
        rasterm_1 = np.array(
            [[values[x + y * rasterm.shape[0]] for y in range(rasterm.shape[1])] for x in range(rasterm.shape[0])])

        trasterm = np.zeros(shape=(rasterm.shape[0] + 2 * w, rasterm.shape[1] + 2 * w))
        trasterm[:] = np.nan
        trasterm[w:w + rasterm.shape[0], w:w + rasterm.shape[1]] = rasterm_1[:]
        
        classes_values = np.array(s.cat.categories)
        nan_to_int64 = np.array([np.nan]).astype('int64')[0]

        for cl in range(w, w + rasterm.shape[1]):
            for rw in range(w, w + rasterm.shape[0]):
                print('Computing pixel {:d},{:d}'.format(rw, cl))

                borderCondition = np.sum(
                    np.invert(np.isnan(trasterm[rw - w:rw + w + 1, cl - w:cl + w + 1]))) < np.power(window, 2) - (
                                          (np.power(window, 2)) * na_tolerance)
                
                # Using JIT-compiled compute_isNaN
                tooManyIntNan = np.sum(compute_isNaN(trasterm[rw - w:rw + w + 1, cl - w:cl + w + 1])) < np.power(window,
                                                                                                                 2) - (
                                        (np.power(window, 2)) * na_tolerance)
                
                if (intNaNPresence and tooManyIntNan) or borderCondition:
                    pass
                else:  
                    tw = pd.Series(trasterm[rw - w:rw + w + 1, cl - w:cl + w + 1].flatten('F'),
                                   dtype="category").value_counts(ascending=True)
                    tw_values = np.array(tw.values)
                    tw_labels = np.array(tw.index)

                    if np.sum(np.isnan(np.array(tw.index))) > 0:
                        pass  

                    if debugging:
                        print("Working on coords ", rw, ",", cl, ". classes length: ", len(tw_values), ". window size=",
                              window)

                    if (len(tw_values) < 2):
                        pass 
                    else:
                        if (nan_to_int64 in tw_labels):
                            temp_index = np.where(tw_labels == nan_to_int64)
                            # NOTE: Original code had 'index' instead of 'temp_index'
                            p = tw_values / (np.sum(tw_values) - tw_values[temp_index]) 
                            p[temp_index] = 0
                        else:
                            p = tw_values / np.sum(tw_values)
                        p1 = np.zeros(shape=(len(tw_labels), len(tw_labels)))
                        for r in range(len(tw_labels)):
                            for c in range(r + 1, len(tw_labels)):
                                p1[r, c] = p[c] * p[r]

                        indices = np.array([int(el) - 1 for el in tw_labels])

                        d2 = np.zeros(shape=(len(indices), len(indices)))
                        for i, r in enumerate(indices):
                            for j, c in enumerate(indices):
                                if ((classes_values[r] == nan_to_int64) | (
                                        classes_values[c] == nan_to_int64)): 
                                    d2[i, j] = 0
                                else:
                                    d2[i, j] = classes_values[c] - classes_values[r]

                        if isfloat:
                            raoqe[rw - w, cl - w] = np.sum(p1 * d2) / mfactor
                        else:
                            raoqe[rw - w, cl - w] = np.sum(p1 * d2)

        print(("\nCalculation of Rao's index complete.\n"))


    # --- MULTIDIMENSION MODE (Numba Optimization Applied) ---
    elif mode == 'multidimension':  
        print("Starting multidimensional Rao's calculation (Numba optimized)...")
        info = data_input[0].profile
        
        # Convertion to Numpy
        numpy_data = []
        for input_raster in data_input:
            # We assume tiff_to_numpy returns (R, C) 2D array for a single band
            numpy_data.append(tiff_to_numpy(input_raster))

        if not numpy_data:
            raise Exception("No valid input data found.")

        rows, cols = numpy_data[0].shape
        num_bands = len(numpy_data)
        
        # Initialize output array
        raoqe = np.zeros(shape=(rows, cols), dtype=np.float64)
        raoqe[:] = np.nan
        
        # 1. Create a single 3D array (Bands, Padded_R, Padded_C) for Numba
        padded_rows = rows + 2 * w
        padded_cols = cols + 2 * w
        
        # Initialize 3D padded array
        trastersm_3d = np.full((num_bands, padded_rows, padded_cols), np.nan, dtype=np.float64)
        
        # Fill the 3D padded array band by band
        for b, mat in enumerate(numpy_data):
            trastersm_3d[b, w:w + rows, w:w + cols] = mat

        # 2. Call the Numba JIT-compiled core function
        _numba_rao_multidim_core(trastersm_3d, raoqe, w, window, na_tolerance, distance_m_flag)

        print("\nNumba calculation complete.\n")
        
    output = [raoqe]
    export_tiff(info, raoqe, output_path)
    return output

### reading in rasters
TGPPTNpct2025_rast = rio.open(os.path.join(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Nitrogen/MeanRaster/TGPP2025TNpct_Mask.tif"))
TGPPPpct2025_rast = rio.open(os.path.join(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Phosphorus/MeanRaster/TGPP2025Ppct_Mask.tif"))
TGPPKpct2025_rast = rio.open(os.path.join(data_dir, "FunctionalTraits/TGPP/Outputs/Rasters/Potassium/MeanRaster/TGPP2025Kpct_Mask.tif"))

### putting individual functional trait rasters in a list for calling spectralrao function
inputs = [TGPPTNpct2025_rast, TGPPPpct2025_rast, TGPPKpct2025_rast]

### calling spectralrao function
output = spectralrao(inputs, os.path.join(data_dir, 'FunctionalDiversity/TGPP/TGPP2025RaoQ.tif'), mode = "multidimension", window = 5, na_tolerance = 0)
