#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
calc_pvalue.py

Description: Calculate p-values (paired t-test) directly from WRF-Chem 
             output files using netCDF4 (avoids dask issues).
             
Author: Nittaya (Knit-nit)
"""

import os
import sys
import glob
import numpy as np
import pandas as pd
from scipy import stats
from netCDF4 import Dataset

# ================================================================
# 1. Input files
# ================================================================

data_dir = "./"
base_files = sorted(glob.glob(os.path.join(data_dir, "wrfout_d02_base_2020-02-*_00:00:00")))
red25_files = sorted(glob.glob(os.path.join(data_dir, "wrfout_d02_red25_2020-02-*_00:00:00")))

print(f"✅ Found Base files: {len(base_files)}")
print(f"✅ Found Red25 files: {len(red25_files)}")

if len(base_files) == 0 or len(red25_files) == 0:
    print("❌ No files found! Please check the directory.")
    sys.exit(1)

# ================================================================
# 2. PCD station coordinates
# ================================================================

stations = {
    'ST01': {'lat': 18.84073, 'lon': 98.96978},
    'ST02': {'lat': 18.79093, 'lon': 98.99000},
    'ST03': {'lat': 18.27833, 'lon': 99.50647},
    'ST04': {'lat': 18.25082, 'lon': 99.76395},
    'ST05': {'lat': 18.42678, 'lon': 99.75763},
    'ST06': {'lat': 18.28263, 'lon': 99.65982},
    'ST07': {'lat': 19.90922, 'lon': 99.82347},
    'ST08': {'lat': 19.30455, 'lon': 97.97165},
    'ST09': {'lat': 18.78888, 'lon': 100.77630},
    'ST10': {'lat': 18.56719, 'lon': 99.03864},
    'ST11': {'lat': 18.12837, 'lon': 100.16240},
    'ST12': {'lat': 19.20023, 'lon': 99.89305},
    'ST13': {'lat': 20.42733, 'lon': 99.88384},
    'ST14': {'lat': 19.57597, 'lon': 101.08150},
    'ST15': {'lat': 16.73441, 'lon': 98.56696},
}

# ================================================================
# 3. Read coordinates from the first file (using netCDF4)
# ================================================================

nc = Dataset(base_files[0], 'r')
lat2d = nc.variables['XLAT'][0, :, :]
lon2d = nc.variables['XLONG'][0, :, :]
nc.close()

nlat, nlon = lat2d.shape
print(f"✅ Grid size: {nlat} x {nlon}")

# ================================================================
# 4. Find nearest grid point
# ================================================================

def find_nearest(lat_target, lon_target):
    """Find the nearest grid indices for a given latitude and longitude."""
    lat_idx = np.argmin(np.abs(lat2d[:, 0] - lat_target))
    lon_idx = np.argmin(np.abs(lon2d[0, :] - lon_target))
    return lat_idx, lon_idx

# ================================================================
# 5. Extract data and calculate p-values
# ================================================================

results = []

print("\n" + "="*60)
print("PAIRED T-TEST: Base vs Red25")
print("="*60)

for name, coord in stations.items():
    lat_idx, lon_idx = find_nearest(coord['lat'], coord['lon'])
    
    # Extract time series from all files
    base_ts = []
    red25_ts = []
    
    for d in range(len(base_files)):
        nc_base = Dataset(base_files[d], 'r')
        nc_red25 = Dataset(red25_files[d], 'r')
        
        # Extract PM2.5_DRY at surface level (bottom_top=0)
        base_val = nc_base.variables['PM2_5_DRY'][0, 0, lat_idx, lon_idx]
        red25_val = nc_red25.variables['PM2_5_DRY'][0, 0, lat_idx, lon_idx]
        
        base_ts.append(base_val)
        red25_ts.append(red25_val)
        
        nc_base.close()
        nc_red25.close()
    
    base_ts = np.array(base_ts)
    red25_ts = np.array(red25_ts)
    
    # Calculate means
    mean_base = np.mean(base_ts)
    mean_red25 = np.mean(red25_ts)
    mean_delta = mean_red25 - mean_base
    pct_change = (mean_delta / mean_base) * 100
    
    # Paired t-test
    t_stat, p_value = stats.ttest_rel(base_ts, red25_ts)
    significant = p_value < 0.05
    
    results.append({
        'Station': name,
        'Mean_Base': mean_base,
        'Mean_Red25': mean_red25,
        'Mean_Delta': mean_delta,
        'Pct_Change': pct_change,
        't_stat': t_stat,
        'p_value': p_value,
        'Significant': significant
    })
    
    status = "✅" if significant else "❌"
    print(f"{status} {name}: p = {p_value:.4f}  (t = {t_stat:.3f})")

# ================================================================
# 6. Save to CSV
# ================================================================

df = pd.DataFrame(results)

# Round decimals
df['Mean_Base'] = df['Mean_Base'].round(2)
df['Mean_Red25'] = df['Mean_Red25'].round(2)
df['Mean_Delta'] = df['Mean_Delta'].round(2)
df['Pct_Change'] = df['Pct_Change'].round(2)
df['t_stat'] = df['t_stat'].round(3)
df['p_value'] = df['p_value'].round(4)

# Save to file
df.to_csv('pvalue_results.csv', index=False)
print("\n✅ Saved: pvalue_results.csv")

# ================================================================
# 7. Display summary table
# ================================================================

print("\n" + "="*60)
print("SUMMARY TABLE")
print("="*60)
print(df.to_string(index=False))

sig_count = sum(df['Significant'])
print(f"\n✅ Stations with significant change (p < 0.05): {sig_count} out of {len(df)}")
if sig_count == len(df):
    print("🎉 All stations show statistically significant changes!")
print("="*60)