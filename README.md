# WRF-Chem PM2.5 Mitigation Northern Thailand

Custom Python and NCL scripts for post-processing and analysis of WRF-Chem simulations for PM2.5 mitigation in Northern Thailand.

## Contents
- `plot_domains_SEA.py` — Figure 1 left panel
- `plot_delta_pm25pct_d02.ncl` — Figure 6 ΔPM2.5 spatial pattern
- `calc_pvalue.py` — station-level statistics and paired t-tests (Tables S4, S6)
- `calc_pvalue.ncl` — NCL version of p-value calculation
- `scale_ebu_25pct.ncl` — 25% biomass-burning emission scaling
- `namelist.wps`, `namelist.input` — WRF-Chem configuration files

## Requirements
- Python 3.x
- numpy, pandas, scipy, matplotlib, xarray, netCDF4
- NCL 6.6.2 or later
- WRF-Chem v4.4 (NCAR)

## Data
WRF-Chem output files are not included due to size. Scripts assume model output is available locally.

## Citation
If you use this code, please cite the manuscript and Zenodo archive.
