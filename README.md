# WRF-Chem PM2.5 Mitigation Northern Thailand

Custom Python and NCL scripts for post-processing and analysis of WRF-Chem simulations for PM2.5 mitigation in Northern Thailand.

## Contents
- `plot_domains_SEA.py` — Figure 1 left panel
- `plot_delta_pm25pct_d02.ncl` — Figure 6 ΔPM2.5 spatial pattern
- `calc_pvalue.py` — station-level statistics and paired t-tests (Tables S4, S6)
- `calc_pvalue.ncl` — NCL version of p-value calculation
- `scale_ebu_25pct.ncl` — 25% biomass-burning emission scaling
- `namelist.wps`, `namelist.input` — WRF-Chem configuration files

## Software Requirements & Citations

The custom scripts in this repository rely on standard scientific Python and NCL libraries. Specifically, the statistical analysis script (`calc_pvalue.py`) uses the `scipy.stats` module to perform paired t-tests. If you use these scripts in your research, please also cite the relevant software dependencies:

- **SciPy:** Virtanen, P., Gommers, R., Oliphant, T. E., et al. SciPy 1.0: fundamental algorithms for scientific computing in Python. *Nature Methods*, 17(3), 261–272 (2020). https://doi.org/10.1038/s41592-019-0686-2
- **NumPy:** Harris, C. R., Millman, K. J., van der Walt, S. J., et al. Array programming with NumPy. *Nature*, 585, 357–362 (2020). https://doi.org/10.1038/s41586-020-2649-2

## Data
WRF-Chem output files are not included due to size. Scripts assume model output is available locally.

## Citation
If you use this code, please cite the manuscript and Zenodo archive.
DOI: [10.5281/zenodo.XXXXXXXX] (จะใส่ DOI จริงหลังจากสร้าง Release เสร็จแล้ว)
