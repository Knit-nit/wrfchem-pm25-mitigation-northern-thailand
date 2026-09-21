; ================================================================
; calc_pvalue.ncl
; Description: Calculate p-values (paired t-test) from WRF-Chem output
;              and write the results to a CSV file.
; ================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

begin

  ndays = 15
  nt    = 0

  ;-----------------------------------------
  ; 1. Input files
  ;-----------------------------------------
  base_files = new(ndays, string)
  red25_files = new(ndays, string)

  do d = 0, ndays-1
    day = sprinti("%0.2i", d+1)
    base_files(d)  = "wrfout_d02_base_2020-02-" + day + "_00:00:00"
    red25_files(d) = "wrfout_d02_red25_2020-02-" + day + "_00:00:00"
  end do

  ;-----------------------------------------
  ; 2. PCD station coordinates (from Paper Table 1)
  ;-----------------------------------------
  nstations = 15
  station_names = (/"ST01", "ST02", "ST03", "ST04", "ST05", \
                    "ST06", "ST07", "ST08", "ST09", "ST10", \
                    "ST11", "ST12", "ST13", "ST14", "ST15"/)
  
  station_lat = (/18.84073, 18.79093, 18.27833, 18.25082, 18.42678, \
                  18.28263, 19.90922, 19.30455, 18.78888, 18.56719, \
                  18.12837, 19.20023, 20.42733, 19.57597, 16.73441/)
  
  station_lon = (/98.96978, 98.99000, 99.50647, 99.76395, 99.75763, \
                  99.65982, 99.82347, 97.97165, 100.77630, 99.03864, \
                  100.16240, 99.89305, 99.88384, 101.08150, 98.56696/)

  ;-----------------------------------------
  ; 3. Open first file to read coordinates
  ;-----------------------------------------
  f0 = addfile(base_files(0), "r")
  lat2d = f0->XLAT(0,:,:)
  lon2d = f0->XLONG(0,:,:)
  
  dims = dimsizes(lat2d)
  nlat = dims(0)
  nlon = dims(1)

  ;-----------------------------------------
  ; 4. Create arrays to store data
  ;-----------------------------------------
  base_ts  = new((/nstations, ndays/), float)
  red25_ts = new((/nstations, ndays/), float)

  ;-----------------------------------------
  ; 5. Find nearest grid point
  ;-----------------------------------------
  do s = 0, nstations-1
    diff = abs(lat2d - station_lat(s)) + abs(lon2d - station_lon(s))
    diff_1d = ndtooned(diff)
    ij_1d = minind(diff_1d)
    i = ij_1d / nlon
    j = ij_1d % nlon
    
    do d = 0, ndays-1
      f_base = addfile(base_files(d), "r")
      f_red25 = addfile(red25_files(d), "r")
      base_ts(s,d)  = f_base->PM2_5_DRY(nt,0,i,j)
      red25_ts(s,d) = f_red25->PM2_5_DRY(nt,0,i,j)
    end do
  end do

  ;-----------------------------------------
  ; 6. Calculate statistics
  ;-----------------------------------------
  print("========================================")
  print("CALCULATING PAIRED T-TEST...")
  print("========================================")

  ; Create arrays for storing results
  mean_base_arr  = new(nstations, float)
  mean_red25_arr = new(nstations, float)
  mean_delta_arr = new(nstations, float)
  pct_change_arr = new(nstations, float)
  t_stat_arr     = new(nstations, float)
  p_value_arr    = new(nstations, float)

  do s = 0, nstations-1
    mean_base  = avg(base_ts(s,:))
    mean_red25 = avg(red25_ts(s,:))
    mean_delta = mean_red25 - mean_base
    pct_change = (mean_delta / mean_base) * 100
    
    d = red25_ts(s,:) - base_ts(s,:)
    mean_d = avg(d)
    std_d = stddev(d)
    
    if (std_d .eq. 0) then
      t_stat = 0.0
      p_value = 1.0
    else
      t_stat = mean_d / (std_d / sqrt(ndays))
      df = ndays - 1
      p_value = 2.0 * (1.0 - cdf_t(abs(t_stat), df))
    end if
    
    mean_base_arr(s)  = mean_base
    mean_red25_arr(s) = mean_red25
    mean_delta_arr(s) = mean_delta
    pct_change_arr(s) = pct_change
    t_stat_arr(s)     = t_stat
    p_value_arr(s)    = p_value
  end do

  ;-----------------------------------------
  ; 7. Display results on screen
  ;-----------------------------------------
  print("")
  print("Results:")
  print("--------")
  do s = 0, nstations-1
    print(station_names(s) + ": p-value = " + p_value_arr(s))
  end do

  ;-----------------------------------------
  ; 8. Write CSV file
  ;-----------------------------------------
  print("")
  print("Writing CSV file...")
  
  filename = "pvalue_results.csv"
  system("rm -f " + filename)  ; Delete old file if exists
  
  ; Create header
  header = "Station,Mean_Base,Mean_Red25,Mean_Delta,Pct_Change,t_stat,p_value,Significant"
  write_table(filename, "w", [/"string"/], "%s")
  write_table(filename, "a", [/" "/], "%s")
  write_table(filename, "a", [/" "/], "%s")
  
  ; Write each row of data
  do s = 0, nstations-1
    if (p_value_arr(s) .lt. 0.05) then
      sig = "YES"
    else
      sig = "NO"
    end if
    
    line = station_names(s) + "," + \
           mean_base_arr(s) + "," + \
           mean_red25_arr(s) + "," + \
           mean_delta_arr(s) + "," + \
           pct_change_arr(s) + "," + \
           t_stat_arr(s) + "," + \
           p_value_arr(s) + "," + sig
           
    write_table(filename, "a", [/" "/], "%s")
  end do

  print("✅ Saved: " + filename)
  print("========================================")

end