; ================================================================
; plot_delta_pm25pct_d02.ncl
; Description: Plot the absolute and relative changes in PM2.5
;              between the Red25 and Base simulations.
;              Uses a shapefile for Thai province boundaries.
; ================================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/wrf/WRF_contributed.ncl"

begin
  ndays = 15
  nt    = 0

  ;-----------------------------------------
  ; 1. Input files setup
  ;-----------------------------------------
  base_files = new(ndays, string)
  red25_files = new(ndays, string)

  do d = 0, ndays-1
    day = sprinti("%0.2i", d+1)
    base_files(d)  = "wrfout_d02_base_2020-02-" + day + "_00:00:00"
    red25_files(d) = "wrfout_d02_red25_2020-02-" + day + "_00:00:00"
  end do

  f0 = addfile(base_files(0), "r")
  x0 = f0->PM2_5_DRY(nt,0,:,:)

  base_sum  = new(dimsizes(x0), float)
  red25_sum = new(dimsizes(x0), float)
  base_sum  = 0.0
  red25_sum = 0.0

  ;-----------------------------------------
  ; 2. Calculate mean values
  ;-----------------------------------------
  do d = 0, ndays-1
    f = addfile(base_files(d), "r")
    x = f->PM2_5_DRY(nt,0,:,:)
    base_sum = base_sum + x

    f2 = addfile(red25_files(d), "r")
    x2 = f2->PM2_5_DRY(nt,0,:,:)
    red25_sum = red25_sum + x2
  end do

  base_mean  = base_sum / ndays
  red25_mean = red25_sum / ndays

  delta_abs = red25_mean - base_mean
  delta_pct = (delta_abs / base_mean) * 100
  
  ; Cap percentage change to [-100, 100] for better visualization
  delta_pct = where(delta_pct .lt. -100, -100, delta_pct)
  delta_pct = where(delta_pct .gt. 100, 100, delta_pct)

  print("========================================")
  print("Delta PM2.5 Statistics (Red25 - Base):")
  print("Min absolute change: " + min(delta_abs) + " ug/m3")
  print("Max absolute change: " + max(delta_abs) + " ug/m3")
  print("Mean absolute change: " + avg(delta_abs) + " ug/m3")
  print("Min relative change: " + min(delta_pct) + " %")
  print("Max relative change: " + max(delta_pct) + " %")
  print("Mean relative change: " + avg(delta_pct) + " %")
  print("========================================")

  ;-----------------------------------------
  ; 3. Plotting setup
  ;-----------------------------------------
  wks = gsn_open_wks("png", "Fig6_delta_PM25_Feb2020")
  gsn_define_colormap(wks, "BlAqGrYeOrReVi200")

  res = True
  res@gsnMaximize         = True
  res@gsnSpreadColors     = True
  res@cnFillOn            = True
  res@cnLinesOn           = False
  res@cnLineLabelsOn      = False
  res@tfDoNDCOverlay      = True
  res@pmTickMarkDisplayMode = "Always"
  res@lbLabelBarOn        = True
  res@lbLabelFontHeightF  = 0.012

  levs_abs  = (/-85, -75, -65, -55, -45, -35, -25, -15, -5, 0, 5/)
  fills_abs = (/200, 190, 180, 160, 140, 120, 100, 80, 60, 40, 20, 2/)

  levs_pct  = (/-100, -90, -80, -70, -60, -50, -40, -30, -20, -10, 0/)
  fills_pct = (/200, 190, 180, 160, 140, 120, 100, 80, 60, 40, 20, 2/)

  ; Map limits (Indochina peninsula)
  res@mpMinLatF = 14.3383
  res@mpMaxLatF = 22.9328
  res@mpMinLonF = 95.1258
  res@mpMaxLonF = 104.206
  res@mpOutlineOn = True
  res@mpDataBaseVersion = "MediumRes"
  res@mpNationalLineColor = "black"
  res@mpNationalLineThicknessF = 2.0
  res@mpGridLineDashPattern = 2
  res@mpGridLatSpacingF = 2.0
  res@mpGridLonSpacingF = 2.0
  res@tmXTOn = False
  res@tmYROn = False
  res@tmXBLabelFontHeightF = 0.010
  res@tmYLLabelFontHeightF = 0.010

  ; Shapefile for Thai provinces
  shp_path = "/home/nittaya/Documents/MapsDomains/Thailand_provinces/Thailand_provinces.shp"
  
  lnres = True
  lnres@gsLineColor = "black"
  lnres@gsLineThicknessF = 1.0

  ;-----------------------------------------
  ; 4. Panel (a): Absolute Change
  ;-----------------------------------------
  res_abs = res
  res_abs@cnLevelSelectionMode = "ExplicitLevels"
  res_abs@cnLevels = levs_abs
  res_abs@cnFillColors = fills_abs
  res_abs@gsnLeftString = "(a) Absolute Change"
  res_abs@gsnLeftStringFontHeightF = 0.016
  res_abs@gsnRightString = ""
  res_abs@lbTitleString = "ug/m3"
  res_abs@lbTitleFontHeightF = 0.014

  plot_abs = gsn_csm_contour_map(wks, delta_abs, res_abs)
  shp_abs = gsn_add_shapefile_polylines(wks, plot_abs, shp_path, lnres)

  ;-----------------------------------------
  ; 5. Panel (b): Relative Change (%)
  ;-----------------------------------------
  res_pct = res
  res_pct@cnLevelSelectionMode = "ExplicitLevels"
  res_pct@cnLevels = levs_pct
  res_pct@cnFillColors = fills_pct
  res_pct@gsnLeftString = "(b) Relative Change"
  res_pct@gsnLeftStringFontHeightF = 0.016
  res_pct@gsnRightString = ""
  res_pct@lbTitleString = "%"
  res_pct@lbTitleFontHeightF = 0.014

  plot_pct = gsn_csm_contour_map(wks, delta_pct, res_pct)
  shp_pct = gsn_add_shapefile_polylines(wks, plot_pct, shp_path, lnres)

  ;-----------------------------------------
  ; 6. Panel assembly and Title
  ;-----------------------------------------
  pres = True
  pres@gsnMaximize = True
  pres@gsnPanelLabelBar = False
  pres@gsnPanelTop = 0.92
  pres@gsnPanelBottom = 0.08
  pres@gsnPanelFigureStringsFontHeightF = 0.014
  pres@gsnPanelFigureStringsJust = "TopLeft"
  pres@gsnPanelFigureStringsXOffset = 0.02
  pres@gsnPanelFigureStringsYOffset = 0.02

  gsn_panel(wks, (/plot_abs, plot_pct/), (/1, 2/), pres)

  txres = True
  txres@txFontHeightF = 0.020
  txres@txFontThicknessF = 2.0
  txres@txJust = "CenterCenter"

  gsn_text_ndc(wks, "Delta PM2.5 (Red25 - Base) - February 2020", 0.5, 0.96, txres)

  frame(wks)

  print("========================================")
  print("Saved figure: Fig6_delta_PM25_Feb2020.png")
  print("========================================")

end