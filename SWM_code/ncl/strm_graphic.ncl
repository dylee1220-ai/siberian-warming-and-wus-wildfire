load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

; =========================================================
; OPEN FILE
; =========================================================

  WksType = "x11"
  filo = "../output/graph/strm_test"

  latS = -90
  latN = 90.
  lonL = 0.
  lonR = 360.

  f = addfile("../output/strm.heat.nc","r")
  x1 = f->strm(:,:,:,:)   ; 3 =0.86605 / 7 =0.46 / 9 =0.2569
  x1 = x1/10^6

; =========================================================
; PLOT
; =========================================================

  wks = gsn_open_wks(WksType,filo)
  gsn_define_colormap(wks,"BlueWhiteOrangeRed")

  res                      = True
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@gsnMaximize          = True
  res@gsnSpreadColors      = True         ; spread out color table
  res@vpHeightF            = 5.0

  res@mpFillOn             = False        ; turn off map fill
;  res@gsnPolar             = "SH"
;  res@mpMaxLatF            = -20.
;  res@mpMinLatF            = -90.
  res@gsnAddCyclic         = False
  res@mpCenterLonF         = 180.
  res@mpGeophysicalLineThicknessF = 1
  res@mpGridLineThicknessF  = 1.
  res@mpPerimLineThicknessF = 2.
  res@mpDataBaseVersion     = "MediumRes"

  res@cnFillOn             = True         ; turn on color fill
  res@cnLinesOn            = False        ; True is default
  res@cnLineLabelsOn       = False        ; True is default
  res@lbLabelBarOn         = False        ; turn off individual lb's
                                          ; set symmetric plot min/max
  res@cnConstFEnableFill   = True
  res@cnInfoLabelOn        = False
  res@cnLevelSelectionMode = "ManualLevels"  ; set manual contour levels
  res@cnMinLevelValF       = -.3             ; set min contour level
  res@cnMaxLevelValF       = .3              ; set max contour level
  res@cnLevelSpacingF      = .05             ; *special* match CPC

  resP                     = True         ; modify the panel plot
  resP@gsnMaximize         = True         ; large format
  resP@gsnPanelLabelBar    = True         ; add common colorbar
  resP@gsnPaperOrientation = "portrait"   ; force portraitorbar
  resP@lbLabelAutoStride   = True         ; auto stride on labels
  resP@lbBoxLinesOn        = True
  resP@lbPerimOn           = False
  resP@lbLabelFontHeightF  = 0.013
  resP@pmLabelBarWidthF    = 0.6
  resP@lbBoxLineColor      = "black"
  resP@lbBoxLineThicknessF = 0.5
  resP@gsnPanelCenter      = False
  resP@txString            = ""

; --------------------------------------------------------------------

 plot = new(3,graphic)

 res@gsnLeftStringFontHeightF = 0.02
 res@gsnRightStringFontHeightF = 0.02

 ; 3 =0.86605 / 7 =0.46 / 9 =0.2569
 res@gsnLeftString   = "(a)"
 res@gsnRightString  = "0.86605"
 plot(0) = gsn_csm_contour_map(wks,x1(29,3,:,:),res)

 res@gsnLeftString   = "(b)"
 res@gsnRightString  = "0.46"
 plot(1) = gsn_csm_contour_map(wks,x1(29,7,:,:),res)

 res@gsnLeftString   = "(c)"
 res@gsnRightString  = "0.2569"
 plot(2) = gsn_csm_contour_map(wks,x1(29,9,:,:),res)

; ---------------------------------------------------------------

  gsn_panel(wks,plot,(/1,3/),resP)     ; now draw as one plot

  delete(wks)
  system("convert -density 200x200 -trim "+filo+".eps "+filo+".jpg")
  system("rm -rf "+filo+".eps")


end


