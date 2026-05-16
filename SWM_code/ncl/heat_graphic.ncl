load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

; =========================================================
; OPEN FILE
; =========================================================

  WksType = "x11"
  filo = "../output/graph/heat_test"

  latS = -90
  latN = 90.
  lonL = 0.
  lonR = 360.

;  yrStrt = 1950
;  yrLast = 2018
;  ymstrt = 195001
;  ymlast = 201812
;  yr = yrLast-yrStrt+1

  f = addfile("../output/strm.heat.nc","r")
  x1 = f->strm
  x2 = f->heat
  x2 = x2*10^5

; =========================================================
; PLOT
; =========================================================

  wks = gsn_open_wks(WksType,filo)
;  gsn_define_colormap(wks,"BlueWhiteOrangeRed")
  gsn_define_colormap(wks,"sunshine_9lev")

  res                      = True
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@gsnMaximize          = True
  res@gsnSpreadColors      = True         ; spread out color table
  res@gsnMaximize          = True

  res@mpFillOn             = True        ; turn off map fill
;  res@gsnPolar             = "SH"
;  res@mpMaxLatF            = -40.
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
  res@cnMinLevelValF       = 0.             ; set min contour level
  res@cnMaxLevelValF       = 2.              ; set max contour level
  res@cnLevelSpacingF      = .1             ; *special* match CPC

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

 plot = new(14,graphic)

 res@gsnLeftString   = ""
 res@gsnRightString  = ""
 res@gsnLeftStringFontHeightF = 0.02
 res@gsnRightStringFontHeightF = 0.02

 do n = 0,13
   plot(n)  = gsn_csm_contour_map(wks,x2(n,:,:),res)
 end do

; ---------------------------------------------------------------

 gsn_panel(wks,plot,(/3,5/),resP)     ; now draw as one plot

  delete(wks)
  system("convert -density 200x200 -trim "+filo+".eps "+filo+".jpg")
  system("rm -rf "+filo+".eps")


end


