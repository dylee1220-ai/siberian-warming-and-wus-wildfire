; ==============================================================
; Create 3-D (x,y,z) U,V,T and 2-D surface pressure file
; Output is Binary
; For MERRA 3D dataset
; ==============================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

; ==============================================================

 yrStrt = 1980
 yrLast = 2017
 ymstrt = 198001
 ymlast = 201712

; ==============================================================
; 19 Level

 lev = (/ 1000., 950., 900., 850., 800., 750., 700., 600., 500., 400., 300., 200., 100., 50., 30., 10., 5., 3., 1. /) 
 dims =  dimsizes(lev)

; ==============================================================
; Open File
; ==============================================================

  fpath  = systemfunc("ls /data4/MERRA2/instM_3d_asm_Np/*.nc")
  f      = addfiles (fpath, "r")
  uu = f[:]->U
  vv = f[:]->V
  tt = f[:]->T
  pp = f[:]->PS

  yyyymm_tem = yyyymm_time(yrStrt, yrLast,"integer")
  yyyymm = yyyymm_tem({ymstrt:ymlast})
  uu&time = yyyymm
  vv&time = yyyymm
  tt&time = yyyymm
  pp&time = yyyymm

  uu@_FillValue = 1e+15
  vv@_FillValue = 1e+15
  tt@_FillValue = 1e+15
  pp@_FillValue = 1e+15


print("Make Mean State")
;================================================================
; Make Mean State
;================================================================

  u = dim_avg_n_Wrap(month_to_season(uu(:,{lev},:,:), "JJA"),0)
  v = dim_avg_n_Wrap(month_to_season(vv(:,{lev},:,:), "JJA"),0)
  t = dim_avg_n_Wrap(month_to_season(tt(:,{lev},:,:), "JJA"),0)
  ps = dim_avg_n_Wrap(month_to_season(pp(:,:,:), "JJA"),0)

;==================================================================
; OUTPUT
;==================================================================
 
 path = "./input/"

system("rm -rf "+path+"U-clm.bdat")
system("rm -rf "+path+"V-clm.bdat")
system("rm -rf "+path+"T-clm.bdat")
system("rm -rf "+path+"PS-clm.bdat")

 fbindirwrite(path+"U-clm.bdat",u)
 fbindirwrite(path+"V-clm.bdat",v)
 fbindirwrite(path+"T-clm.bdat",t)
 fbindirwrite(path+"PS-clm.bdat",ps)

end

