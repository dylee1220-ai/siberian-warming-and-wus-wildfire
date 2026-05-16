; ==============================================================
; Convert SWM output to nc file
; Variable is stream function
; ==============================================================

load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

; ==============================================================
; Open Output
; ==============================================================

; Forcing
; forc = "tran"   ; Vorticity Forcing
 forc = "heat"   ; Heat Forcing

 path = "./output/"
 data = fbindirread(path+"strm."+forc+".tmp",0,(/59,14,80,96/),"float")
 data = data(:,::-1,:,:) ; reverse levels

; heat = new((/1,14,80,96/),"float")
 heat = fbindirread("./temp/heat.forcing",0,(/14,80,96/),"float")
 heat = heat(::-1,:,:)
; ==============================================================
; Make Coordination
; ==============================================================
; Time - not real time but for GrADS
 time = ispan(1,59,1)
 time!0 = "time"
 time@units = "hours since 2000-01-01 00:00:00"
 time@long_name = "Time"
 time@description = "Steps"

; Level
 levels = (/.9966499, .9786499, .9352999, .86605, .77695, .6755, .5681999, .46, .3549, .2569, .17065, .1008999, 5.03500E-02, 1.50000E-02 /)
 levels!0 = "levels"
 levels@units = "hPa"
 levels@long_name = "Levels"
 levels@description = "Sigma level"
 levels@type = "plev"
 levels@positive = "down"

; Latitude
 latitude = fspan(-88.875,88.875,80)
 latitude!0 = "latitude"
 latitude@units = "degrees_north"
 latitude@long_name = "Latitude"

; Longitude
 longitude = fspan(0.0,360,96)
 longitude!0 = "longitude"
 longitude@units = "degrees_east"
 longitude@long_name = "Longitude"

 data!0 = "time"
 data!1 = "levels"
 data!2 = "latitude"
 data!3 = "longitude"
 data&time = time
 data&levels = levels
 data&latitude = latitude
 data&longitude =longitude

; heat!0 = "t"
 heat!0 = "levels"
 heat!1 = "latitude"
 heat!2 = "longitude"
; heat&time = time(0)
 heat&levels = levels
 heat&latitude = latitude
 heat&longitude = longitude
 
printVarSummary(data)

 system("rm -rf "+path+"strm."+forc+".nc")
 fout = addfile(path+"strm."+forc+".nc","c")
 fout->strm = data
 fout->heat = heat

end
