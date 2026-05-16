'reinit'
*
*NDVI
'set rgb 50 139   0   0'
'set rgb 51 160  82  45'
'set rgb 52 244 164  96'
'set rgb 53 255 205   0'

'set rgb 54 149 193 107'
'set rgb 55 128 128   0'
'set rgb 56 127 255 212'
'set rgb 57  32 178 170'
'set rgb 58   0 139 139'
'set rgb 59   0 255   0'
'set rgb 60  60 179 113'
'set rgb 61   0 100   0'
*
*NDVI
'set rgb 79 139  0  0'
'set rgb 78 165 42 42'
'set rgb 77 178 34 34'
'set rgb 76 139 69 19'
'set rgb 75 160 82 45'
'set rgb 74 210 105 30'
'set rgb 73 205 133 63'
'set rgb 72 244 164 96'
'set rgb 71 210 180 140'
'set rgb 70 255 228 181'

'set rgb 80 149 200 107'
'set rgb 81 128 168 100'
'set rgb 82 127 255 212'
'set rgb 83  32 178 170'
'set rgb 84   0 139 139'
'set rgb 85   0 255   0'
'set rgb 86  60 179 113'
'set rgb 87   0 100   0'
**
* Blue to Dark Red
'set rgb 21   0   0  60'
'set rgb 22   0   0 120'
'set rgb 23   0   0 180'
'set rgb 24   0  30 255'
'set rgb 25  40  90 255'
'set rgb 26  80 120 255'
'set rgb 27 120 180 255'
'set rgb 28 160 210 255'
'set rgb 29 200 240 255'
*
'set rgb 40 255 240 200'
'set rgb 41 255 210 160'
'set rgb 42 255 180 120'
'set rgb 43 255 120  80'
'set rgb 44 255  90  40'
'set rgb 45 255  30   0'
'set rgb 46 180   0   0'
'set rgb 47 120   0   0'
'set rgb 48  60   0   0'
*
* Dark gray rgb
'set rgb 99 68 68 68'
*********************************************************** 
* Box map
*********************************************************** 
* 80-280E, 25-80N
'open /Users/dylee/Work-HYU/workstation/stu3/fig/intdecad/paper.v1/map.80-280E_25-80N.ctl'
'set t 1'
'set x 1 360'
'set y 1 181'
'boxnh=mm'
'close 1'

* Siberian region
'open /Users/dylee/Work-HYU/workstation/stu3/fig/intdecad/paper.v1/map.ctl'
'set t 1'
'set x 1 360'
'set y 1 181'
'box=mm'
'close 1'
***********************************************************
* Maskout
***********************************************************
* ERA5-360x181
'sdfopen /Users/dylee/Work-HYU/workstation/stu1/mask/maskera.0-360.360x181.nc'
'set t 1'
'set x 1 360'
'set y 1 181'
'mask=maskera'
'close 1'
***********************************************************
* Z500 JJA P2(2008-2023)-P1(1982-2007)
***********************************************************
'sdfopen /Volumes/PNU_DATA/workstation/RAWDATA_ERA5/ERA5/Z500/mhap/8223/z500.mon.8223.JJA.seamean.time.nc'
'set x 1 360'
'set y 1 181'
'p1 = ave(z500,t=1,t=26)'
'p2 = ave(z500,t=27,t=42)'
'df = p2-p1'
'dfed = (p2-p1)-ave(p2-p1,x=1,x=360)'
*
'clm = ave(z500,t=1,t=42)'
'std = sqrt(ave(pow(z500-clm,2),t=1,t=42))'

'tp1 = sqrt(26)*(p1)/std'
'tp2 = sqrt(16)*(p2)/std'
'tdf = sqrt(42)*(df)/std'
'tdfed = sqrt(42)*(dfed)/std'
'close 1'
*
***********************************************************
* SWM (Stationary Wave Model,360x181) sf Cont 
***********************************************************
'sdfopen /Users/dylee/Work-HYU/workstation/stu3/SWM/SWM_dylee/output/heat_100-160E_60-70N_era5_8223/strm.heat.360x181.nc'
'set x 1 360'
'set y 1 181'
'set z 1'
'set z 8'
'sf1 = ave(strm/1e6,t=30,t=59)'
****************************************
* PCC
'set lat 25 80'
'set lon 80 280'
'set z 1'
'set t 1'
'kk = scorr(dfed,sf1,lon=80,lon=280,lat=25,lat=80)'
'd kk'
xx=subwrd(result,4)
say 'OBS_Z500-SWM_SF-PCC = 'xx
*
***********************************************************
* SWM (Stationary Wave Model) sf Cont
***********************************************************
*A. Z500
*P2-P1
'set vpage 0 8.5 0 11'
'set parea 1.35 6.35 6.2 10.1'
'set map 1 1 1'
'set mpdset mres'
'set grads off'
'set grid on'
'set mproj lambert'
'set lat 25 80'
'set lon 80 280'
'set xlint 20'
'set ylint 10'
'set xlopts 1 3 0.11'
'set ylopts 1 3 0.11'
'set t 1'
'set gxout shaded'
'set clevs     -16   -12   -8   -4    0    4    8   12   16' 
'set ccols   24    26    28   29   0    0    40   42   43   45   46   47   48'
'd dfed'
**************************************** Box squre
'set gxout contour'
'set clevs     -16   -12   -8   -4    0    4    8   12   16' 
'set ccolor 1'
'set cthick 5'
'set clopts 1 1 0.08'
'set clab masked'
'd dfed'
*
'set gxout grid'
'set grid off'
'set gridln off'
'set digsiz 0.010'
'set ccolor 1'
'd skip(const(maskout(tdfed,abs(tdfed)-1.645),0),4,4)'
'/Volumes/PNU_DATA/workstation/glib/xcbar.gs 2.4 6.2 7.40 7.5 -d h -fs 1 -line on -fh 0.09 -fw 0.09'
*
**************************************** Box squre
'set gxout contour'
'set clevs 0.9999'
'set ccols  3'
'set cthick 6'
'set clab off'
'set grid off'
'd box'
*
'set gxout contour'
'set clevs 0.9999'
'set ccols  1'
'set cthick 4'
'set clab off'
'set grid off'
'd boxnh'
'close 1'
***********************************************************
* SWM Surface Heat Forcing
***********************************************************
* EA Box region
'open /Users/dylee/Work-HYU/workstation/stu3/fig/intdecad/paper.v1/map.80-180E_20-81N.ctl'
'set t 1'
'set x 1 360'
'set y 1 181'
'boxtot=mm'
'close 1'
*
'sdfopen /Users/dylee/Work-HYU/workstation/stu3/SWM/SWM_dylee_lon60/output/strm.heat.area.nc'
'set z 1'
'set t 1'
'ht=heat*86400*50'
***********************************************************
*B. Diabatic Heating Area
'set vpage 0 8.5 0 11'
'set parea 0.0 3.0 3.4 6.6'
'set mpdset mres'
'set grads off'
'set grid on'
'set mproj lambert'
'set lat 20 90'
'set lon 80 180'
'set xlint 20'
'set ylint 20'
'set gxout shaded'
'set xlopts 1 3 0.08'
'set ylopts 1 3 0.08'
'set clevs     0   0.5   1   1.5   2   2.5   3   3.5'
'set ccols   0  40    41  42    43  44    45  46    47'
'set t 1'
'd ht'
'/Volumes/PNU_DATA/workstation/glib/xcbar.gs 0.5 3.0 4.5 4.6 -d h -fs 1 -line on -fh 0.09 -fw 0.09'
'set gxout contour'
'set clevs     0   0.5   1   1.5   2   2.5   3   3.5'
'set ccolor 2'
'd ht'
***************************************** Box squre
'set gxout contour'
'set clevs 0.9999'
'set ccols  1'
'set cthick 4'
'set clab off'
'set grid off'
'd boxtot'
*
'close 1'
***********************************************************
* SWM Surface Heat Forcing Vertical Profile
***********************************************************
'sdfopen /Users/dylee/Work-HYU/workstation/stu3/SWM/SWM_dylee/output/heat_100-160E_60-70N_era5_8223//strm.heat.nc'
'set z 1 14'
'set x 1'
'set y 1'
'htz = aave(heat*86400,lon=100,lon=160,lat=60,lat=70)'
*
'set z 1 14'
'set x 1'
'set lat 50 80'
'htzm = ave(heat*86400,lon=100,lon=160)'
***********************************************************
*C. Vertical Profile
'set vpage 0 8.5 0 11'
'set parea 3.5 5.5 4.5 6.6'
'set map 1 1 1'
'set mproj latlon'
'set grads off'
'set grid on'
'set x 1'
'set y 1'
'set z 1 14'
'set gxout line'
'set ylopts 1 3 0.08'
'set xlint 0.2'
'set ylint 0.1'
'set cmark 2'
'set ccolor 1'
'set digsiz 0.04'
'd htz'
'set strsiz 0.10'
'set string 1 tc 3 90'
'draw string 3.08  5.7 Sigma levels'
************************************
*D. Vertical & Merid. Structure
'set vpage 0 8.5 0 11'
'set parea 5.8 7.8 4.5 6.6'
'set map 1 1 1'
'set mproj latlon'
'set grads off'
'set x 1'
'set lat 50 80'
'set xlint 5'
'set z 1 14'
*
'set gxout shaded'
'set clevs    0.05   0.3   0.6   0.9   1.2   1.5   1.8   2.1'
'set ccols   0    40    41    42    43    44    45    46    47'
'd htzm'
'/Volumes/PNU_DATA/workstation/glib/xcbar.gs 7.9 8.0 4.5 6.6 -d v -fs 1 -line on -fh 0.085 -fw 0.085'

'set gxout contour'
'set ylopts 1 3 0.08'
'set clevs    0.05   0.3   0.6   0.9   1.2   1.5   1.8   2.1'
'set ccolor 1'
'set grid off'
'set ylab off'
'set xlab off'
'set frame off'
'd htzm'
*
'close 1'
***********************************************************
* SWM (Stationary Wave Model) sf Cont
***********************************************************
'sdfopen /Users/dylee/Work-HYU/workstation/stu3/SWM/SWM_dylee/output/heat_100-160E_60-70N_era5_8223/strm.heat.nc'
'set x 1 96'
'set y 1 80'
'set z 1'
'ht = heat'
'set z 8'
'sf = ave(strm/1e6,t=30,t=59)'
'stdswm = sqrt(ave(pow(strm/1e6-sf,2),t=30,t=59))'
'tvswm = sqrt(30)*(sf)/stdswm'
******************************************
*E. SWM SF
'set vpage 0 8.5 0 11'
'set parea 1.35 6.35 0.0 3.9'
'set map 1 1 1'
'set mpdset mres'
'set grads off'
'set grid on'
'set mproj lambert'
'set lat 25 80'
'set lon 80 280'
'set xlint 20'
'set ylint 10'
'set xlopts 1 3 0.11'
'set ylopts 1 3 0.11'
'set t 1'
'set gxout shaded'
'set clevs     -1.2  -0.8  -0.6  -0.4  -0.2  -0.1   0   0.1  0.2  0.4  0.6  0.8  1.2  1.6'
'set ccols   24    25    26    27    28    29     0   0    40   41   42   43   44   45   46   47   48'
'd sf'
*
'set gxout contour'
'set clevs     -1.2  -0.8  -0.6  -0.4  -0.2  -0.1   0   0.1  0.2  0.4  0.6  0.8  1.2  1.6'
'set ccolor 1'
'set cthick 5'
'set clopts 1 1 0.08'
'set clab masked'
'd sf'
**************************************** Box squre
'set gxout contour'
'set clevs 0.9999'
'set ccols  3'
'set cthick 6'
'set clab off'
'set grid off'
'd box'
*
'set gxout contour'
'set clevs 0.9999'
'set ccols  1'
'set cthick 4'
'set clab off'
'set grid off'
'd boxnh'
****************************************
'/Volumes/PNU_DATA/workstation/glib/xcbar.gs 2.4 6.2 1.2 1.3 -d h -fs 2 -line on -fh 0.09 -fw 0.09'
*
******************************************
'set strsiz 0.11'
'set string 1 tl 3 0'
'draw string 2.25  10.2 (a) ERA5, Z500'
'draw string 0.5   6.8 (b) Diabatic heating area'
'draw string 3.5   6.8 (c) Vertical profile'
'draw string 5.7   6.8 (d) Vert. & Merid. Structure'
'draw string 2.25   4.0 (e) SWM, SF (at 0.46 sigma level)'
*******************************************
'gxprint fig.5.org.pdf white'
'!convert -quality 100 -density 300 -trim fig.5.org.pdf fig.5.org.png'
******************************************
function arrow(x,y,len,scale)
'set line 1 1 4'
'draw line 'x-len/2.' 'y' 'x+len/2.' 'y
'draw line 'x+len/2.-0.05' 'y+0.025' 'x+len/2.' 'y
'draw line 'x+len/2.-0.05' 'y-0.025' 'x+len/2.' 'y
'set string 1 c'
'set strsiz 0.1'
'draw string 'x' 'y-0.1' 'scale
return
