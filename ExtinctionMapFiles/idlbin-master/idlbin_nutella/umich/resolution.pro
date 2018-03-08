PRO resolution
device, true=24
device, decomposed=0
;colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva7/jkrick/postdoc/resolution.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

zarr = findgen(30)/10 
print, zarr
larr = fltarr(n_elements(zarr))
parr = fltarr(n_elements(zarr))

rad = !PI/180.
print, rad
plot, zarr, lumdist(zarr)*tan(rad)*1000/3600, xtitle = "redshift", ytitle = "resoultion in kpc", linestyle = 2, xthick = 3, ythick = 3,charthick = 4, thick = 3
oplot, zarr, lumdist(zarr)*tan(rad)*1000/3600 *0.1, thick = 3

xarr = findgen(10)/ 10 + 2.8
yarr = fltarr(10) + 15
print, xarr
oplot, xarr, yarr, thick = 3
xyouts, 2.65, 15, "MW", charthick = 3

xarr = findgen(10)/ 10 + 2.8
yarr = fltarr(10) + 40
print, xarr
oplot, xarr, yarr, thick = 3
xyouts, 2.6, 40, "BCG", charthick = 3

device, /close
set_plot, mydevice

END
