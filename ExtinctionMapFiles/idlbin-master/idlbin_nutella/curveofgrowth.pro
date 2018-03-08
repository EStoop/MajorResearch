pro curveofgrowth

;taken from mags.pro

;device, true=24
;device, decomposed=0
;colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'


device, filename = '/Users/jkrick/mmt/apertures.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

xc=[3124.13,1707.92]
yc=[2941.94,3600.71] 


apsize=[1,2,2.5,3,3.5,4,4.5,5,5.5,6,6.5,7]


fwhm=2.5
apsize = apsize*fwhm

image=readfits("/Users/jkrick/mmt/IRACCF.z.coadd.fits",hd)
aper,image,xc,yc,flux,eflux,sky,skyerr,1,apsize,[20,50],[-200,200000],/nan,/exact,/flux,setsky=0.0,/silent
nepgmag=31.5 -2.5*alog10(flux)

print, flux

plot, apsize/fwhm, flux/flux(11), psym=2, yrange=[0.7,1.1], ystyle = 1,$
      xthick=3,ythick=3,charthick=3,xtitle="x fwhm", ytitle="fraction of total flux",$
      title="where to measure magnitudes, z band", thick=3


device, /close
set_plot, mydevice


end
