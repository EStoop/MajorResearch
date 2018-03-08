pro aper_corr

;device, true=24
;device, decomposed=0
;colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'


device, filename = '/Users/jkrick/nutella/palomar/LFC/aper_corr.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

xc=[4567.52, 5225.94]
yc=[5146.85,  5937.14] 

yo = fltarr(50) + 1
apsize=[1,1.2,1.4,1.6,1.8,2,2.2,2.4,2.6,2.8,3.0,3.2,3.4,3.6,3.8,4.0,4.2,4.4,4.6,4.8,5.0]
print, 'apsize', apsize(5), apsize(18)
;---------------------------------------------
image=readfits("/Users/jkrick/nutella/palomar/LFC/coadd_g.fits",hd)
fwhm=7.2
apsizeg = apsize*fwhm
aper,image,xc,yc,flux,eflux,sky,skyerr,1,apsizeg,[40,60],[-200,2000],/nan,/exact,/flux,setsky=0.0,/silent
nepgmag=26.58 -2.5*alog10(flux)
print, 'g', nepgmag(5), nepgmag(18), nepgmag(5) - nepgmag(18)

plot, apsizeg/fwhm, flux/flux(n_elements(apsize)-1), psym=2, yrange=[0.7,1.1], ystyle = 1,$
      xthick=3,ythick=3,charthick=3,xtitle="x fwhm", ytitle="fraction of total flux",$
      title="where to measure magnitudes, g band"
oplot, apsize, yo
xyad, hd, xc, yc, ra, dec

;---------------------------------------------

image=readfits("/Users/jkrick/nutella/palomar/LFC/coadd_r.fits",hd)
adxy, hd, ra, dec, x, y
;find, image, x, y,flux,sharp, round, 5.0,4.8,[-1.0,1.0],[0.2,1.0],/silent
fwhm=5.0
apsizer = apsize*fwhm
gcntrd,image,x,y,xc,yc,6,/silent
aper,image,xc,yc,flux,eflux,sky,skyerr,1,apsizer,[30,50],[-200,20000],/nan,/exact,/flux,setsky=0.0,/silent
neprmag=26.9 -2.5*alog10(flux)
print, 'r', neprmag(5), neprmag(18), neprmag(5) - neprmag(18)

plot, apsizer/fwhm, flux/flux(n_elements(apsize)-1), psym=2, yrange=[0.7,1.1], ystyle = 1,$
      xthick=3,ythick=3,charthick=3,xtitle="x fwhm", ytitle="fraction of total flux",$
      title="where to measure magnitudes, r band"
oplot, apsize, yo

;---------------------------------------------

image=readfits("/Users/jkrick/nutella/palomar/LFC/coadd_i.fits",hd)
adxy, hd, ra, dec, x, y
fwhm=6.7
apsizei = apsize*fwhm
;find, image, x, y,flux,sharp, round, 5.0,4.8,[-1.0,1.0],[0.2,1.0],/silent
gcntrd,image,x,y,xc,yc,6,/silent
aper,image,xc,yc,flux,eflux,sky,skyerr,1,apsizei,[35,55],[-200,20000],/nan,/exact,/flux,setsky=0.0,/silent
nepimag=26.35 -2.5*alog10(flux)
print, 'i', nepimag(5), nepimag(18), nepimag(5) - nepimag(18)

plot, apsizei/fwhm, flux/flux(n_elements(apsize)-1), psym=2, yrange=[0.7,1.1], ystyle = 1,$
      xthick=3,ythick=3,charthick=3,xtitle="x fwhm", ytitle="fraction of total flux",$
      title="where to measure magnitudes, i band"
oplot, apsize, yo

;---------------------------------------------

image=readfits("/Users/jkrick/nutella/palomar/LFC/coadd_u.fits",hd)
adxy, hd, ra, dec, x, y
fwhm=5.5
apsizeu = apsize*fwhm
;find, image, x, y,flux,sharp, round, 5.0,4.8,[-1.0,1.0],[0.2,1.0],/silent
gcntrd,image,x,y,xc,yc,6,/silent
aper,image,xc,yc,flux,eflux,sky,skyerr,1,apsizeu,[30,50],[-200,20000],/nan,/exact,/flux,setsky=0.0,/silent
nepumag=24.3 -2.5*alog10(flux)
print, 'u', nepumag(5), nepumag(18), nepumag(5) - nepumag(18)

plot, apsizeu/fwhm, flux/flux(n_elements(apsize)-1), psym=2, yrange=[0.7,1.1], ystyle = 1,$
      xthick=3,ythick=3,charthick=3,xtitle="x fwhm", ytitle="fraction of total flux",$
      title="where to measure magnitudes, u band"
oplot, apsize, yo

;---------------------------------------------

device, /close
set_plot, mydevice


end
