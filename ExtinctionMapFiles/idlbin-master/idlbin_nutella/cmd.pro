pro cmd

;device, true=24
;device, decomposed=0
;colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'


device, filename = '/Users/jkrick/palomar/LFC/cmd.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color


;make a cmd of the SDSS stars
numostars = 342
sdssstar = replicate({star,umag:0D, gmag:0D,rmag:0D, imag:0D},numostars)
n = 0
openr, lun, "/Users/jkrick/palomar/LFC/SDSS/sdss.csv", /GET_LUN
WHILE (NOT eof(lun)) DO BEGIN
    READF, lun, run,rerun,camcol,field,obj,type,ra,dec,u,g,r,i,z,Err_u,Err_g,Err_r,Err_i,Err_z
    sdssstar[n] = {star, u,g, r, i}
    n = n + 1
ENDWHILE
close, lun
free_lun, lun

plot, sdssstar.gmag, sdssstar.gmag - sdssstar.rmag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="g-r", title = "SDSS CMD", psym = 2, xrange=[12,22], yrange=[-2,4]

plot, sdssstar.gmag, sdssstar.gmag - sdssstar.imag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="g-i", title = "SDSS CMD", psym = 2, xrange=[12,22], yrange=[-2,4]

plot, sdssstar.rmag, sdssstar.rmag - sdssstar.imag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="r magnitude", ytitle="r-i", title = "SDSS CMD", psym = 2, xrange=[12,22], yrange=[-2,4]

plot, sdssstar.gmag, sdssstar.umag - sdssstar.gmag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="u-g", title = "SDSS CMD", psym = 2, xrange=[12,22], yrange=[-2,4]

;make a cmd of the NEP stars

image=readfits("/Users/jkrick/palomar/LFC/coadd_g.fits",hd)
find, image, x, y,flux,sharp, round, 5.0,8.0,[-1.0,1.0],[0.2,1.0],/silent
gcntrd,image,x,y,xc,yc,9,/silent
aper,image,xc,yc,flux,eflux,sky,skyerr,1,[21.25],[20,50],[-200,2000],/nan,/exact,/flux,setsky=0.0,/silent
nepgmag=26.48 -2.5*alog10(flux)

xyad, hd, xc, yc, ra, dec


image=readfits("/Users/jkrick/palomar/LFC/coadd_r.fits",hd)
adxy, hd, ra, dec, x, y
;find, image, x, y,flux,sharp, round, 5.0,4.8,[-1.0,1.0],[0.2,1.0],/silent
gcntrd,image,x,y,xc,yc,6,/silent
aper,image,xc,yc,flux,eflux,sky,skyerr,1,[15],[30,50],[-200,20000],/nan,/exact,/flux,setsky=0.0,/silent
;neprmag=26.67 -2.5*alog10(flux)
neprmag=26.55 -2.5*alog10(flux)

image=readfits("/Users/jkrick/palomar/LFC/coadd_i.fits",hd)
adxy, hd, ra, dec, x, y
;find, image, x, y,flux,sharp, round, 5.0,4.8,[-1.0,1.0],[0.2,1.0],/silent
gcntrd,image,x,y,xc,yc,6,/silent
aper,image,xc,yc,flux,eflux,sky,skyerr,1,[15],[30,50],[-200,20000],/nan,/exact,/flux,setsky=0.0,/silent
;nepimag=26.35 -2.5*alog10(flux)
nepimag=26.0 -2.5*alog10(flux)

image=readfits("/Users/jkrick/palomar/LFC/coadd_u.fits",hd)
adxy, hd, ra, dec, x, y
;find, image, x, y,flux,sharp, round, 5.0,4.8,[-1.0,1.0],[0.2,1.0],/silent
gcntrd,image,x,y,xc,yc,6,/silent
aper,image,xc,yc,flux,eflux,sky,skyerr,1,[15],[30,50],[-200,20000],/nan,/exact,/flux,setsky=0.0,/silent
;nepimag=26.35 -2.5*alog10(flux)
nepumag=24.0-2.5*alog10(flux)

plot, nepgmag, nepgmag-neprmag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="g-r", title = "NEP CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
plot, nepgmag, nepgmag-nepimag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="g-i", title = "NEP CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
plot, neprmag, neprmag-nepimag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="r magnitude", ytitle="r-i", title = "NEP CMD", psym = 2, xrange=[12,22], yrange=[-2,4]
plot, nepgmag, nepumag-nepgmag, thick=3, xthick=3,ythick=3, charthick=3, $
      xtitle="g magnitude", ytitle="u-g", title = "NEP CMD", psym = 2, xrange=[12,22], yrange=[-2,4]

device, /close
set_plot, mydevice

undefine, sdssstar
end
