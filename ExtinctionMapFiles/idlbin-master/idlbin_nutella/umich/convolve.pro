PRO convolve

close, /all

P = fltarr(2)
P(0) = 0.117941
P(1) = 23.1693
devaucarr = fltarr(100)
moffatarr = fltarr(100)

;need to scale the moffat to the galaxy
;need to scale the moffat to the inner seeing disk of the galaxy, ie
;inner 5 pixels
;so pick a galaxy to start with, scale is going to equal f/17.166
;galaxy = a = 788,1108; flux = 84.5723
;oldscale = 519.423/19.83

;scale =6 ;84.5723 / 17.166
FOR x = 0,99,1 DO BEGIN
    devaucarr(x) = (P(0)) * (exp(-7.67*(((x/P(1))^(1.0/4.0)) - 1.0)))
    moffatarr(x) = scale * (1 + ((x) / 3.218)^2)^(-2.592)
ENDFOR
fits_read, "/n/Truffle1/jkrick/A3888/ellipgal.fits", galdata, galheader
fits_read, "/n/Truffle1/jkrick/A3888/npsf.fits", psfdata, psfheader


;this function acutally makes the homemade galaxy and adds it to a blank
;image
radius = 0.0		
a = 0.0                         ;semi-major axis 
;e = 0.0                         ;ellipticity
phi = 0.0                       ;angle measured from positive x axis
;PA = PA                         ;position angle



fd = fft(galdata, -1)
fm = fft(psfdata, -1)

prod = fd * fm

con = fft(prod, 1)

fits_write, "/n/Truffle1/jkrick/A3888/congal.fits", con, galheader

;get the galaxy profile from iraf/pradprof
;radius = fltarr(2824)
;inten = fltarr(2824)
;r =0.D
;i = 0.D
;j = 0
;openr, lun, '/n/Truffle1/jkrick/A3888/atab', /get_lun, ERROR = err
;IF (err NE 0) then PRINT, "file did not open"

;radius = fltarr(1000)
;inten = fltarr(1000)
;ierror = fltarr(1000)
;ellip= fltarr(1000)
;elliperr= fltarr(1000)
;pos = fltarr(1000)
;poserr = fltarr(1000)
;rad = 0.D
;i = 0.D
;ierr = (0.0)
;el = (0.0)
;eerr = (0.0)
;pa = (0.0)
;paerr = (0.0)
;junk = " "
;length = 0.



;WHILE (NOT EOF(lun)) DO BEGIN
;    READF, lun, rad, i, ierr, el, eerr, pa, paerr
    
;    radius(length) = rad
;    inten(length) = i
;    ierror(length) = ierr
;    ellip(length) = el
;    elliperr(length) = eerr
;    pos(length) = pa
;    poserr(length) = paerr

;    length = length +1
;ENDWHILE

;shorten arrays
;radius = radius(0:length-1)
;inten = inten(0:length-1)
;ierror = ierror(0:length-1)
;ellip = ellip(0:length-1)
;elliperr = elliperr(0:length-1)
;pos = pos(0:length-1)
;poserr = poserr(0:length-1)

;make an array of al the y values to put a vertical line at x = 5
;m = 0.D
;yarr = fltarr(1E6+2)
;FOR i = 1E-4, 1E4, 1 DO BEGIN
;    yarr(m) = i
;    m = m +1
;ENDFOR

colors = GetColor(/Load, Start=1)

mydevice = !D.NAME
SET_PLOT, 'ps'

;device, filename='/n/Truffle1/jkrick/A3888/trash.ps', /portrait,$
;  BITS=8,scale_factor=0.9 , /color
;plot, radius,inten, xrange = [0,30], yrange = [1E-4,1E4],/ylog, ytitle = "counts/s", xtitle = "radius (pixels)", title = "effect of convolving deV with a moffat"
;oplot, findgen(100), devaucarr, color = colors.blue
;oplot, findgen(100), moffatarr, color = colors.red
;oplot, findgen(100), con, color = colors.purple
;oplot, findgen(2E4) -findgen(2E5) + 5, yarr, color = colors.black
;xyouts, 10,1000, "deV fit to galaxy", color = colors.blue
;xyouts, 10, 500, "moffat fit scaled to flux of galaxy", color = colors.red
;xyouts, 10, 250, "convolution of moffat and DeV.", color = colors.purple
;xyouts, 10, 140, "galaxy profile data from ellipse", color = colors.black

device, /close
set_plot, mydevice

END
