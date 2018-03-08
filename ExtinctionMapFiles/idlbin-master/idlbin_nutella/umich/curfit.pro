;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;February 2002
;Jessica Krick
;
;This program takes a radial profile from a star (profile.pro output)
;	and fits it with a moffat function. MPFITFUN is used to do the
;	fitting, and moffat.pro is the actual moffat function, which can
;	be replaced with any function
;The moffat function is taken from iraf help file on imexam
;	can also be found in Trujillo 2001 MNRAS 328..977T
;
;input: data to fit
;	function to fit
;	error values for data, or weights
;	start values
;
;output: fit coefficients
;	fit values
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro curfit

;device, true=24
;device, decomposed=0
colors = GetColor(/Load, Start=1)
close, /all		;close all files = tidiness


;OPENR, lun,'/n/Whalen1/jkrick/ICL/satstar/21aug98/totalprof', /GET_LUN
OPENR, lun,'/n/Godiva1/jkrick/satstar/new/ccd3126.m.prof', /GET_LUN
;read in the radial profile
rows= 4001
rad = FLTARR(rows)
averagearr = FLTARR(rows)
r = 0.0				;radius
a = 0.0				;average counts
FOR j=0,rows-1 DO BEGIN
      READF, lun, r,a
      rad(j) = r
      averagearr(j) = a
ENDFOR
rad = rad[0:j-1]
averagearr = averagearr[0:j-1]
close, lun
free_lun, lun

;;make it a full radial profile (both sides of the curve)
;;and make sure it is in the correct order
average = FLTARR(2*rows)
radius = FLTARR(2*rows)
holder = FLTARR(rows)

j = 0.0
FOR i = rows - 1,0,-1 DO BEGIN
	holder(j) =- rad(i)
	j = j +1
ENDFOR
radius = [holder, rad ]

j = 0.0
FOR i = rows -1,0,-1 DO BEGIN
	holder(j) = averagearr(i)
	j = j +1
ENDFOR
average = [holder, averagearr]

average = average[0:j-1]
radius = radius[0:j-1]

err = dindgen(j-1) - dindgen(j-1) + 1

;this is the curve fitting part
;give it starting values, and then ask it to find the best fit
;start = [120.,4.2,1.0]
start = [250, 4.2]
result = MPFITFUN('moffat', radius, average, err, start)

print, result


;set up for plotting
mydevice = !D.NAME
!p.multi = [0, 0, 3]
SET_PLOT, 'ps'


device, filename = '/n/Godiva1/jkrick/satstar/new/moffatfit.ps', /portrait,$
                BITS=8, scale_factor=0.9 , /color
;,  xoffset=0, yoffset=0



;plot the data = black, and the fit = blue
plot, radius, average,XRANGE = [-200,200], ytitle = 'counts normalized to random number' $
	,title = 'Stellar psf profile with moffat function'

;oplot, radius, result(0)*(1 + ((radius) / result(1))^2)^(-result(2)), $
;	thick = 3,color = colors.blue

oplot, radius,result(2)*(1 + (2^(1/result(0))-1) * (radius/result(1))^2)^ (-result(0)) , $
	thick = 3,color = colors.blue


plot, radius, average, XRANGE = [-10,10],YRANGE = [95,125],ytitle = 'counts normalized to random number' $
	,title = 'Looking at the core'

;oplot, radius, result(0)*(1 + ((radius) / result(1))^2)^(-result(2)), $
;	thick = 3,color = colors.blue

oplot, radius,result(2)*(1 + (2^(1/result(0))-1) * (radius/result(1))^2)^ (-result(0)) , $
	thick = 3,color = colors.blue


plot, radius, average, XRANGE = [0,500],YRANGE = [0,0.001],xtitle = 'radius in pixels'$
	, ytitle = 'counts normalized to random number' $
	,title = 'Looking for the wings', subtitle = 'black = data, blue = fit'

;oplot, radius, result(0)*(1 + ((radius) / result(1))^2)^(-result(2)), $
;	thick = 3,color = colors.blue

oplot, radius,result(2)*(1 + (2^(1/result(0))-1) * (radius/result(1))^2)^ (-result(0)) , $
	thick = 3,color = colors.blue


device, /close
set_plot, mydevice
END



