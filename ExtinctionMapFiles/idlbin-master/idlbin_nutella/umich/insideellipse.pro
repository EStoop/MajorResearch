
pro insideellipse, file, number,xcenter, ycenter
COMMON sharegal, inten, ierror, length
xcenter  = xcenter -1
ycenter = ycenter - 1
maxradius = 200.


;device, true=24
;device, decomposed=0

colors = GetColor(/Load, Start=1)
close, /all
radius = fltarr(100)
inten = fltarr(100)
ierror = fltarr(100)
ellip= fltarr(100)
elliperr= fltarr(100)
pos = fltarr(100)
poserr = fltarr(100)
rad = 0.D
i = 0.D
ierr = (0.0)
el = (0.0)
eerr = (0.0)
pa = (0.0)
paerr = (0.0)
junk = " "
length = 0.

;must get rid of all "indef" in the tab file
OPENR, lun,file, /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "file did not open"


WHILE (NOT EOF(lun)) DO BEGIN
;    READF, lun, rad, i, ierr, el, eerr, pa, paerr
    READF, lun, rad, i, ierr, el, pa
    
    radius(length) = rad
    inten(length) = i
    ierror(length) = ierr
    ellip(length) = el
;    elliperr(length) = eerr
    pos(length) = pa
;    poserr(length) = paerr

    length = length +1
ENDWHILE

;shorten arrays
radius = radius(0:length-1)
inten = inten(0:length-1)
ierror = ierror(0:length-1)
ellip = ellip(0:length-1)
;elliperr = elliperr(0:length-1)
pos = pos(0:length-1)
;poserr = poserr(0:length-1)

;take PA and e to be the fourth from the last one in the arrays
PA = pos(length - 4)
ellipticity = ellip(length - 4)

print, "PA and e for file", file,  " is", PA, ellipticity

;print, "common share variables", inten, ierr, length

;-----------------------------------------------------------------------------
;curve fitting
;------------------------------------------------------------------------------
minfit = 25
maxfit = 36

IF (length-1 LT maxfit) THEN maxfit = length-1

print, file, " maxfit", maxfit, " minfit", minfit
err = dindgen(length) - dindgen(length) + 1 
start = [0.2,7.0]
result = MPFITFUN('devauc',radius(minfit:maxfit),inten(minfit:maxfit), ierror(minfit:maxfit), start,  bestnorm = chi2)

openw, lun1, "/n/Truffle1/jkrick/A3888/aoutput", /GET_LUN, /append 
openw, lun2, "/n/Truffle1/jkrick/A3888/boutput", /GET_LUN, /append 
openw, lun3, "/n/Truffle1/jkrick/A3888/coutput", /GET_LUN, /append 
openw, lun4, "/n/Truffle1/jkrick/A3888/doutput", /GET_LUN, /append 
openw, lun5, "/n/Truffle1/jkrick/A3888/houtput", /GET_LUN, /append 
openw, lun6, "/n/Truffle1/jkrick/A3888/ioutput", /GET_LUN, /append 
openw, lun7, "/n/Truffle1/jkrick/A3888/eoutput", /GET_LUN, /append 
openw, lun8, "/n/Truffle1/jkrick/A3888/foutput", /GET_LUN, /append 
openw, lun9, "/n/Truffle1/jkrick/A3888/goutput", /GET_LUN, /append 
openw, lun10, "/n/Truffle1/jkrick/A3888/joutput", /GET_LUN, /append 
openw, lun11, "/n/Truffle1/jkrick/A3888/koutput", /GET_LUN, /append 
openw, lun12, "/n/Truffle1/jkrick/A3888/loutput", /GET_LUN, /append 

IF(file EQ "/n/Truffle1/jkrick/A3888/atab") THEN printf, lun1, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/btab") THEN printf, lun2, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/ctab") THEN printf, lun3, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/dtab") THEN printf, lun4, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/htab") THEN printf, lun5, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/itab") THEN printf, lun6, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/etab") THEN printf, lun7, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/ftab") THEN printf, lun8, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/gtab") THEN printf, lun9, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/jtab") THEN printf, lun10, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/ktab") THEN printf, lun11, file, chi2, result
IF(file EQ "/n/Truffle1/jkrick/A3888/ltab") THEN printf, lun12, file, chi2, result

;_____________________________________________________________________________________________________
;plot as a check on the profile, only plot last one of iterative cycle
;!p.multi = [number, 3, 3]

IF (number EQ 12) THEN BEGIN

    mydevice = !D.NAME
    SET_PLOT, 'ps'

    ending = string(".ps")
    plotname = file + ending
    device, filename=plotname, /landscape,$
      BITS=8,scale_factor=0.9 , /color

    ;!p.multi = [number, 3, 3]
    
    name = strarr(1)
    name = string(xcenter+1) + STRING(ycenter+1)

    x = findgen(30)
    plot, radius, inten, linestyle = 5, title = name, xrange = [0,30], yrange = [0,8],$
      xtitle = 'pixels', ytitle = 'counts/s'
    oplot, x, (result(0)) * (exp(-7.67*(((x/(result(1)))^(1.0/4.0)) - 1.0))), thick = 3, $
      color = colors.green
    name = string(chi2)+ string(result(0)) + string (result(1))
    xyouts, 10., 6., name, color = colors.green

    device, /close
    set_plot, mydevice

ENDIF

;-----------------------------------------------------------------------------------
;make a model of the galaxy and add it to a blank image
;return the galaxy model in array newgal
;-----------------------------------------------------------------------------------

FITS_READ, "/n/Truffle1/jkrick/A3888/centerforellipse.fits", data, header
FITS_READ, '/n/Truffle1/jkrick/A3888/blank.fits', blankdata, blankheader

newgal = subdev2(xcenter, ycenter, maxradius, blankdata, result, PA, ellipticity)
;newgal = trash3(xcenter, ycenter, length, blankdata, result, radius, ellip, pos)
;blankdata = blankdata + newgal  ;already gets done in trash3

data = data - newgal

FITS_WRITE, "/n/Truffle1/jkrick/A3888/centerforellipse.fits", data, header
FITS_WRITE, "/n/Truffle1/jkrick/A3888/ellipgal.fits", newgal, header

close ,/all



END


;read in first line, so I don't have to deal with INDEF's
;READF, lun, rad, i, junk

;radius(length) = rad
;inten(length) = i
;ierror(length) = 0.
;ellip(length) = 0.
;elliperr(length) = 0.
;pos(length) = 0.
;poserr(length) = 0.

;length = length +1

;hrot accepts pa measured clockwise from the x axis, my pa is measured
;ccw from y axis.  damn
;angle = 90 - PA
;IF (angle LT 0 ) THEN angle = angle + 360
;hrot, newgal, header, -1, -1, angle, xcenter, ycenter,  1, missing = 0, /pivot, errmsg=errmsg

;IF (file EQ "/n/Truffle1/jkrick/A3888/atab") THEN maxfit = 34
;IF (file EQ "/n/Truffle1/jkrick/A3888/ctab") THEN maxfit = 32
;IF (file EQ "/n/Truffle1/jkrick/A3888/etab") THEN maxfit = 33
;IF (file EQ "/n/Truffle1/jkrick/A3888/ftab") THEN maxfit = 33
;IF (file EQ "/n/Truffle1/jkrick/A3888/gtab") THEN maxfit = 33
;IF (file EQ "/n/Truffle1/jkrick/A3888/htab") THEN maxfit = 35
;IF (file EQ "/n/Truffle1/jkrick/A3888/itab") THEN maxfit = 34
;IF (file EQ "/n/Truffle1/jkrick/A3888/dtab") THEN BEGIN
;    minfit = 17
;    maxfit = 26
;ENDIF
