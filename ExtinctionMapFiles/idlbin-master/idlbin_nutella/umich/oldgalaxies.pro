
PRO galaxies

close, /all
;device, true=24
;device, decomposed=0


colors = GetColor(/Load, Start=1)
;declare variables
o = 0.0                         ;object number
xcenter = 0.0   		;xcenter
ycenter = 0.0			;ycenter
e = 0.0                         ;ellipticity
a = 0.0                         ;semi-major axis
b = 0.0                         ;semi-minor axis
m = 0.0                         ;magnitude (un-corrected for zpt, X, ...)
j = 0.0                         ;junk
f = 0.0                         ;flux
fwhm = 0.0                      ;full width at half maximum
isoarea = 0.0                   ;area inside of last isophote

filename = strarr(1)            ;declare filename to be a stringarray
maxradius = 100.0                ; maximum radius in pixels to create a star
intensity = 0.0
bkgd = 0.0

true = 0                        ;manual boolean variable

;size of the image
xmax = 2001
ymax = 2001

;read in a blank image to which I can add my own stars			
FITS_READ, '/n/Truffle1/jkrick/A3888/blank.fits', blankdata, blankheader
filename = '/n/Truffle1/jkrick/A3888/centerV.div'
imagefile = filename + '.fits'
datafile = filename + '.cat'
newimagefile = filename + '.s.fits'	

;run SExtractor on the image to find the objects
;commandline = 'lookforstars2 ' + filename 
;spawn, commandline

;read in the fits image with the galaxies that need to be subtracted
FITS_READ, imagefile, data, header

;make a second copy for ellipse fits
data2 = data
blankdata2 = blankdata

;open the SExtractor output datafile (ccd*.cat)
OPENR, lun, datafile, /GET_LUN

mydevice = !D.NAME
SET_PLOT, 'ps'

device, filename='/n/Truffle1/jkrick/A3888/dvprofile.ps', /landscape,$
  BITS=8,scale_factor=0.9 , /color

number = 1


;fit profiles to the ellipse models
;OPENR, lun1, "/n/Truffle1/jkrick/A3888/ellipfits", /GET_LUN
;WHILE (NOT EOF(lun1)) DO BEGIN
;    READF, lun1, filename
;print, filename
;    return= ellipprof(filename,number, xcenter, ycenter, data2,blankdata2, maxradius)
;number = number+1
;ENDWHILE

;device, /close
;set_plot, mydevice

;device, filename='/n/Truffle1/jkrick/A3888/dvprofile.ps', /portrait,$
;BITS=8,scale_factor=0.9 , /color


number = 1
maxfit = 15
WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, o, xcenter, ycenter, a, b, e, f, m, j, isoarea, fwhm, pa, bkgd
	
    ;fit and subtract the center 9 galaxies

    IF(xcenter GT 823 AND xcenter LT 826 AND ycenter GT 1083 AND ycenter LT 1086) THEN BEGIN
        whatever = myprof(xcenter, ycenter, data, blankdata, maxradius, xmax,ymax,number,maxfit)
;        number = number +1
        blah = ellipprof('/n/Truffle1/jkrick/A3888/btab',number, xcenter, ycenter, data2,blankdata2, maxradius)
        number = number +1

    ENDIF

    IF(xcenter GT 787 AND xcenter LT 790 AND ycenter GT 1107 AND ycenter LT 1109) THEN BEGIN
        whatever = myprof(xcenter, ycenter, data, blankdata, maxradius, xmax,ymax,number,maxfit)
;        number = number +1
        blah= ellipprof('/n/Truffle1/jkrick/A3888/atab',number, xcenter, ycenter, data2,blankdata2, maxradius)
        number = number +1
    ENDIF


    IF(xcenter GT 922 AND xcenter LT 925 AND ycenter GT 1148 AND ycenter LT 1151) THEN BEGIN
        whatever = myprof(xcenter, ycenter, data, blankdata, maxradius, xmax,ymax,number,maxfit)
;        number = number +1
        blah= ellipprof('/n/Truffle1/jkrick/A3888/ctab',number, xcenter, ycenter, data2,blankdata2, maxradius)
        number = number +1

    ENDIF

    IF(xcenter GT 921 AND xcenter LT 924 AND ycenter GT 1130 AND ycenter LT 1133) THEN BEGIN
        whatever = myprof(xcenter, ycenter, data, blankdata, maxradius, xmax,ymax,number,maxfit)
;        number = number +1
        blah= ellipprof('/n/Truffle1/jkrick/A3888/dtab',number,xcenter, ycenter, data2,blankdata2, maxradius)
        number = number +1

    ENDIF

    IF(xcenter GT 1059 AND xcenter LT 1062 AND ycenter GT 855 AND ycenter LT 858) THEN BEGIN
        whatever = myprof(xcenter, ycenter, data, blankdata, maxradius, xmax,ymax,number,maxfit)
;        number = number +1
        blah= ellipprof('/n/Truffle1/jkrick/A3888/itab',number,xcenter, ycenter, data2,blankdata2, maxradius)
        number = number +1

    ENDIF

    IF(xcenter GT 1058 AND xcenter LT 1061 AND ycenter GT 890 AND ycenter LT 893) THEN BEGIN
        whatever = myprof(xcenter, ycenter, data, blankdata, maxradius, xmax,ymax,number,maxfit)
;        number = number +1
        blah= ellipprof('/n/Truffle1/jkrick/A3888/htab',number,xcenter, ycenter, data2,blankdata2, maxradius)
        number = number +1

    ENDIF
    
    IF(xcenter GT 977 AND xcenter LT 980 AND ycenter GT 910 AND ycenter LT 913) THEN BEGIN
        whatever = myprof(xcenter, ycenter, data, blankdata, maxradius, xmax,ymax,number,maxfit)
;        number = number +1
        blah= ellipprof('/n/Truffle1/jkrick/A3888/etab',number,xcenter, ycenter, data2,blankdata2, maxradius)
        number = number +1

    ENDIF

    IF(xcenter GT 978 AND xcenter LT 981 AND ycenter GT 882 AND ycenter LT 885) THEN BEGIN
        maxfit = 4
        whatever = myprof(xcenter, ycenter, data, blankdata, maxradius, xmax,ymax,number,maxfit)
;        number = number +1
        blah= ellipprof('/n/Truffle1/jkrick/A3888/ftab',number,xcenter, ycenter, data2,blankdata2, maxradius)
        number = number +1
        maxfit = 12
    ENDIF

    IF(xcenter GT 1000 AND xcenter LT 1003 AND ycenter GT 883 AND ycenter LT 886) THEN BEGIN
        whatever = myprof(xcenter, ycenter, data, blankdata, maxradius, xmax,ymax,number,maxfit)
 ;       number = number +1
        blah= ellipprof('/n/Truffle1/jkrick/A3888/gtab',number,xcenter, ycenter, data2,blankdata2, maxradius)
        number = number +1

    ENDIF



ENDWHILE


;close the SExtractor file, and free up the lun so it can be used again
close, lun
free_lun, lun

;close, lun1
;free_lun, lun1

;write a new fits file with the masked pixels, and one with the stars
fits_write, '/n/Truffle1/jkrick/A3888/subgalsmyprof.fits',shift(blankdata,-1,-1), blankheader

fits_write, '/n/Truffle1/jkrick/A3888/subgalsellip.fits',shift(blankdata2,-1,-1), blankheader


device, /close
set_plot, mydevice


end

