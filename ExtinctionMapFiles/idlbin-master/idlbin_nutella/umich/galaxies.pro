
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
maxradius = 300.0                ; maximum radius in pixels to create a star
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
maxfit = 12

;fit profiles to the ellipse models
OPENR, lun1, "/n/Truffle1/jkrick/A3888/ellipfits", /GET_LUN
openr, lun2, '/n/Truffle1/jkrick/A3888/gallist', /GET_LUN

WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, filename
    print, filename
    READF, lun2, o, xcenter, ycenter, a, b, e, f, m, j, isoarea, fwhm, pa, bkgd
    print, blankdata[824, 1084], "blankdata"
    print, blankdata2[824, 1084], "blankdata2"
    whatever = myprof(xcenter, ycenter, data, blankdata, maxradius, xmax,ymax,number,maxfit)
    blah = ellipprof(filename,number, xcenter, ycenter, data2,blankdata2, maxradius)
    number = number+1
    print, blankdata[824, 1084], "blankdata"
    print, blankdata2[824, 1084], "blankdata2"
ENDWHILE


;close the SExtractor file, and free up the lun so it can be used again
close, lun
free_lun, lun

close, lun2
free_lun, lun2

;write a new fits file with the masked pixels, and one with the stars
fits_write, '/n/Truffle1/jkrick/A3888/subgalsmyprof.fits',shift(blankdata,-1,-1), blankheader
fits_write, '/n/Truffle1/jkrick/A3888/subgalsellip.fits',shift(blankdata2,-1,-1), blankheader


device, /close
set_plot, mydevice


end

