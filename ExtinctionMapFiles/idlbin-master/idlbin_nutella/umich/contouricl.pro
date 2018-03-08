

pro contouricl

;device, true=24
;device, decomposed=0

close, /all		;close all files = tidiness

colors = GetColor(/Load, Start=1)
TEK_COLOR
filename= strarr(1)            ;declare filename to be a stringarray
filename = "/n/Godiva2/jkrick/A141/original/color.icl"
imagefile = filename + '.fits'

; read in image
FITS_READ, imagefile, data, header
;newimage = congrid(data,507.5,510)  	;makes smaller for display on screen;
;tv, hist_equal(newimage)		;displays


openr, lun, "/n/Godiva2/jkrick/A141/original/testx",/get_lun
i = 0
xarr = fltarr(91)
WHILE ( NOT eof(lun)) DO BEGIN
    print, i
    readf, lun, x
    xarr[i] = x
    i = i + 1
ENDWHILE
close, lun
free_lun, lun

openr, lun, "/n/Godiva2/jkrick/A141/original/testy",/get_lun
i = 0
yarr = fltarr(91)
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, x
    yarr[i] = x
    i = i + 1
ENDWHILE
close, lun
free_lun, lun


set_plot, 'ps'
device, filename='/n/Godiva2/jkrick/A141/contour.ps', /portrait,$
                BITS=8, scale_factor=0.9 , /color
!P.thick = 3

;tv, hist_equal(data)		;displays

;contour, data, max_value=1.0, min_value = 0.0, levels = [0.2,0.4,0.6,0.8], c_color = [2,3,4,5], /isotropic
contour, data, xarr, yarr

device, /close

END

