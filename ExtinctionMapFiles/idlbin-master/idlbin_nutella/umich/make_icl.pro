PRO make_icl
close, /all

;make an icl + noise background that has an exponential profile
;blankdata = noisedata
;blankdata[*] = 0
pa = 0.
ellipticity = 0.
maxradius = 400.    ;half size of image
;variables for exponential function
constants = fltarr(2)
constants = [0.01, 500.]
FOR  xcenter = 0 , 400, 50 DO begIN
    FOR ycenter = 0, 400, 50 DO BEGIN
fits_read, "/n/Godiva7/jkrick/galfit/noise400.fits", noisedata, noiseheader
;xcenter = 150
;ycenter = 150
       ;create the background image
        newback = subdev2(xcenter, ycenter, maxradius, noisedata, constants, pa, ellipticity)

        ;save the image
        fits_write, strcompress("/n/Godiva7/jkrick/galfit/newback" + string(xcenter) +"."$
                                + string(ycenter) + ".fits" ,/remove_all), newback, noiseheader
    ENDFOR
ENDFOR


END
