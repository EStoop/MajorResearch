PRO density

close,/all
fits_read,"/n/Godiva2/jkrick/A141/original/larger.fits", data, header

data[*,*] = 0

openr, lun, "/n/Godiva2/jkrick/A141/SExtractor.r.cat", /get_lun
WHILE (NOT eof(lun)) DO BEGIN
    readf, lun, junk, xcenter, ycenter, a, b, e, fluxbest, mbest, $
          isoarea,fwhm,theta,back,fluxaper,fluxisocor, fluxiso

    IF fluxbest GT 0 AND isoarea GT 6 AND fwhm GT 4 AND e GT 0.2 THEN BEGIN
;        print, xcenter,ycenter
        data[xcenter-1, ycenter-1] = 1
    ENDIF
ENDWHILE

fits_write, "/n/Godiva2/jkrick/A141/original/galnum.fits", data, header

close, lun
free_lun, lun

END
