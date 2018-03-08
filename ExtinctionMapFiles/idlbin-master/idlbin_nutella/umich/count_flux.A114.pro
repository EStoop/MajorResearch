PRO count_flux
close,/all

fits_read, "/n/Godiva4/jkrick/A114/original/larger.block.fits", data, header
inner = 0.
innersum = 0.D
outer = 0.
outersum = 0.
;sum = 0.D
npix = 0.D
npix2 = 0.D
xcenter1 = 1564.
ycenter1 = 543.

FOR x = 1200, 1900,1 DO BEGIN
    FOR y= 300, 800, 1 DO BEGIN
        dist1 = sqrt((xcenter1 - x)^2 + (ycenter1 - y)^2)
        IF (dist1 LE 190)  THEN begin
            npix = npix + 1
            IF (data[x,y] GT -.1) AND (data[x,y] NE 0) THEN BEGIN
                inner = inner + 1
                innersum = innersum + data[x,y]
            ENDIF
        ENDIF

    ENDFOR
ENDFOR
print,"r",inner, innersum, innersum/inner, 24.6 - 2.5*alog10((innersum/inner)/(0.259^2))
print, "should be npix: ", npix, inner/npix,"%"
print, "new sum = " , innersum/(inner/npix) ; + (innersum/inner)*(npix - innersum)



END
