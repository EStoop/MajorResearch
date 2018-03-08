PRO make_star

close, /all

maxradius = 1000.0

mag = 10.8
ftot = 10^((24.6 - mag)/2.5)


D = 0.044*ftot
sb4 = 24.6 - 2.5*alog10((D*mymoffat(15.4))/(0.259^2.))
sb2 = 24.6 - 2.5*alog10((D*mymoffat(7))/(0.259^2.))

xcenter = 1000.
ycenter = 1000.
intx = fix(xcenter)
inty = fix(ycenter)

dx = xcenter-intx
dy = ycenter-inty

radius = 0.0		

xmax = 2000
ymax = 2000

FITS_READ, '/n/Godiva1/jkrick/satstar/new/ccd3126.m.fits', blankdata, blankheader
blankdata[*,*] = 0.0

FOR y = inty-maxradius, inty + maxradius,1 DO BEGIN
    FOR x = intx -maxradius, intx + maxradius,1 DO BEGIN
        ;find what the radius is from the center
        radius = sqrt((intx - x)^2 + (inty - y)^2)
        
;        intensity = D*mymoffat(radius)
        IF radius NE 0 THEN BEGIN
;            intensitysb = 3.9 + 2.5*alog10((38.6/radius)^(-3.)) + sb4
            intensitysb = 6.9 + 2.5*alog10((38.6/radius)^(-3.)) + sb2
            intensity = (0.259^2.)*(10.^((24.6-intensitysb)/(2.5)))
        ENDIF;

        ;make sure not off the edge of the chip
        ;need to deal with going from fractional pixels to whole pixels
        IF (x LT xmax AND y LT ymax AND x GT 1 AND y GT 1 ) THEN BEGIN
            blankdata[x,y] = blankdata[x,y] + intensity
        ENDIF
    ENDFOR
ENDFOR


fits_write, '/n/Godiva2/jkrick/A3984/star.fits',blankdata, blankheader

END



;FOR y = ycenter-maxradius, ycenter + maxradius,1 DO BEGIN
;    FOR x = xcenter -maxradius, xcenter + maxradius,1 DO BEGIN
;        ;find what the radius is from the center
;        radius = sqrt((xcenter - x)^2 + (ycenter - y)^2);;;;
;
;        intensity = D*mymoffat(radius)
;      
;        ;make sure not off the edge of the chip
;        IF (x LT xmax AND y LT ymax AND x GT 1 AND y GT 1 ) THEN BEGIN
;            blankdata[x,y] = blankdata[x,y] + intensity
;        ENDIF
;    ENDFOR
;ENDFOR
