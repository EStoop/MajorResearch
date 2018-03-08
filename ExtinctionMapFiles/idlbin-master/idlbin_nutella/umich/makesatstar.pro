PRO make_star

close, /all

maxradius = 1000.0

;A3880
;xcenter =1757.9
;ycenter=2505.7
;xcenter = 1774.3
;ycenter = 2607.8

;A118
xcenter=1224.6
ycenter =953.8

intx = round(xcenter)
inty = round(ycenter)

;dx = xcenter-intx
;dy = ycenter-inty

radius = 0.0		

;xmax = 3460
;ymax = 4590
xmax=2300
ymax=2520
FITS_READ, '/n/Godiva7/jkrick/A118/original/larger.fits', blankdata, blankheader
blankdata[*,*] = 0.0

openr, lun1, "/n/Godiva3/jkrick/satstar/combined.prof", /get_lun
rad = fltarr(10000)
counts = fltarr(10000)
i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    readf, lun1, r, c
    rad[i] = r
    counts[i] = c
    i = i +1
ENDWHILE
rad = rad[0:i-1]
counts = counts[0:i-1]
close, lun1
free_lun, lun1

print, "rad[190]", rad[190]

FOR y = inty-maxradius, inty + maxradius,1 DO BEGIN
    FOR x = intx -maxradius, intx + maxradius,1 DO BEGIN
        ;find what the radius is from the center
        radius = sqrt((intx - x)^2 + (inty - y)^2)
        
        IF radius NE 0 AND radius LT 999 THEN BEGIN
;            intensitysb = 3.9 + 2.5*alog10((38.6/radius)^(-3.)) + sb4
;           intensitysb = 6.9 + 2.5*alog10((38.6/radius)^(-3.)) + sb2
;            intensity = (0.259^2.)*(10.^((24.6-intensitysb)/(2.5)))
            intensity = counts[fix(radius)]
        IF (x LT xmax AND y LT ymax AND x GT 1 AND y GT 1 ) THEN BEGIN
            blankdata[x,y] = blankdata[x,y] + intensity
        ENDIF;

        ;make sure not off the edge of the chip
        ;need to deal with going from fractional pixels to whole pixels
        ENDIF
    ENDFOR
ENDFOR


fits_write, '/n/Godiva7/jkrick/A118/original/star.r.fits',blankdata, blankheader
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
