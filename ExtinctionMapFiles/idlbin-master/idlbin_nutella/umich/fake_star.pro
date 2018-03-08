PRO fake_star

;this function acutally makes the homemade star and adds it to a blank
;image

close, /all
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'
colors = GetColor(/load, Start=1)

device, filename = '/n/Godiva1/jkrick/satstar/new/satstartest.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

f = 0.
xcenter = 1064;662.393 
ycenter = 1009;783.203
maxradius =900. 
xmax = 2040
ymax = 2040
fits_read, "/n/Godiva1/jkrick/satstar/new/ccd3126.m.fits", blankdata, header
blankdata[*,*] = 0.0

sb4=11.6;12.44;18.88;12.48;8.9

radius = 0.0		
intx = fix(xcenter)
inty = fix(ycenter)
dx = xcenter-intx
dy = ycenter-inty
blank1 = blankdata
blank1(*,*) = 0.0

FOR y = inty-maxradius, inty + maxradius,1 DO BEGIN
    FOR x = intx -maxradius, intx + maxradius,1 DO BEGIN
        ;find what the radius is from the center
        radius = sqrt((intx - x)^2 + (inty - y)^2)
        
        IF (radius NE 0) THEN begin
            intensitysb = 3.9 + 2.5*alog10((38.6/radius)^(-3.)) + sb4
;            intensitysb = -1.6 + 2.5*alog10((38.6/radius)^(-3.)) + 17.66
          ;  intensity = (10.^(-0.4*(intensitysb-24.3)))*(0.259)^2.
            intensity = (0.259^2.)*(10.^((24.6-intensitysb)/(2.5)))
        ;make sure not off the edge of the chip
        ;need to deal with going from fractional pixels to whole pixels
            IF (x LT xmax AND y LT ymax AND x GT 1 AND y GT 1 ) THEN BEGIN
                blank1[x,y] = blank1[x,y] + intensity
           ENDIF
        ENDIF
    ENDFOR
ENDFOR

blank1 = shiftf(blank1, dx,dy)
blankdata = blankdata + blank1
       

fits_write, "/n/Godiva1/jkrick/satstar/new/satstartest.fits", blankdata, header   

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;now make a radial profile of the fake star, and plot it
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

datafile = "/n/Godiva1/jkrick/satstar/new/satstartest.prof"
OPENW, lun, datafile, /GET_LUN

xcenter = ulong(1064)	;493;1064	; x coord of center of the star
ycenter = ulong(1009)	;460;1009	; y coord of center of the star

x = ulong(0)			; x pixel coord
y = ulong(0)			; y pixel coord
r = float(0.0)			; radius from center in pixels
p = float(0.0)			; count at a particular radius
max = 3000.0


radius = FLTARR((4*max)^2)	; array of radius (float)
counts = FLTARR((4*max)^2)	; array of counts (float)

; make an array of radii from the central pixel, and counts for the 
; corresponding radii

j = ulong64(0)					;counter

FOR x = 0, 2039, 1 DO BEGIN
    FOR y = 0, 2039, 1 DO BEGIN
        s = abs(xcenter- x)
        t = abs(ycenter - y)
        r = float(sqrt(s^(2)+(t^2))) ;pythagoras
        p = blankdata[x,y] ;-1.527              ;background subtraction
        IF (p GT -1.0 AND P LT 33000) THEN BEGIN ;block the bad data pts
            counts(j) = p 
            radius(j) = r
            j= j + 1
        ENDIF
    ENDFOR 
ENDFOR

;shorten the arrays to be only as long as need be
;radius = newradius(0:z-1)
;counts = newcounts(0:z-1)
radius = radius(0:j-1)
counts = counts(0:j-1)


;sort the radius array, and then apply that sorting order to the counts
sortindex = Sort(radius)
sortedradius = radius[sortindex]
sortedcounts = counts[sortindex]

;average the data into bins
i = ulong64(0)
sum = 0.0
counter = 1.
bin = 0.
binnedrad = fltarr(j-1)
binnedcounts = fltarr(j-1)
WHILE (i LT j-1) DO BEGIN

    IF sortedradius(i) LT (bin + 1) THEN BEGIN
        sum = sum + sortedcounts(i) 
        counter = counter + 1.
        i = i + 1
;        print, "sortedcounts, counter", sortedcounts(i), counter
    ENDIF ELSE BEGIN
        printf, lun,  bin + 0.5, sum/counter
        binnedrad(bin - 0.5) = bin+0.5
        binnedcounts(bin-0.5) = sum/counter
        counter = 1.
        sum = 0.
;        i = i -1
        bin = bin + 1.

    ENDELSE
    IF sortedradius(i) GT bin + 2 THEN BEGIN
        ;have to skip one or more bins
        bin = fix(sortedradius(i))
    ENDIF

ENDWHILE

close, lun
free_lun, lun

binnedrad = binnedrad[0:sortedradius(i-1)]
binnedcounts = binnedcounts[0:sortedradius(i-1)]
plot, alog10(binnedrad*(0.259)), alog10(binnedcounts), thick = 3, color = colors.blue,yrange = [-2,6]


;;;;;;;;;;;;;;;;;;;;;;;
;create profile of real satstar and plot it
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

openw, outlunreal, "/n/Godiva1/jkrick/satstar/new/ccd3126.m.prof2", /get_lun
fits_read, "/n/Godiva1/jkrick/satstar/new/ccd3126.m.fits",satdata, satheader

;set bad column = -100
satdata[1008:1019,640:2039] = -100
satdata[1005:1015,640:661] = -100
satdata[1:10,1:2039] = -100
satdata[2030:2039,1:2039] = -100
satdata[1:2039,1:10] = -100
satdata[1:2039,2030:2039] = -100

;set diffraction spikes = -100
;for ccd3126
satdata[1034:1086,1009:1522] = -1000
satdata[1038:1098,496:1009] = -1000
satdata[1064:1577,990:1046] = -1000
satdata[551:1064,970:1020] = -1000


xcenter = ulong(1064)	;493;1064	; x coord of center of the star
ycenter= ulong(1009)	;460;1009	; y coord of center of the star

x = ulong(0)			; x pixel coord
y = ulong(0)			; y pixel coord
r = float(0.0)			; radius from center in pixels
p = float(0.0)			; count at a particular radius
max = 3000.0


radius = FLTARR((4*max)^2)	; array of radius (float)
counts = FLTARR((4*max)^2)	; array of counts (float)

; make an array of radii from the central pixel, and counts for the 
; corresponding radii

j = ulong64(0)					;counter

FOR x = 0, 2039, 1 DO BEGIN
    FOR y = 0, 2039, 1 DO BEGIN
        s = abs(xcenter- x)
        t = abs(ycenter - y)
        r = float(sqrt(s^(2)+(t^2))) ;pythagoras
        p = satdata[x,y] -3.12
        IF (p GT -1.0 AND P LT 33000) THEN BEGIN ;block the bad data pts
            counts(j) = p 
            radius(j) = r
            j= j + 1
        ENDIF
    ENDFOR 
ENDFOR

;shorten the arrays to be only as long as need be
radius = radius(0:j-1)
counts = counts(0:j-1)


;sort the radius array, and then apply that sorting order to the counts
sortindex = Sort(radius)
sortedradius = radius[sortindex]
sortedcounts = counts[sortindex]

;average the data into bins
i = ulong64(0)
sum = 0.0
counter = 1.
bin = 0.
binnedrad = fltarr(j-1)
binnedcounts = fltarr(j-1)
WHILE (i LT j-1) DO BEGIN

    IF sortedradius(i) LT (bin + 1) THEN BEGIN
        sum = sum + sortedcounts(i) 
        counter = counter + 1.
        i = i + 1
;        IF bin EQ 100 then print, "sortedcounts, counter", sortedcounts(i), counter
    ENDIF ELSE BEGIN
        printf,outlunreal,  bin + 0.5, sum/counter
        binnedrad(bin - 0.5) = bin+0.5
        binnedcounts(bin-0.5) = sum/counter
        counter = 1.
        sum = 0.
;        i = i -1
        bin = bin + 1.

    ENDELSE
    IF sortedradius(i) GT bin + 2 THEN BEGIN
        ;have to skip one or more bins
        bin = fix(sortedradius(i))
    ENDIF

ENDWHILE

close, lun
free_lun, lun
close, outlunreal
free_lun, outlunreal
binnedrad = binnedrad[0:sortedradius(i-1)]
binnedcounts = binnedcounts[0:sortedradius(i-1)]
oplot, alog10(binnedrad*(0.259)), alog10(binnedcounts), thick = 3, color = colors.red

device, /close
set_plot, mydevice

END

;openr, reallun, "/n/Godiva1/jkrick/satstar/new/ccd3126.m.prof", /get_lun, error = err
;IF (err NE 0) then PRINT,  !ERROR_STATE.MSG

;a = 0
;radarr1 = fltarr(1480)
;countarr1 = fltarr(1480)
;rad = 0.0
;counts = 0.0
;WHILE (NOT EOF(reallun)) DO BEGIN
;    READF, reallun, rad,counts
;    radarr1(a) = rad
;    countarr1(a) = counts
;    a = a + 1
;ENDWHILE


;oplot, alog10(radarr1*0.259), alog10(countarr1), color = colors.red


;close, reallun
;free_lun, reallun
