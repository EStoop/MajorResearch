;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;February 2002
;Jessica Krick
;
;This program creates a radial profile of a star in a fits image.  
; first it blocks out the bad column and diffraction spikes, if given,
; then it displays the image on the screen.  The radial distance from the
; center, and the counts at that distance are stored in an array. counts 
; above and below a threshold are not included in the array.  The array
; is binned into 1 pixel radii, and the counts are averaged within the bin.
; The counts are then normalized and plotted.
; 
;
;input: filename, xcenter, ycenter, location of bad columns, and saturated stars
;output: fielname.prof, idlprof.ps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro satstar
;device, true=24
;device, decomposed=0

close, /all		;close all files = tidiness
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'
colors = GetColor(/load, Start=1)

device, filename = '/n/Godiva1/jkrick/A3888/final/starprof.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

;filename = strarr(1)  	; declare filename to be a stringarray
;imagefile= "/n/Godiva1/jkrick/satstar/new/ccd3126.m.div.sub.fits"
;datafile = "/n/Godiva1/jkrick/satstar/new/ccd3126.m.prof"
;imagefile = "/n/Godiva1/jkrick/satstar/new/ccd5062.s.div.fits"
;datafile = "/n/Godiva1/jkrick/satstar/new/ccd5062.s.prof"

imagefile = "/n/Godiva1/jkrick/A3888/final/fullr2.fits"
datafile = "/n/Godiva1/jkrick/A3888/final/starprof.prof"


OPENW, lun, datafile, /GET_LUN
FITS_READ, imagefile, data, header


xcenter = ulong(2704.57)	;493;1064	; x coord of center of the star
ycenter = ulong(2472.95)	;460;1009	; y coord of center of the star

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

FOR x = 0, 3320, 1 DO BEGIN
    FOR y = 0, 3260, 1 DO BEGIN
        s = abs(xcenter - x)
        t = abs(ycenter - y)
        r = float(sqrt(s^(2)+(t^2))) ;pythagoras
        p = data[x,y] ;-1.527              ;background subtraction
        IF (p GT -1.0 AND P LT 33000) THEN BEGIN ;block the bad data pts
            counts(j) = p - 2.0
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
counter = 0.
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
        counter = 0.
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
plot, alog10(binnedrad*(0.259)), alog10(binnedcounts*0.259), thick = 3;[60:140]
;oplot,  alog10(findgen(2000)),  alog10(findgen(2000) - findgen(2000) + 0.8*0.259)

print, sortedradius(i-1)

start = [200000.0]
err = dindgen(i) - dindgen(i) + 1
result = mpfitfun('rcubed', binnedrad[60:140],binnedcounts[60:140] + 2, err[60:140], start)

oplot,  alog10(binnedrad*(0.259)), alog10( (result(0)*binnedrad^(-3))), color = colors.red

plot, binnedrad[60:140], binnedcounts[60:140]
openw, junklun, "/n/Godiva1/jkrick/A3888/final/junk.prof",/get_lun
FOR wxyz = 0,200, 1 DO BEGIN
    printf, junklun, binnedrad(wxyz), binnedcounts(wxyz)
ENDFOR
close, junklun
free_lun, junklun
oplot, binnedrad[60:140], result(0)*(binnedrad[60:140]^(-3.))
;==================================================================
imagefile = "/n/Godiva1/jkrick/A3888/final/satstartest.fits"
datafile = "/n/Godiva1/jkrick/A3888/final/satstarprof.prof"


OPENW, lun, datafile, /GET_LUN
FITS_READ, imagefile, data, header


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

FOR x = 0, 2800, 1 DO BEGIN
    FOR y = 0, 2700, 1 DO BEGIN
        s = abs(xcenter - x)
        t = abs(ycenter - y)
        r = float(sqrt(s^(2)+(t^2))) ;pythagoras
        p = data[x,y] ;-1.527              ;background subtraction
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
oplot, alog10(binnedrad*(0.259)), alog10(binnedcounts*0.259), thick = 3, color = colors.blue
;oplot,  alog10(findgen(2000)),  alog10(findgen(2000) - findgen(2000) + 0.8*0.259)
;=====================================================================

device, /close
set_plot, mydevice

END
;set bad column = -100
;data[1008:1019,640:2039] = -100
;data[1005:1015,640:661] = -100
;data[1:10,1:2039] = -100
;data[2030:2039,1:2039] = -100
;data[1:2039,1:10] = -100
;data[1:2039,2030:2039] = -100

;set diffraction spikes = -100
;for ccd3126
;data[1034:1086,1022:1522] = -1000
;data[1038:1098,496:996] = -1000
;data[1077:1577,990:1046] = -1000
;data[551:1051,970:1020] = -1000

;for ccd5062
;data[486:502,496:2000] = -100
;data[486:502,222:432] = -100
;data[454:530,0:222] = -100
;data[306:900,0:50] = -100
;data[450:540,50:310] = -100
;data[540:990,366:510] = -100
;data[1:410,438:586] = -100
;data[468:488,364:384] = -100


;k = ulong(0)			;location in sortedarray
;l = 0.0				;bin number
;rad = FLTARR(j) 		;array of radii of center of bin
;averagearr = FLTARR(j)  	;array of average values corresponding to radii
;r = sortedradius[k]		;individual radius
;a = sortedcounts[k]		;individual counts

;rad(l) = r
;averagearr(l) = a

;k = k + 1
;r = sortedradius[k]		;individual radius
;a = sortedcounts[k]		;individual counts
;print, rad(l), averagearr(l)
;FOR i =0,max,1 DO BEGIN
;	l = l+1
;	rad(l) = r
;	averagearr(l) = a
;	k = k+1
;	r = sortedradius(k)
;	a = sortedcounts(k)
;	c = 1
;	IF (r = rad(l)) THEN BEGIN
;		averagearr(l) = averagearr(l) + a
;		c = c + 1
;		k = k + 1
;		r = sortedradius(k)
;		a = sortedcounts(k)
;	ENDIF
;	averagearr(l) = averagearr(l) / c
;print, rad(l), averagearr(l);
;ENDFOR

;normalize the averages
;i = 0
;norm = averagearr(0)
;;print, "norm value: ", norm
;WHILE (i LT max)  DO BEGIN
;	averagearr(i) = averagearr(i) / norm
;	printf, lun,rad(i), averagearr(i)
;	i = i +1
;ENDWHILE

;printf, lun,rad, averagearr


;FOR x = xcenter - max, xcenter + max, 1 DO BEGIN
;	FOR y = ycenter -max, ycenter + max, 1 DO BEGIN
;		;print, x, y
;		s = float(xcenter - x)
;		t = float(ycenter - y)
;		r = float(sqrt(s^(2)+(t^2)))  		;pythagoras
;		p = data[x,y]
;                IF (p GT -1 AND P LT 33000) THEN BEGIN ;block the bad data pts
;                        ;printf,lun,r,p			;output to datafile
;                    counts(j) = p
;                    radius(j) = r
;                    j= j + 1
;                ENDIF
;                ;print, x, y    
;	ENDFOR 
;ENDFOR
;before looping, check to make sure that the edges of the box are
;not off the image, 
; ! add one for greater than 2048
;IF (xcenter - max LT 0) THEN BEGIN
;	max = float(xcenter - 0)
;	print, "radius range is not allowed in x, reseting max = ", max
;ENDIF

;IF (ycenter - max LT 0) THEN BEGIN
;	max = float(ycenter - 0)
;	print, "radius range is not allowed in y, reseting max = ", max
;ENDIF
