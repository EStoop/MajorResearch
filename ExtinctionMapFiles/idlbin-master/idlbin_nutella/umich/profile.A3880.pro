;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;February 2002
;Jessica Krick
;
;This program creates a radial profile of a star(or anythiing really) in a fits image.  
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


pro profile;, filename;, xcenter, ycenter    ;if I want to make them command line arguments

device, true=24
device, decomposed=0

close, /all		;close all files = tidiness

colors = GetColor(/Load, Start=1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;declare variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;xcenter = ulong(2180)           ; x coord of center of the star
;ycenter = ulong(405)		; y coord of center of the star
;filename = strarr(1)     	; declare filename to be a stringarray
x = ulong(0)			; x pixel coord
y = ulong(0)			; y pixel coord
r = float(0.0)			; radius from center in pixels
p = float(0.0)			; count at a particular radius
sum = float(0.0)
arad = float(0.0)

max = 900.
radius = FLTARR((3*max)^2 )	; array of radius (float)
counts = FLTARR((3*max)^2 )	; array of counts (float)

filename = '/Users/jkrick/umich/icl/A3880/fullr'

imagefile = filename + '.fits'
datafile = filename + '.prof'

;A3880
;xcenter =1757.9
;ycenter=2505.7
;xcenter = 1774.3
;ycenter = 2607.8
;xcenter = 2188.75
;ycenter = 1101.0
;xcenter = 573.1
;ycenter = 2521.66
xcenter =573.1
ycenter=2521.66
;xcenter=795.99
;ycenter =3968.96

; read in image
FITS_READ, imagefile, data, header


; make an array of radii from the central pixel, and counts for the 
; corresponding radii
;--------------------------------------------------------------------------------

;error check, check to make sure that the edges of the box are
;not off the image, whoops!
; !!!!!!!!!!!!!!need to add one for greater than 2048

;need to make it smarter than this.  will want it to continue, just
;not off the edge, ie don't change max, change boundaries.
xmin = xcenter - 500;1700;1;1100.
xmax = xcenter + 500;2700;2000;2400;400.
ymin = ycenter - 500;600;1;1800.
ymax = ycenter + 500;1600;1000;3100;400.

;diffraction spike removal
data[570,*] = -100
data[571,*] = -100
data[572,*] = -100
data[573,*] = -100
data[574,*] = -100
data[575,*] = -100
data[576,*] = -100




; Two for loops get all x and y in a square = 2*max around the center
;it then outputs the radius and counts to a file and into two arrays

j = ulong64(0)						;counter
increment = 1						;sampling distance
FOR x = xmin,xmax-1, increment DO BEGIN
	FOR y = ymin,ymax-1, increment DO BEGIN
		s = float(xcenter - x)
		t = float(ycenter - y)
		r = float(sqrt(s^(2)+(t^2)))  		;pythagoras
		p = data[x,y]
		IF (p GT -50. AND p LT 31.)  THEN BEGIN;block the bad data pts
			;print,j," j ",r," r ",p, " p "	;output to datafile
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
;this way both arrays get sorted appropriately
sortindex = Sort(radius)
sortedradius = radius[sortindex]
sortedcounts = counts[sortindex]

;binning
;----------------------------------------------------------------------------------

OPENW, lun, datafile, /GET_LUN

c = 0.0				;counter to keep track of radius (starting value of bin)
dx = 1.0			;bin width
k = ulong(0)			;location in sortedarray
l = 0.0				;bin number
average = 0.0			;average counts per bin
rad = FLTARR(j) 		;array of radii of center of bin
averagearr = FLTARR(j)  	;array of average values corresponding to radii
err = FLTARR(j)                 ;array of scatter within the bin
a = sortedradius[k]		;individual radius
b = sortedcounts[k]		;individual counts


;pick off the first one, and then bin the rest
rad(l) = a
averagearr(l) = b
err(l) = 1.0
l = l+1
k = k + 1
a = sortedradius[k]		;individual radius
b = sortedcounts[k]		;individual counts
c = c + 0.5


;this set of loops bins the data, and then averages the counts in each
;bin
c = fix(a)
WHILE (c LT sortedradius(j-1)) DO BEGIN
    sum = 0.0			;sum of counts in the bin
    i = 0                       ;number of individual measurements in the bin
    arad = 0.0			;average radius
    big = -1.0
    small = 1.0
    WHILE (a LT (c + dx)) AND  (a LT sortedradius(j-1)) DO BEGIN ;while in the bin
        
        sum = sum + b           ;add the value of counts to the sum
        k = k+ 1                ;increment how many counts overall
        arad = arad +a
        
        IF b GT big THEN big = b ;find the scatter to put some sort of err bars
        IF b LT small THEN small = b
        
        a = sortedradius[k]     ;read in new set of raddi and counts
        b = sortedcounts[k]
        i = i+1                 ;increment number of measurements in bin
    ENDWHILE
    
    IF (i EQ 0 ) THEN BEGIN     ; no measurements in the bin = totally blocked
                                ;SKIP IT, make it the previous value
        rad(l) = rad(l-1) + dx  ;radius at the average of the radii within the bin
        averagearr(l) = average ;add the radii and averages to arrays
        
    ENDIF ELSE BEGIN
        average = sum/i         ;calculate average in the bin
        arad = arad /i          ;calculate the average radius within the bin
        
        rad(l) = arad           ;record the radius at the average of the radii within the bin
        averagearr(l) = average ;add the radii and averages to arrays
       
    ENDELSE
    
    scatter = big - small
    err(l) = scatter
;    IF(averagearr(l) GT 32000) THEN BEGIN
;        averagearr(l) = averagearr(l-1)
;        print, "definitely shouldn't be here"
;    ENDIF
    c = c+dx                    ;change to the next bin
    l = l+ 1                    ;
    
ENDWHILE


;shorten the arrays to be only as long as they need to be
rad= rad(0:l-1)
averagearr = averagearr(0:l-1)
err = err(0:l-1)
;print, rad, averagearr
;print, "l", l
counter = ulong64(0)
;print out values to a file
FOR counter = 0, l-2, 1 DO BEGIN
    printf, lun, rad(counter), averagearr(counter)
ENDFOR

;plot, rad, averagearr, xrange=[0,200]

df = fltarr(l-1)
dfsum = fltarr(l-1)
FOR j = 0,l - 2,1 DO BEGIN

    df(j) =  averagearr(j)- averagearr(j+1)
    dfsum(j) =  averagearr(j) - (1./(l-1.-j))*(total(averagearr(j:l-1.))) 
ENDFOR





mydevice = !D.NAME
!p.multi = [0, 0, 2]
SET_PLOT, 'ps'

device, filename='/Users/jkrick/umich/icl/A3880/profile.ps', /portrait,$
  BITS=8, scale_factor=0.9 , /color

a = [0,0]
b = [0,0]
;plot, a, b, xrange = [0,50], yrange = [30,20]
plot, alog10(rad[0:l-5]*0.435), 22.04 - 2.5*alog10((averagearr[0:l-5])/(0.435^2)),thick = 3,xrange=[0,3], yrange=[30,15],charthick = 3,xtitle="log radius (arcsec)",ytitle = "SB",xthick=3,ythick=3
;plot, alog10(rad*0.435), averagearr,thick = 3,xrange=[-1,3], yrange=[0,40],charthick = 3,xtitle="radius (arcsec)",ytitle = "counts/s",xthick=3,ythick=3

close, lun
free_lun, lun
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


max = 900.
radiusp = FLTARR((3*max)^2 )	; array of radius (float)
countsp = FLTARR((3*max)^2 )	; array of counts (float)

filename = '/Users/jkrick/umich/icl/A3880/fullrp'

imagefile = filename + '.fits'
datafile = filename + '.prof'


; read in image
FITS_READ, imagefile, data, header
;diffraction spike removal
data[570,*] = -100
data[571,*] = -100
data[572,*] = -100
data[573,*] = -100
data[574,*] = -100
data[575,*] = -100
data[576,*] = -100


; make an array of radii from the central pixel, and counts for the 
; corresponding radii
;--------------------------------------------------------------------------------

;
;
; two for loops get all x and y in a square = 2*max around the center
;it then outputs the radius and counts to a file and into two arrays

j = ulong64(0)						;counter
increment = 1						;sampling distance
FOR x = xmin,xmax-1, increment DO BEGIN
	FOR y = ymin,ymax-1, increment DO BEGIN
		s = float(xcenter - x)
		t = float(ycenter - y)
		r = float(sqrt(s^(2)+(t^2)))  		;pythagoras
		p = data[x,y]
		IF (p GT -50. AND p LT 31.) THEN BEGIN;block the bad data pts
			;print,j," j ",r," r ",p, " p "	;output to datafile
			countsp(j) = p
			radiusp(j) = r
			j= j + 1
		ENDIF
	ENDFOR 
ENDFOR

;shorten the arrays to be only as long as need be
radiusp = radiusp(0:j-1)
countsp = countsp(0:j-1)


;sort the radius array, and then apply that sorting order to the counts
;this way both arrays get sorted appropriately
sortindexp = Sort(radiusp)
sortedradiusp = radiusp[sortindexp]
sortedcountsp = countsp[sortindexp]

;binning
;----------------------------------------------------------------------------------

OPENW, lun, datafile, /GET_LUN

c = 0.0				;counter to keep track of radius (starting value of bin)
dx = 1.0			;bin width
k = ulong(0)			;location in sortedarray
l = 0.0				;bin number
average = 0.0			;average counts per bin
radp = FLTARR(j) 		;array of radii of center of bin
averagearrp = FLTARR(j)  	;array of average values corresponding to radii
err = FLTARR(j)                 ;array of scatter within the bin
a = sortedradiusp[k]		;individual radius
b = sortedcountsp[k]		;individual counts


;pick off the first one, and then bin the rest
radp(l) = a
averagearrp(l) = b
err(l) = 1.0
l = l+1
k = k + 1
a = sortedradiusp[k]		;individual radius
b = sortedcountsp[k]		;individual counts
c = c + 0.5


;this set of loops bins the data, and then averages the counts in each
;bin
c = fix(a)
WHILE (c LT sortedradiusp(j-1)) DO BEGIN
    sum = 0.0			;sum of counts in the bin
    i = 0                       ;number of individual measurements in the bin
    arad = 0.0			;average radius
    big = -1.0
    small = 1.0
    WHILE (a LT (c + dx)) AND  (a LT sortedradiusp(j-1)) DO BEGIN ;while in the bin
        
        sum = sum + b           ;add the value of counts to the sum
        k = k+ 1                ;increment how many counts overall
        arad = arad +a
        
        IF b GT big THEN big = b ;find the scatter to put some sort of err bars
        IF b LT small THEN small = b
        
        a = sortedradiusp[k]     ;read in new set of raddi and counts
        b = sortedcountsp[k]
        i = i+1                 ;increment number of measurements in bin
    ENDWHILE
    
    IF (i EQ 0 ) THEN BEGIN     ; no measurements in the bin = totally blocked
                                ;SKIP IT, make it the previous value
        radp(l) = radp(l-1) + dx  ;radius at the average of the radii within the bin
        averagearrp(l) = average ;add the radii and averages to arrays
        
    ENDIF ELSE BEGIN
        average = sum/i         ;calculate average in the bin
        arad = arad /i          ;calculate the average radius within the bin
        
        radp(l) = arad           ;record the radius at the average of the radii within the bin
        averagearrp(l) = average ;add the radii and averages to arrays
       
    ENDELSE
    
    scatter = big - small
    err(l) = scatter
;    IF(averagearr(l) GT 32000) THEN BEGIN
;        averagearr(l) = averagearr(l-1)
;        print, "definitely shouldn't be here"
;    ENDIF
    c = c+dx                    ;change to the next bin
    l = l+ 1                    ;
    
ENDWHILE


;shorten the arrays to be only as long as they need to be
radp= radp(0:l-1)
averagearrp = averagearrp(0:l-1)
err = err(0:l-1)
;print, rad, averagearr
;print, "l", l
counter = ulong64(0)
;print out values to a file
FOR counter = 0, l-2, 1 DO BEGIN
    printf, lun, radp(counter), averagearrp(counter)
ENDFOR

;plot, rad, averagearr, xrange=[0,200]

df = fltarr(l-1)
dfsum = fltarr(l-1)
FOR j = 0,l - 2,1 DO BEGIN

    df(j) =  averagearrp(j)- averagearrp(j+1)
    dfsum(j) =  averagearrp(j) - (1./(l-1.-j))*(total(averagearrp(j:l-1.))) 
ENDFOR


;plot, a, b, xrange = [0,50], yrange = [30,20]
oplot, alog10(radp[0:l-5]*0.435), 22.04 - 2.5*alog10((averagearrp[0:l-5])/(0.435^2)),thick = 3, color = colors.red

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
print, n_elements(averagearr), n_elements(averagearrp)
delta = averagearr - averagearrp
plot, alog10(radp*0.435), delta,yrange=[-1,1], xrange=[0,3], ytitle = "delta SB"


close, lun
free_lun, lun




;close, outlun
;free_lun, outlun



device, /close
set_plot, mydevice

END


