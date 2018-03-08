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


pro profile2
; xcenter, ycenter    ;if I want to make them command line arguments
;print, xcenter, ycenter
;device, true=24
;device, decomposed=0

close, /all		;close all files = tidiness

;declare variables
xcenter = ulong(1002.5)		; x coord of center of the star
ycenter = ulong(1699.75)		; y coord of center of the star
filename = strarr(1)  	; declare filename to be a stringarray
;filename = "/n/pele/jkrick/ind/testgal"
filename = "/n/Godiva1/jkrick/A3888/final/largeV2"


x = ulong(0)			; x pixel coord
y = ulong(0)			; y pixel coord
r = float(0.0)			; radius from center in pixels
p = float(0.0)			; count at a particular radius
max = 400.0

; create file names
imagefile = filename + '.fits'
datafile = filename + '.prof'
;datafile = "/n/Godiva1/jkrick/satstar/new/ccd3126.m.prof"
; read in image
FITS_READ, imagefile, data, header
;newimage = congrid(data,507.5,510)  	;makes smaller for display on screen
;tv, hist_equal(newimage)		;displays


; make an array of radii from the central pixel, and counts for the 
; corresponding radii

;before looping, check to make sure that the edges of the box are
;not off the image, whoops!
; !!!!!!!!!!!!!!need to add one for greater than 2048
IF (xcenter - max LT 0) THEN BEGIN
	max = float(xcenter - 0)
	print, "radius range is not allowed in x, reseting max = ", max
ENDIF

IF (ycenter - max LT 0) THEN BEGIN
	max = float(ycenter - 0)
	print, "radius range is not allowed in y, reseting max = ", max
ENDIF

radius = FLTARR((4*max)^2)	; array of radius (float)
counts = FLTARR((4*max)^2)	; array of counts (float)

; two for loops get all x and y in a square = 2*max around the center
;it then outputs the radius and counts to a file and into two arrays

j = ulong64(0)					;counter
increment = 1				;sampling distance
FOR x = xcenter - max, xcenter + max, increment DO BEGIN
	FOR y = ycenter -max, ycenter + max, increment DO BEGIN
		;print, x, y
		s = float(xcenter - x)
		t = float(ycenter - y)
		r = float(sqrt(s^(2)+(t^2)))  		;pythagoras
		p = data[x,y]
                IF (p GT -1 AND P LT 33000) THEN BEGIN ;block the bad data pts
                        ;printf,lun,r,p			;output to datafile
                    counts(j) = p
                    radius(j) = r
                    j= j + 1
                ENDIF
                ;print, x, y    
	ENDFOR 
ENDFOR
print, "j = ", j

;shorten the arrays to be only as long as need be
radius = radius(0:j-1)
counts = counts(0:j-1)

;sort the radius array, and then apply that sorting order to the counts
;this way both arrays get sorted appropriately
sortindex = Sort(radius)
sortedradius = radius[sortindex]
sortedcounts = counts[sortindex]



OPENW, lun, datafile, /GET_LUN


k = ulong(0)			;location in sortedarray
l = 0.0				;bin number
rad = FLTARR(j) 		;array of radii of center of bin
averagearr = FLTARR(j)  	;array of average values corresponding to radii
r = sortedradius[k]		;individual radius
a = sortedcounts[k]		;individual counts

rad(l) = r
averagearr(l) = a

k = k + 1
r = sortedradius[k]		;individual radius
a = sortedcounts[k]		;individual counts
;print, rad(l), averagearr(l)
FOR i =0,max,1 DO BEGIN
	l = l+1
	rad(l) = r
	averagearr(l) = a
	k = k+1
	r = sortedradius(k)
	a = sortedcounts(k)
	c = 1
	IF (r = rad(l)) THEN BEGIN
		averagearr(l) = averagearr(l) + a
		c = c + 1
		k = k + 1
		r = sortedradius(k)
		a = sortedcounts(k)
	ENDIF
	averagearr(l) = averagearr(l) / c
;print, rad(l), averagearr(l)
ENDFOR

;print, rad(l), averagearr(l)

;printf, lun,rad, averagearr

;normalize the averages
i = 0
norm = averagearr(0)
;print, "norm value: ", norm
WHILE (i LT max)  DO BEGIN
	averagearr(i) = averagearr(i) / norm
	printf, lun,rad(i), averagearr(i)
	i = i +1
ENDWHILE

print, rad
;printf, lun,rad, averagearr
close, lun
free_lun, lun

;set up for plotting, then plot
!p.multi = [0, 1, 1]
SET_PLOT, 'ps'
device, filename='/n/Godiva1/jkrick/A3888/final/galprofile.ps', /portrait,$
                BITS=8, scale_factor=0.9 , /color

plot, rad, averagearr, title = $
	'radial profile of a galaxy made completely with idl', PSYM=2, $
	XTITLE='radius', YTITLE='normalized counts',XRANGE = [0,100], YRANGE = [0,1]



device, /close
set_plot, 'x'
END

;print, sortedradius[0], sortedcounts[0]
;print, sortedradius[1], sortedcounts[1]
;print, sortedradius[2], sortedcounts[2]
;print, sortedradius[3], sortedcounts[3]
;print, sortedradius[4], sortedcounts[4]
;print, sortedradius[5], sortedcounts[5]
;print, sortedradius[6], sortedcounts[6]
;print, sortedradius[7], sortedcounts[7]
;print, sortedradius[8], sortedcounts[8]
;print, sortedradius[9], sortedcounts[9]
;print, sortedradius[10], sortedcounts[10]
;print, sortedradius[11], sortedcounts[11]
;print, sortedradius[12], sortedcounts[12]
;print, sortedradius[13], sortedcounts[13]
;print, sortedradius[14], sortedcounts[14]
;print, sortedradius[15], sortedcounts[15]
;print, sortedradius[16], sortedcounts[16]
;print, sortedradius[17], sortedcounts[17]
;print, sortedradius[18], sortedcounts[18]
;print, sortedradius[19], sortedcounts[19]

;set bad column = -100
;data[1012:1015,640:2030] = -100
;data[1005:1015,640:661] = -100

;set diffraction spikes = -100
;for ccd3126
;data[1052:1073,1022:1522] = -100
;data[1052:1073,496:996] = -100
;data[1077:1577,1008:1021] = -100
;data[551:1051,998:1009] = -100

;for ccd5061
;data[1546:1555,1580:1913] = -100
;data[1546:1555,1271:1562] = -100
;data[1563:1674,1566:1578] = -100
;data[1431:1543,1570:1579] = -100

;for ccd5062
;data[486:502,496:2000] = -100
;data[486:502,222:432] = -100
;data[454:530,0:222] = -100


;for ccd3120
;data[248:255,359:390] = -100
;data[244:252,410:441] = -100

;for ccd3117
;data[965:972,976:1040] = -100
;data[969:976,892:956] = -100

;display masked image on screen
;twoimage = congrid(data,507.5,510) 
;tv, hist_equal(twoimage)


