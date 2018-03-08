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

;device, true=24
;device, decomposed=0

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

; create file names
;filename = "/n/Godiva1/jkrick/A3888/centerV.div"
;filename = "/n/Godiva1/jkrick/A3888/ellipgal"
;filename = "/n/Godiva3/jkrick/satstar/ccd3041.s"
;filename = "/n/Godiva6/jkrick/A3880/original/fullB.t"
;filename = "/n/Godiva7/jkrick/A118/original/larger.s"
filename = '/Users/jkrick/umich/icl/A3880/fullr'

imagefile = filename + '.fits'
datafile = filename + '.prof'
;datafile = "/n/Godiva1/jkrick/satstar/new/ccd3126.m.prof"

;xcenter = 1014.
;ycenter = 1573.
;xcenter = 1585.8  ;unsat 5090
;ycenter = 1422.96
;xcenter = 433.2  ;very mildly sat 5090
;ycenter = 2402.0
;xcenter = 1507.4
;ycenter = 621.0 ;unsat 5071
;xcenter = 73.5
;ycenter = 1678.
;xcenter=1033.
;ycenter=1717.5

;A3880
;xcenter =1757.9
;ycenter=2505.7
;xcenter = 1774.3
;ycenter = 2607.8
xcenter = 2188.75
ycenter = 1101.0

;A118
;xcenter=1224.6
;ycenter =953.8
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
xmin = 1700;1;1100.
xmax = 2700;2000;2400;400.
ymin = 600;1;1800.
ymax = 1600;1000;3100;400.

;IF (xcenter - max LT 0) THEN BEGIN
;	max = float(xcenter - 0)
;	print, "radius range is not allowed in x, reseting max = ", max
;ENDIF;

;IF (ycenter - max LT 0) THEN BEGIN
;	max = float(ycenter - 0)
;	print, "radius range is not allowed in y, reseting max = ", max
;ENDIF



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
		IF (p GT -50. AND p LT 41.) THEN BEGIN;block the bad data pts
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
print, "l", l
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



;only really want error bars on every tenth point so I can see them on
;the graph, so make smaller arrays with every tenth data point
;radnew = fltarr(l)
;averagenew = fltarr(l)
;errnew = fltarr(l)

;i = 0
;FOR z =0, l-1,10 DO BEGIN
;    radnew(i) = rad(z)
;    averagenew(i) = averagearr(z)
;    errnew(i) = err(z)
;    i = i+1
;ENDFOR;

;radnew = radnew(0:i-1)
;averagenew = averagenew(0:i-1)
;errnew = errnew(0:i-1)


;curve fitting
;------------------------------------------------------------------------------
;err = dindgen(l) - dindgen(l) + 1         ;make a simple error array, all even errors

;weight = findgen(max)                    ;change the weights by eye
;weight(0:200) = 1
;weight(200:500) = 2000
;weight(500:fix(max)-1) = 1
;weight(0:10) = 1
;weight(10:250) = 2000
;weight(250:399) = 1
;pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},2)
;pi(1).limited(0) = 1
;pi(1).limits(0) =0.1
;averagearr = averagearr ;- 2
;sortedcounts = sortedcounts ;- 2

;print, "sortedradius", sortedradius
;print, "sortedcounts", sortedcounts
;
;FOR t = 0, j - 1, 1 DO BEGIN
;    IF (sortedcounts[t] LT 0) THEN sortedcounts[t] = sortedcounts[t - 1]
;endfor

;;newcounts = rebin(sortedcounts, 100)
;;newrad = rebin(sortedradius, 100)
;;start = [0.01,500.0]
;result = MPFITFUN('exponential',rad(4:l-1),averagearr(4:l-1), err, start);,PARINFO=pi);, weights=weight)
;result = MPFITFUN('exponential',rad,averagearr, err, start, weights=weight);,PARINFO=pi);, weights=weight)

;start = [29, 189.0]
;result = MPFITFUN('exponential',rad, averagearr, err, start, weights=weight)

;start = [0.027, 3.0]
;result = MPFITFUN('moffat',rad, averagearr, err, start);, weights=weight)

;print, "exponential result", result

;openw, outlun, "/n/Godiva7/jkrick/galfit/expback/listresults",  /GET_LUN, /append
;printf, outlun, xcenter+1 ,ycenter+1, result
;close, outlun
;free_lun, outlun

;plotting
;---------------------------------------------------------------------------

;value = fltarr(800)

;FOR m=0,j-1,1 DO BEGIN
;    value(m) = rad(m) / 0.435
;endfor

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename='/Users/jkrick/umich/icl/A3880/profile.ps', /portrait,$
  BITS=8, scale_factor=0.9 , /color

a = [0,0]
b = [0,0]
;plot, a, b, xrange = [0,50], yrange = [30,20]
plot, alog10(rad[0:l-5]*0.435), 22.04 - 2.5*alog10((averagearr[0:l-5])/(0.435^2)),thick = 3,xrange=[0,3], yrange=[30,5],charthick = 3,xtitle="log radius (arcsec)",ytitle = "SB",xthick=3,ythick=3
;plot, alog10(rad*0.435), averagearr,thick = 3,xrange=[-1,3], yrange=[0,40],charthick = 3,xtitle="radius (arcsec)",ytitle = "counts/s",xthick=3,ythick=3

print, "after first plot"
;oplot, rad, result(0)*exp(-rad/result(1)), color= colors.blue

;p = fltarr(800)
;FOR m = 0, 799,1 DO BEGIN
;    p(m) = 29.5
;ENDFOR
;n = fltarr(800)
;FOR m = 0, 799,1 DO BEGIN
;    n(m) = 30.6
;ENDFOR

;scale = 55.;165
;scale = 330;900.
scale = 200
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;openw, outlun, "/n/Godiva3/jkrick/satstar/combined.prof", /get_lun

;rad = fltarr(10000)
;counts = fltarr(10000)
;i = 0
;openr, lun1, "/n/Godiva3/jkrick/satstar/ccd5090.s.prof",/get_lun
;WHILE (NOT EOF(lun1)) DO BEGIN
;    readf, lun1, r, c
;    rad[i] = r
;    counts[i] = (100.*(c / 8.))/scale
;    i = i +1
;ENDWHILE
;rad = rad[0:i-1]
;counts = counts[0:i-1]
;oplot, alog10(rad[0:60]*0.435), 22.04 - 2.5*alog10(counts[0:60]/(0.435^2)), thick = 3, linestyle = 5;;color = colors.red
;
;FOR j=0,59,1 DO BEGIN
;printf, outlun, rad[j],counts[j]
;ENDFOR
;
;close, lun1
;free_lun, lun1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;rad = fltarr(10000)
;counts = fltarr(10000)
;i = 0
;openr, lun1, "/n/Godiva3/jkrick/satstar/ccd5071.s.prof",/get_lun
;WHILE (NOT EOF(lun1)) DO BEGIN
;    readf, lun1, r, c
;    rad[i] = r
;    counts[i] = (c / 20.)/scale
;    i = i +1
;ENDWHILE
;rad = rad[0:i-1]
;counts = counts[0:i-1]
;oplot, alog10(rad[60:830]*0.435), 22.04 - 2.5*alog10(counts[60:830]/(0.435^2)), thick = 3, linestyle; = 2;color = colors.cyan;
;
;FOR j=60,1000,1 DO BEGIN
;printf, outlun, rad[j],counts[j]
;ENDFOR
;
;close, lun1
;free_lun, lun1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;rad = fltarr(10000)
;counts = fltarr(10000)
;i = 0
;openr, lun1, "/n/Godiva3/jkrick/satstar/ccd3040.s.prof",/get_lun
;WHILE (NOT EOF(lun1)) DO BEGIN
;    readf, lun1, r, c
;    rad[i] = r
;    counts[i] = (c /57.)/scale
;    i = i +1
;ENDWHILE
;rad = rad[0:i-1]
;counts = counts[0:i-1]
;oplot,alog10(rad*0.435), 22.2 - 2.5*alog10(counts/(0.435^2)), thick = 3, linestyle = 2;color = colors.green;

;FOR j=51,1000,1 DO BEGIN
;printf, outlun, rad[j],counts[j]
;ENDFOR

;close, lun1
;free_lun, lun1
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


close, lun
free_lun, lun
;close, outlun
;free_lun, outlun



device, /close
set_plot, mydevice

END


