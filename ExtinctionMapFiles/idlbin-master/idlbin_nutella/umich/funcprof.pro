;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;February 2002
;Jessica Krick
;
;This program creates a radial profile of a star(or anythiing really) in a fits image.  
; first it blocks out the bad column and diffraction spikes, if given,
; then it displays the image on the screen.  The radial distance from the
; center, and the counts at that distance are stored in an array. counts 
; above and below a threshold are not included in the array.  The array
; is binned into 1 pixel radii, and the counts are averaged within the bin.`
; The counts are then normalized and plotted.
; 
;
;input: filename, xcenter, ycenter, location of bad columns, and saturated stars
;output: fielname.prof, idlprof.ps
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION funcprof, filename, xcenter, ycenter    ;if I want to make them command line arguments

xcenter= xcenter-1.
ycenter=ycenter-1.
device, true=24
device, decomposed=0

colors = GetColor(/Load, Start=1)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;declare variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
x = ulong(0)			; x pixel coord
y = ulong(0)			; y pixel coord
r = float(0.0)			; radius from center in pixels
p = float(0.0)			; count at a particular radius
sum = float(0.0)
arad = float(0.0)
;max = 2000.
max = 1000.
radius = FLTARR((3*max)^2 )	; array of radius (float)
counts = FLTARR((3*max)^2 )	; array of counts (float)

; create file names
imagefile = filename + '.fits'
datafile = filename + '.prof'

; read in image
FITS_READ, imagefile, data, header

xmin = 0.
xmax = 1000.
;xmax = 2000 
ymin = 0.
ymax = 1000.
;ymax = 2000
; two for loops get all x and y in a square = 2*max around the center
;it then outputs the radius and counts to a file and into two arrays

j = ulong64(0)						;counter
increment = 1						;sampling distance
;for all pixels in the image
FOR x = xmin,xmax-1, increment DO BEGIN
	FOR y = ymin,ymax-1, increment DO BEGIN
		s = float(xcenter - x)
		t = float(ycenter - y)
		r = float(sqrt(s^(2)+(t^2)))  		;pythagoras
		p = data[x,y]
                                           IF NOT finite(p) THEN BEGIN
                                                print, "p NaN, returning now"
                                                return, [0,0]
                                            endiF

		IF (p GT -0.1 AND p LT 5000.) THEN BEGIN	;block the bad data pts
  
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

c = 0.0				;counter to keep track of radius (starting value of bin)
dx = 1.0			                      ;bin width
k = ulong(0)			;location in sortedarray
l = 0.0				;bin number
average = 0.0			;average counts per bin
rad = FLTARR(j) 		;array of radii of center of bin
averagearr = FLTARR(j)  	;array of average values corresponding to radii
err = FLTARR(j)                 ;array of scatter within the bin


;pick off the first one, and then bin the rest
a = sortedradius[k]		;individual radius
b = sortedcounts[k]		;individual counts
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

    WHILE (a LT (c + dx)) AND (a LT sortedradius(j-1))DO BEGIN ;while in the bin
        sum = sum + b           ;add the value of counts to the sum
        arad = arad +a
        k = k+ 1                ;increment how many counts overall
      
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
    IF(averagearr(l) GT 28000) THEN BEGIN
        averagearr(l) = averagearr(l-1)
        print, "definitely shouldn't be here"
    ENDIF
    c = c+dx                    ;change to the next bin
    l = l+ 1                    ;
    
ENDWHILE


;shorten the arrays to be only as long as they need to be
rad= rad(0:l-1)
averagearr = averagearr(0:l-1)
err = err(0:l-1)

;output them to a file
;openw, nomask, "/n/Godiva1/jkrick/A3888/galfit/prof.out", /GET_LUN
      
;FOR p = 0, l-1, 1 DO BEGIN
;  printf, nomask, rad(p), averagearr(p)
;ENDFOR
;close, nomask
;free_lun, nomask

dflux = fltarr(l-1)
dfsum = fltarr(l-1)
FOR j = 0,l-2,1 DO BEGIN

    dflux(j) =  averagearr(j)- averagearr(j+1)
;    dfsum(j) =  averagearr(j) - (1./(l-1-j))*(total(averagearr(j:l-1))) 
    dfsum(j) = averagearr(j) - mean(averagearr(j:l-1))
    IF dfsum(j) LE 0.0 THEN dfsum(j) = 0.0001
ENDFOR




;curve fitting
;------------------------------------------------------------------------------
err = dindgen(l) - dindgen(l) + 1         ;make a simple error array, all even errors

weight = findgen(l)                    ;change the weights by eye, good = 50 - 120 "
weight(0:200) = 1
weight(200:500) = 2000
weight(500:l-1) = 1

averagearr = averagearr - 2
;averagearr = averagearr + 0.001
start = [0.01,160.0]
result = MPFITFUN('exponential',rad,averagearr, err, start, weights=weight);,PARINFO=pi);, weights=weight)

;print, "exponential result", result

;-----------------------------------------------------------------------------
openr, nomask, "/n/Godiva1/jkrick/A3888/galfit/prof.out", /GET_LUN
nradius = fltarr(3000)
naverage = fltarr(3000)
i = 0
WHILE ( NOT EOF(nomask)) DO BEGIN
    readf, nomask, nrad, nav
    nradius(i) = nrad
    naverage(i) = nav
    i = i + 1
ENDWHILE
close, nomask
free_lun, nomask

nradius = nradius(0:i-1)
naverage = naverage(0:i-1)
naverage = naverage - 2

print, nradius(100), naverage(100)
print, rad(100), averagearr(100)

;plotting
;---------------------------------------------------------------------------

value = fltarr(800)


;mydevice = !D.NAME
!p.multi = [0, 0, 2]
;SET_PLOT, 'ps'

;device, filename='/n/Godiva1/jkrick/A3888/galfit/expprof.ps', /portrait,$
;  BITS=8, scale_factor=0.9 , /color

ps_open, file = "/n/Godiva7/jkrick/grants/NSF2003/profile.ps", /portrait, /color, ysize = 8
;plot, rad, averagearr, thick = 3, charthick = 3
;, /ylog;, yrange = [1E-6 ,1E-1],thick = 3, charthick = 3

plot, rad*0.259, 24.3 + (-2.5)*(AlOG10(averagearr/(.259^2))),$
  title = 'Radial Profile of Possible ICL Component',$
  YTITLE='!4 l !3 (mag arcsec !E -2!N !3)',$
  thick = 3, color = colors.black,YRANGE=[31, 23], charthick = 3,$
  xthick = 3, ythick = 3;, xtickunits=['numeric','numeric'],$
;  xtickv=['rad','value']  XTITLE='radius(arcseconds)'
;   /ylog, YRANGE=[1E-4, 1E-1]


oplot, nradius*0.259, 24.3 + (-2.5)*(AlOG10(naverage/(.259^2))),$
  thick = 5, color = colors.red

;oplot,  rad*0.259, 24.3 + (-2.5)*(ALOG10((result(0)*exp(- rad*0.259/(result(1))))/(.259^2))),$
;  thick = 3,color = colors.blue
;flux = (result(0)*exp(-rad/result(1)))/ (0.259^2)
;oplot, rad*0.259, 24.3 + (-2.5) * (Alog10(flux)), thick = 3, color = colors.blue

;fluxall = (0.159101*exp(-rad/163.061  ))/ (0.259^2)
;oplot, rad*0.259, 24.3 + (-2.5) * (Alog10(fluxall)), thick = 3, color = colors.red

;fluxblock = ( 0.0434396 *exp(-rad/ 130.638))/ (0.259^2)
;oplot, rad*0.259, 24.3 + (-2.5) * (Alog10(fluxblock)), thick = 3, color = colors.violet

;xyouts, 200, 31, 'fit w/ galaxies subtracted ', color = colors.blue, charthick =3
;xyouts, 90, 24.5, 'fit w/ all abjects subtracted & blocked', color = colors.blue, charthick =3, charsize = 0.8
;xyouts, 90, 25.2, 'fit w/ all objects blocked, not subtracted', color = colors.violet, charthick =3, charsize = 0.8
xyouts, 110, 23.5, 'all object flux', color = colors.red, charthick =3, charsize = 0.9
xyouts, 110, 23.0, 'masks out to FF limit', color = colors.black, charthick =3, charsize = 0.9
;flat fielding error
oplot, rad*0.259, findgen(l) - findgen(l) + 29.6, color = colors.black, thick = 3
xyouts, 10, 29.5, "flat-fielding accuracy on 1' scales", color = colors.black, charthick = 3
;oplot, rad, result(0)*exp(-rad/result(1)), color= colors.blue, thick = 3
;oplot, rad, findgen(800) - findgen(800), color = colors.red, thick=3

;dfplot

;plot, rad*0.259, dflux, thick = 3

plot, rad*0.259, Alog10(dfsum),$
  XTITLE='Radius(arcseconds)', YTITLE='!3Log !4 D!3F Cumulative',$
  thick = 3, color = colors.black, charthick = 3,$
  xthick = 3, ythick = 3, yrange = [-5,-1];, yticks = 4;

;fluxdiff = ( 0.0416589 *exp(-rad/ 129.497))/ (0.259^2)
;print, fluxdiff
;dfsum = fltarr(l-1)
;FOR j = 0,l-2,1 DO BEGIN
;    dfsum(j) = fluxdiff(j) - mean(fluxdiff(j:l-1))
;    IF dfsum(j) LE 0.0 THEN dfsum(j) = 0.0001
;ENDFOR
;print, alog10(dfsum)

;oplot, rad*0.259,alog10(dfsum), thick = 3, color = colors.blue

ps_close, /noprint, /noid
;device, /close
;set_plot, mydevice

return, result

END




;plot, rad*0.259, 25 + (-2.5)*(AlOG10(averagearr/(.259^2))),$
;  title = 'Radial Profile of Possible ICL Component',$
;  XTITLE='radius(pixels)', YTITLE='Surface Brightness (mag/sqarcsec)',$
;  thick = 3, color = colors.black,yrange=[30,20];,  YRANGE=[32, 26];, xtickunits=['numeric','numeric'],$
;  xtickv=['rad','value']
;   /ylog, YRANGE=[1E-4, 1E-1]

;FOR x=0,30,1 DO BEGIN
;    b = 21. + 1.086*(7.6697*( (x/200.)^(1./4.)- 1))
;print, x, b
;ENDFOR

;b = result(0)* exp((-7.67) *(((rad / (result(1)))^(1./4.)) -1.))
;plot, rad, averagearr;, yrange=[0,10]
;oplot, rad, (result(0)) * (exp(-7.67*(((rad/(result(1)))^(1.0/4.0)) - 1.0))), thick = 3, $
;  color = colors.blue
;oplot, rad*0.259,  result(0) + 1.086*(7.6697*( (rad/result(1))^(1./4.)- 1)), color = colors.blue
;oplot, rad*0.259,  28. + 1.086*(7.6697*( (rad/200)^(1./4.)- 1)), color = colors.green

;result(0) = 0.004
;result(1) = 200
;b = (result(0) *exp((-7.67)*((rad/result(1))^(1./4.)-1)))
;oplot, rad, b, color = colors.red


;oplot, rad, df, thick = 3,color = colors.orange
;oplot, rad, dfsum, thick = 3,color = colors.yellow

;errplot, radnew, averagenew - errnew/2, averagenew + errnew/2

;oplot, rad,result(0) +result(1)*rad + result(2)*rad^2, $

;exponential_______________________________________________________________________
;oplot, findgen(800), 25 + (-2.5)*(ALOG10((result(0)*exp(-findgen(800)/(result(1))))/(.259^2))),$
;  thick = 3,color = colors.blue
;oplot, findgen(800),25 + (-2.5)*(ALOG10((0.143368*exp(-findgen(800)/187.231))/(0.259^2))),$
;  thick = 3,color = colors.red
;oplot, findgen(800),25 + (-2.5)*(ALOG10(0.0300519*exp(-findgen(800)/145.164))),$
;  thick = 3,color = colors.green
;oplot, findgen(800),25 + (-2.5)*(ALOG10(0.0287254*exp(-findgen(800)/147.585))),$
;  thick = 3,color = colors.orange
;oplot, rad, result(0)*(1 + ((rad) / result(1))^2)^(-2.592), color = colors.blue

;xyouts, 200, 31, 'fit w/ galaxies subtracted ', color = colors.blue, charthick =3
;xyouts, 100, 35, 'fit w/ galaxies subtracted & blocked', color = colors.orange, charthick =3
;xyouts, 100, 34.5, 'fit w/ galaxies blocked', color = colors.green, charthick =3
;xyouts, 200, 31.5, 'fit including all galaxy flux', color = colors.red, charthick =3

;oplot, findgen(800), p, color = colors.black, linestyle = 1, thick = 3
;oplot, findgen(800), n, color = colors.black, linestyle = 3, thick = 3



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


;normalize the averages
;i = 0
;norm = averagearr(0)
;norm = 1
;print, "norm value: ", norm
;WHILE (i LT c)  DO BEGIN
;	averagearr(i) = averagearr(i) / norm
;	printf, lun,rad(i), averagearr(i)
;	i = i +1
;ENDWHILE
