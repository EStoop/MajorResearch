FUNCTION myprof, xcenter, ycenter,data, blankdata, maxradius, xmax,ymax, number, maxfit
colors = GetColor(/Load, Start=1)
;make a profile of the galaxy
xcenter= xcenter-1.
ycenter= ycenter-1.

radius = FLTARR((3*maxradius)^2 )     ; array of radius (float)
counts = FLTARR((3*maxradius)^2 )     ; array of counts (float)

;error check, check to make sure that the edges of the box are
;not off the image, whoops!
; !!!!!!!!!!!!!!need to add one for greater than 2048
IF (xcenter - maxradius LT 0) THEN BEGIN
    maxradius = float(xcenter - 0)
    print, "radius range is not allowed in x, reseting max = ", max
ENDIF

IF (ycenter - maxradius LT 0) THEN BEGIN
    maxradius = float(ycenter - 0)
    print, "radius range is not allowed in y, reseting max = ", max
ENDIF



; two for loops get all x and y in a square = 2*max around the center
;it then outputs the radius and counts to a file and into two arrays

j = ulong64(0)						;counter
increment = 1						;sampling distance
FOR x = xcenter - maxradius, xcenter + maxradius, increment DO BEGIN
	FOR y = ycenter -maxradius, ycenter + maxradius, increment DO BEGIN
		s = float(xcenter - x)
		t = float(ycenter - y)
		r = float(sqrt(s^(2)+(t^2)))  		;pythagoras
		p = data[x,y]
		IF (p GT -0.1 AND p LT 5000) THEN BEGIN	;block the bad data pts
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

;OPENW, lun, datafile, /GET_LUN

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
WHILE (c LT maxradius) DO BEGIN
	sum = 0.0			;sum of counts in the bin
	i = 0				;number of individual measurements in the bin
	arad = 0.0			;average radius
        big = -1.0
        small = 1.0
	WHILE (a LT (c + dx)) DO BEGIN		;while in the bin

		sum = sum + b			;add the value of counts to the sum
		k = k+ 1			;increment how many counts overall
		arad = arad +a
           
                IF b GT big THEN big = b        ;find the scatter to put some sort of err bars
                IF b LT small THEN small = b
                        
		a = sortedradius[k]		;read in new set of raddi and counts
		b = sortedcounts[k]
		i = i+1				;increment number of measurements in bin
	ENDWHILE

        IF (i EQ 0 ) THEN BEGIN         ; no measurements in the bin = totally blocked
            ;SKIP IT, make it the previous value
            rad(l) = rad(l-1) + dx     ;radius at the average of the radii within the bin
            averagearr(l) = average ;add the radii and averages to arrays

        ENDIF ELSE BEGIN
            average = sum/i     ;calculate average in the bin
            arad = arad /i      ;calculate the average radius within the bin
 
            rad(l) = arad       ;record the radius at the average of the radii within the bin
            averagearr(l) = average ;add the radii and averages to arrays

        ENDELSE

        scatter = big - small
        err(l) = scatter
        IF(averagearr(l) GT 28000) THEN BEGIN
            averagearr(l) = averagearr(l-1)
            print, "definitely shouldn't be here"
        ENDIF
        c = c+dx                ;change to the next bin
        l = l+ 1                ;
        
ENDWHILE


;shorten the arrays to be only as long as they need to be
rad= rad(0:l-1)
averagearr = averagearr(0:l-1)
err = err(0:l-1)


df = fltarr(l-1)
dfsum = fltarr(l-1)
FOR j = 0,l - 2,1 DO BEGIN

    df(j) =  averagearr(j)- averagearr(j+1)
    dfsum(j) =  averagearr(j) - (1./(l-1.-j))*(total(averagearr(j:l-1.))) 
ENDFOR



;only really want error bars on every tenth point so I can see them on
;the graph, so make smaller arrays with every tenth data point
radnew = fltarr(l)
averagenew = fltarr(l)
errnew = fltarr(l)

i = 0
FOR z =0, l-1,10 DO BEGIN
    radnew(i) = rad(z)
    averagenew(i) = averagearr(z)
    errnew(i) = err(z)
    i = i+1
ENDFOR

radnew = radnew(0:i-1)
averagenew = averagenew(0:i-1)
errnew = errnew(0:i-1)


;curve fitting
;------------------------------------------------------------------------------
err = dindgen(l) - dindgen(l) + 1         ;make a simple error array, all even errors

start = [0.2,7.0]
result = MPFITFUN('devauc',rad(3:maxfit),averagearr(3:maxfit), err, start, bestnorm = chi2)

print, "deVauc result", result, "for the galaxy", xcenter, ycenter
;plot as a check on the profile
!p.multi = [number, 3, 3]
name = strarr(1)
name = string(xcenter+1) + STRING(ycenter+1)

plot, rad, averagearr, title = name, xrange = [0,30], yrange = [0,8],xtitle = 'pixels', ytitle = 'counts/s'
oplot, rad, (result(0)) * (exp(-7.67*(((rad/(result(1)))^(1.0/4.0)) - 1.0))), thick = 3, color = colors.red
name = string(chi2) + string(result(0)) + string (result(1))
xyouts, 0., 7., name, charsize = 0.5, color = colors.red

;make a model of the galaxy and add it to a blank image
;return the galaxy model in array newgal

newgal = subdev(xcenter, ycenter, maxradius, blankdata, result)

data = data - newgal

end
