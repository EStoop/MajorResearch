function subdev2, xcenter, ycenter, maxradius, blankdata, constants, PA, ellipticity
;this function acutally makes the homemade galaxy and adds it to a blank image

;position angle in radians
r = 0.0		
a = 0.0                         ;semi-major axis 
;e = 0.0                         ;ellipticity
phi = 0.0                       ;angle measured from positive x axis
;PA = PA                         ;position angle
intx = fix(xcenter)
inty = fix(ycenter)
temppa = 0.0

;for each (x,y) within some radius of the center...
FOR y = inty-maxradius, inty + maxradius,1 DO BEGIN
    FOR x = intx -maxradius, intx + maxradius,1 DO BEGIN

        ;make sure I am not off the edge of the image
        ;need to change this to be more versatile
        IF (x LT 0) OR (x GE 100) OR (y LT 0) OR (y GE 100) THEN BEGIN

        ENDIF ELSE BEGIN
            ;pythagoras
            xprime = abs(x - xcenter)
            yprime = abs(y - ycenter)
            r = sqrt(xprime^2+yprime^2)
            IF (xprime EQ 0) THEN alpha = 0
            IF (xprime NE 0) THEN alpha = atan(yprime / xprime) ;polar angle
            

            ;fold in the position angle
            ;want angle as measured ccw from the semi-major axis
            phiprime = alpha - pa
            
            a = sqrt(r^2.*(sin(phiprime)^2. + (1.-ellipticity)^2.*cos(phiprime)^2.)/(1.-ellipticity)^2.)
            
            print, "x, y, xprime, yprime,alpha, phiprime, a, constants " ,x, y, xprime, yprime,alpha,phiprime, a, constants
            
           intensity = devauc(a, constants)
    ;       intensity = exponential(a, constants)
     ;       intensity = 10
            blankdata[x,y] = blankdata[x,y] + intensity
        ENDELSE

    ENDFOR
ENDFOR
;blankdata[intx,inty] = blankdata[intx+1,inty] ;try getting rid of
;central super-peak
return, blankdata  

END

        ;want the angle as measured ccw from the pos y axis
        ;IF (xprime(i) GT 0) THEN theta = alpha + 3*!DPi/2
        ;IF (xprime(i) LT 0) THEN theta = alpha + !DPi/2

  
