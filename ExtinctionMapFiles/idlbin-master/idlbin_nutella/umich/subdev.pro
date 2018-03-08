function subdev, xcenter, ycenter, maxradius, blankdata, constants, PA, e
;this function acutally makes the homemade galaxy and adds it to a blank
;image
radius = 0.0		

intx = fix(xcenter)
inty = fix(ycenter)

dx = xcenter-intx
dy = ycenter-inty
a = 0.0                         ;semi-major axis 
;e = 0.0                         ;ellipticity
phi = 0.0                       ;angle measured from positive x axis
;PA = PA                         ;position angle
print, "xcenter, ycenter", xcenter, ycenter
print, "inside subdev, pa = ", pa

;whooops, only PA is measured from positive y-axis,,,so
PA = PA - 90
;make position angle positive
IF (PA LT 0) THEN PA = PA + 360

;error check
IF (PA GT 360) THEN print, "shouldn't be happening"

;convert position angle from degrees to radians
PA = PA / !RADEG


IF (PA GE !DPi/2 AND PA LT 3*!DPi/2) THEN PA = PA - !DPI

print, "after manipulation pa = ", pa
blank1 = blankdata
blank1(*,*) = 0.0

;for each (x,y) within some radius of the center...
FOR y = inty-maxradius, inty + maxradius,1 DO BEGIN
    FOR x = intx -maxradius, intx + maxradius,1 DO BEGIN
       
        xdist = x - xcenter  ; length of side on x-axis
        ydist = y - ycenter  ; length of side on y-axis

        ;find phi = angle as measured CCW from the positive x axis
        alpha = atan(ydist / xdist) ;= angle measured from nearest x-axis
       
        IF (xdist LT 0.) THEN phi = !DPI + alpha 
        
        IF (xdist GT 0.) THEN  BEGIN
            IF (alpha GE 0.) THEN phi = alpha
            IF (alpha LT 0.) THEN phi = 2*!DPI + alpha 
        ENDIF
        
        IF (xdist EQ 0.) THEN BEGIN
            IF (ydist GT 0.) THEN phi = !DPI/2
            IF (ydist LT 0.) THEN phi = 3*!DPI/2
            IF (ydist EQ 0.) THEN phi = 0. ; central pixel,don't care
        ENDIF
        
        IF (phi LT 0) THEN BEGIN
            Phi = Phi + 2*!DPI
            print, "interesting"
        ENDIF

        ;find phiprime = angle as measured CCW from the PA
        IF (phi LT PA) THEN phiprime = 2*!DPI + phi - PA
        IF (phi GE PA) THEN phiprime = phi - PA

        ;for reference parametric equations of an ellipse w/ rotation
        ;x = xcenter + a*cos(phiprime)*cos(PA) - a*(1-e)*sin(phiprime)*sin(PA) 
        ;y = ycenter + a*cos(phiprime)*sin(PA) + a*(1-e)*sin(phiprime)*cos(PA) 
        ;definition of ellipticity => e = 1 - b/a

        ;solve the above equations for a
        ax = x-xcenter/(cos(phiprime)*cos(PA) - sin(phiprime)*sin(PA) + sin(phiprime)*sin(PA)*e)
        ay = y-ycenter/(cos(phiprime)*sin(PA) + sin(phiprime)*cos(PA) - sin(phiprime)*cos(PA)*e)
        a = ax
        IF(x EQ 0) THEN a = ay
        IF(ax NE ay) THEN print, "ax ne ay", ax, ay
        

        print,PA, e, x, y, phiprime, a
;        intensity = devauc(a, constants)
            
;        blank1[x,y] = blank1[x,y] + intensity

    ENDFOR
ENDFOR
;blank1 = shiftf(blank1, dx,dy)
;blankdata = blankdata + blank1

return, blankdata             
END


;FOR y = inty-maxradius, inty + maxradius,1 DO BEGIN
;    FOR x = intx -maxradius, intx + maxradius,1 DO BEGIN
;        ;find what the radius is from the center
;radius = sqrt((intx - x)^2 + (inty - y)^2)
        
;        intensity = devauc(radius, constants)
;        blank1[x,y] = blank1[x,y] + intensity
        
;    ENDFOR
;ENDFOR

;two equations, two unknowns, but it is not working

        ;solving these for a, and e, but only care about 'a' here
;        IF (phiprime NE !DPI/2 AND phiprime NE 3*!DPi/2) THEN begin
;            a = -(cos(PA)*x-cos(PA)*xcenter+y*sin(PA)-ycenter*sin(PA))*cos(phiprime)/(sin(phiprime)^2-1)

;            e = (sin(phiprime)*cos(PA)*x-sin(phiprime)*cos(PA)*xcenter+sin(phiprime)*sin(PA)*y-sin(phiprime)*sin(PA)*ycenter+cos(phiprime)*sin(PA)*x-cos(phiprime)*sin(PA)*xcenter-cos(phiprime)*y*cos(PA)+cos(phiprime)*ycenter*cos(PA))/sin(phiprime)/(cos(PA)*x-cos(PA)*xcenter+y*sin(PA)-ycenter*sin(PA))
        ;need to deal with sin(phiprimeprime) = 1 =>phiprime = Pi/2 multiples
        ;maybe other equation?
        
;            print,e, x, y, phiprime, a
;            intensity = devauc(a, constants)
            
;           blank1[x,y] = blank1[x,y] + intensity
;        ENDIF
