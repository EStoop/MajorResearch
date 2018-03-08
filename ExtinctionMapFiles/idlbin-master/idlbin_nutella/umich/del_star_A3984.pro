;need to send pa, increment, a, b, xcenter, ycenter, data

function del_star_A3984, theta, f , a , b, X, Y, data


;figure out what the increment should be 
;based on flux -> size of the object
IF (f GE 30000.0) THEN BEGIN    ;to deal with saturated stars
    increment = (2*!PI)/1000
ENDIF
IF (f LT 30000.0 AND f GT 8000.0) THEN BEGIN
    increment = (2*!PI)/700
ENDIF
IF (f LE 8000.0) THEN BEGIN
    increment =  (2*!PI)/500
ENDIF

;convert theta from degres into radians
theta = theta*(!PI/180)


xmax = 2900;1601;2049;1900;2866;1900;2866;2048;2000;2800;1000
ymax = 3140;1801;3148;2100;3486;2100;3486;3148;2000;2700;1000

;calculate which pixels need to be masked
;to do this change to a new coordinate system (ang)
; 	then work arond the ellipse, each time
;	incrementing angle by some small amount
;	smaller increments for dimmer (smaller) objs.
FOR ang = 0.0, (2*!PI),increment DO BEGIN
    xedge = a*cos(theta)*cos(ang) - b*sin(theta)*sin(ang)
    yedge = a*sin(theta)*cos(ang) + b*cos(theta)*sin(ang)
;   print, "x, y,xedge, yedge", x,y,xedge,yedge

    ; first get the center pixel
    IF (x GT 0 AND x LT xmax) THEN BEGIN
        IF (y GT 0 AND y LT ymax) THEN BEGIN
            data[x,y] = -1E6
                                ;make sure the other pixels are not off the edge
                                ;then mask the little suckers
            IF (xedge+x GT 0 AND xedge+x LT xmax) THEN BEGIN
                IF (yedge+y GT 0 AND yedge+y LT ymax) THEN BEGIN
                                ; have to do four cases to get everything 
                                ;between x and xegde and y and yedge
                    IF ((xedge+x) LT x) THEN BEGIN
                        If (yedge+y LT y) THEN BEGIN
                            data[xedge+x:x,yedge+y:y] = -1000
                        ENDIF ELSE BEGIN
                            data[xedge+x:x,y:yedge+y] = -1000
                        ENDELSE
                    ENDIF ELSE BEGIN
                        If (yedge+y LT y) THEN BEGIN
                            data[x:xedge+x,yedge+y:y] = -1000
                        ENDIF ELSE BEGIN
                            data[x:xedge+x,y:yedge+y] = -1000
                        ENDELSE
                    ENDELSE
                ENDIF
            ENDIF
        ENDIF
    ENDIF
    
ENDFOR


return, data
END













