PRO count_flux
close,/all

fits_read, "/n/Godiva4/jkrick/A114/original/larger.bblock.fits", data, header
inner = 0.
innersum = 0.D
outer = 0.
outersum = 0.
outer2 = 0.
outersum2 = 0.
sum = 0.D
npixinner = 0.D
npixouter = 0.D
npixouter2 = 0.D
xcenter1 = 995.
ycenter1 = 1110.
FOR x = 400, 1600,1 DO BEGIN
    FOR y= 500, 17000, 1 DO BEGIN
        dist1 = sqrt((xcenter1 - x)^2 + (ycenter1 - y)^2)
;        dist2 = sqrt((xcenter2 - x)^2 + (ycenter2 - y)^2)
        IF (dist1 LE 150)  THEN BEGIN
            npixinner = npixinner + 1
            IF (data[x,y] GT -1) AND (data[x,y] NE 0) THEN begin        
                inner = inner + 1
                innersum = innersum + data[x,y] 
            ENDIF
        ENDIF
        
        IF (dist1 GT 150) AND (dist1 LT 300) THEN BEGIN
            npixouter = npixouter + 1
            IF (data[x,y] GT -1) AND (data[x,y] NE 0) THEN BEGIN
                outer = outer +1
                outersum = outersum + data[x,y] 
            ENDIF
        ENDIF

        IF (dist1 GT 300) AND (dist1 LT 500) THEN BEGIN
            npixouter2 = npixouter2 + 1
            IF (data[x,y] GT -1) AND (data[x,y] NE 0) THEN BEGIN
                outer2 = outer2 +1
                outersum2 = outersum2 + data[x,y] 
            ENDIF
        ENDIF

;        IF data[x,y] GT 0 THEN BEGIN
;            sum = sum +data[x,y] -2.
;            npix = npix + 1
         ;   print, data[x,y]
;        ENDIF

    ENDFOR
ENDFOR
print,"r inner",inner, innersum, innersum/inner, 24.6 - 2.5*alog10((innersum/inner)/(0.259^2))
print, "should be npix: ", npixinner, inner/npixinner,"%"
print, "new sum = " , innersum/(inner/npixinner);*(npix - innersum)

print,"r outer",outer, outersum, outersum/outer, 24.6 - 2.5*alog10((outersum/outer)/(0.259^2))
print, "should be npix: ", npixouter, outer/npixouter,"%"
print, "new sum = " , outersum/(outer/npixouter);*(npix - outersum)


print,"r outer2",outer2, outersum2, outersum2/outer2, 24.6 - 2.5*alog10((outersum2/outer2)/(0.259^2))
print, "should be npix: ", npixouter2, outer2/npixouter2,"%"
print, "new sum = " , outersum2/(outer2/npixouter2);*(npix - outersum)

END
