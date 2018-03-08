PRO count_pix

close, /all

fits_read, "/n/Godiva2/jkrick/A141/original/eaperture3.V.fits", data, header
count = 0.D
sum = 0.
histarr = fltarr(160000)

FOR x = 0, 2900, 1 DO BEGIN
    FOR y = 0, 2900, 1 DO BEGIN
        IF data[x,y] GT -100 AND data[x,y] NE 0 THEN BEGIN
            count = count+ 1
            sum = sum + data[x,y]
        ENDIF

    ENDFOR
ENDFOR

print, count, sum, sum/count, 24.3 - 2.5*alog10((sum/count)/(0.259^2))

END
