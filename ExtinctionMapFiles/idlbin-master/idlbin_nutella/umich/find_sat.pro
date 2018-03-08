pro find_sat


close, /all


FITS_READ, '/n/Godiva3/jkrick/A2556/original/fullB.fits', data, header

;OPENR, lun, '/n/Whalen1/jkrick/ICL/A3888/wholes.cat', /GET_LUN
;read the data file into the variables until hit EOF
openw, lun, '/n/Godiva3/jkrick/A2556/original/satlistB',/get_lun
count = 0 
;WHILE (NOT EOF(lun)) DO BEGIN
;   	READF, lun, o, x, y, a, b, e, f, m, j, fwhm, pa, bkgd

FOR x = 0, 2865, 1 DO BEGIN
    FOR y = 0, 3485, 1 DO BEGIN
	IF (data[x,y] GE 23.0) THEN BEGIN
;		IF (x GE 1004 AND x LE 1014 AND y GE 1085 AND y LE 2494) THEN i = 0 ELSE $; do nothing
;		IF (x GE 1450 AND x LE 1458 AND y GE 2017 AND y LE 3421) THEN i = 0 ELSE $
;		IF (x GE 1944 AND x LE 1952 AND y GE 638  AND y LE 2048) THEN i = 0 ELSE $
;		IF (x GE 2383 AND x LE 2392 AND y GE 1570 AND y LE 2965) THEN i = 0 ELSE $
;		IF (x GE 2845 AND x LE 2856 AND y GE 1355 AND y LE 2750) THEN i = 0 ELSE print, "saturated star @", x, y
            printf,lun, x, y
	ENDIF
;count = count +1	
;ENDWHILE
    ENDFOR
ENDFOR

close, lun
free_lun, lun
END
