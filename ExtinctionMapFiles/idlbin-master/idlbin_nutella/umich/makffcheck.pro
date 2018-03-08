PRO makeffcheck

close, /all
openw, outlun, "/n/Godiva1/jkrick/A3888/backlist49.large", /get_lun
;;;;;;;;;;;;;;;;;;;;;;;;
;1
;;;;;;;;;;;;;;;;;;;;;;;
FOR x = 10, 290, 50 DO BEGIN
    FOR y = 1051, 2560, 50 DO BEGIN

        printf, outlun, x, x+49, y, y+49
    ENDFOR
ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;
;2
;;;;;;;;;;;;;;;;;;;;;;;
FOR x = 310, 2500, 50 DO BEGIN
    FOR y = 2400, 2600, 50 DO BEGIN

        printf, outlun, x, x+49, y, y+49
    ENDFOR
ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;
;3
;;;;;;;;;;;;;;;;;;;;;;;
FOR x = 1300, 2700, 50 DO BEGIN
    FOR y = 10, 290, 50 DO BEGIN

        printf, outlun, x, x+49, y, y+49
    ENDFOR
ENDFOR
;;;;;;;;;;;;;;;;;;;;;;;;
;4
;;;;;;;;;;;;;;;;;;;;;;;
FOR x = 2400, 2700, 50 DO BEGIN
    FOR y = 400, 2400, 50 DO BEGIN

        printf, outlun, x, x+49, y, y+49
    ENDFOR
ENDFOR

close, outlun
free_lun,outlun


END

