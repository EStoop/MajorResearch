PRO counter

openw, lun, "/n/Godiva6/jkrick/A3880/backchecksmooth",/get_lun

FOR i = 570,1000, 13 DO BEGIN
    FOR j = 580,3760, 13 DO BEGIN
        printf, lun, strcompress("fullr.bblock.med75["+string(i)+":"+string(i+4)+","  + string(j)+":"+string(j+4)+"]")
    ENDFOR
ENDFOR

FOR i = 530,3110, 13 DO BEGIN
    FOR j = 3760,4200, 13 DO BEGIN
        printf, lun, strcompress("fullr.bblock.med75["+string(i)+":"+string(i+4)+","  + string(j)+":"+string(j+4)+"]")
    ENDFOR
ENDFOR

FOR i = 2660,3110, 13 DO BEGIN       ;2400, 2700
    FOR j = 580, 3760, 13 DO BEGIN            ;1900,2700
        printf, lun, strcompress("fullr.bblock.med75["+string(i)+":"+string(i+4)+","  + string(j)+":"+string(j+4)+"]")
    ENDFOR
ENDFOR

FOR i = 1000,2613, 13 DO BEGIN       ;2400, 2700
    FOR j = 580,880, 13 DO BEGIN            ;1900,2700
        printf, lun, strcompress("fullr.bblock.med75["+string(i)+":"+string(i+4)+","  + string(j)+":"+string(j+4)+"]")
    ENDFOR
ENDFOR


close, lun
free_lun, lun


END
