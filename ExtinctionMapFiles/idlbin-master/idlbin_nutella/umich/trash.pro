PRO trash

openr, lun, "/n/Godiva4/jkrick/oct98/oct21/addhead.cl", /get_lun
hedit = "hello"
openw, outlun, "/n/Godiva4/jkrick/oct98/oct21/test2",/get_lun

FOR i = 0, 83, 1 DO BEGIN
    READF, lun, hedit
    hr = strmid(hedit, 34,2)
    min = strmid(hedit, 37,2)
    sec = strmid(hedit, 40, 2)
    print, hr, min, sec
 ;   min = fix(min) + 3
 ;   sec = fix(sec) + 43

    IF sec GT 60 THEN BEGIN
        min = min + 1
        sec = sec - 60
    ENDIF

    IF min GT 60 THEN BEGIN
        hr = fix(hr) + 1
        min = min - 60
    ENDIF
    hr = hr - 4
    IF hr LT 0 THEN hr = 24 + hr

    print, hr, min, sec
    newtime =string(hr)+":"+string(min)+":"+string(sec)
    newtime = strcompress(newtime,/remove_all)

    newhedit = strmid(hedit,0,34) + newtime + strmid(hedit, 42,65)

    printf,outlun, newhedit

ENDFOR



close, lun
close, outlun
free_lun, lun
free_lun, outlun

END
