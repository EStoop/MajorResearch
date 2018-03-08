PRO addhead
close, /all
hedit = "hello"

openr, lun, "/n/Godiva4/jkrick/oct98/test", /get_lun
openw, outlun, "/n/Godiva4/jkrick/oct98/test2",/get_lun
WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, hedit
    ra = strmid(hedit, 23,8)
    dec = 23.9
    radec, ra, dec, ihr, imin, xsec, ideg, imn, xsc 
    IF ra GT 0 THEN print, ra, ihr, imin,xsec
    newra = string(ihr)+":"+string(imin)+":"+string(xsec)
    newra = strcompress(newra,/remove_all)
    IF ra GT 0 THEN print, newra

    newhedit = strmid(hedit,0,23) + newra + strmid(hedit, 31,62)
    IF ra GT 0 THEN printf, outlun, newhedit
ENDWHILE


close, lun
close, outlun
free_lun, lun
free_lun, outlun
END
