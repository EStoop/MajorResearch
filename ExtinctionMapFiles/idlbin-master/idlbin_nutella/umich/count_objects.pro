PRO count_objects
openr, lun, "/n/Godiva1/jkrick/A3888/final/SExtractor2.r.cat", /get_lun
openw, lunw, "/n/Godiva1/jkrick/A3888/final/tvmarkobj", /get_lun
i = 0
WHILE (NOT EOF(lun)) DO BEGIN
    ;!!!! change this if you change the daofind.param , now usign godiva7/galfit/dao..
    READF, lun, o, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux,isocorflux,isoflux
 
    dist = sqrt((404. - xcenter)^2 + (426.-ycenter)^2)
    IF (dist LT 224) AND (24.6 - 2.5*alog10(isocorflux) LE 24.5) AND fwhm GT 5.1 AND isoarea GT 6 THEN BEGIN
        printf,lunw, xcenter, ycenter, m
        i = i + 1
    ENDIF
    
    
ENDWHILE
print, "number of objects", i

close, lun
close, lunw
END
