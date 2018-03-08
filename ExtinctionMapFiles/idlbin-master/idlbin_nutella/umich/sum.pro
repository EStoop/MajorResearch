PRO sum
sum = 0
openr, rlun, "/n/Godiva1/jkrick/A3888/cmd/block_obj.V.input", /get_lun
WHILE (NOT EOF(rlun)) DO BEGIN
     READF, rlun, o, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux
     IF (apflux GT 0) AND (fwhm GT 0) AND ( f GT 0) THEN BEGIN
         IF (fwhm LT 2.8) OR (fwhm GT 4.0) THEN BEGIN
             IF (xcenter LT 682) AND (xcenter GT 648) AND (ycenter LT 933) AND (ycenter GT 851) THEN BEGIN
                 print, "whoops", xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux
             ENDIF

             sum = sum + f
 ;            print, xcenter, ycenter, fwhm, f
         ENDIF
     ENDIF
 ENDWHILE

print, sum
close, rlun
free_lun, rlun
END
