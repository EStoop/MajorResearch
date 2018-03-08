PRO xyad

fits_read, "/n/Godiva1/jkrick/A4059/original/fullr.wcs.fits", data, header

openr, lun3, "/n/Godiva1/jkrick/A4059/satstar.txt", /GET_LUN
openw, lun2, "/n/Godiva1/jkrick/A4059/extramask.txt",/get_lun

WHILE (NOT eof(lun3)) DO BEGIN
    READF, lun3,rh,rm,rs,dd,dm,ds,rmag,bmag,junk,junk
    r = (rh + ((rm + (rs / 60.)) / 60.))*15.
    d = dd - ((dm + (ds / 60.)) / 60.)
    adxy, header, r, d, xcen, ycen
    IF xcen GT -300 AND xcen LT 3130 AND ycen GT -300 AND ycen LT 3570 THEN begin
        printf,lun2,rmag,bmag,xcen,ycen
    endif
ENDWHILE
close, lun3
free_lun, lun3
close, lun2
free_lun, lun2

END
