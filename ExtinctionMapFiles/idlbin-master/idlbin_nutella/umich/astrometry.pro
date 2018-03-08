PRO astrometry
close,/all

fits_read, "/n/Godiva4/jkrick/A118/original/largeV.wcs.fits", data, header

openr, lun, "/n/Godiva4/jkrick/A118/junk", /get_lun
;openw, outlun, "/n/Godiva4/jkrick/A118/2mass.xy.tbl", /get_lun
ra = 0.D
dec = 0.D

WHILE (NOT EOF(lun)) DO BEGIN
;    READF, lun, ra, dec, junk
;    ADXY, header, ra,dec, x, y
;    printf, outlun, ra, dec,x+1,y+1
    READF, lun,rh,rm,rs,dd,dm,ds
;    print, rh, rm,rs,dd,dm,ds,z,z,z
    r = (rh + ((rm + (rs / 60.)) / 60.))*15.
    d = dd - ((dm + (ds / 60.)) / 60.)
    adxy, header, r, d, xcen, ycen
    print, r, d, xcen+1, ycen+1
ENDWHILE


close, lun
free_lun, lun
;close, outlun
;free_lun, outlun
END




;fits_read,"/n/Godiva2/jkrick/A141/original/A141.fits", data, header

;refernece stars are on A3888V.wcs, x and y come from largeV.fits
;ra = [16.4683,16.33277,16.40428]
;dec = [-24.70576,-24.64125,-24.60788]
;x = [1187.36,2906.98,2016.67]
;y = [2796.93,1939.26,1458.95]
;ra = [ 16.38131 ,16.47215 ,16.33277 ]
;dec = [-24.67585,-24.6354,-24.64125]
;x = [2289.66 ,1154.41 ,2907.03 ]
;y = [2405.37 ,1824.05,1939.22]
;starast, ra,dec,x,y,cd, hdr = header

;print, cd
;fits_write, "/n/Godiva2/jkrick/A141/original/A141.wcs.fits", data, header

