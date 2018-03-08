PRO membership
close, /all
;device, true=24
;device, decomposed=0

colors = GetColor(/load, Start=1)
;fits_read, "/n/Godiva1/jkrick/A3888/centerV2000.wcs.div.fits", data, header
fits_read, "/n/Godiva1/jkrick/A3888/A3888V.wcs.div.fits", data, header
fits_read, "/n/Godiva5/jkrick/xmm.fits", xdata, xheader
numoobjects = 219


;make an array of larcs galaxy objects
galaxylarcs = replicate({object, mag:0, magerr:0D, ra:0D, dec:0D, cz:0D, czerr:0D,$
    x:0D, y:0D},numoobjects)

;read in larcs galaxy parameters into galaxy objects
OPENR, lun1, "/n/Godiva1/jkrick/A3888/pimbblet/larcs.cz.txt", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, rmag, rmagerr, r,d,z,zerr
    galaxylarcs[i] = {object, rmag, rmagerr, r, d, z, zerr, 0.0, 0.0}
    i = i +1
ENDWHILE

close, lun1
free_lun, lun1

;want to shorten the array of galaxies 
numogalaxieslarcs = i
galaxylarcs = galaxylarcs[0:numogalaxieslarcs -1]

;turn ra and dec into x and y given the astrometry in the header
FOR i = 0, numogalaxieslarcs - 1, 1 DO BEGIN
    adxy, header, galaxylarcs[i].ra, galaxylarcs[i].dec, xcen, ycen
    galaxylarcs[i].x = xcen
    galaxylarcs[i].y = ycen
ENDFOR

;----------------------------------------------------------------------
;make an array of larcs galaxy objects
galaxylarcs2 = replicate({object, mag:0, magerr:0D, ra:0D, dec:0D, cz:0D, czerr:0D,$
    x:0D, y:0D},numoobjects)

;read in larcs galaxy parameters into galaxy objects
OPENR, lun1, "/n/Godiva1/jkrick/A3888/pimbblet/nonmembers.txt", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, rmag, rmagerr, r,d,z,zerr
    galaxylarcs2[i] = {object, rmag, rmagerr, r, d, z, zerr, 0.0, 0.0}
    i = i +1
ENDWHILE

close, lun1
free_lun, lun1

;want to shorten the array of galaxies 
numogalaxieslarcs2 = i
galaxylarcs2 = galaxylarcs2[0:numogalaxieslarcs2 -1]

;turn ra and dec into x and y given the astrometry in the header
FOR i = 0, numogalaxieslarcs2 - 1, 1 DO BEGIN
    adxy, header, galaxylarcs2[i].ra, galaxylarcs2[i].dec, xcen, ycen
    galaxylarcs2[i].x = xcen
    galaxylarcs2[i].y = ycen
ENDFOR

;----------------------------------------------------------------------
;make an array of teague galaxy objects
galaxyteague = replicate({object, mag:0, magerr:0D, ra:0D, dec:0D, cz:0D, czerr:0D,$
    x:0D, y:0D},numoobjects)

;read in teague galaxy parameters into galaxy objects
OPENR, lun1, "/n/Godiva1/jkrick/A3888/astrometry/galsincluster.deg", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1,r, d
    galaxyteague[i] = {object,0., 0., r, d, 0., 0., 0.0, 0.0}
    i = i +1
ENDWHILE
close, lun1
free_lun, lun1

;want to shorten the array of galaxies 
numogalaxiesteague = i
galaxyteague = galaxyteague[0:numogalaxiesteague -1]

;turn ra and dec into x and y given the astrometry in the header
FOR i = 0, numogalaxiesteague - 1, 1 DO BEGIN
    adxy, header, galaxyteague[i].ra, galaxyteague[i].dec, xcen, ycen
    galaxyteague[i].x = xcen
    galaxyteague[i].y = ycen
ENDFOR

;--------------------------------------------------------------------------
;make an array of teague galaxy objects
galaxyteague2 = replicate({object, mag:0, magerr:0D, ra:0D, dec:0D, cz:0D, czerr:0D,$
    x:0D, y:0D},numoobjects)

;read in teague galaxy parameters into galaxy objects
OPENR, lun1, "/n/Godiva1/jkrick/A3888/astrometry/nonmembers.deg", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1,r, d
    galaxyteague2[i] = {object,0., 0., r, d, 0., 0., 0.0, 0.0}
    i = i +1
ENDWHILE
close, lun1
free_lun, lun1

;want to shorten the array of galaxies 
numogalaxiesteague2 = i
galaxyteague2 = galaxyteague2[0:numogalaxiesteague2 -1]

;turn ra and dec into x and y given the astrometry in the header
FOR i = 0, numogalaxiesteague2 - 1, 1 DO BEGIN
    adxy, header, galaxyteague2[i].ra, galaxyteague2[i].dec, xcen, ycen
    galaxyteague2[i].x = xcen
    galaxyteague2[i].y = ycen
ENDFOR


;-----------------------------------------------------------------------
ps_open, file = "/n/Godiva1/jkrick/A3888/pimbblet/member.ps", /portrait, /color, ysize = 7

;RESULT = BYTSCL(data, min = 1.9, max = 2.1)
RESULT = BYTSCL(data, min = -0.5, max = 0.5)
plotimage, result


a = makex(0, 360, 2)/!radeg
FOR j = 0, numogalaxieslarcs- 1, 1 DO BEGIN
    x = galaxylarcs[j].x
    y = galaxylarcs[j].y 
    r = 20.
    FOR i = 0, n_elements(r)-1 do begin
        xx = x + r(i)*cos(a)
        yy = y + r(i)*sin(a)
        plots,xx,yy, color = colors.red
    ENDFOR
    
ENDFOR
a = makex(0, 360, 2)/!radeg
FOR j = 0, numogalaxieslarcs2- 1, 1 DO BEGIN
    x = galaxylarcs2[j].x
    y = galaxylarcs2[j].y 
    r = 20.
    FOR i = 0, n_elements(r)-1 do begin
        xx = x + r(i)*cos(a)
        yy = y + r(i)*sin(a)
        plots,xx,yy, color = colors.navy
    ENDFOR
    
ENDFOR
FOR j = 0, numogalaxiesteague- 1, 1 DO BEGIN
    x = galaxyteague[j].x
    y = galaxyteague[j].y 
    r = 14.
    FOR i = 0, n_elements(r)-1 do begin
        xx = x + r(i)*cos(a)
        yy = y + r(i)*sin(a)
        plots,xx,yy, color = colors.orange
    ENDFOR
    
ENDFOR
FOR j = 0, numogalaxiesteague2- 1, 1 DO BEGIN
    x = galaxyteague2[j].x
    y = galaxyteague2[j].y 
    r = 14.
    FOR i = 0, n_elements(r)-1 do begin
        xx = x + r(i)*cos(a)
        yy = y + r(i)*sin(a)
        plots,xx,yy, color = colors.blue
    ENDFOR
    
ENDFOR

;center according to r\Reiprich & Bohringer
adxy, header, 338.6255, -37.7343, x, y
print, x, y
r = 30
FOR i = 0, n_elements(r)-1 do begin
    xx = x + r(i)*cos(a)
    yy = y + r(i)*sin(a)
    plots,xx,yy, color = colors.yellow
ENDFOR


;contour, xdata, nlevels = 10,/oplotimage

;xyouts, 1300,-100, "LARCS member galaxy", color = colors.red
;xyouts, 1300,-150, "LARCS non- member ", color = colors.navy
;xyouts, 1300,-200, "Teague member galaxy", color = colors.orange
;xyouts, 1300,-250, "Teague non- member ", color = colors.blue
xyouts, 2500,-100, "LARCS member galaxy", color = colors.red
xyouts, 2500,-250, "LARCS non- member ", color = colors.navy
xyouts, 2500,-400, "Teague member galaxy", color = colors.orange
xyouts, 2500,-550, "Teague non- member ", color = colors.blue
ps_close, /noprint, /noid
END

