PRO membership
close, /all
device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)

fits_read, "/n/Godiva7/jkrick/A114/original/larger.wcs.fits", data, header
;fits_read, "/n/Godiva5/jkrick/xmm.fits", xdata, xheader


numoobjects = 450
;make an array of couch members

galaxy = replicate({object, ra:0D, dec:0D, cz:0D, x:0D, y:0D},numoobjects)
;numoobjects = 450
;make an array of couch members

;galaxy = replicate({object, ra:0D, dec:0D, cz:0D, x:0D, y:0D},numoobjects)
ngalaxy = replicate({object, ra:0D, dec:0D, cz:0D, x:0D, y:0D},numoobjects)
zgalaxy = replicate({object, ra:0D, dec:0D, cz:0D, x:0D, y:0D},numoobjects)

;read in larcs galaxy parameters into galaxy objects
OPENR, lun1, "/n/Godiva7/jkrick/A114/cmd/2df.A114.txt", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

i = 0
n = 0
c=0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1,rh,rm,rs,dd,dm,ds,junk,z,junk
    r = (rh + ((rm + (rs / 60.)) / 60.))*15.
    d = dd - ((dm + (ds / 60.)) / 60.)
    adxy, header, r, d, xcen, ycen
    s = sqrt((xcen - 860.)^2 + (ycen-2021.)^2)
    IF s LT 100 THEN print, i, s,z
    IF (z GT 0.3) AND (z LT 0.35) THEN BEGIN ;member
        galaxy[i] = {object, r, d, z, xcen, ycen}
        i = i +1
    ENDIF ELSE BEGIN            ;nonmember
        IF (z NE 0.0) AND (z NE -9.0) THEN BEGIN
            ngalaxy[n] = {object, r, d, z, xcen, ycen}
            n = n + 1
            
        ENDIF ELSE BEGIN
            zgalaxy[c] = {object, r, d, z, xcen, ycen}
            c = c+1
        ENDELSE


    ENDELSE

ENDWHILE


;read in larcs galaxy parameters into galaxy objects
OPENR, lun1, "/n/Godiva7/jkrick/A114/astrometry/couch.member", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

;i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, rh, rm, rs, dd, dm, ds, z
    ;turn hours min sec, deg min sec into degrees, degrees
    r = (rh + ((rm + (rs / 60.)) / 60.))*15.
    d = dd - ((dm + (ds / 60.)) / 60.)
    adxy, header, r, d, xcen, ycen
    print, r, d
    galaxy[i] = {object, r, d, z, xcen, ycen}
    i = i +1
ENDWHILE

close, lun1
free_lun, lun1

;want to shorten the array of galaxies 

;turn ra and dec into x and y given the astrometry in the header
;FOR i = 0, numogalaxies - 1, 1 DO BEGIN
;    adxy, header, galaxy[i].ra, galaxy[i].dec, xcen, ycen
;    galaxy[i].x = xcen
;    galaxy[i].y = ycen
;ENDFOR

;----------------------------------------------------------------------
;make an array of couch nonmembers
;ngalaxy = replicate({object, ra:0D, dec:0D, cz:0D, x:0D, y:0D},numoobjects)

;read in larcs galaxy parameters into galaxy objects
OPENR, lun1, "/n/Godiva7/jkrick/A114/astrometry/couch.nonmember", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

;i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, rh, rm, rs, dd, dm, ds, z
    ;turn hours min sec, deg min sec into degrees, degrees
    r = (rh + ((rm + (rs / 60.)) / 60.))*15.
    d = dd - ((dm + (ds / 60.)) / 60.)
    adxy, header, r, d, xcen, ycen
    print, r, d
    ngalaxy[n] = {object, r, d, z, xcen, ycen}
    n = n +1
ENDWHILE

close, lun1
free_lun, lun1

;want to shorten the array of galaxies 
numogalaxies2 = n
ngalaxy = ngalaxy[0:numogalaxies2 -1]
numogalaxies = i
galaxy = galaxy[0:numogalaxies -1]
zgalaxy=zgalaxy[0:c-1]

;turn ra and dec into x and y given the astrometry in the header
;FOR i = 0, numogalaxies2 - 1, 1 DO BEGIN
;    adxy, header, ngalaxy[i].ra, ngalaxy[i].dec, xcen, ycen
;    ngalaxy[i].x = xcen
;    ngalaxy[i].y = ycen
;ENDFOR




;----------------------------------------------------------------------

;-----------------------------------------------------------------------
ps_open, file = "/n/Godiva7/jkrick/A114/membershiptest.ps", /portrait, /color, ysize = 7

;RESULT = BYTSCL(data, min = 1.9, max = 2.1)
RESULT = BYTSCL(data, min = -0.5, max = 0.5)
plotimage, result


a = makex(0, 360, 2)/!radeg
FOR j = 0, numogalaxies- 1, 1 DO BEGIN
    x = galaxy[j].x
    y = galaxy[j].y 
    r = 20.
    FOR i = 0, n_elements(r)-1 do begin
        xx = x + r(i)*cos(a)
        yy = y + r(i)*sin(a)
        plots,xx,yy, thick = 3, color = colors.red
    ENDFOR
    
ENDFOR
a = makex(0, 360, 2)/!radeg
FOR j = 0, numogalaxies2- 1, 1 DO BEGIN
    x = ngalaxy[j].x
    y = ngalaxy[j].y 
    r = 25.
    FOR i = 0, n_elements(r)-1 do begin
        xx = x + r(i)*cos(a)
        yy = y + r(i)*sin(a)
        plots,xx,yy, thick = 3, color = colors.navy
    ENDFOR
    
ENDFOR
a = makex(0, 360, 2)/!radeg
FOR j = 0, c- 1, 1 DO BEGIN
    x = zgalaxy[j].x
    y = zgalaxy[j].y 
    r = 30.
    FOR i = 0, n_elements(r)-1 do begin
        xx = x + r(i)*cos(a)
        yy = y + r(i)*sin(a)
        plots,xx,yy, thick = 3, color = colors.green
    ENDFOR
    
ENDFOR

;center according to r\Reiprich & Bohringer
;adxy, header, 338.6255, -37.7343, x, y
;print, x, y
;r = 30
;FOR i = 0, n_elements(r)-1 do begin
;    xx = x + r(i)*cos(a)
;    yy = y + r(i)*sin(a)
;    plots,xx,yy, color = colors.yellow
;ENDFOR


;contour, xdata, nlevels = 10,/oplotimage

;xyouts, 1300,-100, "LARCS member galaxy", color = colors.red
;xyouts, 1300,-150, "LARCS non- member ", color = colors.navy
;xyouts, 1300,-200, "Teague member galaxy", color = colors.orange
;xyouts, 1300,-250, "Teague non- member ", color = colors.blue
;xyouts, 2500,-100, "LARCS member galaxy", color = colors.red
;xyouts, 2500,-250, "LARCS non- member ", color = colors.navy
ps_close, /noprint, /noid
undefine, galaxy
undefine, ngalaxy
undefine, zgalaxy
END

