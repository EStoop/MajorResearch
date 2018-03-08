PRO membership
close, /all
device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)

fits_read, "/n/Godiva6/jkrick/A3880/original/fullr.fits", data, header
;fits_read, "/n/Godiva5/jkrick/xmm.fits", xdata, xheader


numoobjects = 450
;make an array of couch members

galaxy = replicate({object, ra:0D, dec:0D, cz:0D, x:0D, y:0D},numoobjects)
ngalaxy = replicate({object, ra:0D, dec:0D, cz:0D, x:0D, y:0D},numoobjects)
zgalaxy = replicate({object, ra:0D, dec:0D, cz:0D, x:0D, y:0D},numoobjects)

;read in larcs galaxy parameters into galaxy objects
OPENR, lun1, "/n/Godiva6/jkrick/A3880/cmd/2df.A3880.txt", /GET_LUN, ERROR = err
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
    IF (z GT 0.05) AND (z LT 0.065) THEN BEGIN ;member
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

close, lun1
free_lun, lun1

openr, lun3, "/n/Godiva6/jkrick/A3880/cmd/smith.3880.txt", /GET_LUN
WHILE (NOT eof(lun3)) DO BEGIN
    READF, lun3,r,d,junk, cz,junk,junk
    adxy, header, r, d, xcen, ycen
    IF (cz GT 15000) AND (cz LT 19500) THEN BEGIN ;member
        galaxy[i] = {object,r,d,cz,xcen,ycen}
        i = i+ 1
    ENDIF ELSE BEGIN            ;nonmember
        IF cz NE 0.0 THEN BEGIN
            ngalaxy[n] = {object,r,d,cz,xcen,ycen}
            n = n + 1
            
        ENDIF
    ENDELSE
ENDWHILE
close, lun3
free_lun, lun3

openr, lun3, "/n/Godiva6/jkrick/A3880/cmd/collins.3880.member.txt", /GET_LUN
WHILE (NOT eof(lun3)) DO BEGIN
    READF, lun3,junk, rh,rm,rs,dd,dm,ds,junk, cz,junk
    r = (rh + ((rm + (rs / 60.)) / 60.))*15.
    d = dd - ((dm + (ds / 60.)) / 60.)
    adxy, header, r, d, xcen, ycen
    IF (cz GT 15000) AND (cz LT 19500) THEN BEGIN ;member
        galaxy[i] = {object,r,d,cz,xcen,ycen}
        i = i+ 1
    ENDIF ELSE BEGIN            ;nonmember
        IF cz NE 0.0 THEN BEGIN
            ngalaxy[n] = {object,r,d,cz,xcen,ycen}
            n = n + 1
            
        ENDIF
    ENDELSE
ENDWHILE
close, lun3
free_lun, lun3



;want to shorten the array of galaxies 
numogalaxies = i
numogalaxies2 = n
galaxy = galaxy[0:numogalaxies -1]
ngalaxy = ngalaxy[0:n-1]
zgalaxy = zgalaxy[0:c-1]

print, numogalaxies, numogalaxies2

;----------------------------------------------------------------------
ps_open, file = "/n/Godiva6/jkrick/A3880/membership.ps", /portrait, /color, ysize = 7

;RESULT = BYTSCL(data, min = 1.9, max = 2.1)
RESULT = BYTSCL(data, min = -0.5, max = 0.5)
plotimage, result,/PRESERVE_ASPECT


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
END

