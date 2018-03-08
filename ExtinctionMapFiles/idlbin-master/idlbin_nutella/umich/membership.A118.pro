PRO membership
close, /all
device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)

fits_read, "/n/Godiva4/jkrick/A118/original/largeV.wcs.fits", data, header
;fits_read, "/n/Godiva5/jkrick/xmm.fits", xdata, xheader
;data = reverse(data,2)

numoobjects = 328
;make an array of couch members

galaxy = replicate({object, ra:0D, dec:0D, cz:0D, x:0D, y:0D},numoobjects)

;read in larcs galaxy parameters into galaxy objects
OPENR, lun1, "/n/Godiva4/jkrick/A118/astrometry/member.txt", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, rh, rm, rs, dd, dm, ds, z,z,z
    ;turn hours min sec, deg min sec into degrees, degrees
    r = (rh + ((rm + (rs / 60.)) / 60.))*15.
    d = dd - ((dm + (ds / 60.)) / 60.)
    print, r, d
    galaxy[i] = {object, r, d, z, 0.0, 0.0}
    i = i +1
ENDWHILE

close, lun1
free_lun, lun1

;want to shorten the array of galaxies 
numogalaxies = i
galaxy = galaxy[0:numogalaxies -1]

;turn ra and dec into x and y given the astrometry in the header
FOR i = 0, numogalaxies - 1, 1 DO BEGIN
    adxy, header, galaxy[i].ra, galaxy[i].dec, xcen, ycen
    galaxy[i].x = xcen
    galaxy[i].y = ycen
ENDFOR

;----------------------------------------------------------------------
;make an array of couch nonmembers
ngalaxy = replicate({object, ra:0D, dec:0D, cz:0D, x:0D, y:0D},numoobjects)

;read in larcs galaxy parameters into galaxy objects
OPENR, lun1, "/n/Godiva4/jkrick/A118/astrometry/nonmember.txt", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, rh, rm, rs, dd, dm, ds, z,z,z
    ;turn hours min sec, deg min sec into degrees, degrees
    r = (rh + ((rm + (rs / 60.)) / 60.))*15.
    d = dd - ((dm + (ds / 60.)) / 60.)
    print, r, d
    ngalaxy[i] = {object, r, d, z, 0.0, 0.0}
    i = i +1
ENDWHILE

close, lun1
free_lun, lun1

;want to shorten the array of galaxies 
numogalaxies2 = i
ngalaxy = ngalaxy[0:numogalaxies2 -1]

;turn ra and dec into x and y given the astrometry in the header
FOR i = 0, numogalaxies2 - 1, 1 DO BEGIN
    adxy, header, ngalaxy[i].ra, ngalaxy[i].dec, xcen, ycen
    ngalaxy[i].x = xcen
    ngalaxy[i].y = ycen
ENDFOR

;----------------------------------------------------------------------

;-----------------------------------------------------------------------
ps_open, file = "/n/Godiva4/jkrick/A118/membership.ps", /portrait, /color, ysize = 7

;RESULT = BYTSCL(data, min = 1.9, max = 2.1)
RESULT = BYTSCL(data[300:1900,400:2000], min = -0.1, max = 0.1)
plotimage, result;, imgyrange = [1600,1]


a = makex(0, 360, 2)/!radeg
FOR j = 0, numogalaxies- 1, 1 DO BEGIN
    x = galaxy[j].x -300.
    y = galaxy[j].y -400.
    r = 10.
    FOR i = 0, n_elements(r)-1 do begin
        xx = x + r(i)*cos(a)
        yy = y + r(i)*sin(a)
        plots,xx ,yy, thick = 3, color = colors.red
    ENDFOR
    
ENDFOR
a = makex(0, 360, 2)/!radeg
FOR j = 0, numogalaxies2- 1, 1 DO BEGIN
    x = ngalaxy[j].x - 300.
    y = ngalaxy[j].y  -400.
    r = 12.
    FOR i = 0, n_elements(r)-1 do begin
        xx = x + r(i)*cos(a)
        yy = y + r(i)*sin(a)
        plots,xx,yy, thick = 3, color = colors.navy
    ENDFOR
    
ENDFOR

ps_close, /noprint, /noid

undefine, galaxy
undefine, ngalaxy
END

