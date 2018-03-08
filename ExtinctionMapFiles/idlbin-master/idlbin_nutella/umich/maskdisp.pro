PRO maskdisp

close,/all
mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva2/jkrick/A3984/maskdisp.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

xarr = fltarr(200)
yarr = fltarr(200)
m = 0
;fits_read, "/n/Godiva2/jkrick/A141/A141.wcs.fits", data, header
fits_read, "/n/Godiva2/jkrick/A3984/A3984r.large.fits", data, header
openw, outlun, "/n/Godiva2/jkrick/A3984/gallist",/get_lun
openr,lun, "/n/Godiva7/jkrick/IMACS/masktest", /get_lun
whole = "string"
plotimage, bytscl(data,min=-54,max=60), xrange=[-2200,4900],yrange=[-500,3500]
xyouts, 1200, 3300, "A3984_6"
;colors = GetColor(/load, Start=1)

;A141
;centerra = 16.39804167       
;centerdec = -24.64852222 

;A2556
centerra=348.27500000      
centerdec = -21.56666667

adxy, header, centerra, centerdec, centerx, centery
normal = 0
NS=0
j = 0
WHILE (NOT EOF(lun)) DO BEGIN
    readf,lun,whole
    rh = float(strmid(whole, 0,2))
    rm = float(strmid(whole, 3,2))
    rs = float(strmid(whole, 6,6))

    dd=float(strmid(whole, 13,3))
    dm=float(strmid(whole, 17,2))
    ds=float(strmid(whole, 20,5))

    ra = ((((rs/60.) + rm)/60.) + rh)*15.
    dec = -(((ds/60.) + dm)/60.) + dd
    ADXY, header, ra, dec, x, y

    width = float(strmid(whole, 42, 5))
     IF width GT 1.0 AND width LT 3. THEN w = "s" ELSE  w = " "
    printf, outlun, x,y,"     ",w
    xyouts, x+6, y+6, w, charthick = 3, color = 1
    xarr[m] = x
    yarr[m] = y
    m = m + 1


    ;in a normal orientation, 
    ;how many times does the slit cross a gap?
    FOR k = (fix(y-975.)), fix(y+975.), 1 DO BEGIN
        IF (k LT centery+25.) AND (k GT centery-25) THEN BEGIN
            normal = normal + 1
            k = k + 50
        ENDIF
    ENDFOR

    ;in NS orientation,
    ;how many times does the slit cross a gap?
    FOR k = (fix(y-975.)), fix(y+975.), 1 DO BEGIN
        IF (k LT centery -25-2000) AND (k GT centery - 25 - 2000 -50) THEN BEGIN
            NS = NS + 1
            k = k + 50
        ENDIF

        IF (k LT centery +25) AND (k GT centery - 25) THEN BEGIN
            NS = NS + 1
            k = k + 50
        ENDIF

        IF (k LT centery +25) AND (k GT centery - 25) THEN BEGIN
            NS = NS + 1
            k = k + 50
        ENDIF

    ENDFOR


    j = j + 1  ;total number of slits

   ;length of the slit, irrelevant for now
    l = float(strmid(whole, 48,5))

ENDWHILE

xarr = xarr[0:m-1]
yarr = yarr[0:m-1]

oplot, xarr, yarr, psym = 6, thick = 3, color = 1
    avnormal =float(normal)/ float(j)
    print, "average gap crossing in normal config = ", avnormal

    avNS = float(NS) / float(j)
    print, "average gap crossing in NS config = ", avNS


close, lun
free_lun, lun
close, outlun
free_lun, outlun
device, /close
set_plot, mydevice

END
