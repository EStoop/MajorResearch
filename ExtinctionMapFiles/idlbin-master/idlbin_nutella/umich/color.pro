PRO color

close, /all

numoobjects = 45000
rgalaxy = replicate({object, xcenter:0D,ycenter:0D,fluxa:0D,maga:0D,fluxb:0D,magb:0D,fwhm:0D,isoarea:0D},numoobjects)
mgalaxy = replicate({what, xcen:0D,ycen:0D, vr:0D,rmag:0D, member:0D},numoobjects)

i = 0

openr, lun, "/n/Godiva1/jkrick/A3888/SExtractor.r.cat", /GET_LUN
WHILE (NOT eof(lun)) DO BEGIN
    READF,lun,junk, rx, ry, fluxaper,magaper,erraper,fluxbest,magbest,errbest,fwhm,iso
    IF (fwhm GT 5.1) AND (fluxbest GT 0.0) AND (iso GT 6.0) THEN BEGIN   ;no bother including the stars and junk
        rgalaxy[i] = {object,  rx, ry, fluxaper,magaper,fluxbest,magbest,fwhm,iso}
        i = i + 1
    ENDIF
ENDWHILE

rgalaxy = rgalaxy[0:i - 1]
print, "there are i,",i," r galaxies"
close, lun
free_lun, lun
rflux = 0.D
vflux = 0.D
xdiff = 0D
ydiff = 0.D
vx = 0D
vy = 0D
bgal = 0
openw, nolun, "/n/Godiva1/jkrick/A3888/color.txt",/get_lun
openr, lun2, "/n/Godiva1/jkrick/A3888/SExtractor.V.cat", /GET_LUN
WHILE (NOT eof(lun2)) DO BEGIN
    READF, lun2,junk, vx, vy, vfluxaper,vmagaper,verraper,vfluxbest,vmagbest,verrbest,vfwhm, vra,vdec,viso

    IF (vfwhm GT 5.1) AND (vfluxbest GT 0.0) AND (viso GT 6.0)THEN BEGIN     ;not star or junk

        ;find this object in the r image
        FOR j = 0, i-1, 1 DO begin
            IF (vx LT rgalaxy[j].xcenter + 0.6) AND  (vx GT rgalaxy[j].xcenter - 0.6) AND $
                 (vy LT rgalaxy[j].ycenter + 0.6) AND (vy GT rgalaxy[j].ycenter - 0.6)  THEN BEGIN
                limithigh = (-0.0248)*rgalaxy[j].magb + 0.7067+ 0.12 
                limitlow = (-0.0248)*rgalaxy[j].magb + 0.7067- 0.12 
  
                vr = vmagbest - rgalaxy[j].magb
                IF (vr LT limithigh )AND (vr GT limitlow) THEN BEGIN
                                ;galaxy "a" makes the color cut
                   printf, nolun, vx, vy, vmagbest, rgalaxy[j].magb
                    
                   ;keep track of how much flux is in the background annulus
                   xdiff = (abs(vx - 1936))^2
                   ydiff = (abs(vy - 2224)) ^2
                   IF (sqrt(xdiff + ydiff) LT 1800.) AND (sqrt(xdiff + ydiff) GT 1550.) THEN BEGIN
                       print, vx, vy, vfluxbest, rgalaxy[j].fluxb
                       rflux = rflux + rgalaxy[j].fluxb
                       vflux = vflux + vfluxbest
                       bgal = bgal + 1
                   ENDIF


                ENDIF ELSE BEGIN
                    
                                ;galaxy"a" does not make the color cut
                                ;will want to mask these galaxies to
                                ;determine a cluster light profile
                                ;       printf, nolun, mgalaxy[a].xcen,mgalaxy[a].ycen, mgalaxy[a].rmag
                ENDELSE
                

                
            ENDIF
        ENDFOR
    ENDIF
ENDWHILE

close,nolun
free_lun, nolun
close,lun2
free_lun, lun2

print, "flux in member galaxies in r", rflux, "num of gals ", bgal
print, "flux in member galaxies in V", vflux


;now find out how many pixels are within the background annulus 
;fits_read, "/n/Godiva1/jkrick/A3888/V.div.fits", Vdata, Vheader
;s = 0.D
;r = 0.D
;nx = 0.D
;ny = 0.D
;FOR x = 0.D, 4094.D, 1.D DO BEGIN
;    FOR y = 0.D, 4094.D, 1.D DO BEGIN
;        nx = abs(x-1936.)
;        ny = abs(y-2224.)
;        r = sqrt((nx)^2 + (ny)^2)
;        IF (r GT 1550.0D) AND (r LT 1800.0D)  AND (Vdata[x,y] NE 0.0D) THEN s = s + 1.D
;   ENDFOR
;ENDFOR

;fits_read, "/n/Godiva1/jkrick/A3888/r.div.fits", Vdata, Vheader
;t = 0.D
;FOR x = 0.D, 4094, 1 DO BEGIN
;    FOR y = 0.D, 4094, 1 DO BEGIN
;        nx = abs(x-1936.)
;        ny = abs(y-2224.)
;        r = sqrt((nx)^2 +(ny)^2)
;        IF (r GT 1550.) AND (r LT 1800.) AND (Vdata[x,y] NE 0.0) THEN t = t+ 1
;    ENDFOR
;ENDFOR

;print, "V: ", s, " r: ", t


END

