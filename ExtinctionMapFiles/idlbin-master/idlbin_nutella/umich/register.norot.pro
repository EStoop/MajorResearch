;+
; NAME: Jessica Krick
;
; PURPOSE: this script creates  matched star lists as input to geomap
;
; CALLING SEQUENCE: registernorot  (norotation)
;
; INPUTS: 1) list of files to be registered (should already be set into 4000x4000 images)
;                   2) fwhmcutoffs as a measure of the seeing so I know where to find the stars
;                        I could get around this by choosing a constant, but this is not as efficient
;                       since I would have to choose a large range
;
; OUTPUTS:   ccd????.big.match  which goes right into geomap
;
; PROCEDURE:  takes ~1-2 minutes, longer if there are more stars to check against
;
; MODIFICATION HISTORY:  July 2004
;
;-

FUNCTION equalto, param1, param2, comfortlevel

IF (param1 LT param2 + comfortlevel) AND (param1 GT param2 -comfortlevel) THEN BEGIN
    return, 1
ENDIF ELSE BEGIN
    return, 0
ENDELSE

END
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO registernorot

close, /all

;openr, lun, "/n/Godiva2/jkrick/A3984/phot_ccd8050.dat", /get_lun
;this is the file with the template star locations.

;for now I have chosen three stars in the very center by hand
;for ccd2075;A141
;ax =   2122.03 
;ay = 2093.68
;bx = 2222.89 
;by = 1861.07
;cx =  2045.88 
;cy = 1685.05

;FOR ccd8058;A141
;ax=1277.3
;ay=629.2
;bx=1100.35 
;by= 453.40 
;cx=1176.53  
;cy=861.70

;for 8050;A3984
;ax=1131.194 
;ay=1127.268
;bx=1314.008   
;by=1172.899
;cx= 1246.851  ;  1409.085
;cy=549.207 ;734.332

;A3984 second try
;ax=916.987 
;ay=1323.896
;bx=786.931 
;by=1371.524
;cx=931.064 
;cy=1198.204

;A2556
;ax = 1248.856 
;ay = 1582.925
;bx = 1312.284 
;by=1282.639
;cx=928.858 
;cy=1148.556

;A114
;ax= 1330.052  
;ay=  1193.022 
;bx= 837.663   
;by= 900.002
;cx= 1062.886    
;cy= 871.773

;A114 take2
;ax= 1911.332 
;ay= 650.154
;bx= 837.663   
;by= 900.002
;cx= 1237.153 
;cy= 185.254

;A118
;ax=678.456   
;ay=1144.162
;bx=732.148   
;by=1168.065
;cx=565.818    
;cy=704.057

;A118
;ax=  1534.342 
;ay=284.249
;bx= 1783.600 
;by=193.433
;cx= 1934.057 
;cy=281.774 

;A4010.  974.293   1578.331 
;bx= 974.293   
;by=1578.331
;cx= 1339.575   
;cy=1929.892
;ax=1277.410   
;ay=1356.054 

;ax=1448.56 
;ay=2627.64 
;bx=1529.78 
;by=1903.52 
;cx=1682.38 
;cy=1313.05 

;A2734
;ax=1342.73 
;ay=2106.94 
;bx=1384.04 
;by=1926.44 
;cx= 909.12 
;cy=1168.72 

;A4059
;ax= 1143.10 
;ay= 1556.56
;bx= 798.94 
;by=1513.06
;cx=1119.32 
;cy=1065.82
;bx=1549.27 
;by=1643.93 

;A3880
ax= 688.35 
ay=1049.59
bx= 820.61 
by=1933.15
;cx=1383.15 
;cy=1995.99
ax=930.90 
ay=1221.57
;ax=1317.84 
;ay=2129.69
bx=1340.15 
by=1528.23
cx=1121.96 
cy=1315.10

A = sqrt((abs(ax-bx))^2 + (abs(ay-by))^2)
B = sqrt((abs(ax-cx))^2 + (abs(ay-cy))^2)
C = sqrt((abs(cx-bx))^2 + (abs(cy-by))^2)

openr, filelistlun, "/n/Godiva3/jkrick/A3880/listobj", /get_lun   ;listobj
openr, fwhmlun, "/n/Godiva3/jkrick/A3880/fwhmcutoff2",/get_lun        ;fwhmcutoff
filename = "hello"

WHILE (NOT(EOF(filelistlun)) )DO BEGIN
    xshift = 0.
    yshift = 0.
    readf, filelistlun, filename
    
    infile = strmid(filename,0,40) + ".fits"
    outfile =  strmid(filename,0,40)  + ".big.match"
    print, infile
    print, A, B, C
    commandline = '/n/Godiva7/jkrick/Sex/sex ' + infile + " -c /n/Godiva7/jkrick/A114/default.sex"
    spawn, commandline
    
    sexname = "/n/Godiva7/jkrick/A114/SExtractor.cat"
    openr, sexlun, sexname, /get_lun
    
    ;read in fwhmcuttof values for the bottom of the star distribution
    readf, fwhmlun, junk, fwhmlow
    fwhmhigh = fwhmlow + 2.0
    print, "fwhm from file", fwhmlow, fwhmhigh

    numoobjects = 45000
    star=replicate({object, obj:0D, xcentemp:0D,ycentemp:0D,etemp:0D,fluxtemp:0D,magtemp:0D,fwhmtemp:0D,isoareatemp:0D},numoobjects)

    i = 0
    WHILE (NOT EOF(sexlun)) DO BEGIN  
        READF, sexlun, o, xcenter, ycenter,aimage,bimage,e,fluxbest,magbest,magerr,fluxiso,isoarea,fwhm,theta,back
        IF (e LT 0.3) AND (fluxbest GT 0) AND (isoarea GT 6) AND (magbest LT 22.0) AND $
             (fwhm GT fwhmlow) AND (fwhm LT fwhmhigh) THEN BEGIN
            star[i]={object, o,xcenter, ycenter, e, fluxbest, magbest, fwhm, isoarea}
;            xarr[i] = xcenter 
;            yarr[i] = ycenter 
            i = i + 1
        ENDIF
        
    ENDWHILE
    close, sexlun
    free_lun, sexlun

;    xarr = xarr[0:i-1]
;    yarr = yarr[0:i-1]
    star = star[0:i-1]
    print, "there are ", i, "stars"

;GOTO, printstatement

    FOR j = 0, i-1, 1 DO BEGIN
        FOR k = j+1, i-1, 1 DO BEGIN
            FOR l = k+1, i-1, 1 DO BEGIN
                A2 = sqrt((abs(star[j].xcentemp-star[k].xcentemp))^2 + (abs(star[j].ycentemp-star[k].ycentemp))^2)
                B2 = sqrt((abs(star[j].xcentemp-star[l].xcentemp))^2 + (abs(star[j].ycentemp-star[l].ycentemp))^2)
                C2 = sqrt((abs(star[k].xcentemp-star[l].xcentemp))^2 + (abs(star[k].ycentemp-star[l].ycentemp))^2)

                IF star[j].obj EQ 1455 AND star[k].obj EQ 1620 AND star[l].obj EQ 2227 THEN BEGIN 
                    print, "f1455",A2,B2,C2
                endif

                IF (equalto(A2,A,2.5) EQ 1) AND (equalto(B2,B,2.5) EQ 1) AND (equalto(C2,C,2.5) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C
                    xshift = star[j].xcentemp - ax  
                    yshift = star[j].ycentemp - ay 
                    print, "without rot, shift is", xshift, yshift   
                    GOTO, printstatement
                ENDIF
                IF (equalto(A2,A,2.5) EQ 1) AND (equalto(B2,C,2.5) EQ 1) AND (equalto(C2,B,2.5) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C
                    xshift = star[j].xcentemp - bx 
                    yshift = star[j].ycentemp - by
                    print, "without rot, shift is", star[j].xcentemp - bx, star[j].ycentemp - by
                    GOTO, printstatement
                ENDIF
                IF (equalto(A2,B,2.5) EQ 1) AND (equalto(B2,A,2.5) EQ 1) AND (equalto(C2,C,2.5) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C
                    xshift = star[j].xcentemp - ax 
                    yshift = star[j].ycentemp - ay 
                    print, "without rot, shift is", star[j].xcentemp - ax, star[j].ycentemp - ay
                    GOTO, printstatement
                ENDIF
                IF (equalto(A2,B,2.5) EQ 1) AND (equalto(B2,C,2.5) EQ 1) AND (equalto(C2,A,2.5) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C
                    xshift = star[j].xcentemp - cx 
                    yshift = star[j].ycentemp - cy 
                    print, "without rot, shift is", star[j].xcentemp - cx, star[j].ycentemp - cy
                    GOTO, printstatement
                ENDIF
                IF (equalto(A2,C,2.5) EQ 1) AND (equalto(B2,A,2.5) EQ 1) AND (equalto(C2,B,2.5) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C
                    xshift = star[k].xcentemp - cx 
                    yshift = star[k].ycentemp - cy 
                    print, "without rot, shift is", star[k].xcentemp - cx, star[k].ycentemp - cy
                    GOTO, printstatement
                ENDIF
                IF (equalto(A2,C,2.5) EQ 1) AND (equalto(B2,B,2.5) EQ 1) AND (equalto(C2,A,2.5) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C
                    xshift = star[j].xcentemp - cx 
                    yshift = star[j].ycentemp - cy 
                    print, "without rot, shift is", star[j].xcentemp - cx, star[j].ycentemp - cy
                    GOTO, printstatement
                ENDIF
                
                
            ENDFOR
        ENDFOR
    ENDFOR
    
    ;now try making the output file with the template and new coords
    printstatement: print, "moving on now to the files"

    

;xshift = -894.6
;yshift = 197.3

    IF (xshift NE 0.) AND (yshift NE 0.) THEN BEGIN
    openw, lun4, outfile, /get_lun
    openr, lunstars, "/n/Godiva3/jkrick/A3880/4059.stars", /get_lun

        WHILE (NOT EOF(lunstars)) DO BEGIN 
            READF, lunstars, junk, xcenter66, ycenter66,junk,junk,e66,junk,mag66,junk,junk,junk,fwhm66,junk,junk
            
;            openr, sexlun, sexname, /get_lun
;            WHILE (NOT eof(sexlun)) DO BEGIN
;                READF, sexlun, o2, xcenter2, ycenter2,a2,b2,e2,f2,mag2,magerr2,fluxiso2,isoarea2,fwhm2,theta2,back2 ;,junk
            FOR m=0,i-1,1 DO BEGIN
                
                IF (equalto(xcenter66 +xshift, star[m].xcentemp, 5.0) EQ 1) THEN BEGIN

                    IF (equalto(ycenter66 +yshift, star[m].ycentemp, 5.0) EQ 1)  THEN BEGIN

                       printf, lun4, xcenter66+1000., ycenter66+787., star[m].xcentemp + 1000., star[m].ycentemp + 787.
                       print, "matched star", xcenter66, ycenter66, star[m].xcentemp, star[m].ycentemp
                    ENDIF
                ENDIF
                
            ENDFOR
               
;            ENDWHILE
;            close, sexlun
;            free_lun,sexlun
        ENDWHILE
    close, lunstars
;    close, sexlun
    close, lun4
    free_lun, lunstars
;    free_lun,sexlun
    free_lun,lun4
        
    endif
    
    undefine, star
ENDWHILE

close, fwhmlun
close, filelistlun
free_lun, fwhmlun
free_lun, filelistlun

close, /all
END

            
