;+
; NAME: Jessica Krick
;
; PURPOSE: this script creates  matched star lists as input to geomap
;                       for rotated images
;
; CALLING SEQUENCE: register
;
; INPUTS: 1) list of files to be registered (should already be set into 4000x4000 images)
;                       should already be rotated
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

PRO register

close, /all

;openr, lun, "/n/Godiva2/jkrick/A141/phot_ccd2075.dat", /get_lun
;this is the file with the template star locations.

;for now I have chosen three stars in the very center by hand

;A2734
;ax=1342.73 
;ay=2106.94 
;ax= 657.35 
;ay=1831.21
;bx=849.29 
;by=2135.95
;bx=1384.04 
;by=1926.44 
;cx= 909.12 
;cy=1168.72 


;A4010
;bx= 974.293   
;by=1578.331
;cx= 1339.575   
;cy=1929.892
;ax=1277.410   
;ay=1356.054 

;ax= 627.60 
;ay=1965.80 
;bx= 974.37 
;by=1578.24 
;cx=1277.44 
;cy=1356.02 


;ax= 627.60 
;ay=1965.80
;bx=1529.78 
;by=1903.52
;cx=1277.44 
;cy=1356.02

;A2556
;ax = 1248.856 
;ay = 1582.925
;bx = 1312.284 
;by=1282.639
;cx=928.858 
;cy=1148.556

;A141
;ax =   2122.03 
;ay = 2093.68
;bx = 2222.89 
;by = 1861.07
;cx =  2045.88 
;cy = 1685.05

;for 8050;A3984
;ax=1131.194 
;ay=1127.268
;bx=1314.008   
;by=1172.899
;cx=1409.085    
;cy=734.332 

;A3984
;ax =  177.331 
;ay = 352.118
;bx = 376.410 
;by = 144.161
;cx =  680.838 
;cy = 235.116 

;second A3984
;ax=916.987 
;ay=1323.896
;bx=786.931 
;by=1371.524
;cx=931.064 
;cy=1198.204


;A3880
;ax= 688.35 
;ay=1049.59
;cx= 820.61 
;cy=1933.15
;cx=1383.15 
;cy=1995.99
bx=930.90 
by=1221.57
;ax=1317.84 
;ay=2129.69
;bx=1340.15 
;by=1528.23
;cx=1121.96 
;cy=1315.10
ax= 1227.87 
ay=1296.12
;ax= 853.29 
;ay=1196.28 
 cx=1069.67 
cy=1088.23
 

A = sqrt((abs(ax-bx))^2 + (abs(ay-by))^2)
B = sqrt((abs(ax-cx))^2 + (abs(ay-cy))^2)
C = sqrt((abs(cx-bx))^2 + (abs(cy-by))^2)

openr, filelistlun, "/n/Godiva3/jkrick/A3880/listobj3", /get_lun
openr, fwhmlun, "/n/Godiva3/jkrick/A3880/fwhmcutoff4",/get_lun

filename = "hello"
WHILE (NOT(EOF(filelistlun)) )DO BEGIN
    xshift = 0.
    yshift = 0.
    readf, filelistlun, filename

    infile = strmid(filename,0,40) + ".rotate.fits"
    infile2 = strmid(filename,0,40)+".fits"
    outfile =  strmid(filename,0,40)  + ".big.match"
    print, infile, outfile

    

    commandline = '/n/Godiva7/jkrick/Sex/sex ' + infile + " -c /n/Godiva4/jkrick/A3984/default.sex"
    spawn, commandline
    commandline = '/n/Godiva7/jkrick/Sex/sex ' + infile2 + " -c /n/Godiva4/jkrick/A3984/default2.sex"
    spawn, commandline
    

    tempname = "/n/Godiva4/jkrick/A3984/SExtractor.cat"
    openr, lun2, tempname, /get_lun
    
    ;read in fwhmcuttof values for the bottom of the star distribution
    readf, fwhmlun, junk, fwhmlow
    fwhmhigh = fwhmlow + 2.0
    print, "fwhm from file", fwhmlow, fwhmhigh

    xarr = fltarr(10000)
    yarr = fltarr(10000)
    i = 0
    WHILE (NOT EOF(lun2)) DO BEGIN  
        READF, lun2, o, xcenter, ycenter,aimage,bimage,e,fluxbest,magbest,magerr,fluxiso,isoarea,fwhm,theta,back
        IF (e LT 0.3) AND (fluxbest GT 0) AND (isoarea GT 6) AND (magbest LT 22) AND $
             (fwhm GT fwhmlow) AND (fwhm LT fwhmhigh) THEN BEGIN
            xarr[i] = xcenter 
            yarr[i] = ycenter 
            i = i + 1
        ENDIF
        
    ENDWHILE
    close, lun2
    free_lun, lun2

    xarr = xarr[0:i-1]
    yarr = yarr[0:i-1]
    print, "there are ", i, " values in the center arrays"
    FOR j = 0, i-1, 1 DO BEGIN
        FOR k = j+1, i-1, 1 DO BEGIN
            FOR l = k+1, i-1, 1 DO BEGIN
                A2 = sqrt((abs(xarr[j]-xarr[k]))^2 + (abs(yarr[j]-yarr[k]))^2)
                B2 = sqrt((abs(xarr[j]-xarr[l]))^2 + (abs(yarr[j]-yarr[l]))^2)
                C2 = sqrt((abs(xarr[k]-xarr[l]))^2 + (abs(yarr[k]-yarr[l]))^2)
;               IF yarr[j] EQ 1462.678 AND yarr[k] EQ 2218.865 AND yarr[l] EQ 2399.026 THEN print, "A2,B2,C2",A2,B2,C2

                 
                IF (equalto(A2,A,3.0) EQ 1) AND (equalto(B2,B,3.0) EQ 1) AND (equalto(C2,C,3.0) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A,B,C,xarr[j],yarr[j]
                    xshift = xarr[j] - ax  
                    yshift = yarr[j] - ay 
                    print, "without rot, shift is", xshift, yshift   
                    GOTO, printstatement
                ENDIF
                IF (equalto(A2,A,3.0) EQ 1) AND (equalto(B2,C,3.0) EQ 1) AND (equalto(C2,B,3.0) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C,xarr[j],yarr[j]
                    xshift = xarr[j] - bx 
                    yshift = yarr[j] - by
                    print, "without rot, shift is", xarr[j] - bx, yarr[j] - by
                    GOTO, printstatement
                ENDIF
                IF (equalto(A2,B,3.0) EQ 1) AND (equalto(B2,A,3.0) EQ 1) AND (equalto(C2,C,3.0) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C,xarr[j],yarr[j]
                    xshift = xarr[j] - ax 
                    yshift = yarr[j] - ay 
                    print, "without rot, shift is", xarr[j] - ax, yarr[j] - ay
                    GOTO, printstatement
                ENDIF
                IF (equalto(A2,B,3.0) EQ 1) AND (equalto(B2,C,3.0) EQ 1) AND (equalto(C2,A,3.0) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C,xarr[j],yarr[j]
                    xshift = xarr[j] - cx 
                    yshift = yarr[j] - cy 
                    print, "without rot, shift is", xarr[j] - cx, yarr[j] - cy
                    GOTO, printstatement
                ENDIF
                IF (equalto(A2,C,3.0) EQ 1) AND (equalto(B2,A,3.0) EQ 1) AND (equalto(C2,B,3.0) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C,xarr[j],yarr[j]
                    xshift = xarr[k] - cx 
                    yshift = yarr[k] - cy 
                    print, "without rot, shift is", xarr[k] - cx, yarr[k] - cy
                    GOTO, printstatement
                ENDIF
                IF (equalto(A2,C,3.0) EQ 1) AND (equalto(B2,B,3.0) EQ 1) AND (equalto(C2,A,3.0) EQ 1) THEN BEGIN
                    print, "found something", A2, B2, C2, A, B,C,xarr[j],yarr[j]
                    xshift = xarr[j] - cx 
                    yshift = yarr[j] - cy 
                    print, "without rot, shift is", xarr[j] - cx, yarr[j] - cy
                    GOTO, printstatement
                ENDIF
                
                
            ENDFOR
        ENDFOR
    ENDFOR
    
    ;now try making the output file with the template and new coords
    printstatement: print, "moving on now to the files"
;xshift = -634.
;yshift=-210
    openw, lun4, outfile, /get_lun
    openr, lun, "/n/Godiva3/jkrick/A3880/4059.stars", /get_lun
IF (xshift NE 0.) AND (yshift NE 0.) THEN begin

    WHILE (NOT EOF(lun)) DO BEGIN   ;reference stars
        READF, lun, o, xcenter, ycenter,junk,junk,e,junk,mag,junk,junk,junk,junk,junk,junk
 
        openr, rotlun, "/n/Godiva4/jkrick/A3984/SExtractor.cat", /get_lun ;is rotated
        WHILE (NOT eof(rotlun)) DO BEGIN

            READF, rotlun, orot, xcenterr, ycenterr,aimager,bimager,er,fluxbestr,magbestr,magerrr,$
                           fluxisor,isoarear,fwhmr,thetar,backr
;            IF  fix(xcenterr) EQ 1136  THEN BEGIN
;                print, "I know there is a match here",xcenterr, ycenterr, xcenter + xshift, ycenter + yshift
;            ENDIF

            IF (equalto(xcenter +xshift, xcenterr, 3.0) EQ 1) AND (er LT 0.5) THEN BEGIN

                IF (equalto(ycenter +yshift, ycenterr, 3.0) EQ 1)  THEN BEGIN

                    ;have found the rotated match, but
                    print, "there is a match in rotated version",xcenter, ycenter, xcenterr, ycenterr

                    ;really want to print the un-rotated version
                    openr, sex3lun, "/n/Godiva4/jkrick/A3984/SExtractor2.cat", /get_lun   ;is not-rotated

                    WHILE (NOT eof(sex3lun)) DO BEGIN
                        READF, sex3lun, o3, xcenter3, ycenter3,aimage3,bimage3,e3,fluxbest3,magbest3,$
                                        magerr3,isoarea3,fwhm3,theta3,back3
                        IF  (equalto(magbestr,magbest3, 1.0) EQ 1) AND (equalto(backr,back3, 1.0) EQ 1) AND $
                              (equalto(magerrr,magerr3,0.1) EQ 1) AND $
                              (equalto(xcenterr, xcenter3, 10) EQ 1 ) AND (equalto(ycenterr, ycenter3,10) EQ 1) THEN BEGIN
                            
                            printf, lun4, xcenter+1000., ycenter+787., xcenter3+ 1000., ycenter3+ 787.
                            print, "matched star", xcenter, ycenter,xcenterr+ 1000., ycenterr+ 787., xcenter3+ 1000., $
                                        ycenter3+ 787.
                        ENDIF

                    ENDWHILE
                    close, sex3lun
                    free_lun, sex3lun
                ENDIF
            ENDIF
            
            
        ENDWHILE
        close, rotlun
        free_lun,rotlun
    ENDWHILE
endif
    close, lun
    close, lun2
    close, lun4
    free_lun, lun
    free_lun,lun2
    free_lun,lun4
    
ENDWHILE

close, fwhmlun
close, filelistlun
free_lun, fwhmlun
free_lun, filelistlun

close, /all
END

            
