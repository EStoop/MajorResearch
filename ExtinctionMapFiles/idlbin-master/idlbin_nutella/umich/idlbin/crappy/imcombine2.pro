pro imcombine2
l = 0
sum = 0
numocounts=0
nreject = 0
numovals = 0

xmax = 100                      ;xsize of the image
ymax = 100                      ;ysize of the image
numoimages = 5                  ;number of images to combine
hsigma = 5.5                    ;for rejection, in sigma
lsigma = 5.5                    ;for rejection, in sigma
R = 7.                          ;read noise of the tek5
g = 3.                          ;gain 
cntarr = fltarr(numoimages)
newarr = fltarr(numoimages)

;read in blank image
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combblank.fits',newdata,newheader

;read in registered exposure time corrected images to be combined
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest1.fits',data0 ,header0
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest2.fits',data1 ,header1
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest3.fits',data2 ,header2
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest4.fits',data3 ,header3
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest5.fits',data4 ,header4

FOR x = 0, xmax, 1 DO BEGIN
    FOR y = 0, ymax, 1 DO BEGIN
        cntarr(0) = data0(x,y)
        cntarr(1) = data1(x,y)
        cntarr(2) = data2(x,y)
        cntarr(3) = data3(x,y)
        cntarr(4) = data4(x,y)
        
        FOR i = 0, numoimages - 1, 1 DO BEGIN
            IF (cntarr(i) NE 0) THEN BEGIN
                sum = sum + cntarr(i)
                numocounts = numocounts + 1
            ENDIF
        ENDFOR
        
        IF (numocounts EQ 0) THEN newdata(x,y) = 0
        IF (numocounts EQ 1) THEN newdata(x,y) = sum
        IF (numocounts EQ 2) THEN newdata(x,y) = sum/2
        IF (numocounts GE 3) THEN BEGIN
            FOR j = 0, numoimages - 1, 1 DO BEGIN
                IF (cntarr(j) NE 0) THEN BEGIN
                    newarr(l) = cntarr(j)
                    l = l+1
                ENDIF
            ENDFOR
            
            newarr2 = fltarr(l)
            FOR h = 0, l- 1, 1 DO newarr2(h) = newarr(h)
                                

            med = MEDIAN (newarr2, /EVEN)

            IF ((R/g)^2 + med/g LT 0) THEN BEGIN
                print, "sigma neg", x, y, med
                sigma = sqrt(-((R/g)^2 + med/g))
            ENDIF ELSE BEGIN
                sigma = sqrt((R/g)^2 + med/g)
            ENDELSE

            FOR k = 0, numoimages-1, 1 DO BEGIN
                IF(cntarr(k) LT (med - lsigma*sigma) OR cntarr(k) $
                   GT (med + hsigma*sigma)) THEN BEGIN
                    nreject = nreject +1
                    ;print,cntarr(k), x,y,med,sigma,cntarr(0), cntarr(1),cntarr(2),cntarr(3)
                    IF (k EQ 0) THEN BEGIN
                        data0(x,y) = 0
                        IF (x-2 GE 0 AND x+2 LE xmax) THEN BEGIN
                            IF(y-2 GE 0 AND y+2 LE ymax) THEN BEGIN
                                data0(x-2:x+2,y-2:y+2) = 0
                            ENDIF 
                        ENDIF
                    ENDIF
                    IF (k EQ 1) THEN BEGIN
                        data1(x,y) = 0
                        IF (x-2 GE 0 AND x+2 LE xmax) THEN BEGIN
                            IF(y-2 GE 0 AND y+2 LE ymax) THEN BEGIN
                                data1(x-2:x+2,y-2:y+2) = 0
                            ENDIF 
                        ENDIF
                    ENDIF
                    IF (k EQ 2) THEN BEGIN
                        data2(x,y) = 0
                        IF (x-2 GE 0 AND x+2 LE xmax) THEN BEGIN
                            IF(y-2 GE 0 AND y+2 LE ymax) THEN BEGIN
                                data2(x-2:x+2,y-2:y+2) = 0
                            ENDIF 
                        ENDIF
                    ENDIF
                    IF (k EQ 3) THEN BEGIN
                        data3(x,y) = 0
                        IF (x-2 GE 0 AND x+2 LE xmax) THEN BEGIN
                            IF(y-2 GE 0 AND y+2 LE ymax) THEN BEGIN
                                data3(x-2:x+2,y-2:y+2) = 0
                            ENDIF 
                        ENDIF
                    ENDIF

                    IF (k EQ 4) THEN BEGIN
                        data4(x,y) = 0
                        IF (x-2 GE 0 AND x+2 LE xmax) THEN BEGIN
                            IF(y-2 GE 0 AND y+2 LE ymax) THEN BEGIN
                                data4(x-2:x+2,y-2:y+2) = 0
                            ENDIF 
                        ENDIF
                    ENDIF

                ENDIF ELSE BEGIN
                    IF ( cntarr(k) NE 0) THEN BEGIN
                        sum = sum + cntarr(k)
                        numovals = numovals +1
                    ENDIF
                ENDELSE
                
            ENDFOR

;if all the values are rejected,
;likely in the center of a star, should change this to not include the
;negative values
            IF (numovals EQ 0) THEN BEGIN
                print,"numovals = 0", x, y
                sum = cntarr(0)+cntarr(1)+cntarr(2)+cntarr(3)
                numovals = 4
            ENDIF

            newdata(x,y) = sum/numovals
            
        ENDIF
        sum = 0
        numocounts = 0
        l = 0
        nreject = 0
        numovals = 0
    ENDFOR
ENDFOR

FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/testcombine.fits', newdata,newheader
;FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest1.fits',data0, header0
;FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest2.fits',data1 ,header1
;FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest3.fits',data2 ,header2
;FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest4.fits',data3 ,header3
;FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest5.fits',data4 ,header4

END
