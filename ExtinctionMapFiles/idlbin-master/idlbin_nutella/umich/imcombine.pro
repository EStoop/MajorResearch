;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;May 2002
;
; This program is a crappy form of IRAF's imcombine
; all versions of this program are not good.  they take forever,
; and so have not been sufficiently debugged
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro imcombine

k = 0
sum = float(0)
numnotzero = 0
ncombine = 0
xmax = 100                      ;xsize of the image
ymax = 100                      ;ysize of the image
numoimages = 5                  ;number of images to combine
hsigma = 5.5                    ;for rejection, in sigma
lsigma = 5.5                    ;for rejection, in sigma
R = 7.                          ;read noise of the tek5
g = 3.                          ;gain 
nreject = 0                     ;number of pixels rejected
im = fltarr(numoimages)
;namearr = strarr("data0","data1","data2","data3","data4")

;read in blank image
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combblank.fits',newdata,newheader

;read in registered exposure time corrected images to be combined
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest1.fits',data0 ,header0
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest2.fits',data1 ,header1
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest3.fits',data2 ,header2
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest4.fits',data3 ,header3
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest5.fits',data4 ,header4

;make an array of pointers to the data arrays
imarr = ptrarr(numoimages, /allocate_heap)
*imarr[0] = data0
*imarr[1] = data1
*imarr[2] = data2
*imarr[3] = data3
*imarr[4] = data4

;at each pixel...
FOR x = 0, xmax, 1 DO BEGIN
    FOR y = 0, ymax, 1 DO BEGIN
        ;throw away the images with no value at that particular pixel
        ; keep the ones with real data in an array of floats
        FOR i =0, numoimages - 1, 1 DO BEGIN
            dataarr = *imarr[i]
            IF (dataarr(x,y) NE 0) THEN BEGIN
                im[numnotzero] = dataarr(x,y)
                ;namearr(numnotzero) = imarr(i)
                numnotzero = numnotzero +1
            ENDIF
        ENDFOR
        
        ;truncate the arrays to be only as long as they need to be
        nonzeroarr = ptrarr(numnotzero, /allocate_heap)
        FOR k = 0, numnotzero - 1, 1 DO nonzeroarr(k) = im(k)
        ;namenonzero = ptrarr(numnotzero, /allocate_heap)
        ;FOR k = 0, numnotzero - 1, 1 DO namenonzero(k) = namearr(k)
        
        ;if none of the pixels to be combined have data
        IF (numnotzero EQ 0) THEN newdata(x,y) = 0

        ;if there is only one good data point
        IF (numnotzero EQ 1) THEN newdata(x,y) = nonzeroarr[0]

        ;if there are two good data points
        IF (numnotzero EQ 2) THEN newdata(x,y) = (nonzeroarr[0] + nonzeroarr[1])/2

        ;if there are three or more good data points
        IF (numnotzero GE 3) THEN BEGIN  		
            ;could do this iteratively,find med,reject,find new med until not 
            ;rejecting anymore
            med = MEDIAN(nonzeroarr, /EVEN)
            ;result = MOMENT(nonzeroarr,SDEV = sigmanormal)
            
            ;using sigma based on assumption of just shot noise
            ;need to deal with neg values in the sqrt
            IF ((R/g)^2 + med/g LT 0) THEN BEGIN
                print, "sigma neg", x, y, med
                sigma = sqrt(-((R/g)^2 + med/g))
            ENDIF ELSE BEGIN
                sigma = sqrt((R/g)^2 + med/g)
            ENDELSE

            FOR j = 0, numnotzero-1, 1 DO BEGIN
                IF ((nonzeroarr(j) LT (med + (hsigma * sigma))) AND $
                    (nonzeroarr(j) GT (med - (lsigma * sigma)))) THEN BEGIN
                    IF(x EQ 31 AND y EQ 67) THEN BEGIN
                        print, x, y, j, nonzeroarr(j), sigma,med+(hsigma*sigma), $
                          med - (lsigma*sigma)
                    ENDIF
                    sum = sum + nonzeroarr[j]
                    ncombine = ncombine +1
                ENDIF ELSE BEGIN
                    nreject = nreject +1
                    ;grow out the bad value 
                    ;need to make sure not off the edge of the array
                    ;this will definitely get the center, but won't grow if on the edge
                    
                    IF (x-2 GT 0 AND x+2 LT xmax) THEN BEGIN
                        data = *nonzeroarr(j)
                        data(x-2:x+2,y) = 0
                        IF (y-2 GT 0 AND y+2 LT ymax) THEN BEGIN
                            data(x-2:x+2,y+1) = 0
                            data(x-2:x+2,y+2) = 0
                            data(x-2:x+2,y-1) = 0
                            data(x-2:x+2,y-2) = 0
                        ENDIF
                    ENDIF
                ENDELSE
            ENDFOR
            IF (ncombine EQ 0) THEN BEGIN
                ;plan for if it rejects all values
                ;print, med, nonzeroarr(0),nonzeroarr(1), nonzeroarr(2), nonzeroarr(3)
                sum = nonzeroarr(0)+nonzeroarr(1)+nonzeroarr(2)+nonzeroarr(3)
                newdata(x,y) = sum/4
            ENDIF ELSE BEGIN
                newdata(x,y) = sum/ncombine
            ENDELSE
            ;rezero numbers for the next pixel
            sum = 0.
            ncombine = 0
        ENDIF 
        numnotzero = 0
    ENDFOR
ENDFOR

;newdata0 = *namenonzero(0)
;if (newdata0(0,0) EQ data0(0,0) then newdata0 = data0
;fitswrite,newdata0,header0

FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/testcombine.fits', newdata,newheader
FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest1.fits',data0, header0
FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest2.fits',data1 ,header1
FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest3.fits',data2 ,header2
FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest4.fits',data3 ,header3
FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest5.fits',data4 ,header4


print, "nreject %", nreject/40000.
END






;uneeded
;IF (data0(x,y) NE 0.) THEN BEGIN
;	im[numnotzero] = data0(x,y)
;	numnotzero = numnotzero + 1
;ENDIF
;IF (data1(x,y) NE 0.) THEN BEGIN
;	im[numnotzero] = data1(x,y)
;	numnotzero = numnotzero +1
;ENDIF
;IF (data2(x,y) NE 0.) THEN BEGIN
;	im[numnotzero] = data2(x,y)
;	numnotzero = numnotzero +1
;ENDIF
;IF (data3(x,y) NE 0.) THEN BEGIN
;	im[numnotzero] = data3(x,y)
;	numnotzero = numnotzero +1
;ENDIF
;IF (data4(x,y) NE 0.) THEN BEGIN
;	im[numnotzero] = data4(x,y)
;	numnotzero = numnotzero +1
;ENDIF
