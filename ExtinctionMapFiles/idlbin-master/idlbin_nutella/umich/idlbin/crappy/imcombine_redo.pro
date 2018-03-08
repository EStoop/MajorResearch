;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;May 2002
;
; This program is a form of IRAF's imcombine
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro imcombine_redo

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
namearr = ptrarr(numoimages,/allocate_heap)

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

nonzeroflt_a = fltarr(numoimages)
nonzeroptr_a = ptrarr(numoimages,/allocate_heap)

;at each pixel...
FOR x = 0, xmax, 1 DO BEGIN
    FOR y = 0, ymax, 1 DO BEGIN
        ;throw away the images with no value at that particular pixel
        ; keep the ones with real data in an array of floats
        FOR i =0, numoimages - 1, 1 DO BEGIN
            IF ((*imarr(i))(x,y) NE 0) THEN BEGIN
                nonzeroflt_a(numnotzero) = (*imarr(i))(x,y)
                nonzeroptr_a(numnotzero) = imarr(i)
                numnotzero = numnotzero +1
            ENDIF
        ENDFOR
        
        ;truncate the arrays to be only as long as they need to be
        nonzeroflt = fltarr(numnotzero)
        nonzeroptr = ptrarr(numnotzero, /allocate_heap)
        FOR j=0,numnotzero-1,1 DO BEGIN
            nonzeroflt(j) = nonzeroflt_a(j)
            nonzeroptr(j) = nonzeroptr_a(j)
        ENDFOR 
        
        ;if none of the pixels to be combined have data
        IF (numnotzero EQ 0) THEN newdata(x,y) = 0

        ;if there is only one good data point
        IF (numnotzero EQ 1) THEN newdata(x,y) = nonzeroarr[0]

        ;if there are two good data points
        IF (numnotzero EQ 2) THEN newdata(x,y) = (nonzeroarr[0] + nonzeroarr[1])/2

        ;if there are three or more good data points
        IF (numnotzero GE 3) THEN BEGIN  		
            med = MEDIAN(nonzeroflt, /EVEN)
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
                IF ((nonzeroflt(j) LT (med + (hsigma * sigma))) AND $
                    (nonzeroflt(j) GT (med - (lsigma * sigma)))) THEN BEGIN
                    sum = sum + nonzeroflt(j)
                    ncombine = ncombine +1
                ENDIF ELSE BEGIN
                    nreject = nreject +1
                    data = **nonzeroptr(j)
                    data(x,y) = 0
                    ;(*nonzeroptr(j))(x,y) = 0
                    ;print, "rejecting",x,y
                ENDELSE
            ENDFOR
            IF (ncombine EQ 0) THEN BEGIN
                ;plan for if it rejects all values
                ;print, med,nonzeroflt(0),nonzeroflt(1), nonzeroflt(2), nonzeroflt(3)
                FOR m = 0,numnotzero-1 DO BEGIN
                    sum = sum+nonzeroflt(m)
                    ncombine = ncombine +1
                ENDFOR
                newdata(x,y) = sum/ncombine
                ;sum = nonzeroflt(0)+nonzeroflt(1)+nonzeroflt(2)+nonzeroflt(3)
                ;newdata(x,y) = sum/4
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


FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/testcombine.fits', newdata,newheader
FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest1.fits',data0, header0
FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest2.fits',data1 ,header1
FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest3.fits',data2 ,header2
FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest4.fits',data3 ,header3
FITS_WRITE, '/n/Sheriff1/jkrick/sep00/A3888/imcombine/combtest5.fits',data4 ,header4


print, "nreject %", (nreject/40000.)*100
END



                                ;grow out the bad value 
                                ;need to make sure not off the edge of the array
                                ;this will definitely get the center, but won't grow if on the edge
                    
;IF (x-2 GT 0 AND x+2 LT xmax) THEN BEGIN
;    data = *namenonzer0(j)
;    data(x-2:x+2,y) = 0
;    IF (y-2 GT 0 AND y+2 LT ymax) THEN BEGIN
;        data(x-2:x+2,y+1) = 0
;        data(x-2:x+2,y+2) = 0
;        data(x-2:x+2,y-1) = 0
;        data(x-2:x+2,y-2) = 0
;    ENDIF
;ENDIF


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
