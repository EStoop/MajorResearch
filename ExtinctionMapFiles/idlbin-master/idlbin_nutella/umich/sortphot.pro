PRO sortphot

;winter 2003/2004
;this code is written to compare my galaxy photometry to LARCS galaxy
;photometry

close, /all
numoobjects = 4500
larcsgalaxy = replicate({object, lra:0D, ldec:0D, lRmag:0D, lRmagerr:0D,lgunnrmag:0D},numoobjects)

;--------------------------------------------------------------------------------
;read in the larcs magnitudes, and store them in an array
openr, lun1, "/n/Godiva1/jkrick/A3888/pimbblet/larcs.mag.txt", /GET_LUN
i = 0
WHILE (NOT eof(lun1)) DO BEGIN
    READF, lun1,  larcsra, larcsdec,larcsRmag,larcsRmagerr, larvsVmag, larcsr
;    larcsr = larcsRmag + 0.366
;    readf, lun1, larcsra, larcsdec, larcsRmag, larcsVmag, larcsRmagerr, fwhm
    larcsgalaxy[i] = {object, larcsra, larcsdec, larcsRmag, larcsRmagerr,larcsr}
 
    i = i +1 
ENDWHILE
close, lun1
free_lun, lun1

larcsgalaxy = larcsgalaxy[0:i-1]

;--------------------------------------------------------------------------------
;read in my magnitudes and see if they match the Larcs magnitudes

diffarr = fltarr(numoobjects)
mag = fltarr(numoobjects)
m = 0
fits_read,"/n/Godiva1/jkrick/A3888/V.wcs.fits",data,header
openr, huhlun2, "/n/Godiva1/jkrick/A3888/pimbblet/SExtractor.r.cat", /GET_LUN
WHILE (NOT eof(huhlun2)) DO BEGIN
;    READF, huhlun2, myra, mydec, myVflux, myVmag, myVmagerr, fwhm
    READF, huhlun2, num, myx, myy, myrflux, myrmag, myrmagerr, fwhm
     ;get the astrometry from the better astrometry image (same exact area
     ;on sky)
    xyad,header,myx,myy,myra,mydec 

    FOR j = 0, i - 1, 1 DO BEGIN
        IF (myra LT larcsgalaxy[j].lra + 0.0009) AND  (myra GT larcsgalaxy[j].lra - 0.0009) AND $
             (mydec LT larcsgalaxy[j].ldec + 0.0009) AND (mydec GT larcsgalaxy[j].ldec - 0.0009) $
             AND (fwhm LT 30) AND (fwhm GT 2) THEN BEGIN
;            mymag = 24.29 -(2.5*ALOG10(myVflux)) - (0.2*(1.24))
            diffarr(m) = myrmag - larcsgalaxy[j].lgunnrmag 
            mag(m) = myrmag
;            print, mymag, myVmag
;            print, m, larcsgalaxy[j].lra, larcsgalaxy[j].ldec, myra, mydec
;            IF diffarr(m) GT 6.0 THEN print, myra, mydec, myVmag, diffarr(m)
            m = m + 1
        ENDIF
    ENDFOR
ENDWHILE

close, huhlun2
free_lun, huhlun2
diffarr = diffarr[0:m-1]
mag = mag[0:m-1]

;---------------------------------------------------------------------------

;chop them into regions of mag < 20 and mag > 20
;sort the magnitudes into numerical order
sortindex = sort(mag)
sortedmag = mag[sortindex]
FOR b = 0, m - 1 DO BEGIN
    IF sortedmag(b) GT 19.0 THEN BEGIN
        print, "val at which mag goes above 19: ",b
        BREAK
    ENDIF
ENDFOR
;sort the difference the same way
sorteddiff = diffarr[sortindex]

;make a cut into 2 different diffarrays
sorteddiff1 = sorteddiff[0:b-1]
sorteddiff2 = sorteddiff[b-1:m-1]

;find mean, want to throw in a 3 sigma clip, so not fitting the mag limit part.

;--------------------------------
resistant_mean,sorteddiff1, 3, mean1, sigmean1
print, mean1, sigmean1, " mean1"
resistant_mean,sorteddiff2, 3, mean2, sigmean2
print, mean2, sigmean2, " mean2"

;-------------------------------

;for the mags less than 20;..............................
;sort diffarr

sortindex= sort(sorteddiff1)
sorteddiff = sorteddiff1[sortindex]

newarr = fltarr(b)
nrejlow = 5
nrejhigh = 5
reject = 0
print, "before everything there are ",b," vals in the array"
meanval = mean(sorteddiff(0 + nrejlow : b-nrejhigh - 2))
sigmaval = sqrt(abs(meanval))

WHILE (reject EQ 0) DO BEGIN
    print, "mean, sigma", meanval, sigmaval
    n = 0
    reject = 1
    FOR k = 0, b-(1+nrejlow+nrejhigh), 1 DO BEGIN
        IF(sorteddiff(k) LT meanval + 1* sqrt(abs(meanval))) $
               AND (sorteddiff(k) GT meanval - 1*sqrt(abs(meanval))) THEN BEGIN
            newarr(n) = sorteddiff(k)
            n = n+1
        ENDIF ELSE BEGIN
            IF sorteddiff(k) LT meanval THEN nrejlow = nrejlow + 1
            iF sorteddiff(k) GT meanval THEN nrejhigh= nrejhigh+ 1
            reject = 0
        endeLSE
    ENDFOR
    sorteddiff = newarr(0:n-1)
    meanval = mean(sorteddiff)
    sigamval = sqrt(abs(meanval))
    print, "nrejlow, nrejhigh ", nrejlow, nrejhigh
ENDWHILE

print, "after while loop", meanval, sigmaval
meanvallow = meanval
sigmavallow = sigmaval

;for the mags greater than 20;..............................
;sort diffarr
sortindex= sort(sorteddiff2)
sorteddiff = sorteddiff2[sortindex]

newarr = fltarr(m -b)
nrejlow = 4
nrejhigh = 4
reject = 0
print, "before everything there are ",m-b," vals in the array"
meanval = mean(sorteddiff(0 + nrejlow : m-b-nrejhigh - 2))
sigmaval = sqrt(abs(meanval))

WHILE (reject EQ 0) DO BEGIN
    print, "mean, sigma", meanval, sigmaval
    n = 0
    reject = 1
    FOR k = 0, m-b-(1+nrejlow+nrejhigh), 1 DO BEGIN
        IF(sorteddiff(k) LT meanval + 3* sqrt(abs(meanval))) $
               AND (sorteddiff(k) GT meanval - 3*sqrt(abs(meanval))) THEN BEGIN
            newarr(n) = sorteddiff(k)
            n = n+1
        ENDIF ELSE BEGIN
            IF sorteddiff(k) LT meanval THEN nrejlow = nrejlow + 1
            iF sorteddiff(k) GT meanval THEN nrejhigh= nrejhigh+ 1
            reject = 0
        endeLSE
    ENDFOR
    sorteddiff = newarr(0:n-1)
    meanval = mean(sorteddiff)
    sigamval = sqrt(abs(meanval))
    print, "nrejlow, nrejhigh ", nrejlow, nrejhigh
ENDWHILE

print, "after while loop", meanval, sigmaval
meanvalhigh = meanval
sigmavalhigh = sigmaval



;------------------------------------------------------------------------------------
;plot it up


device, true=24
device, decomposed=0
colors = GetColor(/load, Start=1)

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva1/jkrick/A3888/pimbblet/mag.r.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plot, mag, diffarr, thick = 3, psym = 2, xthick = 3, ythick = 3,charthick = 3,$
                        xrange = [16,24], yrange = [-2,2], ystyle = 1,$
	  ytitle = "This survey's r mag - LARCS r mag", title = 'Photometry vs. The Literature', $
                        xtitle = "This survey's r magnitude"


;oplot, findgen(20),0.0*mag + mean1, color = colors.green, thick = 3
;oplot, findgen(20), 0.0*mag + mean1 +sigmean1, thick = 3, color = colors.honeydew
;oplot, findgen(20), 0.0*mag + mean1 -sigmean1, thick = 3, color = colors.honeydew

;oplot, findgen(21) + 19,0.0*mag + mean2, color = colors.green, thick = 3
;oplot, findgen(21)+ 19, 0.0*mag + mean2 +sigmean2, thick = 3, color = colors.honeydew
;oplot, findgen(21)+ 19, 0.0*mag + mean2 -sigmean2, thick = 3, color = colors.honeydew

oplot, findgen(21), 0.0*mag + meanvallow, thick = 3, color = colors.blue
oplot, findgen(21), 0.0*mag + meanvallow + sigmavallow, thick = 3, color = colors.skyblue
oplot, findgen(21), 0.0*mag +meanvallow - sigmavallow, thick = 3, color = colors.skyblue

oplot, findgen(20) +20, 0.0*mag + meanvalhigh, thick = 3, color = colors.blue
oplot, findgen(20)+ 20, 0.0*mag + meanvalhigh + sigmavalhigh, thick = 3, color = colors.skyblue
oplot, findgen(20)+ 20, 0.0*mag +meanvalhigh - sigmavalhigh, thick = 3, color = colors.skyblue


;for the magnitude limit
;x = [18.84,21.84,22.84,23.84,24.84,25.84]    ;for the V images
x = [18.66,21.66,22.66,23.66,24.66,25.66]
y = [-3.0,0.0,1.0,2.0,3.0,4.0]
oplot, x, y, thick = 3, color = colors.red

xyouts, 21.5, -1.2, "Mean value", charsize = 0.9,charthick = 2,color = colors.blue
xyouts, 21.5, -1.35, "One sigma", charsize = 0.9,charthick = 2,color = colors.skyblue
xyouts, 21.5, -1.5, "LARCS magnitude limit", charsize = 0.9,charthick = 2,color = colors.red

device, /close
set_plot, mydevice



undefine, larcsgalaxy
END
