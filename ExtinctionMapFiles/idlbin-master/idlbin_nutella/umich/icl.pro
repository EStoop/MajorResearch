PRO icl
close, /all

fits_read, "/n/Godiva1/jkrick/A3888/centerV2000.div.s.fits", galdata, galheader
fits_read, "/n/Godiva1/jkrick/A3888/centerV2000.block.fits", blockdata, blockheader
blockdata = blockdata - 2.
galdata = galdata - 2.
xcenter = 1030.
ycenter = 1030.
inner = 885.
outer = 1000.

total = double(0.)
totalicl =double(0.)
ann = double(0)
i = long(0)
j = long(0)
k = long(0)
FOR x = 0, 2000, 1 DO BEGIN
    FOR y = 0, 2000, 1 DO BEGIN
        IF galdata[x,y] GT -10. THEN BEGIN
            total = total + galdata[x,y] ;total flux in the image, minus bad pixels etc
            j = j + 1     ;number of good pixels
            r = sqrt((abs(x - xcenter)^2) + (abs(y - ycenter)^2))
            IF r GT inner AND r LT outer THEN BEGIN
                 ann = ann + galdata[x,y] ;flux in the background annulus
                i = i + 1       ;number of pixels in the background annulus
            ENDIF
        ENDIF

        IF blockdata[x,y] GT -70. THEN BEGIN
            totalicl = totalicl + blockdata[x,y]
            k = k + 1
        ENDIF

    ENDFOR
ENDFOR

;get an estimage of background from large image instead
;!!!need to block stars
bigxcenter = 1908.
bigycenter = 2195.
inner = 1550.
outer = 1800.
fits_read, "/n/Godiva1/jkrick/A3888/A3888V.div.mask.fits", bigdata, bigheader
FOR x = 0, 4094, 1 DO BEGIN
    FOR y = 0, 4094, 1 DO BEGIN
        IF bigdata[x,y] GT -1. THEN BEGIN
            r = sqrt((abs(x - bigxcenter)^2) + (abs(y - bigycenter)^2))
            IF r GT inner AND r LT outer THEN BEGIN
                ann = ann + bigdata[x,y] ;flux in the background annulus
                i = i + 1       ;number of pixels in the background annulus
            ENDIF
        ENDIF
    ENDFOR
ENDFOR


;need to scale the background to the total good area in the image
multiplier = float(j)/float(i)
ann = ann * multiplier
print, "total, ann, j,i,totalicl", total, ann,j,i, totalicl

;background subtracted total light in the galaxy
totalb = total - ann

print, "fraction of light in the ICL w/ 'background subtraction' ", (totalicl / totalb) * 100, "%"
print, "fraction of light in the ICL w/o 'background subtraction' ", (totalicl / total) * 100, "%"

END


