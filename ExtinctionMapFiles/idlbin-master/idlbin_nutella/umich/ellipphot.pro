PRO ellipphot
close,/all

fits_read, "/n/Godiva1/jkrick/A3888/larger.fits", rdata, rheader 
fits_read, "/n/Godiva1/jkrick/A3888/largeV.fits", Vdata, Vheader 

rdata = rdata - 2
Vdata = Vdata - 2

;feature 1
;r
rsuma = total(rdata[1372:1378,1555:1574]) 
rnpixa = (1378-1372) * (1574-1555)

rmeanb = mean(rdata[1379:1384,1553:1572]) 
rmeanb2= mean(rdata[1367:1371,1568:1577]) 
rmeanback = (rmeanb+rmeanb2) / 2
rback = rmeanback * rnpixa

rsuma = rsuma - rback

;V
Vsuma = total(Vdata[1372:1378,1555:1574]) 
Vnpixa = (1378-1372) * (1574-1555)

Vmeanb = mean(Vdata[1379:1384,1553:1572]) 
Vmeanb2= mean(Vdata[1367:1371,1568:1577]) 
Vmeanback = (Vmeanb+Vmeanb2) / 2
Vback = Vmeanback * Vnpixa

Vsuma = Vsuma - Vback

print, "feature 1 rsum", rsuma, "which is ", 24.6 + (-2.5)*alog10(rsuma), " rmag"
print, "feature 1 Vsum", Vsuma, "which is ", 24.3 + (-2.5)*alog10(Vsuma), " Vmag"

;feature 2
fits_read, "/n/Godiva1/jkrick/A3888/back1V.fits", Vdata2, Vheader2 
npix = 0
sum = 0.
FOR x = 1000, 2000, 1 DO BEGIN
    FOR y = 1000, 2000, 1 DO BEGIN
        IF (Vdata2[x,y] NE 0.) THEN BEGIN
;             print, x, y, Vdata2[x,y]
            npix = npix + 1
            sum = sum + Vdata2[x,y] - 2
        ENDIF
    ENDFOR
ENDFOR

print, npix, sum

;feature 3
;r
rsuma = total(rdata[1282:1439,711:739]) 
rnpixa = (1439-1282) * (739-711)

rmeanb = mean(rdata[1274:1371,688:704]) 
rmeanb2= mean(rdata[1374:1416,743:753]) 
rmeanback = (rmeanb+rmeanb2) / 2
rback = rmeanback * rnpixa

rsuma = rsuma - rback

;V
Vsuma = total(Vdata[1282:1439,711:739]) 
Vnpixa = (1439-1282) * (739-711)

Vmeanb = mean(Vdata[1274:1371,688:704]) 
Vmeanb2= mean(Vdata[1374:1416,743:753]) 
Vmeanback = (Vmeanb+Vmeanb2) / 2
Vback = Vmeanback * Vnpixa
print, Vsuma, Vback
Vsuma = Vsuma - Vback

print, "feature 3 rsum", rsuma, "which is ", 24.6 + (-2.5)*alog10(rsuma), " rmag"
print, "feature 3 Vsum", Vsuma, "which is ", 24.3 + (-2.5)*alog10(Vsuma), " Vmag"

END
