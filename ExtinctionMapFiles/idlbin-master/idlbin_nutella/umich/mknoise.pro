;**********************************************************************
;summer 2003
;
;this program makes an image of the poisson noise of a combined
;background noise dominated image
;
;input is the the images that were combined and their background
;values
;
;I find this program necessary since I subtract the background from
;all of my sience images, mknoise, or other noise generators do not
;know what the correct backgrounds are and therefor do not find the
;correct noise.
;
;*******************************************************************

PRO mknoise

close, /all

junk = " "
numexposures = 16
gain = 3.

;open blank image
fits_read, "/n/Godiva1/jkrick/A3888/original/blank.shift.fits", blankdata, blankheader

OPENR, lun, '/home/grad/jkrick/idlbin/imagenamesV', /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "file did not open"

names = strarr(16)
i = 0
WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, junk
    g = strsplit(junk, /extract)
    filename = g(0)
    maskname = g(1)
    dataname = g(2)
    headername = g(3)
    bkg = g(4)
    airmass = g(5)
    
    names(i) = dataname

    fits_read, filename, dataname, headername
    fits_read, maskname, maskdata, maskheader
    print, "working on ", filename
  
    FOR x = 0, 4094, 1 DO BEGIN
        FOR y = 0, 4094, 1 DO BEGIN
            IF (dataname(x,y) NE 0.) THEN dataname(x,y) = bkg
            IF (maskdata(x,y) EQ 0.0) THEN maskdata(x,y) = 1.0
            IF (maskdata(x,y) EQ 1.0E6) THEN maskdata(x,y) = 0.0
        ENDFOR
    ENDFOR

    ;correct for airmass of the images
    dataname = dataname / airmass

    ;apply the mask, so as to not add rejected data points
    dataname = dataname * maskdata

    blankdata = blankdata + dataname
    print, "blankdata(2000,2000)", blankdata(2000,2000)
    i = i + 1
ENDWHILE

;write out the background image
bkgdata = blankdata
fits_write, "/n/Godiva1/jkrick/A3888/bkg.fits", bkgdata, blankheader

print, "background", bkgdata(2000,2000)

;want to take the poisson deviates of the bkg in electrons, so need to
;include the gain, and then take it back out later
bkggain = bkgdata * gain
print, "background* gain", bkggain(2000,2000)

;add the readnoise, although this is negligible other than on the very
;edges where I am not fitting galaxies anyway
bkggain = bkggain + 49.0

;find the poisson deviates with mean = background data
imnoise = POIDEV(bkggain)
print, "poisson(background*gain)", imnoise(2000,2000)

;take the gain back out
imnoise = imnoise / gain
print, "poisson(background*gain) / gain", imnoise(2000,2000)

;need to subtract back out the background values
imnoise = imnoise - bkgdata
print, "poisson(background*gain) / gain - bkg", imnoise(2000,2000)

;find out the exposure times involved
fits_read, "/n/Godiva1/jkrick/A3888/3888Vpl.fits", pldata, plheader
exptime = (numexposures - pldata)*900

;correct the noise estimate for the exptimes
;imnoise = imnoise / exptime
imnoise = imnoise / (numexposures - pldata)
print, "poisson(background*gain) / gain - bkg / exptime", imnoise(2000,2000)


fits_write, "/n/Godiva1/jkrick/A3888/noisea.fits", imnoise, blankheader

close ,/all
END

;need to do r and V seperately, seperate name lists or some other way
;of separating them
