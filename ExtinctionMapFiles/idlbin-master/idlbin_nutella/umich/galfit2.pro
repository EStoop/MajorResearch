PRO galfit2, imagename, noiseimagename, psfname, segmentationname, directoryname

;example way of running this program-------------------------------
;galfit, "/n/Godiva7/jkrick/galfit/thumb.1.fits",
;"/n/Godiva7/jkrick/galfit/thumbnoise.1.fits",
;"/n/Godiva7/jkrick/galfit/psf.small.fits",
;"/n/Godiva7/jkrick/galfit/segmentation.fits", "/n/Godiva7/jkrick/galfit/"
;-----------------------------------------------------------------------------

close, /all
cd, directoryname, current = idldir

outputname = "model.fits"
fits_read, imagename, data, header

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;run SExtractor on the image to find the objects
;need to have name/location of segmentation image as well as input
;parameters set
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

commandline = '/n/Godiva7/jkrick/Sex/sex ' + imagename + " -c  default.sex"
spawn, commandline


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read in the galaxy information for all galaxies
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

fits_read, segmentationname, segdata, segheader
;figure out the size of the image from the header information
junk = strsplit( segheader[3], /extract)
nxaxis = fix(junk[2])
junk = strsplit( segheader[4], /extract)
nyaxis = fix(junk[2])

print, "nxaxis,nyaxis", nxaxis,nyaxis

;figure out how many SExtractor objects there are
;must be a better way of doing this
OPENR, lun1, "SExtractor.cat", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"
numoobjects = 0
junk = " "
WHILE (NOT eof(lun1)) DO BEGIN
    READF, lun1, junk
    numoobjects= numoobjects + 1
ENDWHILE
close, lun1

;make an array of galaxy objects
galaxy = replicate({object, num:0, xcenter:0D, ycenter:0D, a:0D, b:0D, ellip:0D, flux:0D, $
     mag:0D, isoarea:0D, fwh:0D, theta:0D, bkgd:0D, good:0},numoobjects)

;read in galaxy parameters into galaxy objects
OPENR, lun1, "SExtractor.cat", /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, number, xcen, ycen, alength, blength, ellipticity, fluxbest, magnitude, $
         isoareagal, fwhmin, thetaimage, background
   
    ;if I think it is a galaxy, add it to the list of galaxies
    ;also weed out small junky things, and things on edges
    IF(fwhmin GT 4.0) AND (xcen GT 557) AND (xcen LT nxaxis - 720) AND (ycen GT 635) $
               AND (ycen LT nyaxis - 720) AND (magnitude LT 22.0) AND (isoareagal GT 50) $  ;mag 22.6
               THEN BEGIN 
;    IF(fwhmin GT 4.0) AND (xcen GT 650) AND (xcen LT nxaxis - 880) AND (ycen GT 700) $
;              AND (ycen LT nyaxis - 720) AND (magnitude LT 22.0) AND (isoareagal GT 50) $  ;mag 22.6
;               THEN BEGIN 
 

       galaxy[i] = {object, number, xcen, ycen, alength, blength, ellipticity, fluxbest, $
             magnitude, isoareagal, fwhmin, thetaimage - 90, background, 0}
        i = i +1
        print, xcen, ycen, number
    ENDIF

ENDWHILE

;want to shorten the array of galaxies to only be as long as it needs
;to be
numogalaxies = i
galaxy = galaxy[0:numogalaxies -1]

;-------------------------------------------------------------------------------
;want to iterate through the fitting procedure to lessen the effects
;of overlaping galaxies
;-------------------------------------------------------------------------------
numoiters= 5;number of iterations
galaxy[3].good = 1
FOR m = 1, numoiters,1 DO BEGIN

    IF (m EQ 1) THEN BEGIN

         ;--------------------------------------------------------------------------------
         ;fit each galaxy, one at a time, using the SExtractor segmentation
         ;image as a mask image
         ;--------------------------------------------------------------------------------
        FOR j = 0, numogalaxies-1, 1 DO BEGIN
        ;j= 7
            
            print, "working on galaxy", galaxy[j].num, "which is number ", j," out of ", numogalaxies
            
            ;manually make the segmentation image
            fits_read, segmentationname, segdata, segheader
            FOR x = 1, nxaxis - 1, 1 DO BEGIN   
                FOR y = 1,nyaxis - 1, 1 DO BEGIN
                    IF (segdata(x,y) EQ galaxy[j].num) THEN segdata(x,y) = 0
                    IF (segdata(x,y) GT 0) THEN segdata(x,y) = 1
                ENDFOR
            ENDFOR
            maskname = "mask.fits"
            fits_write, maskname, segdata, segheader
            
            ;call funtion to print out the input file for galfit
            shifts = print_feedme2(imagename, outputname, noiseimagename, psfname, maskname,$
                                  galaxy[j].isoarea, galaxy[j].theta, galaxy[j].a, galaxy[j].b, galaxy[j].xcenter, $
                                  galaxy[j].ycenter, galaxy[j].mag, galaxy[j].bkgd, 0, 0)
            
            ;run galfit
            spawn, '/n/Godiva7/jkrick/galfit/galfit/galfit  galfit.pro.feedme ', exit_status = y
            IF (y GT 0) THEN BREAK ;galfit has crashed-----------------!
            
            ;***************************************************************************************************
            ;now that I have a model, create a model image of a larger size to be 
            ;subtracted from the original image
            ;***************************************************************************************************
            
            ;call funtion to print out the input file for galfit
            ;this time using the parameters from the best fit model made above
            
            shifts = print_feedme2("none", "bigmodel.fits", noiseimagename, $
                                  psfname,maskname,galaxy[j].isoarea, galaxy[j].theta, galaxy[j].a, galaxy[j].b, $
                                  galaxy[j].xcenter, galaxy[j].ycenter, galaxy[j].mag, galaxy[j].bkgd, nxaxis, nyaxis)
            
             ;run galfit
            spawn, '/n/Godiva7/jkrick/galfit/galfit/galfit  galfit.pro.feedme ', exit_status = y
            IF (y GT 0) THEN BREAK ;galfit has crashed-----------------!
            
            ;shift image, add it to bigger image, output the best model image
            fits_read, "bigmodel.fits", modeldata, modelheader
            data[*,*] = 0
            data[shifts[0],shifts[1]] = modeldata
            fits_write, strcompress("best" + string(j) + ".fits", $
                                    /remove_all), data, modelheader
            
             ;if the fit was crappy flip the switch on the "good" parameter
            IF (shifts[2] GT 100) OR (shifts[3] LT 0.09)THEN BEGIN
                openw, dlun2, "diagnostics.log", /GET_LUN, /append
                printf, dlun2, "NOT subtracting galaxy ",  galaxy[j].xcenter, "  ", galaxy[j].ycenter, " because re or b/a are off"
                close, dlun2
                galaxy[j].good = 1
            ENDIF


            ;get rid of the galfit output file, otherwise it will increment the number       
            spawn, "rm galfit.01"
            spawn, "rm galfit.02"
            
        ENDFOR
        
        close, /all
        
    ENDIF ELSE BEGIN
  
       print, "doing iteration # ", m
        FOR j = 0, numogalaxies-1, 1 DO BEGIN
        ;j= 7
            print, "working on galaxy", galaxy[j].num, "which is number ", j," out of ", numogalaxies
      
            ;make a new input image with all
            ;galaxies but the one being worked on subtracted away
            fits_read, imagename, imagedata, imageheader
            newdata = imagedata
            FOR q = 0,numogalaxies - 1, 1 DO BEGIN
                IF (q NE j AND galaxy[q].good EQ 0) THEN BEGIN    ;if not the galaxy being worked on...
                     galname = strcompress("best" + string(q) + ".fits", /remove_all)
                    fits_read, galname, galdata, galheader
                    newdata = newdata - galdata
                ENDIF
            ENDFOR
            fits_write, "subgals.fits", newdata, imageheader
            
            ;manually make the segmentation image
            ;would be nice if i remasked to something smaller, say by rerunning SExtractor?
            fits_read, segmentationname, segdata, segheader
            FOR x = 1, nxaxis - 1, 1 DO BEGIN   
                FOR y = 1,nyaxis - 1, 1 DO BEGIN
                    IF (segdata(x,y) EQ galaxy[j].num) THEN segdata(x,y) = 0
                    IF (segdata(x,y) GT 0) THEN segdata(x,y) = 1
                ENDFOR
            ENDFOR
            maskname = "mask.fits"
            fits_write, maskname, segdata, segheader

            shifts = print_feedme2("subgals.fits", outputname, $
                                   noiseimagename, psfname, maskname, galaxy[j].isoarea, $
                                   galaxy[j].theta, galaxy[j].a, galaxy[j].b, galaxy[j].xcenter, $
                                  galaxy[j].ycenter, galaxy[j].mag, galaxy[j].bkgd, 0, 0)
            
            ;run galfit
            spawn, '/n/Godiva7/jkrick/galfit/galfit/galfit  galfit.pro.feedme ', exit_status = y
            IF (y GT 0) THEN BREAK ;galfit has crashed-----------------!
            
            ;***************************************************************************************************
            ;now that I have a model, create a model image of a larger size to be 
            ;subtracted from the original image
            ;***************************************************************************************************
            
            ;call funtion to print out the input file for galfit
            ;this time using the parameters from the best fit model made above
            
            shifts = print_feedme2("none", "bigmodel.fits", noiseimagename, $
                                  psfname,maskname,galaxy[j].isoarea, galaxy[j].theta, galaxy[j].a, galaxy[j].b, $
                                  galaxy[j].xcenter, galaxy[j].ycenter, galaxy[j].mag, galaxy[j].bkgd, nxaxis, nyaxis)
            
             ;run galfit
            spawn, '/n/Godiva7/jkrick/galfit/galfit/galfit  galfit.pro.feedme ', exit_status = y
            IF (y GT 0) THEN BREAK ;galfit has crashed-----------------!
            
            ;shift image, add it to bigger image, output the best model image
            fits_read, "bigmodel.fits", modeldata, modelheader
            data[*,*] = 0
            data[shifts[0],shifts[1]] = modeldata
            fits_write, strcompress("best" + string(j) + ".fits", $
                                    /remove_all), data, modelheader
            
            ;reset the good parameter in case the fit is really good this time
            galaxy[j].good = 0

            ;if the fit was crappy flip the switch on the "good" parameter
            IF (shifts[2] GT 100) OR (shifts[3] LT 0.09)THEN BEGIN
                openw, dlun2, "diagnostics.log", /GET_LUN, /append
                printf, dlun2, "NOT subtracting galaxy ",  galaxy[j].xcenter, "  ", galaxy[j].ycenter, " because re or b/a are off"
                close, dlun2
                galaxy[j].good = 1
            ENDIF
            
             ;get rid of the galfit output file, otherwise it will increment the number       
            spawn, "rm galfit.01"
            spawn, "rm galfit.02"
            
        ENDFOR

        
    ENDELSE
    
ENDFOR

;--------------------------------------------------
;output an image which is the model galaxies
;and an image which is the original - galaxies
;--------------------------------------------------
allgals = data
allgals[*,*] = 0
FOR k = 0, numogalaxies - 1, 1 DO BEGIN
    fits_read, strcompress("best" + string(k) + ".fits", /remove_all), datagal, headergal

   allgals = allgals + datagal
ENDFOR

fits_write, "allgals.fits", allgals, headergal
fits_read, imagename, data, header
fits_write, "nogals.fits", data - allgals, header

close, /all
undefine, galaxy

END





