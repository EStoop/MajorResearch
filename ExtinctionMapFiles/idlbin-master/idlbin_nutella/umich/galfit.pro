PRO galfit, imagename, noiseimagename, psfname, segmentationname, directoryname

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
free_lun, lun1

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
   neigh = 0
    ;if I think it is a galaxy, add it to the list of galaxies
    ;also weed out small junky things, and things on edges
    IF(fwhmin GT 4.0) AND (xcen GT 200) AND (xcen LT nxaxis - 200) AND (ycen GT 200) $
               AND (ycen LT nyaxis - 200) AND (magnitude LT 22.) AND (isoareagal GT 50) $  ;mag 22.6
               THEN BEGIN 
;    IF(fwhmin GT 4.0) AND (xcen GT 557) AND (xcen LT nxaxis - 720) AND (ycen GT 635) $
;               AND (ycen LT nyaxis - 720) AND (magnitude LT 21.4) AND (isoareagal GT 50) $  ;mag 22.6
;               THEN BEGIN 
;    IF(fwhmin GT 4.0) AND (xcen GT 650) AND (xcen LT nxaxis - 880) AND (ycen GT 700) $
;              AND (ycen LT nyaxis - 720) AND (magnitude LT 22.0) AND (isoareagal GT 50) $  ;mag 22.6
;               THEN BEGIN 

;don't include galaxies in the big bright clumps
;        IF  (xcen LT 406. + 30.)  AND (xcen GT 406. - 30.) AND (ycen LT 598. + 30.) AND (ycen GT 598. - 30.)$
;         THEN neigh = 1
;        IF  (xcen LT 523. + 30.)  AND (xcen GT 523. - 30.) AND (ycen LT 638. + 30.) AND (ycen GT 638. - 30.)$
;          THEN neigh = 1
;        IF (xcen LT 586. + 30.)  AND (xcen GT 586. - 30.) AND (ycen LT 398. + 30.) AND (ycen GT 398. - 30.)$
;         THEN neigh = 1
;        IF (xcen LT 654. + 30.)  AND (xcen GT 654. - 30.) AND (ycen LT 372. + 30.) AND (ycen GT 372. - 30.)$
;         THEN neigh = 1
; don't include the star in the center
         IF  (xcen LT 580. + 1.)  AND (xcen GT 580. - 1.) AND (ycen LT 526. + 1.) AND (ycen GT 526. - 1.)$
         THEN neigh = 1
   
        IF neigh EQ 0 THEN BEGIN

            galaxy[i] = {object, number, xcen, ycen, alength, blength, ellipticity, fluxbest, $
              magnitude, isoareagal, fwhmin, thetaimage - 90, background, 0}
            i = i +1
        ENDIF ELSE print, xcen, ycen, "is not in the array"
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

FOR m = 1, numoiters,1 DO BEGIN

    IF (m EQ 1) THEN BEGIN

         ;--------------------------------------------------------------------------------
         ;fit each galaxy, one at a time, using the SExtractor segmentation
         ;image as a mask image
         ;--------------------------------------------------------------------------------
        FOR j = 0, numogalaxies-1, 1 DO BEGIN
        ;j= 7
            
            print, "working on galaxy", galaxy[j].num, "which is number ", j," out of ", numogalaxies
            
            ;manually make the segmentation image (imreplace)
  ;          fits_read, segmentationname, segdata, segheader
  ;          FOR x = 1, nxaxis - 1, 1 DO BEGIN   
  ;              FOR y = 1,nyaxis - 1, 1 DO BEGIN
  ;                  IF (segdata(x,y) EQ galaxy[j].num) THEN segdata(x,y) = 0
  ;                  IF (segdata(x,y) GT 0) THEN segdata(x,y) = 1
  ;              ENDFOR
  ;          ENDFOR
  ;          maskname = "mask.fits"
  ;          fits_write, maskname, segdata, segheader
   
                    ;make a mask image
            ;add the galaxy numbers to the bumpy image
            fits_read, "/n/Godiva7/jkrick/galfit/examples/mask.fits", bumpy, header
            fits_read, segmentationname, segdata, segheader
            FOR x = 1, nxaxis - 1,1 DO BEGIN
                FOR y =1,nyaxis - 1, 1 DO BEGIN
                     IF (bumpy(x,y) GT 1.01) THEN bumpy(x,y) = segdata(x,y) ELSE bumpy(x,y) = 0.
                ENDFOR
            ENDFOR

            ;get rid of the galaxy being worked on
            FOR x = 1, nxaxis - 1, 1 DO BEGIN   
                FOR y = 1,nyaxis - 1, 1 DO BEGIN
                    IF (bumpy(x,y) EQ galaxy[j].num) THEN bumpy(x,y) = 0
                    IF (bumpy(x,y) GT 0) THEN bumpy(x,y) = 1
                ENDFOR
            ENDFOR

            maskname = "galmask.fits"
            fits_write, maskname, bumpy, segheader

            ;call funtion to print out the input file for galfit
            shifts = print_feedme(imagename, outputname, noiseimagename, psfname, maskname,$
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
            
            shifts = print_feedme("none", "bigmodel.fits", noiseimagename, $
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
            
            ;save a few of the individual images
            mvcomm = strcompress("mv model.fits " + strcompress("model.72."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num EQ 72) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.312."+string(m-1)+".fits",/remove_all))
            IF ( galaxy[j].num EQ 312) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.397."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num EQ 397) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.852."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num EQ 852) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.1."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num EQ 1) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.427."+string(m-1)+".fits",/remove_all))
            IF ( galaxy[j].num EQ 427) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.284."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num  EQ 284) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.477."+string(m-1)+".fits",/remove_all))
            IF ( galaxy[j].num EQ 477) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.197."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num EQ 197) THEN spawn,  mvcomm, exit_status = y



             ;if the fit was crappy flip the switch on the "good" parameter
            IF (shifts[2] GT 99) OR (shifts[3] LT 0.09)THEN BEGIN
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
            mname = strcompress("subgals."+string(m-1) + ".fits")
            fits_write,mname, newdata, imageheader
 
            ;make a mask image
            ;add the galaxy numbers to the bumpy image
            fits_read, "/n/Godiva7/jkrick/galfit/examples/mask.fits", bumpy, header
            fits_read, segmentationname, segdata, segheader
            FOR x = 1, nxaxis - 1,1 DO BEGIN
                FOR y =1,nyaxis - 1, 1 DO BEGIN
                     IF (bumpy(x,y) GT 1.01) THEN bumpy(x,y) = segdata(x,y) ELSE bumpy(x,y) = 0.
                ENDFOR
            ENDFOR

            ;get rid of the galaxy being worked on
            FOR x = 1, nxaxis - 1, 1 DO BEGIN   
                FOR y = 1,nyaxis - 1, 1 DO BEGIN
                    IF (bumpy(x,y) EQ galaxy[j].num) THEN bumpy(x,y) = 0
                    IF (bumpy(x,y) GT 0) THEN bumpy(x,y) = 1
                ENDFOR
            ENDFOR

            maskname = "galmask.fits"
            fits_write, maskname, bumpy, segheader


            shifts = print_feedme("subgals.fits", outputname, $
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
            
            shifts = print_feedme("none", "bigmodel.fits", noiseimagename, $
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

            ;save a few of the individual images
            mvcomm = strcompress("mv model.fits " + strcompress("model.72."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num EQ 72) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.312."+string(m-1)+".fits",/remove_all))
            IF ( galaxy[j].num EQ 312) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.397."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num EQ 397) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.852."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num EQ 852) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.1."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num EQ 1) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.427."+string(m-1)+".fits",/remove_all))
            IF ( galaxy[j].num EQ 427) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.284."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num  EQ 284) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.477."+string(m-1)+".fits",/remove_all))
            IF ( galaxy[j].num EQ 477) THEN spawn,  mvcomm, exit_status = y
            mvcomm = strcompress("mv model.fits " + strcompress("model.197."+string(m-1)+".fits",/remove_all))
            IF (galaxy[j].num EQ 197) THEN spawn,  mvcomm, exit_status = y

            ;reset the good parameter in case the fit is really good this time
            galaxy[j].good = 0

            ;if the fit was crappy flip the switch on the "good" parameter
            IF (shifts[2] GT 99) OR (shifts[3] LT 0.09)THEN BEGIN
                openw, dlun2, "diagnostics.log", /GET_LUN, /append
                printf, dlun2, "NOT subtracting galaxy ",  galaxy[j].xcenter, "  ", galaxy[j].ycenter, $
                     " because re or b/a are off"
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
openw, gallun, "gallist", /GET_LUN, /append
FOR k = 0, numogalaxies - 1, 1 DO BEGIN
    fits_read, strcompress("best" + string(k) + ".fits", /remove_all), datagal, headergal
   ;if there are NAN's in the image (fit crappy) then zero the whole thing
    IF total(finite(datagal)) NE N_ELEMENTS(datagal) THEN BEGIN
        datagal[*,*] = 0.
        printf, gallun, galaxy[k].xcenter, galaxy[k].ycenter, k
                
    ENDIF

   allgals = allgals + datagal
ENDFOR
 close, gallun
free_lun, gallun
fits_write, "allgals.fits", allgals, headergal
fits_read, imagename, data, header
fits_write, "nogals.fits", data - allgals, header

close, /all
undefine, galaxy

END



            ;manually make a mask image using the galaxy properties from SExtractor
            ;maybe nice to rerun SExtractor?
;            fits_read, segmentationname, segdata, segheader

            ;blank the image
;           segdata(*,*) = 0.0
;           maskimage = segdata
           ;add masks that are the shapes of the galaxies (SEX)
           ;since the shapes from galfit could be wrong at this stage
           
           ;for all galaxies except the one working on
;           FOR r = 0,numogalaxies - 1, 1 DO BEGIN
;               IF (r NE j AND galaxy[r].good EQ 0) THEN BEGIN
                   ;make masks the shapes of the galaxies
                   ;with half the semi-major and semi-minor axes
;                   newimg = del_star(galaxy[r].theta, 40000.0, 0.7 *galaxy[r].a, 0.7 *galaxy[r].b, galaxy[r].xcenter, $
;                                  galaxy[r].ycenter, segdata)
                   
                   ;add the individual mask image to the others
;                   maskimage = maskimage + newimg
;               ENDIF
;           ENDFOR

           ;re-assaign values to be input for galfit
;           FOR x = 1, nxaxis - 1, 1 DO BEGIN   
;               FOR y = 1,nyaxis - 1, 1 DO BEGIN
;                   IF (maskimage(x,y) LT 0 ) THEN maskimage(x,y) = 1
                  
;               ENDFOR
;           ENDFOR

;            maskname = "mask.fits"
;            fits_write, maskname, maskimage, segheader


;try not including those galaxies with near neighbors
;newgalaxy = galaxy
;k = 0
;FOR j = 0, i-1 DO BEGIN
;        neigh = 0
;    FOR other = 0, i-1 DO BEGIN
                                ;the ones that are near each other...
;        IF (j NE other)  AND (galaxy[j].xcenter GT galaxy[other].xcenter - 23.) $
;           AND (galaxy[j].xcenter LT galaxy[other].xcenter + 23.) $
;           AND (galaxy[j].ycenter GT galaxy[other].ycenter - 23.) $
;           AND (galaxy[j].ycenter LT galaxy[other].ycenter + 23.) THEN BEGIN
;            print, "galaxy ", j, " has a neighbor, ignoring it" 
;            print, galaxy[j].xcenter, galaxy[j].ycenter, galaxy[other].xcenter, galaxy[other].ycenter
;            neigh = 1
;        ENDIF


;    ENDFOR
;    IF (neigh EQ 0) THEN BEGIN
;        newgalaxy[k] = galaxy[j]
;        k = k + 1
;    ENDIF

;ENDFOR
;numogalaxies = k
;galaxy = newgalaxy[0:numogalaxies - 1]
