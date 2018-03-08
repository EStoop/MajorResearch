PRO simulatemulti

close, /all
;want a statistical sample

numogals = 4
fakegalaxy = replicate({fobject, fnum:0, fxcenter:0D, fycenter:0D, fmag:0D, fre:0D, faxisratio:0D, fpa:0D}, numogals)

fmagnitude = findgen(30)/5. + 18.
fregal = findgen(37)/2. + 0.1
fpagal = findgen(90)*2.
faxisrat = findgen(19)/20 + 0.1
fdskyx = findgen(20)
fdskyx[0:9] = fdskyx[0:9] *1E-4+ 1E-4
fdskyx[10:19] = -(findgen(10) * 1E-4  + 1E-4) 
fdskyy= fdskyx

openr, namelun, "/n/Godiva7/jkrick/galfit/listbkgnames", /get_lun
numonames = 0
junk = " "
WHILE (NOT eof(namelun)) DO BEGIN
    READF, namelun, junk
    numonames = numonames+ 1
ENDWHILE
close, namelun

fbkg = strarr(numonames)

openr, namelun, "/n/Godiva7/jkrick/galfit/listbkgnames", /get_lun
junk = " "
i = 0
WHILE (NOT eof(namelun)) DO BEGIN
    READF, namelun, junk
    fbkg(i) = junk
    i = i + 1
ENDWHILE
close, namelun


seed = 20329

FOR dist = 10, 36, 2 DO BEGIN
    
   FOR k = 3, 20, 1 DO BEGIN
;dist = 0
;k = 0 
       spawn, "rm /n/Godiva7/jkrick/galfit/multiexp/diagnostics.log"
        
        FOR a = 0, numogals - 1, 1 DO BEGIN
            randomnumber = randomu(seed, 7)
            fakegalaxy[a]= {fobject, a, 0,0,fmagnitude(fix(randomnumber(0)*30)), fregal(fix(randomnumber(1)*37)), $
             faxisrat(fix(randomnumber(3)*19)), fpagal(fix(randomnumber(2)*90))}
            
        ENDFOR
        fakegalaxy[0].fxcenter = 160 + dist
        fakegalaxy[0].fycenter = 160 + dist
        fakegalaxy[1].fxcenter = 160 + dist
        fakegalaxy[1].fycenter = 240 - dist
        fakegalaxy[2].fxcenter = 240 - dist
        fakegalaxy[2].fycenter = 160 + dist
        fakegalaxy[3].fxcenter = 240 - dist
        fakegalaxy[3].fycenter = 240 - dist
        
        xmin = 1
        ymin = 1
        xmax = 400
        ymax = 400
        
        convolx = 60
        convoly = 60
        
        bkgd = 2.0
        fdx = fdskyx(fix(randomnumber(4)*20))
        fdy = fdskyy(fix(randomnumber(5)*20))

;_________________________________________________________________
;input file to make fake galaxy image
;-----------------------------------------------------------------------------------        
        openw, outlun, "/n/Godiva7/jkrick/galfit/multiexp/simulate.pro.feedme", /GET_LUN
        
        printf, outlun, "#IMAGE PARAMETERS"
        printf, outlun, "A) ", "none" ;, "   # Input data image (FITS file)"
        printf, outlun, "B) ", "/n/Godiva7/jkrick/galfit/multiexp/fake.fits" ;, "     # Name for the output image"
        printf, outlun, "C) ", "/n/Godiva7/jkrick/galfit/multiexp/noise400.fits" ;, "   # Noise image name(made from data if blank or none) "
        printf, outlun, "D) ", "/n/Godiva7/jkrick/galfit/multiexp/psf.small.fits" ;,"            # Input PSF image for convolution (FITS file)"
        printf, outlun, "E) ", "none" ;, "           # Pixel mask (ASCII file or FITS file with non-0 values)"
        printf, outlun, "F) /n/Godiva7/jkrick/galfit/multiexp/const.test	       # Parameter constraint file (ASCII)"
        printf, outlun, "G) ",xmin ,xmax ,ymin ,ymax ;,"  # Image region to fit (xmin xmax ymin ymax)"
        printf, outlun, "I) ", convolx, convoly , "            # Size of convolution box (x y)"
        printf, outlun, "J) 25.0              # Magnitude photometric zeropoint "
        printf, outlun, "K) 0.259   0.259         # Plate scale (dx dy)  [arcsec/pix.  Only for Nuker.]"
        printf, outlun, "O) both                # Display type (regular, curses, both)"
        printf, outlun, "P) 0                   # Create output image only? (1=yes; 0=optimize) "
        printf, outlun, "S) 0		       # Modify/create objects interactively?"
        
        printf, outlun, "###########################################################################"
        
        printf, outlun, " 0) sky"
        printf, outlun, " 1)", bkgd, "       1       # sky background       [ADU counts]"
        printf, outlun, " 2) 0.0    1       # dsky/dx (sky gradient in x) "
        printf, outlun, " 3) 0.0    1       # dsky/dy (sky gradient in y) "
        printf, outlun, " Z) 0                  # output image"
        
        
        printf, outlun, "###########################################################################"

        openw, dlun3, "/n/Godiva7/jkrick/galfit/multiexp/diagnostics.log", /GET_LUN, /append
        
        FOR b = 0, numogals - 1, 1 DO BEGIN
            printf, dlun3, "real galaxy",  fakegalaxy[b], " ", fdx, " " , fdy
            
            
            printf, outlun, " 0) devauc             # Object type"
            printf, outlun, " 1) " , fakegalaxy[b].fxcenter,  "  ", fakegalaxy[b].fycenter, " 1 1    # position x, y        [pixel]"
            printf, outlun, " 3) " ,fakegalaxy[b].fmag,"        1       # total magnitude    "
            printf, outlun, " 4) ", fakegalaxy[b].fre, "        1       #     R_e              [Pixels]" ;fudge factor
            printf, outlun, " 8) ", fakegalaxy[b].faxisratio, "       1       # axis ratio (b/a)   "
            printf, outlun, " 9) ", fakegalaxy[b].fpa, " 1       # position angle (PA)  [Degrees: Up=0, Left=90]"
            printf, outlun," 10) 0          0       # diskiness (< 0) or boxiness (> 0)"
            printf, outlun, " Z) 0                  # output image (see above)"
            printf, outlun, "###########################################################################"
            
        ENDFOR
        free_lun, outlun
        close, outlun
         spawn, "/n/Godiva7/jkrick/galfit/galfit/galfit  /n/Godiva7/jkrick/galfit/multiexp/simulate.pro.feedme"
        
        spawn, "rm galfit.??"
        
        fits_read, "/n/Godiva7/jkrick/galfit/multiexp/fake.fits", fakedata, fakeheader
        fits_read, "/n/Godiva7/jkrick/galfit/multiexp/noise400.fits", noisedata, noiseheader
        name = fbkg(fix(randomnumber(6)*numonames))
        fits_read, name, noisedata, noiseheader
        printf, dlun3, name

        fakedata = fakedata + noisedata
        
        fits_write, "/n/Godiva7/jkrick/galfit/multiexp/fake.fits", fakedata, fakeheader
       free_lun, dlun3
        close, dlun3
        
;        galfit , "/n/Godiva7/jkrick/galfit/fake.fits", "/n/Godiva7/jkrick/galfit/noise400.fits", "/n/Godiva7/jkrick/galfit/psf.small.fits", "/n/Godiva7/jkrick/galfit/segmentation.fits"
 
;__________________________________________________________________________________
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;run SExtractor on the image to find the objects
;need to have name/location of segmentation image as well as input
;parameters set
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

        commandline = '/n/Godiva7/jkrick/Sex/sex  /n/Godiva7/jkrick/galfit/multiexp/fake.fits  -c /n/Godiva7/jkrick/galfit/multiexp/default.sex'
        spawn, commandline
        

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;read in the galaxy information for all galaxies
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        outputname = "/n/Godiva7/jkrick/galfit/multiexp/model.fits"

        fits_read, "/n/Godiva7/jkrick/galfit/multiexp/segmentation.fits", segdata, segheader
;figure out the size of the image from the header information
        junk = strsplit( segheader[3], /extract)
        nxaxis = fix(junk[2])
        junk = strsplit( segheader[4], /extract)
        nyaxis = fix(junk[2])
        
        
;figure out how many SExtractor objects there are
;must be a better way of doing this
        OPENR, lun1, '/n/Godiva7/jkrick/galfit/multiexp/SExtractor.cat', /GET_LUN, ERROR = err
        IF (err NE 0) then PRINT, "catalog file did not open"
        numoobjects = 0
        junk = " "
        WHILE (NOT eof(lun1)) DO BEGIN
            READF, lun1, junk
            numoobjects= numoobjects + 1
        ENDWHILE
        free_lun, lun1

        close, lun1
        
;make an array of galaxy objects
        galaxy2 = replicate({object, num:0, xcenter:0D, ycenter:0D, a:0D, b:0D, ellip:0D, flux:0D, $
     mag:0D, isoarea:0D, fwhm:0D, theta:0D, bkgd:0D, good:0},numoobjects)
        
        
;read in galaxy parameters into galaxy objects
        OPENR, lun1, '/n/Godiva7/jkrick/galfit/multiexp/SExtractor.cat', /GET_LUN, ERROR = err
        IF (err NE 0) then PRINT, "catalog file did not open"
        
        i = 0
        WHILE (NOT EOF(lun1)) DO BEGIN
            READF, lun1, number, xcen, ycen, alength, blength, ellipticity, fluxbest, magnitude, $
         isoareagal, fwhmin, thetaimage, background
            
                                ;if I think it is a galaxy, add it to the list of galaxies
                                ;also weed out small junky things, and things on edges
            IF(fwhmin GT 4.0) AND (xcen GT 100) AND (xcen LT nxaxis - 100) AND (ycen GT 100) $
               AND (ycen LT nyaxis - 100) AND (magnitude LT 23.9) AND (isoareagal GT 40) $  ;mag 22.6
               THEN BEGIN 
;               print, "all the params before going into galaxies ", number,
;                 xcen, ycen, alength, blength, ellipticity, fluxbest,  magnitude, $
;                  isoareagal, fwhmin, thetaimage, background
                galaxy2[i] = {object, float(number), float(xcen), float(ycen), float(alength), float(blength), $
                              float(ellipticity), float(fluxbest), float(magnitude), float(isoareagal), float(fwhmin), $
                              float(thetaimage - 90.), float(background), 0}
                i = i +1
            ENDIF
        ENDWHILE
        free_lun, lun1
        close, lun1

        ;want to shorten the array of galaxies to only be as long as it needs to be
        numogalaxies = i
        galaxy2 = galaxy2[0:numogalaxies -1]
        
;________________________________________________________________
;input parameters as beginning values to try to fit fake galaxies
;----------------------------------------------------------------------------------
        
        openw, outlun, "/n/Godiva7/jkrick/galfit/multiexp/galfit2.pro.feedme", /GET_LUN
        
        ;need to change all of these names to real files 
        printf, outlun, "#IMAGE PARAMETERS"
        printf, outlun, "A) ","/n/Godiva7/jkrick/galfit/multiexp/fake.fits" ;, "   # Input data image (FITS file)"
        printf, outlun, "B) ", "/n/Godiva7/jkrick/galfit/multiexp/junk.fits" ;, "     # Name for the output image"
        printf, outlun, "C) ", "/n/Godiva7/jkrick/galfit/multiexp/noise400.fits" ;, "   # Noise image nam
        printf, outlun, "D) ", "/n/Godiva7/jkrick/galfit/multiexp/psf.small.fits" ;,"            # Input PSF image 
;       printf, outlun, "E) ", "/n/Godiva7/jkrick/galfit/multiexp/segmentation.fits";, "     # Pixel mask (
        printf, outlun, "F) /n/Godiva7/jkrick/galfit/multiexp/const.test	       # Parameter constraint file (ASCII)"
        printf, outlun, "G) ",100 ,300 ,100,300 ;,"  # Image region to fit (xmin xmax ymin ymax)"
        printf, outlun, "I) ", 60, 60 , "            # Size of convolution box (x y)"
        printf, outlun, "J) 25.0              # Magnitude photometric zeropoint "
        printf, outlun, "K) 0.259   0.259         # Plate scale (dx dy)  [arcsec/pix.  Only for Nuker.]"
        printf, outlun, "O) both                # Display type (regular, curses, both)"
        printf, outlun, "P) 0                   # Create output image only? (1=yes; 0=optimize) "
        printf, outlun, "S) 0		       # Modify/create objects interactively?"
        
        printf, outlun, "###########################################################################"
        
        FOR p = 0,numogalaxies  - 1, 1 DO BEGIN
            bkgdval = galaxy2[p].bkgd
            
            re = 0.2*sqrt(galaxy2[p].isoarea) ;initial guess
            axisratio = galaxy2[p].b / galaxy2[p].a
            
            
            dx = 0
            dy = 0
            
            
            printf, outlun, " 0) devauc             # Object type"
            printf, outlun, " 1) " , galaxy2[p].xcenter,  "  ", galaxy2[p].ycenter, " 1 1    # position x, y        [pixel]"
            printf, outlun, " 3) " ,galaxy2[p].mag,"        1       # total magnitude    "
            printf, outlun, " 4) ", re, "        1       #     R_e              [Pixels]" ;fudge factor
            printf, outlun, " 8) ", axisratio, "       1       # axis ratio (b/a)   "
            printf, outlun, " 9) ", galaxy2[p].theta , " 1       # position angle (PA)  [Degrees: Up=0, Left=90]"
            printf, outlun," 10) 0          0       # diskiness (< 0) or boxiness (> 0)"
            printf, outlun, " Z) 0                  # output image (see above)"
            
            printf, outlun, "###########################################################################"
            
        ENDFOR
        
        printf, outlun, " 0) sky"
        printf, outlun, " 1) ", bkgdval, "       1       # sky background       [ADU counts]"
        printf, outlun, " 2) 0.0005     1       # dsky/dx (sky gradient in x) "
        printf, outlun, " 3) 0.0005     1       # dsky/dy (sky gradient in y) "
        printf, outlun, " Z) 0                  # output image"
        free_lun, outlun
        close, outlun
        
         ;run galfit
        spawn, '/n/Godiva7/jkrick/galfit/galfit/galfit  /n/Godiva7/jkrick/galfit/multiexp/galfit2.pro.feedme ', exit_status = y
;        IF (y GT 0) THEN BREAK  ;galfit has crashed-----------------!
        
;))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))))
        
;        OPENR, lun2, 'galfit.01', /GET_LUN, ERROR = err
;        IF (err NE 0) then PRINT, "galfit file did not open"
        string = " "
        xcentervals = fltarr(12)
        ycentervals = fltarr(12)
        magvals = fltarr(12)
        revals=fltarr(12)
        axisratiovals = fltarr(12)
        thetavals = fltarr(12)
        dsxvals = fltarr(12)
        dsyvals = fltarr(12)

        spawn, "rm junkm"
        spawn, "tail -20 fit.log > junkm"
        openr, inlun, "junkm", /GET_LUN
        r = 0
 
        WHILE (NOT EOF(inlun)) DO BEGIN
            READF, inlun, string
            filestring = strsplit(string, /extract)
            s= N_elements(filestring)
            IF (s GT 1) THEN BEGIN
                IF (filestring[0] EQ "Chi^2/nu") THEN chi2 = filestring[2]
                                ;not working, but don't care about it anyway
                chi2 = 0.0
                
                IF (filestring[0] EQ "devauc") THEN BEGIN
 
                    blah = strsplit(filestring[2], "(", /extract)
                    blah2 = strsplit(blah[0], ",", /extract)
                    xcentervals[r]= blah2[0]

                    blah3 = strsplit(filestring[3], ")", /extract)
                    ycentervals[r]= blah3[0]
                    magvals[r] = filestring[4]
                    revals[r] = filestring[5]
                    axisratiovals[r]= filestring[6]
                    thetavals[r] = filestring[7]
                    r = r + 1

                ENDIF
                
                IF (filestring[0] EQ "sky") THEN BEGIN
                    bkgd = filestring[4]
                    dsxvals[*] = filestring[5]
                    dsyvals[*] = filestring[6]
                ENDIF
            ENDIF

        ENDWHILE
        free_lun, inlun
        close, inlun
        IF (re EQ 0) THEN re = 0.5
        

  

;-------------------------------------------------------------------------------
;output final parameters to a file for post galfit diagnostics
;-------------------------------------------------------------------------------

        openw, dlun, "/n/Godiva7/jkrick/galfit/multiexp/diagnostics.log", /GET_LUN, /append
        FOR t = 0, r, 1 DO BEGIN
            IF (xcentervals[t] NE 0) THEN BEGIN
                printf, dlun, xcentervals[t], "  ", ycentervals[t], " ", 0.0, " ", magvals[t], " ", revals[t], " ", axisratiovals[t], $
                        " ", thetavals[t], " ", bkgd, " ", dsxvals[t], " " , dsyvals[t]
            ENDIF


        ENDFOR

        free_lun, dlun
        close, dlun
        
        
;---------------------------------------------------------------------------------
;also want to run galfit one more time to make an image without
;galaxies  and without background
;
;so take parameters from best fit, and make a best image without
;background values.  then subtract that image from the original fake
;image
;---------------------------------------------------------------------------------
        
        openw, outlun, "/n/Godiva7/jkrick/galfit/multiexp/junk.feedme", /GET_LUN
        
        ;need to change all of these names to real files 
        printf, outlun, "#IMAGE PARAMETERS"
;        printf, outlun, "A) ","/n/Godiva7/jkrick/galfit/multiexp/fake.fits" ;, "   # Input data image (FITS file)"
        printf, outlun, "B) ", "/n/Godiva7/jkrick/galfit/multiexp/gals.fits" ;, "     # Name for the output image"
        printf, outlun, "C) ", "/n/Godiva7/jkrick/galfit/multiexp/noise400.fits" ;, "   # Noise image nam
        printf, outlun, "D) ", "/n/Godiva7/jkrick/galfit/multiexp/psf.small.fits" ;,"            # Input PSF image 
;       printf, outlun, "E) ", "/n/Godiva7/jkrick/galfit/multiexp/segmentation.fits";, "     # Pixel mask (
        printf, outlun, "F) /n/Godiva7/jkrick/galfit/multiexp/const.test	       # Parameter constraint file (ASCII)"
        printf, outlun, "G) ",1 ,400 ,1,400 ;,"  # Image region to fit (xmin xmax ymin ymax)"
        printf, outlun, "I) ", 60, 60 , "            # Size of convolution box (x y)"
        printf, outlun, "J) 25.0              # Magnitude photometric zeropoint "
        printf, outlun, "K) 0.259   0.259         # Plate scale (dx dy)  [arcsec/pix.  Only for Nuker.]"
        printf, outlun, "O) both                # Display type (regular, curses, both)"
        printf, outlun, "P) 0                   # Create output image only? (1=yes; 0=optimize) "
        printf, outlun, "S) 0		       # Modify/create objects interactively?"
        
        printf, outlun, "###########################################################################"
        
        FOR p = 0,numogalaxies  - 1, 1 DO BEGIN
              
            printf, outlun, " 0) devauc             # Object type"
            printf, outlun, " 1) " , xcentervals[p],  "  ", ycentervals[p], " 1 1    # position x, y        [pixel]"
            printf, outlun, " 3) " ,magvals[p],"        1       # total magnitude    "
            printf, outlun, " 4) ", revals[p], "        1       #     R_e              [Pixels]" ;fudge factor
            printf, outlun, " 8) ", axisratiovals[p], "       1       # axis ratio (b/a)   "
            printf, outlun, " 9) ", thetavals[p] , " 1       # position angle (PA)  [Degrees: Up=0, Left=90]"
            printf, outlun," 10) 0          0       # diskiness (< 0) or boxiness (> 0)"
            printf, outlun, " Z) 0                  # output image (see above)"
            
            printf, outlun, "###########################################################################"
            
        ENDFOR
        
        printf, outlun, " 0) sky"
        printf, outlun, " 1) 0.0      0      # sky background       [ADU counts]"
        printf, outlun, " 2) 0.000     0       # dsky/dx (sky gradient in x) "
        printf, outlun, " 3) 0.000     0       # dsky/dy (sky gradient in y) "
        printf, outlun, " Z) 0                  # output image"
        free_lun, outlun
        close, outlun
        
;-------------------------------------------------------------------------------------------------------------------

        spawn, "/n/Godiva7/jkrick/galfit/galfit/galfit /n/Godiva7/jkrick/galfit/multiexp/junk.feedme"

        fits_read, "/n/Godiva7/jkrick/galfit/multiexp/gals.fits", galsdata, galsheader
;        fits_read, "/n/Godiva7/jkrick/galfit/multiexp/noise2.fits",noisedata, noiseheader
        fits_read, "/n/Godiva7/jkrick/galfit/multiexp/fake.fits", fakedata, fakeheader
        residualdata = fakedata  -  galsdata 

        fits_write, strcompress("/n/Godiva7/jkrick/galfit/multiexp/nogals" + string(k) + "." + string(dist)+ ".fits", $
                                /remove_all), residualdata, fakeheader
        newname = strcompress("/n/Godiva7/jkrick/galfit/multiexp/diagnostics" + string(k)  + "." + string(dist) + ".log", /remove_all)
        newcommand = "mv /n/Godiva7/jkrick/galfit/multiexp/diagnostics.log "+  newname
        ;newimg = strcompress("/n/Godiva7/jkrick/galfit/multiexp/fake" + string(k) + "." + string(dist)+ ".fits", /remove_all)
        ;newcommand2 = "mv /n/Godiva7/jkrick/galfit/multiexp/fake.fits " + newimg
        ;newimg2 = strcompress("/n/Godiva7/jkrick/galfit/multiexp/gals" + string(k) + "." + string(dist)+ ".fits", /remove_all)
        ;newcommand3 = "mv /n/Godiva7/jkrick/galfit/multiexp/junk.fits " + newimg2
        
        spawn, newcommand
                                ;spawn, newcommand2
        ;spawn, newcommand3
        
        spawn, "rm galfit.??"
        
    ENDFOR
ENDFOR


END

;        WHILE (NOT EOF(lun2)) DO BEGIN
;            READF, lun2, string
;            filestring = strsplit(string, /extract)
;            s= N_elements(filestring)
            
;            IF (s GT 1) THEN BEGIN

;                IF (filestring[1] EQ "devauc") THEN BEGIN
                    
;                    readf, lun2, string
;                    vals = strsplit(string, /extract)
;                    xcentervals[r] = vals[1]
;                    ycentervals[r] = vals[2]
                    
;                    readf, lun2, string
;                    vals = strsplit(string, /extract)
;                    magvals[r] = vals[1]
                    
;                    readf, lun2, string
;                    vals = strsplit(string, /extract)
;                    revals[r] = vals[1]
                    
;                    IF (re EQ 0) THEN re = 0.5
;                    
;                    readf, lun2, string
;                    readf, lun2, string
;                    readf, lun2, string
;                    
;                    readf, lun2, string
;                    vals = strsplit(string, /extract)
;                    axisratiovals[r] = vals[1]
;                    
;                    readf, lun2, string
;                    vals = strsplit(string, /extract)
;                    thetavals[r] = vals[1]
;                    r = r + 1
;                ENDIF
;                
;                IF (filestring[1] EQ "sky") THEN BEGIN
;                    readf, lun2, string
;                    vals = strsplit(string, /extract)
;                    bkgd = vals[1]
;                ENDIF
;                
                                ;now that I know what the bkgd value is, don't use it
                                ;bkgd = 0
;            ENDIF

;        ENDWHILE
            
;            free_lun, lun2
;            close, lun2
