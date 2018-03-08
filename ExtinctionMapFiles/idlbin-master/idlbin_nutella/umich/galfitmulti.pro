PRO galfitmulti, imagename, noiseimagename, psfname, segmentationname

;example way of running this program
;galfit, "/n/Godiva7/jkrick/galfit/thumb.1.fits", "/n/Godiva7/jkrick/galfit/thumbnoise.1.fits", "/n/Godiva7/jkrick/galfit/psf.small.fits", "/n/Godiva7/jkrick/galfit/segmentation.fits"

close, /all

outputname = "/n/Godiva7/jkrick/galfit/model.fits"
fits_read, imagename, data, header
;spawn, "rm /n/Godiva7/jkrick/galfit/diagnostics.log"
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;run SExtractor on the image to find the objects
;need to have name/location of segmentation image as well as input
;parameters set
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

commandline = '/n/Godiva7/jkrick/Sex/sex ' + imagename + ' -c /n/Godiva7/jkrick/galfit/default.sex'
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
OPENR, lun1, '/n/Godiva7/jkrick/galfit/SExtractor.cat', /GET_LUN, ERROR = err
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
     mag:0D, isoarea:0D, fwhm:0D, theta:0D, bkgd:0D, good:0},numoobjects)

;read in galaxy parameters into galaxy objects
OPENR, lun1, '/n/Godiva7/jkrick/galfit/SExtractor.cat', /GET_LUN, ERROR = err
IF (err NE 0) then PRINT, "catalog file did not open"

i = 0
WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, number, xcen, ycen, alength, blength, ellipticity, fluxbest, magnitude, $
         isoareagal, fwhmin, thetaimage, background
   
    ;if I think it is a galaxy, add it to the list of galaxies
    ;also weed out small junky things, and things on edges
    IF(fwhmin GT 4.0) AND (xcen GT 100) AND (xcen LT nxaxis - 100) AND (ycen GT 100) $
               AND (ycen LT nyaxis - 100) AND (magnitude LT 24.0) AND (isoareagal GT 40) $  ;mag 22.6
               THEN BEGIN 
        galaxy[i] = {object, number, xcen, ycen, alength, blength, ellipticity, fluxbest, $
             magnitude, isoareagal, fwhmin, thetaimage - 90, background, 0}
        i = i +1
  ENDIF
ENDWHILE

numogalaxies = i
galaxy = galaxy[0:numogalaxies -1]

   ;--------------------------------------------------------------------------------
         ;fit each galaxy, one at a time, using the SExtractor segmentation
         ;image as a mask image
         ;--------------------------------------------------------------------------------
       openw, outlun, "multi.pro.feedme", /GET_LUN
        
        printf, outlun, "#IMAGE PARAMETERS"
        printf, outlun, "A) center.fits" ;, "   # Input data image (FITS file)"
        printf, outlun, "B) centermulti.fits" ;, "     # Name for the output image"
        printf, outlun, "C) centernoise.fits" ;, "   # Noise image name(made from data if blank or none) "
        printf, outlun, "D) psf.small.fits" ;,"            # Input PSF image for convolution (FITS file)"
        printf, outlun, "E)  none" ;, "           # Pixel mask (ASCII file or FITS file with non-0 values)"
        printf, outlun, "F) const.test	       # Parameter constraint file (ASCII)"
        printf, outlun, "G) 100 500 100 500  # Image region to fit (xmin xmax ymin ymax)"
        printf, outlun, "I) 80 80            # Size of convolution box (x y)"
        printf, outlun, "J) 25.0              # Magnitude photometric zeropoint "
        printf, outlun, "K) 0.259   0.259         # Plate scale (dx dy)  [arcsec/pix.  Only for Nuker.]"
        printf, outlun, "O) both                # Display type (regular, curses, both)"
        printf, outlun, "P) 0                   # Create output image only? (1=yes; 0=optimize) "
        printf, outlun, "S) 0		       # Modify/create objects interactively?"
        
        printf, outlun, "###########################################################################"
        
        printf, outlun, " 0) sky"
        printf, outlun, " 1) 2       1       # sky background       [ADU counts]"
        printf, outlun, " 2) 0.0    1       # dsky/dx (sky gradient in x) "
        printf, outlun, " 3) 0.0    1       # dsky/dy (sky gradient in y) "
        printf, outlun, " Z) 0                  # output image"
        
        
        printf, outlun, "###########################################################################"

         
        FOR j= 0, numogalaxies - 1, 1 DO BEGIN
             
            
            printf, outlun, " 0) devauc             # Object type"
            printf, outlun, " 1) " ,galaxy[j].xcenter,  "  ", galaxy[j].ycenter, " 1 1    # position x, y        [pixel]"
            printf, outlun, " 3) " ,galaxy[j].mag,"        1       # total magnitude    "
            printf, outlun, " 4) ",  0.1*sqrt(galaxy[j].isoarea), "        1       #     R_e              [Pixels]" ;fudge factor
            printf, outlun, " 8) ",galaxy[j].b / galaxy[j].a,"       1       # axis ratio (b/a)   "
            printf, outlun, " 9) ", galaxy[j].theta, " 1       # position angle (PA)  [Degrees: Up=0, Left=90]"
            printf, outlun," 10) 0          0       # diskiness (< 0) or boxiness (> 0)"
            printf, outlun, " Z) 0                  # output image (see above)"
            printf, outlun, "###########################################################################"
            
        ENDFOR
        free_lun, outlun
        close, outlun

            
   
    close, /all
undefine, galaxy
END





