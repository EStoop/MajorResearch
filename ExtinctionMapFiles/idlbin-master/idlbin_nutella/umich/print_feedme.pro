FUNCTION print_feedme, imagename, outputname,noiseimagename, psfname, segmentationname, isoarea, theta, a, b, xcenter, ycenter, mag, bkgd, nxaxis, nyaxis

;-------------------------------------------------------------------------------------------------------
;This function prints out an input file for galfit
;it is called by galfit.pro
;
;there are 2 modes
;1)makes  a guess at the parameters given the SExtractor parameters
;and then writes the input file for galfit to fit for the parameters
;
;2)takes the best fit parameters from a previous run of galfit and
;outputs another input file which will just make a larger version of
;the model
;---------------------------------------------------------------------------------------------------------

;is there some smarter way to figure what this should be?
convolx = 50
convoly = 50

IF (outputname EQ "bigmodel.fits") THEN BEGIN
     ;read back in the parameters from the last fit, and make them into a
     ;model over the whole image, don't include the background
;    OPENR, lun2, 'galfit.01', /GET_LUN, ERROR = err
;    IF (err NE 0) then PRINT, "galfit file did not open"
   string = " "
;    
    bkgd = 0.
    dsx = 0.
    dsy  = 0.

;try reading in the parameters from fit.log instead of galfit.01
    spawn, "rm junk"
    spawn, "tail -15 fit.log > junk"
    openr, inlun, "junk", /GET_LUN
    WHILE (NOT EOF(inlun)) DO BEGIN
        READF, inlun, string
        filestring = strsplit(string, /extract)
        s= N_elements(filestring)
        
        IF (s GT 1) THEN BEGIN
            ;IF (filestring[0] EQ "Chi^2/nu") THEN chi2 = filestring[2]
           ;not working, but don't care about it anyway
            chi2 = 0.0
            
            IF (filestring[0] EQ "devauc") THEN BEGIN
                blah = strsplit(filestring[2], "(", /extract)
                blah2 = strsplit(blah[0], ",", /extract)
                xcenter = blah2[0]
                blah3 = strsplit(filestring[3], ")", /extract)
                ycenter = blah3[0]
                mag = filestring[4]
                re = filestring[5]
                axisratio= filestring[6]
                theta = filestring[7]
            ENDIF
            
            IF (filestring[0] EQ "sersic") THEN BEGIN
                blah = strsplit(filestring[2], "(", /extract)
                blah2 = strsplit(blah[0], ",", /extract)
                xcenter = blah2[0]
                blah3 = strsplit(filestring[3], ")", /extract)
                ycenter = blah3[0]
                mag = filestring[4]
                re = filestring[5]
                nindex = filestring[6]
                axisratio= filestring[7]
                theta = filestring[8]
            ENDIF

           IF (filestring[0] EQ "expdisk") THEN BEGIN
                blah = strsplit(filestring[2], "(", /extract)
                blah2 = strsplit(blah[0], ",", /extract)
                xcentere = blah2[0]
                blah3 = strsplit(filestring[3], ")", /extract)
                ycentere  = blah3[0]
                mage = filestring[4]
                rs = filestring[5]
                axisratioe= filestring[6]
                thetae = filestring[7]
            ENDIF

            IF (filestring[0] EQ "sky") THEN BEGIN
                bkgd = filestring[4]
                dsx = filestring[5]
                dsy = filestring[6]
            ENDIF
        ENDIF
    ENDWHILE
    free_lun, inlun
    close, inlun
    print, re, " re"
    IF (re LT 0.5) THEN re = 0.5

                ;figure out how large to make the model galaxy
                ;guessing at 15, should look in the literature for this
    xmax = fix(xcenter + 15.*re)
    ymax = fix(ycenter + 15.*re)
    xmin = fix(xcenter - 15.*re)
    ymin = fix(ycenter -15.*re)
    IF (xmin LT 0) THEN xmin = 0
    IF (ymin LT 0) THEN ymin = 0
    IF (xmax GT nxaxis) THEN xmax = nxaxis - 1
    IF (ymax GT nyaxis) THEN ymax = nyaxis - 1
                                   
                ;then need to figure out shift 
                ;to subtract from larger image
    IF (xmin NE 0) THEN dx = xmin - 1 ELSE dx = 0
    IF (ymin NE 0) THEN dy = ymin - 1 ELSE dy = 0            

            
;    print, xcenter, "  ", ycenter, " ", chi2, " ", mag, " ", re, " ",nindex, " " , axisratio, " ", theta, $
;       "  ", bkgd, "  ",  dsx, "  ", dsy
 ;  print, xcenter, "  ", ycenter, " ", chi2, " ", mag, " ", re, " ", axisratio, " ", theta, $
 ;      "  ", bkgd, "  ",  dsx, "  ", dsy

                
;-------------------------------------------------------------------------------
;output final parameters to a file for post galfit diagnostics
;-------------------------------------------------------------------------------

    openw, dlun, "diagnostics.log", /GET_LUN, /append
;    printf, dlun, xcenter, "  ", ycenter, " ", chi2, " ", mag, " ", re, " ",nindex, " " , axisratio, " ", theta, $
;       "  ", bkgd, "  ",  dsx, "  ", dsy
   printf, dlun, xcenter, "  ", ycenter, " ", chi2, " ", mag, " ", re, " " , axisratio, " ", theta, $
       "  ", bkgd, "  ",  dsx, "  ", dsy
    free_lun, dlun
    close, dlun

    ;now that I know what the bkgd value is, don't use it
    ;I don't want to subtract the background
    bkgd = 0.
    dsx = 0.
    dsy  = 0.

ENDIF ELSE BEGIN
    
   ;figure out how large the fitting region should be, as per Peng's prescription
   ; majx = abs(sqrt(isoarea) * 4. * sin (theta/180.*!PI)) 
   ; majy = abs(sqrt(isoarea) * 4. * cos (theta/180.*!PI)) 
   ; minx = (b /a)* sqrt(isoarea) * 4. * abs(cos ((theta+!PI/2)/180.*!PI)) 
   ; miny = (b /a)* sqrt(isoarea) * 4. * abs(sin ((theta+!PI/2)/180.*!PI)) 
    
   ; xsize = majx                    
   ; if (minx GE majx) THEN xsize = minx                                                 
   ; if (xsize LE 30.) THEN xsize = 30.                     
    
    ;ysize = majy                    
    ;if (miny GE majy) THEN ysize = miny                                 
    ;if (ysize LE 30.) THEN ysize = 30.                    

    ;make the thumbnail just big enough to encompass a and b
    npoints = 120
    phi = 2 * !PI * (Findgen(npoints) / (npoints-1))
    xsize = a*cos(phi)
    ysize = b*sin(phi)
    
    xprime = (xsize*cos(theta/180.*!PI)) - (ysize*sin(theta/180.*!PI))
    yprime = (xsize*sin(theta/180.*!PI)) - (ysize*cos(theta/180.*!PI))
    
    xsortindex = sort(xprime)
    xsort = xprime[xsortindex]
    ysortindex = sort(yprime)
    ysort = yprime[ysortindex]
         
    xrad = xsort(N_ELEMENTS(xprime) - 1)
    yrad = ysort(N_ELEMENTS(yprime) - 1)

    ;impose summer lower limit to the radius of the fitting box
    IF (xrad LT 5) THEN xrad = 7
    IF (yrad LT 5) THEN yrad = 7

    xmin = fix(xcenter - 1.5*xrad)
    xmax = fix(xcenter + 1.5*xrad)
    ymin= fix(ycenter - 1.5*yrad)
    ymax = fix(ycenter + 1.5*yrad)

    re = 0.1*sqrt(isoarea)                   ;initial guesses
    axisratio = b / a
    xcentere = xcenter
    ycentere = ycenter
    mage = mag
    axisratioe = axisratio
    thetae = theta
    rs = re  
    
    dx = 0
    dy = 0

    ;need to make an initial guess at dsky's
    ;is there a better way, SEx?
    dsx = 0.0005
    dsy = 0.0005
    

    ;initial guess at n = average of the distributio
    nindex = 3.

    print, "dx, dy", dx, dy
    print, "xcenter, ycenter, xrad, yrad", xcenter, ycenter, xrad, yrad

ENDELSE

openw, outlun, "galfit.pro.feedme", /GET_LUN

printf, outlun, "#IMAGE PARAMETERS"
printf, outlun, "A) ", imagename;, "   # Input data image (FITS file)"
printf, outlun, "B) ", outputname;, "     # Name for the output image"
printf, outlun, "C) ", noiseimagename;, "   # Noise image name(made from data if blank or none) "
printf, outlun, "D) ", psfname;,"            # Input PSF image for convolution (FITS file)"
printf, outlun, "E) ", segmentationname;, "           # Pixel mask (ASCII file or FITS file with non-0 values)"
printf, outlun, "F) /n/Godiva1/jkrick/A3888/galfit/const.test";	       # Parameter constraint file (ASCII)"
printf, outlun, "G) ",xmin ,xmax ,ymin ,ymax;,"  # Image region to fit (xmin xmax ymin ymax)"
printf, outlun, "I) ", convolx, convoly , "            # Size of convolution box (x y)"
printf, outlun, "J) 24.3           # Magnitude photometric zeropoint "
printf, outlun, "K) 0.259   0.259         # Plate scale (dx dy)  [arcsec/pix.  Only for Nuker.]"
printf, outlun, "O) both                # Display type (regular, curses, both)"
printf, outlun, "P) 0                   # Create output image only? (1=yes; 0=optimize) "
printf, outlun, "S) 0		       # Modify/create objects interactively?"

printf, outlun, "###########################################################################"

printf, outlun, " 0) devauc             # Object type"
;printf, outlun, " 0) sersic            # Object type"
printf, outlun, " 1) " , xcenter,  "  ", ycenter, " 1 1    # position x, y        [pixel]"
printf, outlun, " 3) " ,mag,"        1       # total magnitude    "
printf, outlun, " 4) ", re, "        1       #     R_e              [Pixels]";fudge factor
;printf, outlun, " 5) ", nindex,"       1       # Sersic exponent (deVauc=4, expdisk=1)  "
printf, outlun, " 8) ", axisratio, "       1       # axis ratio (b/a)   "
printf, outlun, " 9) ", theta , " 1       # position angle (PA)  [Degrees: Up=0, Left=90]"
printf, outlun," 10) 0          0       # diskiness (< 0) or boxiness (> 0)"
printf, outlun, " Z) 0                  # output image (see above)"

printf, outlun, "###########################################################################"

;printf, outlun, " 0) expdisk            # Object type"
;printf, outlun, " 1) " , xcentere,  "  ", ycentere, " 1 1    # position x, y        [pixel]"
;printf, outlun, " 3) " ,mage,"        1       # total magnitude    "
;printf, outlun, " 4) ", rs, "        1       #     Rs            [Pixels]";fudge factor
;printf, outlun, " 8) ", axisratioe, "       1       # axis ratio (b/a)   "
;printf, outlun, " 9) ", thetae , " 1       # position angle (PA)  [Degrees: Up=0, Left=90]"
;printf, outlun, " 10) 0          0       # diskiness (< 0) or boxiness (> 0)"
;printf, outlun, " Z)"

;printf, outlun, "###########################################################################"


printf, outlun, " 0) sky"
printf, outlun, " 1)", bkgd, "       1       # sky background       [ADU counts]"
printf, outlun, " 2) ", dsx, "     1      # dsky/dx (sky gradient in x) "
printf, outlun, " 3) ", dsy, "      1       # dsky/dy (sky gradient in y) "
printf, outlun, " Z) 0                  # output image"

close, outlun
close, /all
shifts = fltarr(4)
shifts = [dx,dy, float(re), float(axisratio)]
return, shifts
END



 
;    WHILE (NOT EOF(lun2)) DO BEGIN
;        READF, lun2, string
;        filestring = strsplit(string, /extract)
;        s= N_elements(filestring)
;    
;        IF (s GT 1) THEN BEGIN
;            IF (filestring[1] EQ "Chi^2/nu") THEN chi2 = filestring[3]
;
;            IF (filestring[1] EQ "devauc") THEN BEGIN
;                
;                readf, lun2, string
;                vals = strsplit(string, /extract)
;                xcenter = vals[1]
;                ycenter = vals[2]

;                readf, lun2, string
;                vals = strsplit(string, /extract)
;                mag = vals[1]

;                readf, lun2, string
;                vals = strsplit(string, /extract)
;                re = vals[1]

;               IF (re EQ 0) THEN re = 0.5

                ;figure out how large to make the model galaxy
                ;guessing at 15, should look in the literature for this
;                xmax = fix(xcenter + 15.*re)
;                ymax = fix(ycenter + 15.*re)
;                xmin = fix(xcenter - 15.*re)
;                ymin = fix(ycenter -15.*re)
;                IF (xmin LT 0) THEN xmin = 0
;                IF (ymin LT 0) THEN ymin = 0
;                IF (xmax GT nxaxis) THEN xmax = nxaxis - 1
;                IF (ymax GT nyaxis) THEN ymax = nyaxis - 1
                                   
                ;then need to figure out shift 
                ;to subtract from larger image
;                IF (xmin NE 0) THEN dx = xmin - 1 ELSE dx = 0
;                IF (ymin NE 0) THEN dy = ymin - 1 ELSE dy = 0

;                readf, lun2, string
;                readf, lun2, string
;               readf, lun2, string
                
;                readf, lun2, string
;                vals = strsplit(string, /extract)
;                axisratio = vals[1]

;                readf, lun2, string
;                vals = strsplit(string, /extract)
;                theta = vals[1]
;           ENDIF
            
;            IF (filestring[1] EQ "sky") THEN BEGIN
;                readf, lun2, string
;                vals = strsplit(string, /extract)
;                bkgd = vals[1]
                
;                readf, lun2, string
;                vals = float(strsplit(string, /extract))
;                dsx = vals[1]

;                readf, lun2, string
;                print, string
;                vals = float(strsplit(string, /extract))
;                print, vals[0], " ", vals[1], " ", vals[2]
;                dsy = vals[1]
;            ENDIF

;         ENDIF
        
;    ENDWHILE
