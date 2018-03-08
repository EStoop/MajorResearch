PRO simulate6

close, /all
;want a statistical sample

numogals = 6
fakegalaxy = replicate({fobject, fnum:0, fxcenter:0D, fycenter:0D, fmag:0D, fre:0D, faxisratio:0D, $
     fpa:0D}, numogals)

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

;fakegalaxy[0].fxcenter = 210 
;fakegalaxy[0].fycenter = 133 
;fakegalaxy[1].fxcenter = 210
;fakegalaxy[1].fycenter = 273
;fakegalaxy[2].fxcenter = 141 
;fakegalaxy[2].fycenter = 175 
;fakegalaxy[3].fxcenter = 141 
;fakegalaxy[3].fycenter = 231 
;fakegalaxy[4].fxcenter = 279 
;fakegalaxy[4].fycenter = 175 
;fakegalaxy[5].fxcenter = 279 
;fakegalaxy[5].fycenter = 231 
fakegalaxy[0].fxcenter = 210 
fakegalaxy[0].fycenter = 157
fakegalaxy[1].fxcenter = 210
fakegalaxy[1].fycenter = 249
fakegalaxy[2].fxcenter = 165
fakegalaxy[2].fycenter = 187 
fakegalaxy[3].fxcenter = 165
fakegalaxy[3].fycenter = 219 
fakegalaxy[4].fxcenter = 255
fakegalaxy[4].fycenter = 187 
fakegalaxy[5].fxcenter = 255 
fakegalaxy[5].fycenter = 219 

seed = 20362

FOR dist = 50, 70, 2 DO BEGIN
    fakegalaxy[0].fxcenter = fakegalaxy[0].fxcenter
    fakegalaxy[0].fycenter = fakegalaxy[0].fycenter + 2 
    fakegalaxy[1].fxcenter = fakegalaxy[1].fxcenter
    fakegalaxy[1].fycenter = fakegalaxy[1].fycenter - 2
    fakegalaxy[2].fxcenter = fakegalaxy[2].fxcenter + 2
    fakegalaxy[2].fycenter = fakegalaxy[2].fycenter + 1
    fakegalaxy[3].fxcenter = fakegalaxy[3].fxcenter + 2
    fakegalaxy[3].fycenter = fakegalaxy[3].fycenter - 1
    fakegalaxy[4].fxcenter = fakegalaxy[4].fxcenter -2
    fakegalaxy[4].fycenter = fakegalaxy[4].fycenter +1
    fakegalaxy[5].fxcenter = fakegalaxy[5].fxcenter -2
    fakegalaxy[5].fycenter = fakegalaxy[5].fycenter -1   

    FOR k = 0, 20, 1 DO BEGIN
;dist = 0
;k = 0       
        spawn, "rm /n/Godiva7/jkrick/galfit/sixexp/diagnostics.log"
        
        FOR a = 0, numogals - 1, 1 DO BEGIN
            randomnumber = randomu(seed, 7)
            fakegalaxy[a]= {fobject, a, fakegalaxy[a].fxcenter,fakegalaxy[a].fycenter,$
              fmagnitude(fix(randomnumber(0)*30)), fregal(fix(randomnumber(1)*37)), $
              faxisrat(fix(randomnumber(3)*19)), fpagal(fix(randomnumber(2)*90))}
        ENDFOR

   
        xmin = 1
        ymin = 1
        xmax = 400
        ymax = 400
        
        convolx = 60
        convoly = 60
        
        bkgd = 0.0
        fdx = fdskyx(fix(randomnumber(4)*20))
        fdy = fdskyy(fix(randomnumber(5)*20))

        openw, outlun, "/n/Godiva7/jkrick/galfit/sixexp/simulate.pro.feedme", /GET_LUN
        
        printf, outlun, "#IMAGE PARAMETERS"
        printf, outlun, "A) ", "none" ;, "   # Input data image (FITS file)"
        printf, outlun, "B) ", "/n/Godiva7/jkrick/galfit/sixexp/fake.fits" ;, "     # Name for the output image"
        printf, outlun, "C) ", "/n/Godiva7/jkrick/galfit/sixexp/noise400.fits" ;, "   # Noise image name(made from data if blank or none) "
        printf, outlun, "D) ", "/n/Godiva7/jkrick/galfit/sixexp/psf.small.fits" ;,"            # Input PSF image for convolution (FITS file)"
        printf, outlun, "E) ", "none" ;, "           # Pixel mask (ASCII file or FITS file with non-0 values)"
        printf, outlun, "F) /n/Godiva7/jkrick/galfit/sixexp/const.test	       # Parameter constraint file (ASCII)"
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
        printf, outlun, " 2)  0.0     1      # dsky/dx (sky gradient in x) "
        printf, outlun, " 3) 0.0    1       # dsky/dy (sky gradient in y) "
        printf, outlun, " Z) 0                  # output image"
        
        
        printf, outlun, "###########################################################################"
        openw, dlun3, "/n/Godiva7/jkrick/galfit/sixexp/diagnostics.log", /GET_LUN, /append
        
        FOR b = 0, numogals-1, 1 DO BEGIN
            
            printf, dlun3, "real galaxy",  fakegalaxy[b], fdx, fdy
            
            
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
        
        spawn, "/n/Godiva7/jkrick/galfit/galfit/galfit  /n/Godiva7/jkrick/galfit/sixexp/simulate.pro.feedme"
        
        spawn, "rm galfit.??"
        
        fits_read, "/n/Godiva7/jkrick/galfit/sixexp/fake.fits", fakedata, fakeheader

        ;instead of just plain noise, add a
        ;noise image with an exponential
        ;background to the fake galaxies

        name = fbkg(fix(randomnumber(6)*numonames))
        fits_read, name, noisedata, noiseheader
        printf, dlun3, name
;        print, "using ", name
        fakedata = fakedata + noisedata
        
        fits_write, "/n/Godiva7/jkrick/galfit/sixexp/fake.fits", fakedata, fakeheader

        free_lun, dlun3
        close, dlun3
        
        galfit, "/n/Godiva7/jkrick/galfit/sixexp/fake.fits", "/n/Godiva7/jkrick/galfit/sixexp/noise400.fits", "/n/Godiva7/jkrick/galfit/sixexp/psf.small.fits", "/n/Godiva7/jkrick/galfit/sixexp/segmentation.fits","/n/Godiva7/jkrick/galfit/sixexp/"
        
        newname = strcompress("/n/Godiva7/jkrick/galfit/sixexp/diagnostics" + string(k)  + "." + $
            string(dist) + ".log", /remove_all)
        newcommand = "mv /n/Godiva7/jkrick/galfit/sixexp/diagnostics.log "+  newname

 ;       newimg = strcompress("/n/Godiva7/jkrick/galfit/sixexp/fake" + string(k) + "." + $
 ;                string(dist)+ ".fits", /remove_all)
;       newcommand2 = "mv /n/Godiva7/jkrick/galfit/sixexp/fake.fits " + newimg
        newimg2 = strcompress("/n/Godiva7/jkrick/galfit/sixexp/nogals" + string(k) + "." + $
            string(dist)+ ".fits", /remove_all)
        newcommand3 = "mv /n/Godiva7/jkrick/galfit/sixexp/nogals.fits " + newimg2
        
        spawn, newcommand
;        spawn, sexcommand
;        spawn, newcommand2
        spawn, newcommand3
    
    ENDFOR
ENDFOR

undefine, fakegalaxy
END
