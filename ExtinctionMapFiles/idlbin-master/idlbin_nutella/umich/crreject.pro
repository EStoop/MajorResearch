pro crreject

filename = "hello"  	;declare filename to be a stringarray

close, /all		;close all files = tidiness
crcount = 0

;OPENR, lun1, '/n/Godiva1/jkrick/A3888/listfits', /GET_LUN
;OPENR, lun1, '/n/Godiva1/jkrick/A3888/final/listcr', /GET_LUN
;OPENR, lun2, '/n/Godiva3/jkrick/sep99/A3888/fwhm', /GET_LUN
OPENR, lun1, '/n/Godiva6/jkrick/A3880/listcr', /GET_LUN
OPENR, lun2, '/n/Godiva6/jkrick/A3880/fwhmcutoff', /GET_LUN
;OPENR, lun1, '/n/Godiva2/jkrick/A141/listcr', /GET_LUN
;OPENR, lun2, '/n/Godiva2/jkrick/A141/fwhmcutoff', /GET_LUN

WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, filename

    print, 'working on', filename
           
    imagefile = filename + ".fits"
    outfile = filename + '.cr.fits'
    datafile = filename + '.cat'
;    print, 'working on', imagefile
    
    ;run SExtractor on the image to find the objects
    ;commandline = 'lookforstars2 ' + filename 
    commandline = 'sex ' + imagefile + ' -c  /n/Godiva4/jkrick/A3984/default2.sex'
    spawn, commandline
    
    ;determine the fwhm cuttof for CR!!!!!!!!!!
    READF, lun2, junk, fwhmcutoff
;    fwhmcutoff = fwhmcutoff - .1
    print, "fwhm cutoff ", fwhmcutoff
    
;    fwhmcutoff = 3.8;2.8
    ;read in the fits image with the stars that need to be subtracted
    FITS_READ, imagefile, data, header
        
     ;open the SExtractor output datafile (ccd*.cat)
    OPENR, lun, '/n/Godiva4/jkrick/A3984/SExtractor2.cat', /GET_LUN
    
    WHILE (NOT EOF(lun)) DO BEGIN
        READF, lun, o, xcenter, ycenter, a, b, e, f, m, j, isoarea, fwhm, pa, bkgd
        IF (xcenter GE 800. AND xcenter LE 808. OR xcenter GT 2040) THEN BEGIN
            ;print, "bad column"
        ENDIF ELSE BEGIN

            IF(isoarea LT 6.0 OR fwhm LT fwhmcutoff ) THEN BEGIN
            ;add 2 pixels to both a and b, and then block out to that radius
                a = a + 1
                b = b + 1
                return = del_star(pa, f, a , b, xcenter, ycenter, data)
                
 ;               print, "removing cosmic ray", o, xcenter, ycenter
                crcount = crcount +1
            ENDIF
        ENDELSE
        
    ENDWHILE
    print, "outfile", outfile
    fits_write, outfile, data, header
    print, "number of CR removed", crcount
    crcount = 0

    close, lun
    free_lun, lun

ENDWHILE
;close, lun1
;free_lun, lun1
;close, lun2
;free_lun, lun2

END
