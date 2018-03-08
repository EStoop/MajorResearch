pro crreject2

filename = strarr(1)  	;declare filename to be a stringarray

close, /all		;close all files = tidiness
crcount = 0

OPENR, lun1, '/n/Sheriff1/jkrick/A141/listfits2', /GET_LUN
;OPENR, lun2, '/n/Sheriff3/jkrick/sep99/A3888/fwhm', /GET_LUN

WHILE (NOT EOF(lun1)) DO BEGIN
    READF, lun1, filename
    print, 'working on', filename
           
    imagefile = filename + '.fits'
    outfile = filename + '.cr.fits'
    datafile = filename + '.cat'
    
    ;run SExtractor on the image to find the objects
    commandline = 'lookforstars2 ' + filename 
    spawn, commandline
    
    ;determine the fwhm cuttof for CR
    ;READF, lun2, fwhmcutoff
    ;print, "fwhm cutoff ", fwhmcutoff
    
    fwhmcutoff = 2.8
    ;read in the fits image with the stars that need to be subtracted
    FITS_READ, imagefile, data, header
        
     ;open the SExtractor output datafile (ccd*.cat)
    OPENR, lun, datafile, /GET_LUN
    
    WHILE (NOT EOF(lun)) DO BEGIN
        READF, lun, o, xcenter, ycenter, a, b, e, f, m, j, isoarea, fwhm, pa, bkgd
        IF (xcenter GE 1010 AND xcenter LE 1020 OR xcenter GT 2040) THEN BEGIN
            ;print, "bad column"
        ENDIF ELSE BEGIN

            IF(isoarea LT 6.0 OR fwhm LT fwhmcutoff ) THEN BEGIN
            ;add 2 pixels to both a and b, and then block out to that radius
                a = a + 2
                b = b + 2
                return = del_star(pa, f, a , b, xcenter, ycenter, data)
                
                ;print, "removing cosmic ray", o, xcenter, ycenter
                crcount = crcount +1
            ENDIF
        ENDELSE
        
    ENDWHILE
    print, "outfile", outfile
    fits_write, outfile, data, header
    print, "number of CR removed", crcount
    crcount = 0
ENDWHILE
close, lun
free_lun, lun
close, lun1
free_lun, lun1


END
