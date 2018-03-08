FUNCTION printsex, catname, fwhmaver
openw, outlun, "/n/Godiva3/jkrick/sep9940/new.sex", /get_lun


printf, outlun, "# Default configuration file for SExtractor V1.2"
printf, outlun, "# EB 18/08/97"
printf, outlun, "# (*) indicates parameters which can be omitted from this config file."
printf, outlun, ""
printf, outlun, "#-------------------------------- Catalog ------------------------------------"
printf, outlun, ""
;printf, outlun, "CATALOG_NAME  /n/Godiva4/jkrick/aug98/aug19/ccd1080.cat"
printf, outlun, "CATALOG_NAME  " + catname
printf, outlun, "CATALOG_TYPE	ASCII "	
printf, outlun, ""
printf, outlun, "PARAMETERS_NAME	/n/Godiva3/jkrick/sep9940/phot.param	#file catalog contents"
printf, outlun, ""
printf, outlun, "#------------------------------- Extraction ----------------------------------"
printf, outlun, ""
printf, outlun, "DETECT_TYPE	CCD	"
printf, outlun, "#DETECT_IMAGE	SAME	"
printf, outlun, "#FLAG_IMAGE	flag.fits"
printf, outlun, "DETECT_MINAREA	6	"
printf, outlun, "DETECT_THRESH	1.5	"
printf, outlun, "ANALYSIS_THRESH	1"
printf, outlun, ""
printf, outlun, "FILTER		Y		"
printf, outlun, "FILTER_NAME	/n/Godiva7/jkrick/Sex/sextractor2.1.0/config/gauss_1.5_3x3.conv "
printf, outlun, ""
printf, outlun, "DEBLEND_NTHRESH	40		# Number of deblending sub-thresholds"
printf, outlun, "DEBLEND_MINCONT	0.0001	# Minimum contrast parameter for deblending"
printf, outlun, ""
printf, outlun, "CLEAN		Y		# Clean spurious detections? (Y or N)?"
printf, outlun, "CLEAN_PARAM	1.0		# Cleaning efficiency"
printf, outlun, ""
printf, outlun, "#BLANK		Y		# Blank detected objects (Y or N)?"
printf, outlun, ""
printf, outlun, "#------------------------------ Photometry -----------------------------------"
printf, outlun, ""
;printf, outlun, "PHOT_APERTURES   5              #32"
printf, outlun, "PHOT_APERTURES   " + string(5.*fwhmaver)
printf, outlun, "PHOT_AUTOPARAMS	2.5, 3.0	# MAG_AUTO parameters: <Kron_fact>,<min_radius>"
printf, outlun, ""
printf, outlun, "SATUR_LEVEL	28000.0		# level (in ADUs) at which arises saturation"
printf, outlun, ""
printf, outlun, "MAG_ZEROPOINT	25.0 	# magnitude zero-point for r"
printf, outlun, "MAG_GAMMA	4.0		# gamma of emulsion (for photographic scans)"
printf, outlun, "GAIN		3.0		# detector gain in e-/ADU."
printf, outlun, "PIXEL_SCALE	0.435		# size of pixel in arcsec (0=use FITS WCS info)."
printf, outlun, ""
printf, outlun, "#------------------------- Star/Galaxy Separation ----------------------------"
printf, outlun, ""
;printf, outlun, "SEEING_FWHM		1.0	# stellar FWHM in arcsec"
printf, outlun, "SEEING_FWHM		" + string(fwhmaver *0.259) 
printf, outlun, "STARNNW_NAME	default.nnw	# Neural-Network_Weight table filename"
printf, outlun, ""
printf, outlun, "#------------------------------ Background -----------------------------------"
printf, outlun, ""
printf, outlun, "BACK_SIZE	60		# Background mesh: <size> or <width>,<height>"
printf, outlun, "BACK_FILTERSIZE	4		# Background filter: <size> or <width>,<height>"
printf, outlun, ""
printf, outlun, "BACKPHOTO_TYPE	GLOBAL		# can be "
printf, outlun, "BACKPHOTO_THICK	24		# thickness of the background LOCAL annulus (*)"
printf, outlun, ""
printf, outlun, "#------------------------------ Check Image ----------------------------------"
printf, outlun, ""
printf, outlun, "CHECKIMAGE_TYPE none "
printf, outlun, ""			  
printf, outlun, ""		  
printf, outlun, ""		   
printf, outlun, "CHECKIMAGE_NAME /n/Godiva4/jkrick/aug98/background.fits"
printf, outlun, ""
printf, outlun, "#--------------------- Memory (change with caution!) -------------------------"
printf, outlun, ""
printf, outlun, "MEMORY_OBJSTACK	2000		# number of objects in stack"
printf, outlun, "MEMORY_PIXSTACK	100000		# number of pixels in stack"
printf, outlun, "MEMORY_BUFSIZE	512		# number of lines in buffer"
printf, outlun, ""
printf, outlun, "#---------------- Scanning parameters (change with caution!) -----------------"
printf, outlun, ""
printf, outlun, "#SCAN_ISOAPRATIO	0.6		# maximum isoph. to apert ratio allowed (*)"
printf, outlun, ""
printf, outlun, "#----------------------------- Miscellaneous ---------------------------------"
printf, outlun, ""
printf, outlun, "VERBOSE_TYPE	QUIET		# can be QUIET, NORMAL or FULL (*)"
printf, outlun, ""


close, outlun
free_lun, outlun


return, 0
END



PRO photometry

close, /all

filename = "hello"
openr, luns, "/n/Godiva6/jkrick/oct98/liststandard", /get_lun
WHILE (NOT EOF(luns)) DO BEGIN
    READF, luns, filename
    print, "working on ", filename

    fitsname = filename + ".fits"
    catname = filename + ".cat"
    starname = filename + ".star"
    
;run SExtractor on the standard frame
    commandline = '/n/Godiva7/jkrick/Sex/sex ' + fitsname +$
      " -c /n/Godiva3/jkrick/sep9940/phot.sex"
    spawn, commandline

;figure out what the fwhm of stars is on that frame
    fwhmaver = 0.0
    n = 0
    openr, lun, "/n/Godiva3/jkrick/sep9940/SExtractor.cat", /get_lun
    WHILE (NOT EOF(lun)) DO BEGIN
        READF, lun, o, xcenter, ycenter, a, b, e, fbest, mbest,$
          fiso, miso, faper, maper, mapererr, isoarea, fwhm
        IF (fbest GT 0.0) AND (isoarea GT 12) AND (maper LT 19.)$
          AND (e LT 0.15) AND (xcenter LT 2040) AND $
          (fwhm GT 0.) AND (ycenter LT 3140) THEN BEGIN
        ;a decently bright star
            fwhmaver = fwhmaver + fwhm
            n = n + 1
        ENDIF
    ENDWHILE
    print, "fwhmaver", fwhmaver
    
    fwhmaver = fwhmaver / n
    print, "fwhmaver", fwhmaver
    close, lun
    free_lun, lun
    
;re-run SExtractor with aperture sized 5 x fwhm

    holder = printsex(catname, fwhmaver)

    commandline = '/n/Godiva7/jkrick/Sex/sex ' + fitsname +$
      " -c /n/Godiva3/jkrick/sep9940/new.sex"
    spawn, commandline
    


;open catname, write new file with just the candidate stars.
    openr, lun1, catname, /get_lun
    openw, lun2, starname, /get_lun
    WHILE (NOT EOF(lun)) DO BEGIN
        READF, lun1, o, xcenter, ycenter, a, b, e, fbest, mbest,$
          fiso, miso, faper, maper, mapererr, isoarea, fwhm
        IF (fbest GT 0.0) AND (isoarea GT 10) AND (maper LT 19.0)$
          AND (e LT 0.25) AND (xcenter LT 2040) AND (ycenter LT 3140) THEN BEGIN
                                ;a standard star
            printf, lun2, xcenter, ycenter, maper, mapererr
        ENDIF
    ENDWHILE
    
    close, lun1
    free_lun, lun1
    close, lun2
    free_lun, lun2
    
    
ENDWHILE
close, luns
free_lun, luns        

END
