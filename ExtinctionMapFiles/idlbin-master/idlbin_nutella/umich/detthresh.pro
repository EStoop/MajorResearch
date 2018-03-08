PRO detthresh

close,/all
colors = GetColor(/Load, Start=1)
ps_open, file = "/n/Godiva1/jkrick/A3888/final/detthresh.ps", /portrait, /color, ysize = 7

numoobjects=10
cluster = replicate({object, name:"joe", z:0D, det:0D, sig:0D,m0:0D, pix:0D },numoobjects)

cluster[0] = {object, "A4059", 0.048, 24.64, 0.0,22.04, 0.435}
cluster[1] = {object, "A3880", 0.058, 24.79, 0.0,22.04, 0.435}
cluster[2] = {object, "A2734", 0.062, 24.47, 0.0,22.04, 0.435}
cluster[3] = {object, "A2556", 0.087, 24.25, 0.0,22.04, 0.435}
cluster[4] = {object, "A4010", 0.096, 24.51, 0.0,22.04, 0.435}
cluster[5] = {object, "A3888", 0.150, 25.64, 0.0, 24.6, 0.259}
cluster[6] = {object, "A3984", 0.180, 25.48, 0.0, 24.6, 0.259}
cluster[7] = {object, "A0141", 0.230, 25.16, 0.0, 24.6, 0.259}
cluster[8] = {object, "AC114", 0.310, 24.73, 0.0, 24.6, 0.259}
cluster[9] = {object, "AC118", 0.310, 24.83, 0.0, 24.6, 0.259}

FOR n=0,9, 1 DO BEGIN
    delta = 26.4 - (cluster[5].det -cluster[n].det)
    cluster[n].sig = (1/0.008521)*((0.259)^2)*10^(((delta)-24.6)/(-2.5))
    print, cluster[n].name, delta, cluster[n].sig

openw, outlun, "/n/Godiva1/jkrick/A3888/final/test.sex",/get_lun

printf, outlun, "CATALOG_NAME  /n/Godiva1/jkrick/A3888/final/test.cat"
printf, outlun, "CATALOG_TYPE	ASCII		"

printf, outlun, "PARAMETERS_NAME	/n/Godiva1/jkrick/A3888/daofind.param	#file catalog contents"



printf, outlun, "DETECT_TYPE	CCD		"
printf, outlun, "DETECT_MINAREA	6		# 6 minimum number of pixels above threshold"
printf, outlun, "DETECT_THRESH   ", cluster[n].sig
printf, outlun, "ANALYSIS_THRESH	1 #25.5,24.6	#31.1, 24.6	"

printf, outlun, "FILTER		Y		# apply filter for detection "
printf, outlun, "FILTER_NAME	/n/Godiva7/jkrick/Sex/sextractor2.1.0/config/gauss_1.5_3x3.conv "

printf, outlun, "DEBLEND_NTHRESH	32		# Number of deblending sub-thresholds"
printf, outlun, "DEBLEND_MINCONT	0.000001	# Minimum contrast parameter for deblending"

printf, outlun, "CLEAN		Y		# Clean spurious detections? (Y or N)?"
printf, outlun, "CLEAN_PARAM	1.0		# Cleaning efficiency"


printf, outlun, "PHOT_APERTURES    19 #32"
printf, outlun, "PHOT_AUTOPARAMS	2.5, 3.0	# MAG_AUTO parameters: <Kron_fact>,<min_radius>"

printf, outlun, "SATUR_LEVEL	31.0		# level (in ADUs) at which arises saturation"

printf, outlun, "MAG_ZEROPOINT	24.6  	# magnitude zero-point for V"
printf, outlun, "MAG_GAMMA	4.0		# gamma of emulsion (for photographic scans)"
printf, outlun, "GAIN		3.0		# detector gain in e-/ADU."
printf, outlun, "PIXEL_SCALE	0.259		# size of pixel in arcsec (0=use FITS WCS info)."

printf, outlun, "SEEING_FWHM		1.0	# stellar FWHM in arcsec"
printf, outlun, "STARNNW_NAME	default.nnw	# Neural-Network_Weight table filename"

printf, outlun, "BACK_SIZE	80 #100 for V	# Background mesh: <size> or <width>,<height>"
printf, outlun, "BACK_FILTERSIZE	4  #6 for V	# Background filter: <size> or <width>,<height>"

printf, outlun, "BACKPHOTO_TYPE	GLOBAL		# can be "
printf, outlun, "BACKPHOTO_THICK	24		# thickness of the background LOCAL annulus (*)"

printf, outlun, "CHECKIMAGE_TYPE apertures 	# can be one of "
printf, outlun, "CHECKIMAGE_NAME /n/Godiva1/jkrick/A3888/final/segmentation.fits	"

printf, outlun, "MEMORY_OBJSTACK	2000		# number of objects in stack"
printf, outlun, "MEMORY_PIXSTACK	100000		# number of pixels in stack"
printf, outlun, "MEMORY_BUFSIZE	512		# number of lines in buffer"

printf, outlun, "VERBOSE_TYPE	quiet		# can be "



free_lun, outlun
close, outlun

commandline = '/n/Godiva7/jkrick/Sex/sex ' + "/n/Godiva1/jkrick/A3888/final/fullr2.fits" + " -c /n/Godiva1/jkrick/A3888/final/test.sex"
spawn, commandline


rarr = fltarr(10000)
rprime = fltarr(10000)
i = 0

OPENR, lun2, "/n/Godiva1/jkrick/A3888/final/test.cat", /GET_LUN

WHILE (NOT EOF(lun2)) DO BEGIN
    READF, lun2, o, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux,isocorflux,isoflux
    
    IF (fwhm GT 4.0) AND (isoarea GT 6.0) THEN BEGIN
        OPENR, lun6, "/n/Godiva1/jkrick/A3888/final/fullr2.cat", /GET_LUN
        WHILE (NOT EOF(lun6)) DO BEGIN
            READF, lun6, o, xcen6, ycen6, a6, b6, e6, f6, m6, isoarea6, fwhm6, pa6, bkgd6, junk,junk,junk

            IF (fwhm GT 4.0) AND (isoarea GT 6.0)AND (xcen6 GT xcenter - 1.0)$
              AND (xcen6 LT xcenter + 1.0) AND $
              (ycen6 GT ycenter - 1.0) AND (ycen6 LT ycenter + 1.0) THEN BEGIN
                rarr[i] = sqrt(isoarea6/(!PI*(1-e6)))
                rprime[i] = sqrt(isoarea/(!PI*(1-e)))
                i = i +1
            ENDIF
        ENDWHILE
        close, lun6
        free_lun, lun6
    ENDIF

ENDWHILE
print, i
rarr = rarr[0:i-1]
rprime = rprime[0:i-1]



close, lun2
free_lun, lun2

ytitl = strcompress("sma"+string(cluster[n].det)+"detthresh")

plot, rarr, rprime, psym = 2, xtitle ="sma A3888 25.64 dettthresh", ytitle = ytitl,$
  thick = 3, xthick=3, ythick=3,charthick=3,xrange=[0,80],yrange=[0,80]
err = fltarr(i) + 1
start = [0.66,0.]
result = MPFITFUN('linear',rarr,rprime,err, start)

oplot, findgen(100), result(0)*findgen(100) + result(1), thick = 3, linestyle = 2, color = colors.red

xyouts, 10, 70, cluster[n].name, thick = 3
xyouts, 10,50, string(result(0)) + string(result(1)), thick=3

endfor

ps_close, /noprint, /noid

undefine, cluster
END
