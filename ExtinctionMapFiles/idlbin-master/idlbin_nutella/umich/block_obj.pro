;May 2002
;
;The purpose of this program is to block all the objects in an image.
; The stars are fit with a moffat profile, and subtracted as well as
; blocked in the centers.
; bad columns, cosmic rays, faint objects, and galaxies are also dealt
; with.
; saturated stars are ignored and need to be blocked by hand
;
;input: image with objects to be blocked, and a blank image of the
;       same size
;
;output: an image named ccd****.s.fits with the objects set equal to 
; 	-100, and an image called substars with the homemade stars on a 
;	a zero background ready to be subtracted
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


FUNCTION starmask, f, xcenter, ycenter, maxradius, xmax, ymax,blankdata

;this function acutally makes the homemade star and adds it to a blank
;image

;need to change this for final version. is a way of normalizing the flux
;D = f * 1.36/38.6597
;D = f * 1.36/53.4
;D = f * 1.36/39.6
;D = f* 1.36/35.6
;D = f 
D = 0.044*f 
sb4 = 24.6 - 2.5*alog10((D*mymoffat(15.4))/(0.259^2.))
sb2 = 24.6 - 2.5*alog10((D*mymoffat(7))/(0.259^2.))
;print, "xcen, ycenf, D, sb4",xcenter, ycenter, f, D, sb4, sb2
radius = 0.0		

intx = fix(xcenter)
inty = fix(ycenter)

dx = xcenter-intx
dy = ycenter-inty

;print, "star inside starmask", intx, dx, inty, dy

blank1 = blankdata
blank1(*,*) = 0.0

FOR y = inty-maxradius, inty + maxradius,1 DO BEGIN
    FOR x = intx -maxradius, intx + maxradius,1 DO BEGIN
        ;find what the radius is from the center
        radius = sqrt((intx - x)^2 + (inty - y)^2)
        
;        intensity = D*mymoffat(radius)
        IF radius NE 0 THEN BEGIN
;            intensitysb = 3.9 + 2.5*alog10((38.6/radius)^(-3.)) + sb4
            intensitysb = 6.9 + 2.5*alog10((38.6/radius)^(-3.)) + sb2
            intensity = (0.259^2.)*(10.^((24.6-intensitysb)/(2.5)))
        ENDIF

        ;make sure not off the edge of the chip
        ;need to deal with going from fractional pixels to whole pixels
        IF (x LT xmax AND y LT ymax AND x GT 1 AND y GT 1 ) THEN BEGIN
            blank1[x,y] = blank1[x,y] + intensity
        ENDIF
    ENDFOR
ENDFOR

blank1 = shiftf(blank1, dx,dy)
blankdata = blankdata - blank1
return, 0              
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO block_obj

close, /all

openw, starlun, "/n/Godiva1/jkrick/A3888/final/starlist", /get_lun
openw, gallun, "/n/Godiva1/jkrick/A3888/final/gallist", /get_lun

;declare variables
o = 0.0                         ;object number
xcenter = 0.0   		;xcenter
ycenter = 0.0			;ycenter
e = 0.0                         ;ellipticity
a = 0.0                         ;semi-major axis
b = 0.0                         ;semi-minor axis
m = 0.0                         ;magnitude (un-corrected for zpt, X, ...)
j = 0.0                         ;junk
f = 0.0                         ;flux
fwhm = 0.0                      ;full width at half maximum
isoarea = 0.0                   ;area inside of last isophote

filename = strarr(1)            ;declare filename to be a stringarray
maxradius = 58.0                ; maximum radius in pixels to create a star
intensity = 0.0
bkgd = 0.0

true = 0                        ;manual boolean variable

;size of the image
xmax = 3320;2801;2001;1001
ymax = 3260;2701;2001;1001

;read in a blank image to which I can add my own stars			
;FITS_READ, '/n/Godiva1/jkrick/A3888/largeblank.fits', blankdata, blankheader
;filename = '/n/Godiva1/jkrick/A3888/final/largeV2'

FITS_READ, '/n/Godiva1/jkrick/A3888/final/fullblank.fits', blankdata, blankheader
filename = '/n/Godiva1/jkrick/A3888/final/fullr2'

;FITS_READ, '/n/Godiva1/jkrick/A3888/blank4000.fits', blankdata, blankheader
;filename = '/n/Godiva1/jkrick/A3888/A3888V.div'
;filename = '/n/Godiva1/jkrick/A3888/galfit/centerV.div.2'

imagefile = filename + '.fits'
;datafile = filename + '.cat'
datafile = '/n/Godiva1/jkrick/A3888/final/fullr2.cat'
newimagefile = filename + '.b.fits'	

;run SExtractor on the image to find the objects
;commandline = 'lookforstars2 ' + filename 
;commandline = '/n/Godiva7/jkrick/Sex/sex ' + "/n/Godiva1/jkrick/A3888/galfit/centerV.div.2.fits" + " -c /n/Godiva1/jkrick/A3888/default.sex"

;commandline = '/n/Godiva7/jkrick/Sex/sex ' + "/n/Godiva5/jkrick/aug98/19aug98/ccd1100.div.fits" + " -c /n/Godiva1/jkrick/A3888/final/iclr.sex"

commandline = '/n/Godiva7/jkrick/Sex/sex ' + "/n/Godiva1/jkrick/A3888/final/fullr2.fits" + " -c /n/Godiva1/jkrick/A3888/default.sex"
spawn, commandline

;for the galaxies that did not make the color cut in cmd.pro
;datafile = '/n/Godiva1/jkrick/A3888/final/block_obj.r.input'

;read in the fits image with the stars that need to be subtracted
FITS_READ, imagefile, data, header


;read the data file into the variables until hit EOF
starcount = 0
galcount = 0
crcount = 0
dimcount = 0
satcount = 0
totalflux = 0.0  				;in case I want to know
galbestflux = 0.0
galisoflux = 0.0
subgalscount = 0
;-________________________________________________________________________________
;saturated stars from USNO
satx = fltarr(13)
saty = fltarr(13)
satrmag = fltarr(13)
satVmag = fltarr(13)
satx=[1078.44 ,2094.74, 2438.01 ,2704.27 ,2799.16,187.75 ,782.55 ,1923.36 ,1768.55 , 2114.58 , 3016.94 ,3373.43,3353.47]
saty=[2945.76,2456.28,2587.25,2473.15, 2737.52,2242.52,1040.93,1278.38,443.71,286.04, 988.10,1661.01,721.67]
satrmag=[11.2,12.6,13.7,13.1,13.1,14.4,11.3,14.0,14.8,14.0,12.9,11.6,11.6]
satVmag=[11.6,13.4,14.4,13.7,13.7,14.9,11.6,14.6,14.0,14.8,13.6,11.9,11.8]

FOR satstar = 0, 12, 1 DO BEGIN
;    holder = starmask (10^((satrmag[satstar] - 24.6)/(-2.5)), satx[satstar], saty[satstar],500.,xmax,ymax,blankdata)
;    holder = starmask (10^((satrmag[satstar] - 24.6)/(-2.5)), satx[satstar], saty[satstar],500.,xmax,ymax,data)
;    holder = starmask (10^((satVmag[satstar] - 24.6)/(-2.5)), satx[satstar], saty[satstar],500.,xmax,ymax,blankdata)
ENDFOR



;-------------------------------------------------------------------------------
;open the SExtractor output datafile (ccd*.cat)
OPENR, lun, datafile, /GET_LUN

WHILE (NOT EOF(lun)) DO BEGIN
    ;!!!! change this if you change the daofind.param , now usign godiva7/galfit/dao..
    READF, lun, o, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux,isocorflux,isoflux
    ;print, "working on object", o
    
   
    ;___________________________________________________________________________
    ; ignore the other bad spots with flux less than 0
    IF (f LT 0 ) THEN true = 1
   
    ;__________________________________________________________________________
    ;get rid of dim objects
    ;these objects aren't putting enough into their wings to worry about them
    IF (apflux LE 0.14 ) AND (isoarea GE 6.0) AND (true EQ 0) THEN BEGIN
        ;go out to two times the fwhm then zap 'em
        a = 1.5*fwhm
        b = 1.5*fwhm

        IF (isoarea LT 14) THEN isoarea = 14 ;small galaxies

        sma = sqrt(isoarea/(!PI*(1-e)))
        smb = isoarea/(!PI*sma)
       ;print, "reomving dim object", o,xcenter, ycenter
        return = del_star(pa, f, sma , smb, xcenter, ycenter, data)
        return = del_star(pa, f, sma , smb, xcenter, ycenter, blankdata)
        
        true = 1
        IF (xcenter GT 950) AND (xcenter LT 1800) AND (ycenter GT 1100) AND (ycenter LT 1950) THEN dimcount = dimcount +1

    ENDIF

    ;_________________________________________________________________________
    ;get rid of cosmic rays
    IF(isoarea LT 6.0 OR fwhm LT 3.8 AND true EQ 0) THEN BEGIN  ;2.8
        ;add 2 pixels to both a and b, and then block out to that radius
        ;print, "CR", o,xcenter,ycenter
        a = a + 2
        b = b + 2
       return = del_star(pa, f, a , b, xcenter, ycenter, data)

		
	;print, "removing cosmic ray", o, xcenter, ycenter
        true = 1
        IF (xcenter GT 950) AND (xcenter LT 1800) AND (ycenter GT 1100) AND (ycenter LT 1950) THEN crcount = crcount +1
    ENDIF
		
    ;__________________________________________________________________________
    ;get rid of the stars that are not very saturated
        
;  IF (fwhm GE 2.8 AND fwhm LE 4.0 AND true EQ 0) THEN BEGIN    
;    V fwhm limits
    IF (fwhm GE 3.8 AND fwhm LE 5.4  AND true EQ 0) THEN BEGIN    ;5.1
        openw, dimlun, "/n/Godiva1/jkrick/A3888/tvmark.out", /GET_LUN, /append
        printf, dimlun, xcenter, ycenter, fwhm
        close, dimlun
        free_lun, dimlun
        IF (m LE 12.6) THEN BEGIN
            print,"shouldn't have gotten here", xcenter, ycenter
            true = 1
        ENDIF 

        ;for the dim things that could be
        ;stars, or galaxies, or junk, I won't subtract a profile, instead just
        ;block out to 2"
        IF (m GT 22.1) THEN BEGIN
            return = del_star(pa, f, 8, 8, xcenter, ycenter, data)

            true = 1

            ;holder = starmask (f, xcenter, ycenter, maxradius,xmax,ymax,blankdata)
            IF (xcenter GT 436) AND (xcenter LT 2436) AND (ycenter GT 524) AND (ycenter LT 2524) THEN BEGIN
                starcount = starcount + 1
;                printf, starlun, "star in center with mag ", m
            ENDIF

        ENDIF
        IF (true EQ 0 AND m GT 17) THEN BEGIN
;            return = del_star(pa, f, 22, 22, xcenter, ycenter, data)
            return = del_star(pa, f, 30, 30, xcenter, ycenter, data)
 ;           return = del_star(pa, f, 30, 30, xcenter, ycenter, blankdata)
            true = 1
;           print, "removing star at   ",xcenter, ycenter,m
;;            holder = starmask (f, xcenter, ycenter, maxradius,xmax,ymax,blankdata)
            holder = starmask (f, xcenter, ycenter, maxradius,xmax,ymax,data)

            IF (xcenter GT 436) AND (xcenter LT 2436) AND (ycenter GT 524) AND (ycenter LT 2524) THEN BEGIN
                starcount = starcount + 1
                printf, starlun, "star in center with mag ", m
            ENDIF

        ENDIF
        IF (true EQ 0) THEN BEGIN
	; want to make the flux within a radius of 5" the same as
	; the flux within the inner 5" taken from mymoffat profile
	; that function is sitting in mymoffat.pro 
	; D is the scale factor for making mymoffat have the same flux as the star
	; 38.65 is gotten by making a normalized star with the same specs
	; 54 should be the number according to the math
	;need to change this if Sextractor params change or params in next loop change

;;            holder = starmask (f, xcenter, ycenter, maxradius,xmax,ymax,blankdata)
            holder = starmask (f, xcenter, ycenter, maxradius,xmax,ymax,data)


		
	;and now for the part where I block out the inner 'a' pixels of all of the stars
	;just to make sure that we get it all, and so I don't have to worry about
	;seeing, focus, ...

            b= 30;19               ;block out star up to a radius of B1 pixels
            a= 30;19
            return = del_star(pa, f, a , b, xcenter, ycenter, data)
;;            return = del_star(pa, f, a , b, xcenter, ycenter, blankdata)

            IF (xcenter GT 436) AND (xcenter LT 2436) AND (ycenter GT 524) AND (ycenter LT 2524) THEN BEGIN
                starcount = starcount + 1  
                printf, starlun, "star in center with mag ", m
            ENDIF

            totalflux = totalflux + f
		
;           print, "removing star at   ",xcenter, ycenter
            true = 1

        ENDIF
        IF (true EQ 0 AND m LE 17) THEN BEGIN
            print, "m LT 17 ", xcenter, ycenter
            true = 1
        ENDIF

    ENDIF

    ;________________________________________________________________________________
    ;get rid of the galaxies
;    IF (fwhm GT 4.0 AND true EQ 0) THEN BEGIN
        IF (fwhm GT 5.4 AND true EQ 0) THEN BEGIN;5.1
  
        ;to make the dim galaxy masks big enough
        IF (isoarea LT 14) THEN isoarea = 14 ;small galaxies

        sma = sqrt(isoarea/(!PI*(1-e)))
        smb = isoarea/(!PI*sma)
		
        ;go out further than the isoarea
        IF (sma GT 15.0 OR smb GT 15.0) THEN BEGIN ;big galaxies
            IF(isoarea GT 2000.) THEN BEGIN
 ;               print, "central galaxy bigger than 2000", xcenter, ycenter, isoarea
                sma = 2.0*sma ; 2.6,2.0,1.4                  
                smb = 2.0*smb;  2.6,2.0,1.4      
            ENDIF ELSE BEGIN

                sma = 2.3*sma; 3.0,2.3 ,1.6      
                smb = 2.3*smb; 3.0,2.3,1.6       
            ENDELSE
        
        ENDIF ELSE BEGIN        ;normal galaxies
            sma = 2.0*sma; 2.6,2.0,1.4           
            smb = 2.0*smb; 2.6,2.0,1.4                       
        ENDELSE

         theta = pa
        return = del_star(theta, f, sma , smb, xcenter, ycenter, data)
;        print, "sending data to del_star, pa, sma, smb", theta, sma, smb

        return = del_star(pa, f, sma , smb, xcenter, ycenter, blankdata)
;        print, "sending blankdata to del_star, pa, sma, smb", pa, sma, smb

        openw, blocklun, "/n/Godiva1/jkrick/A3888/final/masksize.log", /GET_LUN, /append
        printf, blocklun, xcenter, ycenter, sma, smb
        close, blocklun
        free_lun, blocklun
       
        true = 1
        IF (xcenter GT 836) AND (xcenter LT 2036) AND (ycenter GT 924) AND (ycenter LT 2124) THEN BEGIN
            galcount = galcount+1
            printf, gallun, xcenter, ycenter, sma, smb
        ENDIF

        dist = sqrt((1436. - xcenter)^2 + (1524.-ycenter)^2)
        IF (dist LT 616 ) THEN begin
            galbestflux = galbestflux + isoflux
            galisoflux = galisoflux + isocorflux
        ENDIF
    ENDIF
    IF (true EQ 0) THEN print, xcenter,ycenter,"did not get blocked"
 
    true = 0

ENDWHILE

print,"number of stars is:", starcount
print,"number of galaxies is:", galcount
print,"number of cosmic rays is:", crcount
print,"number of dim objects is:", dimcount
;print,"number of saturated objects is:", satcount
print,"total flux_iso in galaxies is:",galbestflux
print,"total flux_isocor in galaxies is:",galisoflux

;print,"number of subtracted galaxies is:", subgalscount
;close the SExtractor file, and free up the lun so it can be used again
close, lun
free_lun, lun

;want to put the sextractor paramters in the header to avoid loosing
;them
j = "j"
OPENR, lunicl, "/n/Godiva1/jkrick/A3888/final/iclr.sex", /GET_LUN
WHILE (NOT EOF(lunicl)) DO BEGIN
    readf, lunicl, j
    sxaddpar, header, "COMMENT",j
ENDWHILE
close, lunicl
free_lun, lunicl

;write a new fits file with the masked pixels, and one with the stars
;fits_write, newimagefile, shift(data, -1,-1), header

;;fits_write, '/n/Godiva1/jkrick/A3888/final/samemask/mask.bblock.galr.fits',shift(blankdata,-1,-1), blankheader

;lessstars = data ;- blankdata   ;;;minus
;fits_write, '/n/Godiva1/jkrick/A3888/galfit/nogals.s-stars.fits', shift(lessstars, -1, -1), header

;;fits_write, '/n/Godiva4/jkrick/aug98/aug19/ccd1100.bblock.fits', shift(data, -1, -1), header

;make a profile of what is left
;random center choice
;result = funcprof("/n/Godiva1/jkrick/A3888/galfit/nogals.s-stars", 492., 503.)

;undefine, subgals

close, /all
END
;take out later
;FOR badx=1010, 1020,1 DO BEGIN
;    FOR bady = 1,2047,1 DO BEGIN
;data[badx,bady] = -1E6
;ENDFOR

;-------------------------------------------------------------------------------------------
;read in the list  of galaxies that have been subtracted and hence
;only need to be blocked out to a small radius
;numoobjects = 93
;subgals= replicate({object, num:0, centerx:0D, centery:0D},numoobjects)
;i = 0
;openr, gallun, "/n/Godiva1/jkrick/A3888/galfit/gallist", /get_lun
;WHILE ( NOT EOF(gallun)) DO BEGIN
;    readf, gallun, number, xcen, ycen
;    subgals[i] = {object, number, xcen, ycen}
;    i = i + 1
;ENDWHILE
;close, gallun
;free_lun, gallun
;------------------------------------------------------------------------------
      ;except those galalxies which have already been subtracted
        ;FOR j = 0, numoobjects - 1, 1 DO BEGIN
                                ; print, xcenter, ycenter, subgals[j].centerx, subgals[j].centery
        ;    IF (xcenter LT subgals[j].centerx + 1. AND xcenter GT subgals[j].centerx - 1. AND $
        ;         ycenter LT subgals[j].centery + 1. AND ycenter GT subgals[j].centery -1.) THEN BEGIN
        ;        sma = 30.
        ;        smb = 30.
        ;        subgalscount = subgalscount + 1
                
        ;    ENDIF
            
        ;ENDFOR


;fits_write, "/n/Godiva1/jkrick/A3888/galfit/test.fits", data, header
;get the bad columns out here
;data[3030:3038,2100:3450] = -1E6
;data[3440:3448,2520:2872] = -1E6
;data[3512:3520,480:2528] = -1E6
;data[3456:3464,48:776] = -1E6
;data[3000:3008,48:900] = -1E6
;data[2600:2608,1:1000] = -1E6
;data[3102:3352,2864:2872] = -1E6
;data[3440:3502,2520:2528] = -1E6


;for ccd3126 bleedout columns
;data[1052:1073,1022:1522] = -100
;data[1052:1073,496:996] = -100
	
;display the new image with the masked pixels
;twoimage = congrid(data,507.5,510)  ;makes smaller for display on screen
;tv, twoimage	;displays


;make a mask file out of the ***.s file
;FOR x = 1,xmax-1,1 DO BEGIN
;    FOR y = 1,ymax-1,1 DO BEGIN
;        IF(data[x,y] GE -99.0) THEN data[x,y] = 0
;        IF (data[x,y] LT -99) THEN data[x,y] = 1
;    ENDFOR
;ENDFOR

;fits_write,"/n/Godiva2/jkrick/sep00/A3888/mask.fits", data, header



    ;want to make a cut on total flux as well, using the SExtractor param
    ;flux_max should be always greater than 5-7sigma (<10sigma)
    
    ;____________________________________________________________________________
    ;find the saturated objects
    ;IF (data(xcenter,ycenter) GE 39) THEN BEGIN
    ;    print, "saturated object at", xcenter,ycenter
    ;    satcount = satcount + 1
     ;   true = 1
    ;ENDIF


    ;____________________________________________________________________________
    ;ignore the bad columns
;    IF (xcenter GE 3030 AND xcenter LE 3038 AND ycenter GE 2100 AND ycenter LE 3450) $
;      THEN true = 1


 ;   IF (xcenter GE 3440 AND xcenter LE 3448 AND ycenter GE 2520 AND ycenter LE 2872) $
;THEN true = 1
    
;    IF (xcenter GE 3512 AND xcenter LE 3520 AND ycenter GE 480  AND ycenter LE 2528) $
;      THEN true = 1
	
 ;   IF (xcenter GE 3456 AND xcenter LE 3464 AND ycenter GE 48 AND ycenter LE 776) $
 ;     THEN true = 1
;	
;    IF (xcenter GE 3000 AND xcenter LE 3008 AND ycenter GE 48 AND ycenter LE 900) $
;      THEN true = 1

;    IF (xcenter GE 2600 AND xcenter LE 2608 AND ycenter GE 1 AND ycenter LE 1000) $
;      THEN true = 1

;    IF (xcenter GE 3102 AND xcenter LE 3352 AND ycenter GE 2864 AND ycenter LE 2872) $
;      THEN  true = 1

;    IF (xcenter GE 3440 AND xcenter LE 3520 AND ycenter GE 2520 AND ycenter LE 2528) $
;      THEN  true = 1

 
 ;__________________________________________________________________________

    ;make special blocks for the galaxies that were subtracted
;    IF(xcenter GT 787 AND xcenter LT 790 AND ycenter GT 1107 AND ycenter LT 1109) THEN BEGIN
;        sma = 30
;        smb = 30
;        return = del_star(pa, f, sma , smb, xcenter, ycenter, data)

;        print, "removing special galaxy at ",xcenter, ycenter,isoarea
;        true = 1
;        galcount = galcount+1
;        galflux = galflux + f
;    ENDIF

;    IF(xcenter GT 823 AND xcenter LT 826 AND ycenter GT 1083 AND ycenter LT 1086) THEN BEGIN
;        sma = 30
;        smb = 30
  ;      return = del_star(pa, f, sma , smb, xcenter, ycenter, data)

;        print, "removing special galaxy at ",xcenter, ycenter,isoarea
;        true = 1
;        galcount = galcount+1
;        galflux = galflux + f
;    ENDIF

 ;   IF(xcenter GT 922 AND xcenter LT 925 AND ycenter GT 1148 AND ycenter LT 1151) THEN BEGIN
 ;       sma = 30
 ;       smb = 30
 ;       return = del_star(pa, f, sma , smb, xcenter, ycenter, data)

;	print, "removing special galaxy at ",xcenter, ycenter, isoarea
;        true = 1
;        galcount = galcount+1
;        galflux = galflux + f
;    ENDIF

 ;   IF(xcenter GT 921 AND xcenter LT 924 AND ycenter GT 1130 AND ycenter LT 1133) THEN BEGIN
 ;      sma = 30
 ;      smb = 30
 ;      return = del_star(pa, f, sma , smb, xcenter, ycenter, data)

 ;      print, "removing special galaxy at ",xcenter, ycenter, isoarea
 ;      true = 1
 ;      galcount = galcount+1
 ;      galflux = galflux + f
 ;   ENDIF

 ;   IF(xcenter GT 1059 AND xcenter LT 1062 AND ycenter GT 855 AND ycenter LT 858) THEN BEGIN
 ;       sma = 30
 ;       smb = 30

;        return = del_star(pa, f, sma , smb, xcenter, ycenter, data)

;       print, "removing special galaxy at ",xcenter, ycenter, isoarea
;        true = 1
;        galcount = galcount+1
;        galflux = galflux + f
;    ENDIF

;    IF(xcenter GT 1058 AND xcenter LT 1061 AND ycenter GT 890 AND ycenter LT 893) THEN BEGIN
;        sma = 30
;        smb = 30
;        return = del_star(pa, f, sma , smb, xcenter, ycenter, data)

;	print, "removing special galaxy at ",xcenter, ycenter, isoarea
;        true = 1
;        galcount = galcount+1
;        galflux = galflux + f
;    ENDIF
    
;    IF(xcenter GT 977 AND xcenter LT 980 AND ycenter GT 910 AND ycenter LT 913) THEN BEGIN
;       sma = 30
;        smb = 30

;        return = del_star(pa, f, sma , smb, xcenter, ycenter, data)

;	print, "removing special galaxy at ",xcenter, ycenter, isoarea
;        true = 1
;        galcount = galcount+1
;        galflux = galflux + f
;    ENDIF

;    IF(xcenter GT 978 AND xcenter LT 981 AND ycenter GT 882 AND ycenter LT 885) THEN BEGIN
;        sma = 30
;        smb = 30

;        return = del_star(pa, f, sma , smb, xcenter, ycenter, data)

;	print, "removing special galaxy at ",xcenter, ycenter, isoarea
;        true = 1
;        galcount = galcount+1
;        galflux = galflux + f
;    ENDIF

;   IF(xcenter GT 1000 AND xcenter LT 1003 AND ycenter GT 883 AND ycenter LT 886) THEN BEGIN
;        sma = 30
;        smb = 30
 
;        return = del_star(pa, f, sma , smb, xcenter, ycenter, data)

;	print, "removing special galaxy at ",xcenter, ycenter, isoarea
;        true = 1
;        galcount = galcount+1
;        galflux = galflux + f
;    ENDIF















