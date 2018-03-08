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
D = 0.00044*f    ; for the moffat function

sb4 = 22.04 - 2.5*alog10((D*mymoffat(15.4))/(0.435^2.))   ;sb @ 4arcsec
sb2 = 22.04 - 2.5*alog10((D*mymoffat(7))/(0.435^2.))      ;sb @ 2arcsec
;print, "xcen, ycen,f, D, sb4,sb2",xcenter, ycenter, f, D, sb4, sb2
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
        intensity = 0.
        IF radius LT 45.0 AND radius NE 0.0 THEN BEGIN
            intensitysb = (3.8 + 2.5*alog10((43.5/radius)^(-3.9)) + sb2)
            intensity = (0.435^2.)*(10.^((22.04-intensitysb)/(2.5)))
        ENDIF

        IF (radius GE 45.0) THEN begin
            intensitysb = (3.9 + 2.5*alog10((43.5/radius)^(-2.0)) + sb2)
;            intensitysb2 = (4.3 + 2.5*alog10((43.5/radius)^(-3.9)) + sb2)

            intensity = (0.435^2.)*(10.^((22.04-intensitysb)/(2.5))); + (0.435^2.)*(10.^((22.04-intensitysb2)/(2.5)))
        ENDIF

        ;make sure not off the edge of the chip
        ;need to deal with going from fractional pixels to whole pixels
        IF (x LT xmax AND y LT ymax AND x GT 1 AND y GT 1 ) THEN BEGIN
            blank1[x,y] = blank1[x,y] + intensity
        ENDIF
;    ENDIF

    ENDFOR
ENDFOR

blank1 = shiftf(blank1, dx,dy)
blankdata = blankdata + blank1
return, 0              
END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO block_obj_A2556

close, /all
filename = strarr(1) 
;openr, lun5, "/n/Godiva2/jkrick/A114/original/listfull", /get_lun
;WHILE (NOT EOF(lun5)) DO BEGIN
;    readf, lun5, filename
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

;filename = strarr(1)            ;declare filename to be a stringarray
maxradius = 58.0                ; maximum radius in pixels to create a star
intensity = 0.0
bkgd = 0.0
smaarr = fltarr(20000)
smacount = 0
true = 0                        ;manual boolean variable

;size of the image
xmax =2866;2001;
ymax =3486;2161;2900;

;read in a blank image to which I can add my own stars			

FITS_READ, '/n/Godiva3/jkrick/A2556/original/fullblank.fits', blankdata, blankheader

filename = "/n/Godiva3/jkrick/A2556/original/fullr"

imagefile = filename + '.fits'
;datafile = '/n/Godiva3/jkrick/A2556/SExtractor.r.cat'
datafile = '/n/Godiva3/jkrick/A2556/block_obj.r.input'
newimagefile = filename + '.s.fits'	
outfile = filename + '.bblock.fits'


;run SExtractor on the image to find the objects
commandline = '/n/Godiva7/jkrick/Sex/sex ' + imagefile + " -c /n/Godiva3/jkrick/A2556/iclr.sex"
spawn, commandline


;read in the fits image with the stars that need to be subtracted
FITS_READ, imagefile, data, header   ;changed this to be the  1 second image (not 900)


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

;_______________________________________________________________________________
;saturated stars from USNO
satx = fltarr(38)
saty = fltarr(38)
satrmag = fltarr(38)
  
;r
satx = [355.00, 686.31,1535.07,1849.54,2319.32, 276.40,2213.71, 124.86, 470.32, 786.27,$
1091.60,1890.10,1028.37,1879.20,2039.59,2629.88, 651.50,  72.97, 249.75, 724.95,2377.39, $
703.90,1072.46,1155.64,1240.96,2099.30, 995.32,1308.23,1384.00,1139.94,1105.00,1630.12,$
1801.35,1979.30,1998.53,2287.39,2649.04,2778.51]

saty = [50.40 ,  79.48 , 271.55 , 292.33 ,  35.53 , 421.39 , 970.78 , 980.61 ,1248.81 ,$
1327.97 ,1672.62 ,1734.40 ,1927.50 ,1944.94 ,1889.02 ,1879.51 ,2064.24 ,2284.41 ,2300.62$
 ,2294.34 ,2074.71 ,2620.26 ,2574.94 ,2571.04 ,2525.64 ,2514.65 ,2843.49 ,2810.94 ,$
2950.53 ,3199.58 ,3301.70 ,3184.69 ,3297.39 ,3129.11 ,3022.94 ,2835.92 ,3271.04 ,3048.13]

satrmag=[13.9,13.6,13.9,15.6,12.1,9.4,15.2,13.6,15.6,14.7,14.7,13.1,15.0,12.9,13.9,12.3,13.6,$
        15.2,14.2,15.1,12.7,15.7,15.1,15.9,15.1,15.2,12.3,15.9,14.6,14.3,15.5,15.7,12.1,14.0,$
         12.9,15.8,10.2,14.7]

satBmag=[14.6,13.8,14.1,16.7,12.1,11.3,15.4,14.0,16.0,15.3,14.8,13.8,15.1,13.3,14.5,12.8,$
         14.0,16.1,14.9,15.4,13.0,15.0,15.6,16.4,15.1,15.5,12.6,16.5,15.1,14.8,15.4,15.7,$
         12.3,15.0,14.1,16.4,11.5,15.1]
  
;satx = satx - 500.
;saty = saty - 600.

FOR satstar = 0, 37, 1 DO BEGIN
;    holder = starmask (10^((satrmag[satstar] - 22.04)/(-2.5)), satx[satstar], saty[satstar],500.,xmax,ymax,blankdata)
    holder = starmask (10^((satrmag[satstar] - 22.04)/(-2.5)), satx[satstar], saty[satstar],500.,xmax,ymax,data)
ENDFOR

;______________________________________________________________________________________


;open the SExtractor output datafile (ccd*.cat)
OPENR, lun, datafile, /GET_LUN

WHILE (NOT EOF(lun)) DO BEGIN
    ;!!!! change this if you change the daofind.param , now usign godiva7/galfit/dao..
    READF, lun, o, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux,isocorflux,isoflux
;    print, "working on object", o
    
   
    ;___________________________________________________________________________
    ; ignore the other bad spots with flux less than 0
    IF (f LT 0) AND isoarea GT 50 THEN BEGIN    ;is a real object on top of the bad column
        sma = sqrt(isoarea/(!PI*(1-e)))
        smb = isoarea/(!PI*sma)
        sma = 2.3*sma           ; 3.0,2.3 ,1.6                     3.2;3.0*sma;2.6*sma;2*sma
        smb = 2.3*smb           ; 3.0,2.3,1.6                     3.2;3.0*smb;2.6*smb;2*smb
        return = del_star_A2556(pa, f, sma , smb, xcenter, ycenter, data)
        true = 1
    ENDIF

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
        return = del_star_A2556(pa, f, sma , smb, xcenter, ycenter, data)
        
        true = 1
        IF (xcenter GT 950) AND (xcenter LT 1800) AND (ycenter GT 1100) AND (ycenter LT 1950) THEN dimcount = dimcount +1

    ENDIF

    ;_________________________________________________________________________
    ;get rid of cosmic rays
    IF(isoarea LT 6.0 OR fwhm LT 2.9 AND true EQ 0) THEN BEGIN  ;2.5
        ;add 2 pixels to both a and b, and then block out to that radius
        ;print, "CR", o,xcenter,ycenter
        a = a + 2
        b = b + 2
       return = del_star_A2556(pa, f, a , b, xcenter, ycenter, data)
		
	;print, "removing cosmic ray", o, xcenter, ycenter
        true = 1
        IF (xcenter GT 950) AND (xcenter LT 1800) AND (ycenter GT 1100) AND (ycenter LT 1950) THEN crcount = crcount +1
    ENDIF
		
    ;__________________________________________________________________________
    ;get rid of the stars that are not very saturated
        
;  IF (fwhm GE 2.8 AND fwhm LE 4.2 AND true EQ 0) THEN BEGIN    

    IF (fwhm GE 2.9 AND fwhm LE 4.0 AND true EQ 0) THEN BEGIN   
        openw, dimlun, "/n/Godiva3/jkrick/A2556/tvmark.out", /GET_LUN, /append
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
            return = del_star_A2556(pa, f, 5, 5, xcenter, ycenter, data)
;            return = del_star_A2556(pa, f, 5, 5, xcenter, ycenter, blankdata)
            true = 1

            ;holder = starmask (f, xcenter, ycenter, maxradius,xmax,ymax,blankdata)
            IF (xcenter GT 436) AND (xcenter LT 2436) AND (ycenter GT 524) AND (ycenter LT 2524) THEN BEGIN
                starcount = starcount + 1
;                printf, starlun, "star in center with mag ", m
            ENDIF

        ENDIF
        IF (true EQ 0 AND m GT 17) THEN BEGIN
;            return = del_star_A2556(pa, f, 22, 22, xcenter, ycenter, data)
           return = del_star_A2556(pa, f, 30, 30, xcenter, ycenter, data)
;           return = del_star_A2556(pa, f, 30, 30, xcenter, ycenter, blankdata)
            true = 1
;           print, "removing star at   ",xcenter, ycenter,m
            holder = starmask (f, xcenter, ycenter, maxradius,xmax,ymax,data)
;            holder = starmask (f, xcenter, ycenter, maxradius,xmax,ymax,blankdata)
            IF (xcenter GT 436) AND (xcenter LT 2436) AND (ycenter GT 524) AND (ycenter LT 2524) THEN BEGIN
                starcount = starcount + 1
;                printf, starlun, "star in center with mag ", m
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

            holder = starmask (f, xcenter, ycenter, maxradius,xmax,ymax,data)
;            holder = starmask (f, xcenter, ycenter, maxradius,xmax,ymax,blankdata)

		
	;and now for the part where I block out the inner 'a' pixels of all of the stars
	;just to make sure that we get it all, and so I don't have to worry about
	;seeing, focus, ...

            b= 30;19               ;block out star up to a radius of B1 pixels
            a= 30;19
            return = del_star_A2556(pa, f, a , b, xcenter, ycenter, data)
;            return = del_star_A2556(pa, f, a , b, xcenter, ycenter, blankdata)

            IF (xcenter GT 436) AND (xcenter LT 2436) AND (ycenter GT 524) AND (ycenter LT 2524) THEN BEGIN
                starcount = starcount + 1  
;                printf, starlun, "star in center with mag ", m
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
;    IF (fwhm GT 4.2 AND true EQ 0) THEN BEGIN
    IF (fwhm GT 4.0 AND true EQ 0) THEN BEGIN ;5.1
  
        ;to make the dim galaxy masks big enough
        IF (isoarea LT 9) THEN isoarea = 9 ;small galaxies

        sma = sqrt(isoarea/(!PI*(1-e)))
        smb = isoarea/(!PI*sma)
		
        ;go out further than the isoarea
        IF (sma GT 15.0 OR smb GT 15.0) THEN BEGIN ;big galaxies
            IF(isoarea GT 2000.) THEN BEGIN
 ;               print, "central galaxy bigger than 2000", xcenter, ycenter, isoarea
                sma = ((2.0*sma)+1.44)/0.788   ; 2.6,2.0,1.4                 
                smb = ((2.0*smb)+1.44)/0.788   ;  2.6,2.0,1.4                 
            ENDIF ELSE BEGIN

                sma = ((2.3*sma)+1.44)/0.788   ; 3.0,2.3 ,1.6                 
                smb = ((2.3*smb)+1.44)/0.788   ; 3.0,2.3,1.6                  
            ENDELSE
                
        ENDIF ELSE BEGIN        ;normal galaxies
            sma = ((2.0*sma)+1.44)/0.788       ; 2.6,2.0,1.4                      
            smb = ((2.0*smb )+1.44)/0.788      ; 2.6,2.0,1.4                      
        ENDELSE
        
        theta = pa
        return = del_star_A2556(theta, f, sma , smb, xcenter, ycenter, data)
        return = del_star_A2556(pa, f, sma , smb, xcenter, ycenter, blankdata)

        ;openw, blocklun, "/n/Godiva3/jkrick/A2556/galfit/block.log", /GET_LUN, /append
        ;printf, blocklun, xcenter, ycenter, sma, smb
        ;close, blocklun
        ;free_lun, blocklun
       
        true = 1
        IF (xcenter GT 836) AND (xcenter LT 2036) AND (ycenter GT 924) AND (ycenter LT 2124) THEN BEGIN
            galcount = galcount+1
;            printf, gallun, xcenter, ycenter, sma, smb
        ENDIF

        dist = sqrt((1436. - xcenter)^2 + (1524.-ycenter)^2)
        IF (dist LT 616 ) THEN begin
            galbestflux = galbestflux + isoflux
            galisoflux = galisoflux + isocorflux
        ENDIF

        smaarr[smacount] = sma
        smacount = smacount + 1
    ENDIF

    IF (true EQ 0) THEN print, xcenter,ycenter,"did not get blocked"
 
    true = 0

ENDWHILE

smaarr = smaarr[0:smacount-1]
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
OPENR, lunicl, "/n/Godiva3/jkrick/A2556/iclr.sex", /GET_LUN
WHILE (NOT EOF(lunicl)) DO BEGIN
    readf, lunicl, j
    sxaddpar, header, "COMMENT",j
ENDWHILE
close, lunicl
free_lun, lunicl

;write a new fits file with the masked pixels, and one with the stars

fits_write, newimagefile, shift(data, -1,-1), header

fits_write, '/n/Godiva3/jkrick/A2556/original/mask.bblock2.galr.fits',shift(blankdata,-1,-1), blankheader

lessstars = data - blankdata   ;;;minus
;fits_write, '/n/Godiva3/jkrick/A2556/galfit/nogals.s-stars.fits', shift(lessstars, -1, -1), header

fits_write,'/n/Godiva3/jkrick/A2556/original/galr.bblock2.fits', data, header

;make a profile of what is left
;random center choice
;result = funcprof("/n/Godiva3/jkrick/A2556/galfit/nogals.s-stars", 492., 503.)

;undefine, subgals
;ENDWHILE

;device, true=24
;device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness

mydevice = !D.NAME
!p.multi = [0, 0, 1]
SET_PLOT, 'ps'

device, filename = '/n/Godiva3/jkrick/A2556/masksizes.ps', /portrait, $
  BITS=8, scale_factor=0.9 , /color

plothist, smaarr,xhist, yhist, bin=0.5,/noplot
plot, xhist, yhist, thick=3, charthick=3, xthick=3,ythick=3

device, /close
set_plot, mydevice


close, /all
END
