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


PRO galphottest

FOR run = 1, 1, 1 DO begin

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
xmax =3190;3011;1601;2866;2001;
ymax =3650;3461;1801;3486;2161;2900;

;read in a blank image to which I can add my own stars			

FITS_READ, '/n/Godiva4/jkrick/A2734/original/fullblank.fits', blankdata, blankheader

filename = "/n/Godiva4/jkrick/A2734/original/fullr"

imagefile = filename + '.fits'
datafile = '/n/Godiva3/jkrick/A2556/SExtractor.r.cat'
;datafile = '/n/Godiva4/jkrick/A4010/block_obj.r.input'
newimagefile = filename + '.s.fits'	
outfile = filename + '.bblock.fits'


;run SExtractor on the image to find the objects
commandline = '/n/Godiva7/jkrick/Sex/sex ' + imagefile + " -c /n/Godiva3/jkrick/A2556/iclr.sex"
;spawn, commandline


;read in the fits image with the stars that need to be subtracted
FITS_READ, imagefile, data, header   ;changed this to be the  1 second image (not 900)
datacopy = data

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

openw, outlun, "/n/Godiva4/jkrick/A2734/galphot", /get_lun
fracarr = fltarr(20000)
magarr = fltarr(20000)
;_______________________________________________________________________________

;______________________________________________________________________________________


;open the SExtractor output datafile (ccd*.cat)
OPENR, lun, datafile, /GET_LUN

WHILE (NOT EOF(lun)) DO BEGIN
    ;!!!! change this if you change the daofind.param , now usign godiva7/galfit/dao..
    READF, lun, o, xcenter, ycenter, a, b, e, f, m, isoarea, fwhm, pa, bkgd, apflux,isocorflux,isoflux
;    print, "working on object", o
    dist = sqrt((1530. - xcenter)^2 + (1677.-ycenter)^2)

    IF (f LT 0 ) THEN true = 1
     IF (apflux LE 1.0 ) AND (isoarea GE 6.0) THEN true = 1
      IF(isoarea LT 6.0) THEN true = 1

    ;________________________________________________________________________________
    ;get rid of the galaxies
;    IF (fwhm GT 4.2 AND true EQ 0) THEN BEGIN
    IF (fwhm GT 4.2 AND true EQ 0 AND dist LT 1000) THEN BEGIN ;5.1
  
        ;to make the dim galaxy masks big enough
;        IF (isoarea LT 9) THEN isoarea = 9 ;small galaxies

        sma = sqrt(isoarea/(!PI*(1-e)))
        smb = isoarea/(!PI*sma)
		
        ;go out further than the isoarea
        IF run EQ 0 THEN begin
            IF (sma GT 15.0 OR smb GT 15.0) THEN BEGIN ;big galaxies
                IF(isoarea GT 2000.) THEN BEGIN
 ;               print, "central galaxy bigger than 2000", xcenter, ycenter, isoarea
                    sma = ((1.4*sma)+1.9)/0.879 ; 2.6,2.0,1.4                 
                    smb =  ((1.4*smb)+1.9)/0.879 ;  2.6,2.0,1.4                 
                ENDIF ELSE BEGIN

                    sma =  ((1.6*sma)+1.9)/0.879 ; 3.0,2.3 ,1.6                 
                    smb =  ((1.6*smb)+1.9)/0.879 ; 3.0,2.3,1.6                  
                ENDELSE
                
            ENDIF ELSE BEGIN    ;normal galaxies
                sma =  ((1.4*sma)+1.9)/0.879   ; 2.6,2.0,1.4                      
                smb =  ((1.4*smb)+1.9)/0.879   ; 2.6,2.0,1.4                      
            ENDELSE
            
        ENDIF 
        IF run EQ 1 THEN begin
            IF (sma GT 15.0 OR smb GT 15.0) THEN BEGIN ;big galaxies
                IF(isoarea GT 2000.) THEN BEGIN
 ;               print, "central galaxy bigger than 2000", xcenter, ycenter, isoarea
                    sma =  ((2.0*sma)+1.9)/0.879 ; 2.6,2.0,1.4                 
                    smb =  ((2.0*smb)+1.9)/0.879 ;  2.6,2.0,1.4                 
                ENDIF ELSE BEGIN

                    sma =  ((2.3*sma)+1.9)/0.879 ; 3.0,2.3 ,1.6                 
                    smb =  ((2.3*smb)+1.9)/0.879 ; 3.0,2.3,1.6                  
                ENDELSE
                
            ENDIF ELSE BEGIN    ;normal galaxies
                sma =  ((2.0*sma)+1.9)/0.879   ; 2.6,2.0,1.4                      
                smb =  ((2.0*smb)+1.9)/0.879   ; 2.6,2.0,1.4                      
            ENDELSE
            
        ENDIF 
        IF run EQ 2 THEN begin
            IF (sma GT 15.0 OR smb GT 15.0) THEN BEGIN ;big galaxies
                IF(isoarea GT 2000.) THEN BEGIN
 ;               print, "central galaxy bigger than 2000", xcenter, ycenter, isoarea
                    sma =  ((2.6*sma)+1.9)/0.879 ; 2.6,2.0,1.4                 
                    smb =  ((2.6*smb)+1.9)/0.879 ;  2.6,2.0,1.4                 
                ENDIF ELSE BEGIN

                    sma =  ((3.0*sma)+1.9)/0.879 ; 3.0,2.3 ,1.6                 
                    smb =  ((3.0*smb)+1.9)/0.879 ; 3.0,2.3,1.6                  
                ENDELSE
                
            ENDIF ELSE BEGIN    ;normal galaxies
                sma =  ((2.6*sma)+1.9)/0.879   ; 2.6,2.0,1.4                      
                smb =  ((2.6*smb)+1.9)/0.879   ; 2.6,2.0,1.4                      
            ENDELSE
            
        ENDIF 
        
        theta = pa
        return = del_star_gal(pa, f, sma , smb, xcenter, ycenter, datacopy)
       
;        newdata = data - datacopy
        masksum = 0
        galsum= 0
        IF xcenter - 5*sma GT 0 AND xcenter + 5*sma LT xmax $
          AND ycenter - 5*sma GT 0 AND ycenter + 5*sma LT xmax THEN BEGIN
            FOR x = xcenter - 5*sma, xcenter+5*sma DO BEGIN
                FOR y = ycenter - 5*sma, ycenter+5*sma DO BEGIN

                    masksum = masksum + datacopy[x,y]
                    galsum = galsum + data[x,y]

                ENDFOR
            ENDFOR
        galsum= galsum - masksum
;        print," a galaxy", xcenter, ycenter, galsum, isocorflux, sma, smb
         
        fracarr[j] = isocorflux / galsum
        magarr[j] = 22.04 - 2.5*alog10(isocorflux)
        printf,outlun, xcenter, ycenter, galsum, isocorflux, isoarea

;        fits_write, "/n/Godiva4/jkrick/A2734/mask.fits", datacopy, header
        endiF        

        j = j + 1
        datacopy = data
    ENDIF
    
    true = 0
    
ENDWHILE
;smaarr = smaarr[0:smacount-1]

fracarr = fracarr[0:j - 1]
magarr = magarr[0:j-1]

device, true=24
device, decomposed=0

colors = GetColor(/load, Start=1)
close, /all		;close all files = tidiness


ps_open, file = "/n/Godiva4/jkrick/galphot.ps", /portrait, /color, ysize = 5
plot, magarr, fracarr, psym = 2, thick = 3, xtitle = "isophotal corrected magnitude", $
  ytitle = "isophotal / masked ",yrange = [0,2.0],$
   charthick = 3,xstyle=1, ystyle = 1,xrange=[13,24],$
  xthick = 3, ythick = 3
ps_close, /noprint, /noid

close, lun
free_lun, lun



close, outlun
free_lun, outlun
close, /all

ENDFOR
END
