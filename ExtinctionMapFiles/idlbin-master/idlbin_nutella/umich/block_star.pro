;March 2002
;
;The purpose of this program is to use a moffat psf created by curfit.pro
; and use that profile to subtract away all of the stars in the image.
; after subtracting the stars, the centers of the stars are blocked out
; since it is nearly impossible to get the centers correct due to seeing/focus.
; bad columns, cosmic rays, and faint objects are also dealt with
;
;input: in the same directory need to have a file named
;	ccd****.cat with sextractor output from lookforstars2
;
;output: an image named ccd****.s.fits with the centers set equal to 
; 	-100, and an image called substars with the homemade stars on a 
;	a zero background ready to be subtracted
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

FUNCTION testcall
  print, "called from function testcall"
end

pro block_star

close, /all

;declare variables
o = 0.0          	;object number
xcenter = 0.0   		;xcenter
ycenter = 0.0			;ycenter
e = 0.0   		;ellipticity
a = 0.0   		;semi-major axis
b = 0.0			;semi-minor axis
m = 0.0   		;magnitude (un-corrected for zpt, X, ...)
j = 0.0 		;junk
f = 0.0			;flux
fwhm = 0.0		;full width at half maximum
isoarea = 0.0		;area inside of last isophote

filename = strarr(1)  	;declare filename to be a stringarray
maxradius = 58.0	; maximum radius in pixels to create a star
intensity = 0.0
bkgd = 0.0

true = 0		;manual boolean variable

;size of the image
xmax = 501
ymax = 501

;read in a blank image to which I can add my own stars			
FITS_READ, '/n/Sheriff1/jkrick/sep00/A3888/blocktest/blank1.fits', blankdata, blankheader

filename = '/n/Sheriff1/jkrick/sep00/A3888/blocktest/frame1'
imagefile = filename + '.fits'
datafile = filename + '.cat'
newimagefile = filename + '.s.fits'	

;read in the fits image with the stars that need to be subtracted
FITS_READ, imagefile, data, header

;open the SExtractor output datafile (ccd*.cat)
OPENR, lun, datafile, /GET_LUN

;get the bad columns out here
;data[1004:1014,1085:2494] = -100
;data[1449:1460,2017:3421] = -100
;data[1942:1954,638:2048] = -100
;data[2383:2392,1570:2965] = -100
;data[2845:2856,1355:2750] = -100

;read the data file into the variables until hit EOF
starcount = 0
galcount = 0
crcount = 0
dimcount = 0
 
totalflux = 0.0  				;in case I want to know
WHILE (NOT EOF(lun)) DO BEGIN
    READF, lun, o, xcenter, ycenter, a, b, e, f, m, j, isoarea, fwhm, pa, bkgd
	
;want to make a cut on total flux as well, using the SExtractor param
;flux_max should be always greater than 5-7sigma (<10sigma)



	;____________________________________________________________________________
	;ignore the bad columns
	IF (xcenter GE 1004 AND xcenter LE 1014 AND ycenter GE 1085 AND ycenter LE 2494) THEN true = 1
	
	IF (xcenter GE 1449 AND xcenter LE 1460 AND ycenter GE 2017 AND ycenter LE 3421) THEN true = 1
		
	IF (xcenter GE 1942 AND xcenter LE 1954 AND ycenter GE 638  AND ycenter LE 2048) THEN true = 1
		
	IF (xcenter GE 2383 AND xcenter LE 2392 AND ycenter GE 1570 AND ycenter LE 2965) THEN true = 1
		
	IF (xcenter GE 2845 AND xcenter LE 2856 AND ycenter GE 1355 AND ycenter LE 2750) THEN true = 1
		
	;____________________________________________________________________________
	; ignore the other bad spots with flux less than 0
	IF (f LT 0 ) THEN true = 1

	;__________________________________________________________________________
	;get rid of dim objects
	;these objects aren't putting enough into their wings to worry about them
	IF(f LE 0.1 AND isoarea GE 6.0 AND true EQ 0) THEN BEGIN
		;go out to two times the fwhm then zap 'em
		a = 1.5*fwhm
		b = 1.5*fwhm

		IF (isoarea LT 14) THEN isoarea = 14	   ;small galaxies

		sma = sqrt(isoarea/(!PI*(1-e)))
		smb = isoarea/(!PI*sma)
		return = del_star(pa, f, sma , smb, xcenter, ycenter, data)
		
		print, "reomving dim object", o,xcenter, ycenter
		true = 1
		dimcount = dimcount +1
	ENDIF

	;_________________________________________________________________________
	;get rid of cosmic rays
	IF(isoarea LT 6.0 OR fwhm LT 2.8 AND true EQ 0) THEN BEGIN
		;add 2 pixels to both a and b, and then block out to that radius
		a = a + 2
		b = b + 2
		return = del_star(pa, f, a , b, xcenter, ycenter, data)
		
		print, "removing cosmic ray", o, xcenter, ycenter
		true = 1
		crcount = crcount +1
	ENDIF
	
	
	;__________________________________________________________________________
	;get rid of the stars that are not very saturated
        
	IF (fwhm GE 2.8 AND fwhm LE 4.0 AND true EQ 0) THEN BEGIN
            IF (m LE 12.6) THEN BEGIN
                print,"ignoring supersaturated star at", xcenter, ycenter
                true = 1
            ENDIF 
                ;for the dim things that could be
                ;stars, or galaxies, or junk, I won't subtract a profile, instead just
                ;block out to 2"
            ;15.8 and 2" are both gotten by experimenting
            IF (m GT 15.8) THEN BEGIN
                return = del_star(pa, f, 8, 8, xcenter, ycenter, data)
                true = 1
                starcount = starcount + 1
            ENDIF
            IF (true EQ 0) THEN BEGIN
		; want to make the flux within a radius of 5" the same as
		; the flux within the inner 5" taken from mymoffat profile
		; that function is sitting in mymoffat.pro 
		; D is the scale factor for making mymoffat have the same flux as the star
		; 38.65 is gotten by making a normalized star with the same specs
		; 54 should be the number according to the math
		;need to change this if Sextractor params change or params in next loop change
		D = f * 1/38.6597
		
		X = 0.0					;radius of my homemade star
		WHILE (X LT maxradius) DO BEGIN
			intensity= D* mymoffat(X) 	;value of the moffat function at radius X
			
			IF (X LT 15) THEN i = !PI/80
			IF(X LT 51 AND X GE 15) THEN i = !PI/200
			;this loop makes a star
			FOR theta = -!PI/2, !PI/2, i DO BEGIN
				x1 = X*cos(theta)
				y1 = X*sin(theta)
			
				;make sure not off the edge of the chip
                                IF (xcenter+x1 LT xmax AND ycenter+y1 LT ymax AND xcenter-x1 GT 1 $
                                    AND ycenter-y1 GT 1AND xcenter+x1 GT 1 and ycenter+y1 GT 1 $
                                    AND xcenter-x1 LT xmax AND ycenter-y1 LT ymax) THEN BEGIN
					blankdata[xcenter + x1, ycenter + y1] = intensity
					blankdata[xcenter - x1, ycenter - y1] = intensity
				ENDIF
			ENDFOR
	
			X = X + 0.1
		ENDWHILE

		;and now for the part where I block out the inner 'a' pixels of all of the stars
		;just to make sure that we get it all, and so I don't have to worry about
		;seeing, focus, ...

		b= 19			;block out star up to a radius of B1 pixels
		a= 19
		return = del_star(pa, f, a , b, xcenter, ycenter, data)

		starcount = starcount + 1  
		totalflux = totalflux + f
		
		print, "removing star at   ", o,xcenter, ycenter
		true = 1

            ENDIF


	ENDIF

	;_______________________________________________________________________________________
	;get rid of the galaxies
	IF (fwhm GT 4.0 AND true EQ 0) THEN BEGIN
		;to make the dim galaxy masks big enough
		IF (isoarea LT 14) THEN isoarea = 14	   ;small galaxies

		sma = sqrt(isoarea/(!PI*(1-e)))
		smb = isoarea/(!PI*sma)
		
		;go out further than the isoarea
		IF (sma GT 15.0 OR smb GT 15.0) THEN BEGIN  ;big galaxies
			sma = 1*sma
			smb = 1*smb
		ENDIF ELSE BEGIN			    ;normal galaxies
			sma = 1*sma
			smb = 1*smb
		ENDELSE

		return = del_star(pa, f, sma , smb, xcenter, ycenter, data)

		print, "removing galaxy at ", o,xcenter, ycenter
		true = 1
		galcount = galcount+1
	ENDIF
	true = 0

ENDWHILE

print,"number of stars is:", starcount
print,"number of galaxies is:", galcount
print,"number of cosmic rays is:", crcount
print,"number of dim objects is:", dimcount
print,"total flux in 'stars' is", totalflux
testcall;close the SExtractor file, and free up the lun so it can be used again
close, lun
free_lun, lun

;for ccd3126 bleedout columns
;data[1052:1073,1022:1522] = -100
;data[1052:1073,496:996] = -100
	
;display the new image with the masked pixels
;twoimage = congrid(data,507.5,510)  ;makes smaller for display on screen
;tv, twoimage	;displays

;write a new fits file with the masked pixels, and one with the stars
fits_write, newimagefile, data, header
fits_write, '/n/Sheriff1/jkrick/sep00/A3888/ellipse/substars.fits',blankdata, blankheader


end




















