;April 2002
;
;The purpose of this program is to block out the galaxies, given
; isophote measurements from SExtractor
;
;input: in the same directory need to have a file named
;	ccd****.cat with sextractor output from lookforstars2
;
;output: an image named ccd****.g.fits with the galaxies set equal to 
; 	-100
;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


pro block_gal

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
isoarea = 0.0		;isoarea
filename = strarr(1)  	;declare filename to be a stringarray
maxradius = 58.0
intensity = 0.0
bkgd = 0.0

;size of the image
xmax = 2030
ymax = 2040


filename = '/n/Sheriff1/jkrick/sep00/A3888/ellipse/sextest'
imagefile = filename + '.fits'
datafile = filename + '.cat'
newimagefile = filename + '.g.fits'	

;read in the fits image with the galaxies that need to be subtracted
FITS_READ, imagefile, data, header

;open the SExtractor output datafile (ccd*.cat)
OPENR, lun, datafile, /GET_LUN

;get the bad columns out here

;read the data file into the variables until hit EOF
count = 0 
totalflux = 0.0  				;in case I want to know
WHILE (NOT EOF(lun)) DO BEGIN
   	READF, lun, o, xcenter, ycenter, a, b, e, f, m, j,isoarea, fwhm, pa, bkgd
	
	

	;get rid of the galaxies  /////////////////////////////////////////////////////
	IF (e LT 0.13 AND f GT 0 AND fwhm GT 2.0 AND fwhm LE 5.0) THEN BEGIN
		b= 19			;block out star up to a radius of B1 pixels
		a= 19
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
				IF (xcenter+x1 LT xmax AND ycenter+y1 LT ymax AND xcenter-x1 GT 1 AND ycenter-y1 GT 1AND xcenter+x1 GT 1 and ycenter+y1 GT 1 AND xcenter-x1 LT xmax AND ycenter-y1 LT ymax) THEN BEGIN
					blankdata[xcenter + x1, ycenter + y1] = intensity
					blankdata[xcenter - x1, ycenter - y1] = intensity
				ENDIF
			ENDFOR
	
			X = X + 0.1
		ENDWHILE

		;and now for the part where I block out the inner 'a' pixels of all of the stars
		;just to make sure that we get it all, and so I don't have to worry about
		;seeing, focus, ...

		return = del_star(pa, f, a , b, xcenter, ycenter, data)

		count = count + 1  
		totalflux = totalflux + f
	ENDIF
        
ENDWHILE

print,"number of stars is:", count
print,"total flux in 'stars' is", totalflux

;close the SExtractor file, and free up the lun so it can be used again
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
fits_write, '/n/Whalen1/jkrick/ICL/A3888/substars.fits',blankdata, blankheader


end
