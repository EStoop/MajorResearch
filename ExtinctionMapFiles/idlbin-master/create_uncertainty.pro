pro create_uncertainty, im, h,  WARM=warm
                       
;+
; PROCEDURE:
;    CREATE_UNCERTAINTY
;
; PURPOSE: 
;
; USAGE:
;  fits_read, '/Users/jkrick/irac_warm/pcrs_planets/wasp43/r42615040/ch2/bcd/SPITZER_I2_42615040_0150_0000_1_bcd.fits', im, h
;   create_uncertainty, im, h, /warm
;
; INPUT:
;    im: array containing  data
;      h: header of the data fits file
;
; OUTPUTS:
;
; OPTIONAL FLAGS:
;    WARM: if set, use warm IRAC factors for noise
;
; METHOD:
;
;
; HISTORY:
;pulled relevant parts from get_centroids_for_calstars_jk.pro
;-

; Array of readnoise in electrons, per channel for 0.02, 0.1, 0.4 frametimes
; 2s sub, 0.4, 0.6, 1.2, 2, 6, 12, 30, 100 second full
	if (keyword_set(WARM)) then begin
		readnoise = [[23.9, 23.8, 12.7, 11.9], [16.9, 16.8, 9.0, 8.4], $
	             [11.8, 12.1, 9.1, 7.1], [9.4, 9.4, 8.8, 6.7], $
	             [23.9, 23.8, 12.7, 11.9], [23.9, 23.8, 12.7, 11.9], $
	             [23.9, 23.8, 12.7, 11.9], [11.8, 12.1, 9.1, 7.1], $
	             [9.4, 9.4, 8.8, 6.7], [9.4, 9.4, 8.8, 6.7], $
	             [6.0, 8.4, 4.5, 4.2], [6.0, 8.4, 4.5, 4.2]] * 1.D 
; Bump up readnoise by 5% for warm mission for FPA 1
	    readnoise[0, *] = readnoise[0, *] * 1.05	    
	endif else $
		readnoise = [[23.9, 23.8, 12.7, 11.9], [16.9, 16.8, 9.0, 8.4], $
	             [11.8, 12.1, 9.1, 7.1], [9.4, 9.4, 8.8, 6.7], $
	             [23.9, 23.8, 12.7, 11.9], [23.9, 23.8, 12.7, 11.9], $
	             [23.9, 23.8, 12.7, 11.9], [11.8, 12.1, 9.1, 7.1], $
	             [9.4, 9.4, 8.8, 6.7], [9.4, 9.4, 8.8, 6.7], $
	             [8.5, 8.4, 4.5, 4.2], [8.5, 8.4, 4.5, 4.2]] * 1.D

; Gain (conversion from DN to e-) per channel
	if (keyword_set(WARM)) then gain = [3.7, 3.7, 3.8, 3.8] * 1.D$
	else gain = [3.3, 3.7, 3.8, 3.8] * 1.D
	
; Sigma level for outlier rejection
	sigma = 3.0

; Get channel number of image
	nch = sxpar(h, 'CHNLNUM')

; Set image dimensions
        nx = sxpar(h, 'NAXIS1')
        ny = sxpar(h, 'NAXIS2')
	naxis = sxpar(h, 'NAXIS')
	if (naxis eq 3) then ns = sxpar(h, 'NAXIS3') else ns = 1
	

; Check to see if BCD or not! (1: BCD, 0: Raw, -1: pbcd?)
	areadmod = sxpar(h, 'AREADMOD', COUNT=acount)
	aortime = sxpar(h, 'AORTIME', COUNT=xcount)
	if (xcount gt 0) then bcdflag = 1 else bcdflag = 0
	if (acount eq 0) then bcdflag = -bcdflag

        ft = sxpar(h, 'FRAMTIME')

; If raw image need to convert to cube
	if (bcdflag eq 0 and areadmod eq 1) then begin
		nx = 32
		ny = 32
		ns = 64
		
; transform to cube
		cube = reform(im, nx, ny, ns)
; Perform InSb sign flip
		if (nch lt 3) then cube = fix(65535 - cube)
; Rotate to BCD coordinates
		if (nch lt 3) then cube = reverse(cube, 2)
; Blank bad pixel - only know about bad pixel in channel 1 raws so far
		if (nch eq 1 and not keyword_set(WARM)) then cube[24, 6, *] = 0
	endif else cube = im

; Get sigma image
; Need to determine which readnoise and sat level to use
; 0.02, 0.1, 0.4, 2.0 sub, 0.4, 0.6, 1.2, 2, 6, 12, 30, 100 frametimes
	if (abs(ft-0.02) lt 0.01) then findex = 0 $
	else if (abs(ft-0.1) lt 0.01) then findex = 1 $
	else if (abs(ft-0.4) lt 0.01 and areadmod eq 1) then findex = 2 $
	else if (abs(ft-0.4) lt 0.01 and areadmod eq 0) then findex = 4 $
	else if (abs(ft-0.6) lt 0.01) then findex = 5 $
	else if (abs(ft-1.2) lt 0.01) then findex = 6 $
	else if (abs(ft-2.) lt 0.01 and areadmod eq 1) then findex = 3 $
	else if (abs(ft-2.) lt 0.01 and areadmod eq 0) then findex = 7 $
	else if (abs(ft-6.) lt 0.01) then findex = 8 $
	else if (abs(ft-12.) lt 0.01) then findex = 9 $
	else if (abs(ft-30.) lt 0.01) then findex = 10 $
	else if (abs(ft-100.) lt 0.01) then findex = 11 else findex = 0

; If no sigma image then make one -- BCD case first poisson + readnoise
	if (bcdflag ne 0) then begin
		sbtoe = sxpar(h, 'GAIN') * sxpar(h, 'EXPTIME') / sxpar(h,'FLUXCONV')
		sig = abs(cube) * sbtoe
; Make noise positive definite
;			if (min(sig) lt 0.0) then sig = sig - min(sig)
		sig = sig + readnoise[nch-1, findex] * readnoise[nch-1, findex]
		sig = sqrt(sig)
; Convert back to MJy/sr
		sig = sig / sbtoe
	endif else begin
; Sigma image for raw data, poisson noise plus readnoise
		dntoe = gain[nch-1] 
		sig = cube * dntoe
		sig = sig + readnoise[nch-1, findex] * readnoise[nch-1, findex]
		sig = sqrt(sig)
; Convert back to DN
		sig = sig / dntoe
	endelse

fits_write, '/Users/jkrick/idlbin/test.fits', sig, h		
;
end

