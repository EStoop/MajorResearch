FUNCTION correct_pixel_phase,ch,x,y,flux,SUB=sub,verbose=verbose
;;; CH - channel number (1 or 2)
;;; X  - X centroid of star (middle of bottom left pixel is 0,0)
;;; Y  - Y centroid of star 
;;; FLUX - measured flux of star from aperture photometry
;;; /SUB - the observations were done at the nominal subarray FOV (otherwise
;;;       assumed to be the nominal full array FOV)
;;;
;;; CH, X, Y, FLUX are scalars or arrays with the same number of elements.
;;; 
dataloc = '/Users/jkrick/iwic/'

interpolated_relative_flux = MAKE_ARRAY(size=size(flux),VALUE=0.0)
xphase = x - ROUND(x)
yphase = y - ROUND(y)

IF keyword_set(SUB) EQ 1 THEN label = '_sub' ELSE label = '_full'
ich1 = WHERE(ch EQ 1,nch1)
ich2 = WHERE(ch EQ 2,nch2)
ion = WHERE(ch NE 0,non)
IF nch1 NE 0 THEN BEGIN
   RESTORE,dataloc+'pixel_phase_img_ch1'+label+'.sav'
   interpolated_relative_flux[ich1] = interp2d(relative_flux,x_center,y_center,xphase[ich1],yphase[ich1],/regular)
ENDIF
IF nch2 NE 0 THEN BEGIN
   RESTORE,dataloc+'pixel_phase_img_ch2'+label+'.sav'
   interpolated_relative_flux[ich2] = interp2d(relative_flux,x_center,y_center,xphase[ich2],yphase[ich2],/regular)
ENDIF
corrected_flux = flux/interpolated_relative_flux

IF keyword_set(verbose) EQ 1 THEN BEGIN
   print,'Channel   X    Y     XPHASE  YPHASE   ORIG FLUX   INTERP. REL. FLUX   CORR. FLUX'
   FOR i = 0,non-1 DO PRINT,ch[ion[i]],x[ion[i]],y[ion[i]],xphase[ion[i]],yphase[ion[i]],flux[ion[i]],interpolated_relative_flux[ion[i]],corrected_flux[ion[i]]
ENDIF

RETURN,corrected_flux
END
