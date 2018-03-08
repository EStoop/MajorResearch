FUNCTION fitspec,freq,spec,weight,f_char=f_char,parm0=parm0,active=active, $
         cmb_temp=cmb_temp,bb1_temp=bb1_temp,bb2_temp=bb2_temp, $
         opdep1=opdep1,opdep2=opdep2,emis1=emis1,emis2=emis2, $
         mask=mask,tolerance=tolerance,out_parm=out_parm,parm_sig=parm_sig, $
         residual=residual,nodust=nodust,dust=dust,fail=fail,chi2=chi2, $
         n_iter=n_iter,funits=funits,eplees=eplees,wcm2=wcm2,mjy=mjy, $
	 no_show=no_show
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    FITSPEC is a multiple blackbody fitting facility.
;
;DESCRIPTION:
;    This IDL function fits spectra to a three blackbody model.
;    One blackbody typically corresponds to the cosmic background
;    while the other two can be used to fit the residual dust.  These
;    last two spectra include corresponding optical depth and emissivity
;    parameters.  Therefore there are seven fitting parameter in all.
;    In most cases it is not possible to fit all seven at once,
;    therefore parameters can be designated as fixed.  The initial
;    values of the parameters are supplied by the user (either using a
;    command line keyword or interactively) as well as their
;    variability (fixed or variable).  It is also possible to designate
;    fixed value parameters that are different for each input spectrum.
;    Particular ranges of the spectra can be masked if desired.
;    This is useful for avoiding spectral lines.
;
;    Model: fit = P(T0;freq) + op_depth1 * (freq/f_char)^emiss1 * P(T1;freq)
;                            + op_depth2 * (freq/f_char)^emiss2 * P(T2;freq)
;
;    where P(T;freq) is a Planck (Blackbody) spectrum.
;
;CALLING SEQUENCE:
;    fit = fitspec(freq,spec,[weight],[parm0=parm0],[active=active], $
;      [cmb_temp=cmb_temp],[bb1_temp=bb1_temp],[bb2_temp=bb2_temp], $
;      [opdep1=opdep1],[opdep2=opdep2],[emis1=emis1],[emis2=emis2], $
;      [mask=mask],[tolerance=tolerance],[out_parm=out_parm] $
;      [parm_sig=parm_sig],[residual=residual],[nodust=nodust], $
;      [dust=dust],[fail=fail],[chi2=chi2],[n_iter=n_iter], $
;      [funits=funits],[/eplees],[/wcm2],[/mjy],[/no_show])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    freq        I     flt arr             Spectra Frequency
;    spec        I     flt arr             Spectra
;    [weight]    I     flt arr             Spectral Weights
;    [parm0]     I     flt arr             Initial Parameter Values
;    [active]    I     byte arr            Parameter Variability
;                                          [1 - vary, 0 - fixed]
;    [cmb_temp]  I     flt arr             CMP Temp value array
;    [bb1_temp]  I     flt arr             BB 1 Temp value array
;    [bb2_temp]  I     flt arr             BB 2 Temp value array
;    [opdep1]    I     flt arr             Optical Depth 1 value array
;    [opdep2]    I     flt arr             Optical Depth 2 value array
;    [emis1]     I     flt arr             Emissivity 1 value array
;    [emis2]     I     flt arr             Emissivity 1 value array
;    [mask]      I     int arr             Masked element array
;    [tolerance] I     flt scl             Convergence Tolerance
;    [out_parm]  O     flt arr             Final Fit Parameter values
;    [parm_sig]  O     flt arr             Final Fit Parameter sigmas
;    [residual]  O     flt arr             Fit Residuals
;    [nodust]    O     flt arr             Spectrum w/dust fit subtracted
;    [dust]      O     flt arr             Spectrum w/cmb fit subtracted
;    [fail]      O     byte arr            Fit Failure Flag
;    [chi2]      O     flt arr             Chi_squared/df
;    [n_iter]    O     int arr             # iterations for fit
;    fit         O     flt arr             Model fits to spectra
;    [funits]    I     str                 Frequency Units
;                                          'icm' (default), 'microns',
;                                          'GHz'
;    [/eplees]   I     qual                Results in eplees (default)
;    [/wcm2]     I     qual                Results in W/cm^2/sr/K
;    [/mjy]      I     qual                Results in MJy/sr
;    [/no_show]  I     qual                Suppress initial parameter
;                                          display
;
;    The spectrum and weight arrays must have the same number of elements
;    and the same number of dimensions, which can be one (a single
;    spectrum), two (an array of spectra), or three (an array of spectra
;    from a skymap).  The last dimension should correspond to the
;    frequency array.
;
;    The order of the values in the initial parameter array, 'parm0',
;    is:  Pure Blackbody Temperature (in K), First Residual Blackbody
;    Temperature, First Residual Optical Depth, First Residual
;    Emissivity, Second Residual Blackbody Temperature, Second
;    Residual Optical Depth, Second Residual Emissivity.
;
;    The cmb_temp, bb1_temp, bb2_temp, opdep1, opdep2, emis1, and emis2
;    input arrays are used to define fixed parameter values for the
;    individual spectra.
;
;    The mask vector contains the indicies of those frequencies
;    where no attempt is made to fit the spectra.
;
;    The tolerance value is the minimum relative change in the
;    chi_squared for the fit that must be reached for convergence.
;    The default value is 1e-4.  A maximum of 20 iterations is
;    allowed.
;
;    If the fit fails to converge in 20 iterations or one of the fitted
;    values becomes negative, then the fit is considered to have failed
;    and the failure flag is set to 1, otherwise it is set to 0.  It is,
;    of course, the responsibility of the user to check the final fitted
;    values for credibility.
;    
;WARNINGS:
;    Some initial parameter values may lead to non-convergence.  The
;    user is advised to experiment with various sets of initial values
;    on a representative set of spectra to determine a sufficiently
;    robust set.  We have found in practice, that the optical depth
;    parameters converge better if they are overestimated.
;
;EXAMPLE: 
;
; To fit a set of high frequency FIRAS spectra to a pure blackbody
; (for the CMB) and a dust blackbody with the emissivity fixed at 1.5:
;
; fit_hi = fitspec(freq_hi,spec_hi,weight_hi, $
;                  parm0=[2.72,20,2e-5,1.5,0,0,0],f_char=60.0, $
;                  active=[1,1,1,0,0,0,0],mask=[82,107,108,109], $
;                  out_parm=parm_hi,fail=fail_hi,parm_sig=parm_sig)
;
; The initial values of the parameters are 2.72K for the CMB BB, 20K
; for the dust BB, 2e-5 for the dust optical depth, and 1.5 for the
; dust emissivity.  Setting the second residual blackbody temperature
; to zero eliminates it from the model.  The first three parameters
; are allowed to vary.  The emissivity is fixed (active(3) = 0).
; The spectrum values at the four frequency points (82,107,108,109) are
; ignored (masked) for the fit.  The optical depth coefficients give
; the ratio of the fit to a Planck at the fit temperature at 60 icm.
; The final fit parameters are returned in the parm_rhss array and their
; sigmas in the parm_sig array.
;
; If either the 'parm0' or 'active' keywords were missing, then the
; user would be prompted for their values.
;
; ------
;
; To fit a set of FIRAS low frequency spectra to a CMB BB and a
; residual dust spectrum using the dust temperatures derived
; from the high frequency fit as fixed parameters:
;
; fit_lo = fitspec(freq_lo,spec_lo,weight_lo, $
;                  parm0=[2.72,0,2e-5,1.5,0,0,0],f_char=60.0, $
;                  active=[1,1,0,0,0,0,0], bb1_temp=parm_hi(*,1), $
;                  out_parm=parm_lo,fail=fail_lo)
;
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Check for compatibility of frequency, spectra, and weight arrays.
;    Check for variable (pixel dependent) fixed parameter arrays.
;    If parm0 not defined then query user for inital parameter values.
;    If active not defined then query user for parameter variability.
;    Loop through spectra calling fit_cdd to determine fit.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: fit_cdd
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jan 93
;    Initial Delivery  SPR 10579
;
;.TITLE
;Routine FITSPEC
;-

IF (KEYWORD_SET(funits) EQ 0) THEN funits = 'I' ELSE BEGIN
	i = STRPOS('IMG',STRUPCASE(STRMID(funits,0,1)))
	IF (i EQ -1) THEN BEGIN
		str = 'Improper Frequency Units: ' + funits
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDELSE

res_units = 'E'		; default unit is 'eplees'
IF (KEYWORD_SET(eplees)) THEN res_units = 'E'
IF (KEYWORD_SET(wcm2)) THEN res_units = 'W'
IF (KEYWORD_SET(mjy)) THEN res_units = 'M'


default = [2.7,20,2e-5,1.5,300,1e-9,1]
dflt_str = ['<2.7>','<20>','<2e-5>','<1.5>','<300>','<1e-9>','<1>']
parm_title=['Cosmic Background Temperature ', $
            'First Black Body Temperature ', $
            'First Black Body Optical Depth ', $
            'First Black Body Emissivity ', $
            'Second Black Body Temperature ', $
            'Second Black Body Optical Depth ', $
            'Second Black Body Emissivity ']

inp = ''

sz = SIZE(spec)
sz_weight = SIZE(weight)
n_dim = sz(0)
n_dim_weight = sz_weight(0)
n_freq = N_ELEMENTS(freq)

IF (n_freq NE sz(sz(0))) THEN BEGIN
	str = 'Frequency and Spectrum Arrays not compatible'
	MESSAGE,str,/CONT
	out_fit = -1
	GOTO, exit
ENDIF

IF (sz_weight(0) EQ 0) THEN BEGIN
	weight = 0*spec + 1.0
	sz_weight = SIZE(weight)
	n_dim_weight = sz_weight(0)
ENDIF
	; if weights not present set them to 1


IF (n_dim NE n_dim_weight) THEN BEGIN
	str = 'Spectrum and Weight Arrays must have same number of dimensions'
	MESSAGE,str,/CONT
	out_fit = -1
	GOTO, exit
ENDIF

IF (TOTAL(sz(0:n_dim)) NE TOTAL(sz_weight(0:n_dim))) THEN BEGIN
	str = 'Spectrum and Weight Arrays must have same number of elements'
	MESSAGE,str,/CONT
	out_fit = -1
	GOTO, exit
ENDIF

pix_parm = BYTARR(7)

IF (KEYWORD_SET(cmb_temp)) THEN BEGIN
	pix_parm(0) = 1
	sz_fix = SIZE(cmb_temp)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and CMB_TEMP Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF

IF (KEYWORD_SET(bb1_temp)) THEN BEGIN
	pix_parm(1) = 1
	sz_fix = SIZE(bb1_temp)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and BB1_TEMP Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF

IF (KEYWORD_SET(opdep1)) THEN BEGIN
	pix_parm(2) = 1
	sz_fix = SIZE(opdep1)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and OPDEP1 Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF

IF (KEYWORD_SET(emis1)) THEN BEGIN
	pix_parm(3) = 1
	sz_fix = SIZE(emis1)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and EMIS1 Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF

IF (KEYWORD_SET(bb2_temp)) THEN BEGIN
	pix_parm(4) = 1
	sz_fix = SIZE(bb2_temp)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and BB2_TEMP Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF


IF (KEYWORD_SET(opdep2)) THEN BEGIN
	pix_parm(5) = 1
	sz_fix = SIZE(opdep2)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and OPDEP2 Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF


IF (KEYWORD_SET(emis2)) THEN BEGIN
	pix_parm(6) = 1
	sz_fix = SIZE(emis2)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and EMIS2 Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF

IF (TOTAL(pix_parm) NE 0) THEN pp_flag = 1 ELSE pp_flag = 0

IF (KEYWORD_SET(parm0) EQ 0) THEN BEGIN
	parm0 = FLTARR(7)
	FOR i=0,6 DO BEGIN
		IF (pix_parm(i) EQ 1) THEN GOTO,lbl1
		IF ((i EQ 2 OR i EQ 3) AND (parm0(1) EQ 0) $
		AND (pix_parm(1) EQ 0)) $
		THEN GOTO,lbl1
		IF (i EQ 4 AND parm0(1) EQ 0) THEN GOTO,lbl1
		IF ((i EQ 5 OR i EQ 6) AND (parm0(4) EQ 0) $
		AND (pix_parm(4) EQ 0)) $
		THEN GOTO,lbl1
		READ,parm_title(i)+dflt_str(i)+': ',inp
		IF (inp EQ '') THEN parm0(i) = default(i) $
		ELSE parm0(i) = FLOAT(inp)
		; if parameter read in from array then skip to next
		; if BB temp=0 then don't prompt for op dep and emis

	lbl1:
	ENDFOR
ENDIF


IF (KEYWORD_SET(active) EQ 0) THEN BEGIN
	PRINT,' '
	active = BYTARR(7)
	FOR i=0,6 DO BEGIN
		IF (pix_parm(i) EQ 1) THEN GOTO,lbl2
		IF ((i GE 1 AND i LE 3) AND (parm0(1) EQ 0) $
		AND (pix_parm(1) EQ 0)) $
		THEN GOTO,lbl2
		IF ((i GE 4 AND i LE 6) AND (parm0(4) EQ 0) $
		AND (pix_parm(4) EQ 0)) $
		THEN GOTO,lbl2
		READ,'Vary ' + parm_title(i) + '<Yes>? ',inp
		IF (inp EQ '' OR STRUPCASE(STRMID(inp,0,1)) EQ 'Y') $
		THEN active(i) = 1
	lbl2:
	ENDFOR
ENDIF

IF (TOTAL(pix_parm) NE 0) THEN active(WHERE(pix_parm) GT 0) = 0

IF (KEYWORD_SET(no_show) EQ 0) THEN BEGIN
	PRINT,' '
	PRINT,'Initial Paramater Values'
	PRINT,'------------------------'
	FOR i=0,6 DO BEGIN
		IF (pix_parm(i) EQ 0) THEN str = STRING(parm0(i)) $
		ELSE str = 'Read From Input Array'
		PRINT,FORMAT='(a,t40,a,t55,a)',parm_title(i)+': ',str, $
					       STRING(active(i))
	ENDFOR
	PRINT,' '
ENDIF

IF (KEYWORD_SET(tolerance) EQ 0) THEN tolerance = 1.e-4

CASE n_dim OF

   1:	BEGIN

	out_parm = FLTARR(7)
	fail = BYTARR(1)
	chi2 = FLTARR(1)
	n_iter = INTARR(1)
	out_fit = FLTARR(n_freq)
	parm = parm0
	y_fit = fit_cdd(freq,spec,weight,parm,active,mask, $
		        tolerance,c2,n_i,fail0,nodust=nodust, $
		        dust=dust,funits=funits,res_units=res_units, $
			parm_sig=parm_sig,f_char=f_char)

	out_fit = y_fit
	out_parm = parm
	fail(0) = fail0
	chi2(0) = c2
	n_iter(0) = n_i

	END


   2:	BEGIN

	out_parm = FLTARR(sz(1),7)
	parm_sig = FLTARR(sz(1),7)
	fail = BYTARR(sz(1))
	chi2 = FLTARR(sz(1))
	n_iter = INTARR(sz(1))
	out_fit = FLTARR(sz(1),sz(2))
	nodust = FLTARR(sz(1),sz(2))
	dust = FLTARR(sz(1),sz(2))

	FOR i=0l,sz(1)-1 DO BEGIN
	IF (i/100 eq i/100.) THEN PRINT,'Processing Spectrum:',i
	parm = parm0

	IF (pp_flag NE 0) THEN BEGIN
		IF (pix_parm(0)) THEN parm(0) = CMB_TEMP(i)	
		IF (pix_parm(1)) THEN parm(1) = BB1_TEMP(i)	
		IF (pix_parm(2)) THEN parm(2) = OPDEP1(i)	
		IF (pix_parm(3)) THEN parm(3) = EMIS1(i)
		IF (pix_parm(4)) THEN parm(4) = BB2_TEMP(i)
		IF (pix_parm(5)) THEN parm(5) = OPDEP2(i)
		IF (pix_parm(6)) THEN parm(6) = EMIS2(i)
	ENDIF

	y_fit = fit_cdd(freq,spec(i,*),weight(i,*),parm,active, $
		        mask,tolerance,c2,n_i,fail0, $
                        nodust=nodust0,dust=dust0, $
		        funits=funits,res_units=res_units, $
			parm_sig=parm_sig0,f_char=f_char)

	out_fit(i,*) = y_fit
	out_parm(i,*) = parm
	parm_sig(i,*) = parm_sig0
	nodust(i,*) = nodust0
	dust(i,*) = dust0
	fail(i) = fail0
	chi2(i) = c2
	n_iter(i) = n_i
	ENDFOR

	END


   3:	BEGIN

	out_parm = FLTARR(sz(1),sz(2),7)
	parm_sig = FLTARR(sz(1),sz(2),7)
	fail = BYTARR(sz(1),sz(2))
	chi2 = FLTARR(sz(1),sz(2))
	n_iter = INTARR(sz(1),sz(2))
	out_fit = FLTARR(sz(1),sz(2),n_freq)
	nodust = FLTARR(sz(1),sz(2),n_freq)
	dust = FLTARR(sz(1),sz(2),n_freq)

	FOR j=0l,sz(2)-1 DO BEGIN
	FOR i=0l,sz(1)-1 DO BEGIN
	k = i+(j*sz(1))
	IF (k/100 eq k/100.) THEN PRINT,'Processing Spectrum:',k
	parm = parm0

	IF (pp_flag NE 0) THEN BEGIN
		IF (pix_parm(0)) THEN parm(0) = CMB_TEMP(i,j)	
		IF (pix_parm(1)) THEN parm(1) = BB1_TEMP(i,j)	
		IF (pix_parm(2)) THEN parm(2) = OPDEP1(i,j)	
		IF (pix_parm(3)) THEN parm(3) = EMIS1(i,j)
		IF (pix_parm(4)) THEN parm(4) = BB2_TEMP(i,j)
		IF (pix_parm(5)) THEN parm(5) = OPDEP2(i,j)
		IF (pix_parm(6)) THEN parm(6) = EMIS2(i,j)
	ENDIF

	y_fit = fit_cdd(freq,spec(i,j,*),weight(i,j,*),parm,active, $
		        mask,tolerance,c2,n_i,fail0, $
                        nodust=nodust0,dust=dust0, $
		        funits=funits,res_units=res_units, $
			parm_sig=parm_sig0,f_char=f_char)

	out_fit(i,j,*) = y_fit
	out_parm(i,j,*) = parm
	parm_sig(i,j,*) = parm_sig0
	nodust(i,j,*) = nodust0
	dust(i,j,*) = dust0
	fail(i,j) = fail0
	chi2(i,j) = c2
	n_iter(i,j) = n_i
	ENDFOR
	ENDFOR

	END



ENDCASE

residual = spec - out_fit

exit:

RETURN,out_fit
END
;DISCLAIMER:
;
;This software was written at the Cosmology Data Analysis Center in
;support of the Cosmic Background Explorer (COBE) Project under NASA
;contract number NAS5-30750.
;
;This software may be used, copied, modified or redistributed so long
;as it is not sold and this disclaimer is distributed along with the
;software.  If you modify the software please indicate your
;modifications in a prominent place in the source code.  
;
;All routines are provided "as is" without any express or implied
;warranties whatsoever.  All routines are distributed without guarantee
;of support.  If errors are found in this code it is requested that you
;contact us by sending email to the address below to report the errors
;but we make no claims regarding timely fixes.  This software has been 
;used for analysis of COBE data but has not been validated and has not 
;been used to create validated data sets of any type.
;
;Please send bug reports to CGIS@ZWICKY.GSFC.NASA.GOV.


