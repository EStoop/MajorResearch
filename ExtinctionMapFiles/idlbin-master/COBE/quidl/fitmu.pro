FUNCTION fitmu,freq,spec,weight,parm0=parm0,active=active, $
         cmb_temp=cmb_temp,cmb_mu=cmb_mu,bb_temp=bb_temp, $
         opdep=opdep,emis=emis,mask=mask,out_parm=out_parm, $
	 parm_sig=parm_sig,dust=dust,nodust=nodust,residual=residual, $
         fail=fail,chi2=chi2,n_iter=n_iter,tolerance=tolerance, $
         funits=funits,eplees=eplees,wcm2=wcm2,mjy=mjy
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    FITMU is a Bose_Einstein / Dust Residual fitting facility.
;
;DESCRIPTION:
;    This IDL function fits spectra to a Bose_Einstein spectrum plus a
;    residual Planck (blackbody) spectrum.  The blackbody includes
;    corresponding optical depth and emissivity parameters.  Therefore
;    there are five fitting parameter in all.  In most cases it is not
;    possible to fit all five at once, therefore parameters can be
;    designated as fixed.  The initial values of the parameters are
;    supplied by the user (either using a command line keyword or
;    interactively) as well as their variability (fixed or variable).
;    It is also possible to designate  fixed value parameters that are
;    different for each input spectrum.  Particular ranges of the spectra
;    can be masked if desired.  This is useful for avoiding spectral
;    lines.
;
;    Model: fit = BE(T0,mu;freq) + op_depth * (freq)^emiss * P(T1;freq)
;
;    where BE(T,mu;freq) is a Bose-Einstein spectrum and
;          P(T;freq) is a Planck (Blackbody) spectrum.
;
;CALLING SEQUENCE:
;    fit = fitmu(freq,spec,[weight],[parm0=parm0],[active=active], $
;      [cmb_temp=cmb_temp],[cmb_mu=cmb_mu],[bb_temp=bb_temp], $
;      [opdep=opdep],[emis=emis],[mask=mask],[tolerance=tolerance], $
;      [out_parm=out_parm],[parm_sig=parm_sig],[residual=residual], $
;      [fail=fail],[chi2=chi2],[n_iter=n_iter],[funits=funits], $
;      [/eplees],[/wcm2],[/mjy])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    freq        I     flt arr             Spectra Frequency
;    spec        I     flt arr             Spectra
;    [weight]    I     flt arr             Spectral Weights
;    [parm0]     I     flt arr             Initial Parameter Values
;    [active]    I     byte arr            Parameter Variability
;                                          [1 - vary, 0 - fixed]
;    [cmb_temp]  I     flt arr             CMP Temp value array
;    [cmb_mu]    I     flt arr             CMB chemical potential
;    [bb_temp]   I     flt arr             BB Temp value array
;    [opdep]     I     flt arr             Optical Depth value array
;    [emis]      I     flt arr             Emissivity value array
;    [mask]      I     int arr             Masked element array
;    [tolerance] I     flt scl             Convergence Tolerance
;    [out_parm]  O     flt arr             Final Fit Parameter values
;    [parm_sig]  O     flt arr             Final Fit Parameter sigmas
;    [residual]  O     flt arr             Fit Residuals
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
;
;    The spectrum and weight arrays must have the same number of elements
;    and the same number of dimensions, which can be one (a single
;    spectrum), two (an array of spectra), or three (an array of spectra
;    from a skymap).  The last dimension should correspond to the
;    frequency array.
;
;    The order of the values in the initial parameter array, 'parm0',
;    is:  Bose-Einstein Temperature (in K), Bose-Einstein Chemical
;    Potential, Residual Blackbody Temperature, Residual Optical Depth,
;    Residual Emissivity.
;
;    The cmb_temp, cmb_mu, bb_temp, opdep, and emis1 input arrays are
;    used to define fixed parameter values for the individual spectra.
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
; To fit a set of FIRAS low frequency spectra to a Bose_Einstein
; spectrum plus a residual dust spectrum using the dust temperatures
; derived from a previous fit to the high frequency spectra as fixed
; parameters:
;
; fit_be = fitmu(freq_lo,spec_lo,weight_lo, $
;                parm0=[2.72,0.001,0,2e-7,1.5], $
;                active=[1,1,0,1,0], bb_temp=hi_temp, $
;                out_parm=parm_be,fail=fail_be)
;
; The initial values of the parameters are 2.72K for the BE temperature,
; 0.001 for the BE chemical potential, 2e-7 for the dust optical depth,
; and 1.5 for the dust emissivity.  The BE temperature and chemical
; potential and the dust optical depth are allowed to vary.  The dust
; emissivity is fixed (active(4) = 0).  The dust temperature is read
; in from the 'hi_temp' array.
;
; If either the 'parm0' or 'active' keywords were missing, then the
; user would be prompted for their values.
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
;    Loop through spectra calling fit_bed to determine fit.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    Subroutines called: fit_bed
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jan 93
;    Initial Delivery  SPR #10579
;
;.TITLE
;Routine FITMU
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


default = [2.7,0,20,2e-7,1.5]
dflt_str = ['<2.7>','<0>','<20>','<2e-7>','<1.5>']
parm_title=['Cosmic Background Temperature ', $
	    'Cosmic Background Chemical Potential ', $
            'Resisual Black Body Temperature ', $
            'Residual Black Body Optical Depth ', $
            'Residual Black Body Emissivity ']

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

IF (KEYWORD_SET(cmb_mu)) THEN BEGIN
	pix_parm(1) = 1
	sz_fix = SIZE(cmb_mu)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and CMB_MU Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF

IF (KEYWORD_SET(bb_temp)) THEN BEGIN
	pix_parm(2) = 1
	sz_fix = SIZE(bb_temp)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and BB_TEMP Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF

IF (KEYWORD_SET(opdep)) THEN BEGIN
	pix_parm(3) = 1
	sz_fix = SIZE(opdep)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and OPDEP Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF

IF (KEYWORD_SET(emis)) THEN BEGIN
	pix_parm(4) = 1
	sz_fix = SIZE(emis)
	IF (TOTAL(sz(1:n_dim-1)) NE TOTAL(sz_fix(1:n_dim-1))) THEN BEGIN
		str = 'Spectrum and EMIS Arrays not compatible'
		MESSAGE,str,/CONT
		out_fit = -1
		GOTO, exit
	ENDIF
ENDIF


IF (TOTAL(pix_parm) NE 0) THEN pp_flag = 1 ELSE pp_flag = 0

IF (KEYWORD_SET(parm0) EQ 0) THEN BEGIN
	parm0 = FLTARR(5)
	FOR i=0,4 DO BEGIN
		IF (pix_parm(i) EQ 1) THEN GOTO,lbl1
		IF ((i EQ 3 OR i EQ 4) AND (parm0(2) EQ 0) $
		AND (pix_parm(2) EQ 0)) $
		THEN GOTO,lbl1
		READ,parm_title(i)+dflt_str(i)+': ',inp
		IF (inp EQ '') THEN parm0(i) = default(i) $
		ELSE parm0(i) = FLOAT(inp)
	lbl1:
	ENDFOR
ENDIF

PRINT,' '
IF (KEYWORD_SET(active) EQ 0) THEN BEGIN
	active = BYTARR(5)
	FOR i=0,4 DO BEGIN
		IF (pix_parm(i) EQ 1) THEN GOTO,lbl2
		IF ((i GE 2 AND i LE 4) AND (parm0(2) EQ 0) $
		AND (pix_parm(1) EQ 0)) $
		THEN GOTO,lbl2
		READ,'Vary ' + parm_title(i) + '<Yes>? ',inp
		IF (inp EQ '' OR STRUPCASE(STRMID(inp,0,1)) EQ 'Y') $
		THEN active(i) = 1
	lbl2:
	ENDFOR
ENDIF

PRINT,' '
PRINT,'Initial Paramater Values'
PRINT,'------------------------'
FOR i=0,4 DO BEGIN
	IF (pix_parm(i) EQ 0) THEN str = STRING(parm0(i)) $
	ELSE str = 'Read From Input Array'
	PRINT,FORMAT='(a,t40,a,t55,a)',parm_title(i)+': ',str, $
				       STRING(active(i))
ENDFOR
PRINT,' '

IF (KEYWORD_SET(tolerance) EQ 0) THEN tolerance = 1.e-4

CASE n_dim OF

   1:	BEGIN

	out_parm = FLTARR(5)
	fail = BYTARR(1)
	chi2 = FLTARR(1)
	n_iter = INTARR(1)
	out_fit = FLTARR(n_freq)
	parm = parm0

	y_fit = fit_bed(freq,spec,weight,parm,active,mask, $
		        tolerance,c2,n_i,fail0,nodust=nodust, $
		        dust=dust,funits=funits,res_units=res_units, $
			parm_sig=parm_sig)

	out_fit = y_fit
	out_parm = parm
	fail(0) = fail0
	chi2(0) = c2
	n_iter(0) = n_i

	END


   2:	BEGIN

	out_parm = FLTARR(sz(1),5)
	parm_sig = FLTARR(sz(1),5)
	fail = BYTARR(sz(1))
	chi2 = FLTARR(sz(1))
	n_iter = INTARR(sz(1))
	out_fit = FLTARR(sz(1),sz(2))
	nodust = FLTARR(sz(1),sz(2))
	dust = FLTARR(sz(1),sz(2))

	FOR i=0l,sz(1)-1 DO BEGIN
	IF (i/100 EQ i/100.) THEN PRINT,'Processing Spectrum:',i
	parm = parm0

	IF (pp_flag NE 0) THEN BEGIN
		IF (pix_parm(0)) THEN parm(0) = cmb_temp(i)
		IF (pix_parm(1)) THEN parm(1) = cmb_mu(i)
		IF (pix_parm(2)) THEN parm(2) = bb_temp(i)	
		IF (pix_parm(3)) THEN parm(3) = opdep(i)	
		IF (pix_parm(4)) THEN parm(4) = emis(i)
	ENDIF

	y_fit = fit_bed(freq,spec(i,*),weight(i,*),parm,active, $
		        mask,tolerance,c2,n_i,fail0, $
                        nodust=nodust0,dust=dust0, $
		        funits=funits,res_units=res_units, $
			parm_sig=parm_sig0)

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

	out_parm = FLTARR(sz(1),sz(2),5)
	parm_sig = FLTARR(sz(1),sz(2),5)
	fail = BYTARR(sz(1),sz(2))
	chi2 = FLTARR(sz(1),sz(2))
	n_iter = INTARR(sz(1),sz(2))
	out_fit = FLTARR(sz(1),sz(2),n_freq)
	nodust = FLTARR(sz(1),sz(2),n_freq)
	dust = FLTARR(sz(1),sz(2),n_freq)

	FOR j=0l,sz(2)-1 DO BEGIN
	FOR i=0l,sz(1)-1 DO BEGIN
        k=i+(j*sz(1))
	IF (k/100 EQ k/100.) THEN PRINT,'Processing Spectrum:',k
	parm = parm0

	IF (pp_flag NE 0) THEN BEGIN
		IF (pix_parm(0)) THEN parm(0) = cmb_temp(i,j)
		IF (pix_parm(1)) THEN parm(1) = cmb_mu(i,j)
		IF (pix_parm(2)) THEN parm(2) = bb_temp(i,j)	
		IF (pix_parm(3)) THEN parm(3) = opdep(i,j)	
		IF (pix_parm(4)) THEN parm(4) = emis(i,j)
	ENDIF

	y_fit = fit_bed(freq,spec(i,j,*),weight(i,j,*),parm,active, $
		        mask,tolerance,c2,n_i,fail0, $
                        nodust=nodust0,dust=dust0, $
		        funits=funits,res_units=res_units, $
			parm_sig=parm_sig0)

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


