FUNCTION fit_bed,freq,spec,weight,parm,active,mask,tol,chi2,n_iter, $
                 fail,nodust=nodust,dust=dust,funits=funits, $
                 res_units=res_units,parm_sig=parm_sig
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    FIT_BED performs the least-squared fit for FITSPEC.
;
;DESCRIPTION:
;    This IDL function fits spectra to a Bose-Einstein plus dust model.
;    The dust model includes optical depth and emissivity parameters.  
;
;    Model: fit = BE(T0,mu;freq) + op_depth1 * (freq)^emiss1 * P(T1;freq)
;
;    where BE(T,mu;freq) is a Bose-Einstein spectrum.
;
;CALLING SEQUENCE:
;    fit = fit_bed(freq,spec,weight,parm,active,mask,tol,chi2,n_iter, $
;                fail,nodust=nodust,dust=dust,funits=funits, $
;		 res_units=res_units,parm_sig=parm_sig)
;
;ARGUMENTS (I = input, O = output, [] = optional):
;    freq        I     flt arr             Spectra Frequency (in icm)
;    spec        I     flt arr             Spectra (in eplees)
;    weight      I     flt arr             Spectral Weights
;    parm        I/O   flt arr             Parameter Values
;    active      I     byte arr            Parameter Variability
;                                          [1 - vary, 0 - fixed]
;    mask        I     int arr             Masked element array
;    tol         I     flt scl             Convergence Tolerance
;    chi2        O     flt arr             Chi_squared/df
;    n_iter      O     int arr             # iterations for fit
;    fail        O     byte arr            Fit Failure Flag
;    nodust      O     flt arr             Spectrum w/dust fit subtracted
;    dust        O     flt arr             Spectrum w/cmb fit subtracted
;    funits      I     str                 Frequency Units
;    res_units   I     str                 Results Units
;    parm_sig    O     flt arr             Final Fit Parameter sigmas
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
;#
;COMMON BLOCKS:
;    None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES):
;    Iteratively calculates fit using:
;
;    (y-fit) = (df/da) * del(a) + (df/db) * del(b),
;
;    solving for del(a) and del(b) by matrix inversion.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;    None
;
;MODIFICATION HISTORY
;    Written by J.M. Gales, Applied Research Corp.   Jan 93
;    Initial Delivery  SPR 10579
;
;    SPR 11188   Fix computation of derivative of model
;                w.r.t. spectral indices.
;    08/04/93    J.M. Gales
;
;-

m2ep = 2.99792458d-7	; megajanskies to eplees
act_parm = WHERE(active GT 0)
		; get active (variable) paramater indicies

n_parm = N_ELEMENTS(act_parm)
n_freq = N_ELEMENTS(freq)
ndf = n_freq - n_parm
		; get # active parms, # frequencies, # dof

scl = EXP(TOTAL(ALOG(freq))/n_freq)
s_fq = freq / scl
parm(3) = parm(3) * (scl^parm(4))
sqrt_wt = SQRT(weight)
IF (N_ELEMENTS(mask) GT 0) THEN sqrt_wt(mask) = 0
		; scale frequencies
		; get sqrt(weights)
		; deweight desired points


wt_derv = DBLARR(n_freq,n_parm)
parm_sig = DBLARR(5)
		; allocate weighted derv and parameter sigma arrays

CASE STRUPCASE(STRMID(res_units,0,1)) OF

	'E' : 	BEGIN
		be = bose_ein(parm(0),parm(1),freq,dbe, $
				    units=funits,/mjy) * m2ep
		dbe = dbe * m2ep
		END

	'W' : 	BEGIN
		be = bose_ein(parm(0),parm(1),freq,dbe, $
				    units=funits,/wmc2)
		END

	'M' : 	BEGIN
		be = bose_ein(parm(0),parm(1),freq,dbe, $
				    units=funits,/mjy)
		END

ENDCASE


plnck = 0*freq
IF (parm(2) NE 0.0) THEN plnck = planck(parm(2),freq,dbbdt,units='i',/m)
		; calculate bose-einstein spectrum
		; calcuate planck and dbdt if dust temp not 0


ft = be + plnck * s_fq^parm(4) * parm(3)
		; calculate initial fit

chi2 = TOTAL(((spec-ft)*sqrt_wt)^2)/ndf
		; calculate initial chi_squared


n_iter = 0
REPEAT BEGIN
		; interation fit loop

n_iter = n_iter + 1

FOR i=0,n_parm-1 DO BEGIN
	CASE act_parm(i) OF
		0: derv = dbe(*,0)
		1: derv = dbe(*,1)
		2: derv = dbbdt
		3: derv = s_fq^parm(4) * plnck
		4: derv = parm(3) * s_fq^parm(4) * plnck * alog(s_fq)
	ENDCASE
	wt_derv(*,i) = derv * sqrt_wt
ENDFOR
		; calculate d(ft)/d(parm) for all active parms


wt_spec_diff = REFORM((spec-ft)*sqrt_wt,1,n_freq)
		; calculate weighted difference between spec and fit

beta = REFORM(wt_spec_diff # wt_derv)
alpha = TRANSPOSE(wt_derv) # wt_derv
		; build fit vector (beta) and array (alpha)

IF (n_parm EQ 1) THEN BEGIN
	err_mat = 1 / alpha(0,0)
	del_parm = beta(0) * err_mat
ENDIF ELSE BEGIN
	err_mat = INVERT(alpha)
	del_parm = err_mat # beta
ENDELSE
		; calculate error matrix (=inverse of alpha)
		; calculate change in parameters

FOR j=0,n_parm-1 DO parm(act_parm(j)) = parm(act_parm(j)) + del_parm(j)
		; calculate new parameters


CASE STRUPCASE(STRMID(res_units,0,1)) OF

	'E' : 	BEGIN
		be = bose_ein(parm(0),parm(1),freq,dbe, $
				    units=funits,/mjy) * m2ep
		dbe = dbe * m2ep
		END

	'W' : 	BEGIN
		be = bose_ein(parm(0),parm(1),freq,dbe, $
				    units=funits,/wmc2)
		END

	'M' : 	BEGIN
		be = bose_ein(parm(0),parm(1),freq,dbe, $
				    units=funits,/mjy)
		END

ENDCASE

IF (active(2) NE 0) THEN BEGIN

	CASE STRUPCASE(STRMID(res_units,0,1)) OF

	'E' : 	plnck = planck(parm(2),freq,dbbdt,units='i',/mjy) * m2ep
	'W' : 	plnck = planck(parm(2),freq,dbbdt,units='i',/wmc2)
	'M' : 	plnck = planck(parm(2),freq,dbbdt,units='i',/mjy)

	ENDCASE

ENDIF
		; calculate bose-einstein spectrum
		; calcuate planck and dbdt if dust temp active

ft = be + plnck * s_fq^parm(4) * parm(3)
		; calculate new fit

last_chi2 = chi2
chi2 = TOTAL(((spec-ft)*sqrt_wt)^2)/ndf
		; calculate new chi_squared

ENDREP UNTIL (ABS((chi2-last_chi2)/last_chi2) LT tol OR n_iter GE 20)
		; repeat if relative change in chi2 gt tolerance

nodust = spec - (plnck * s_fq^parm(4) * parm(3))
		; return no dust residual

dust = spec - be
		; return no cmb residual


i = WHERE(parm * [1,1,1,1,0] LT 0)
IF (i(0) NE -1 OR n_iter GE 20) THEN fail = 1B ELSE fail = 0B


FOR j=0,n_parm-1 DO parm_sig(act_parm(j)) = SQRT(err_mat(j,j))
		; calculate parameter uncertainties

parm(3) = parm(3) / (scl^parm(4))
parm_sig(3) = parm_sig(3) / (scl^parm(4))

RETURN,ft
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


