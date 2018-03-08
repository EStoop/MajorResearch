FUNCTION fit_cdd,freq,spec,weight,parm,active,mask,tol,chi2,n_iter, $
                 fail,nodust=nodust,dust=dust,funits=funits, $
		 res_units=res_units,parm_sig=parm_sig,f_char=f_char
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;    FIT_CDD performs the least-squared fit for FITSPEC.
;
;DESCRIPTION:
;    This IDL function fits spectra to a three blackbody model.
;    One blackbody typically corresponds to the cosmic background
;    while the other two can be used to fit the residual dust.  These
;    last two spectra include corresponding optical depth and emissivity
;    parameters.  
;
;    Model: fit = P(T0;freq) + op_depth1 * (freq/f_char)^emiss1 * P(T1;freq)
;                            + op_depth2 * (freq/f_char)^emiss2 * P(T2;freq)
;
;    where P(T;freq) is a Planck (Blackbody) spectrum and f_char is the
;    "characteristic" frequency that sets the frequency scale.
;
;CALLING SEQUENCE:
;    fit = fit_cdd(freq,spec,weight,parm,active,mask,tol,chi2,n_iter, $
;                fail,nodust=nodust,dust=dust,funits=funits, $
;		 res_units=res_units,parm_sig=parm_sig,f_char=f_char)
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
;    f_char      I     flt scl             Characteristic frequency
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

s_fq = freq / f_char

sqrt_wt = SQRT(weight)
IF (N_ELEMENTS(mask) GT 0) THEN sqrt_wt(mask) = 0
		; scale frequencies
		; get sqrt(weights)
		; deweight desired points


plnck = DBLARR(n_freq,3)
dbdt = DBLARR(n_freq,3)
wt_derv = DBLARR(n_freq,n_parm)
parm_sig = DBLARR(7)
		; allocate stored planck, dbdt, weighted derv 
		; and parameter uncertainties arrays

FOR i=0,2 DO BEGIN
	IF (parm(i*i) NE 0.0) THEN BEGIN
		CASE STRUPCASE(STRMID(res_units,0,1)) OF

		'E' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/mjy) * m2ep
			dbdt(*,i) = db * m2ep
			END

		'W' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/wcm2)
			dbdt(*,i) = db
			END

		'M' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/mjy)
			dbdt(*,i) = db
			END

		ENDCASE

	ENDIF
ENDFOR
		; for temp parms if temp ne 0 then calculate and
		; store planck and dbdt


ft = plnck(*,0) + plnck(*,1) * s_fq^parm(3) * parm(2) + $
                  plnck(*,2) * s_fq^parm(6) * parm(5)
		; calculate initial fit

chi2 = TOTAL(((spec-ft)*sqrt_wt)^2)/ndf
		; calculate initial chi_squared


n_iter = 0
REPEAT BEGIN
		; interation fit loop

n_iter = n_iter + 1

FOR i=0,n_parm-1 DO BEGIN
	CASE act_parm(i) OF
		0: derv = dbdt(*,0)
		1: derv = parm(2) * s_fq^parm(3) * dbdt(*,1)
		2: derv = s_fq^parm(3) * plnck(*,1)
		3: derv = parm(2) * s_fq^parm(3) * plnck(*,1) * alog(s_fq)
		4: derv = parm(5) * s_fq^parm(6) * dbdt(*,2)
		5: derv = s_fq^parm(6) * plnck(*,2)
		6: derv = parm(5) * s_fq^parm(6) * plnck(*,2) * alog(s_fq)
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

FOR i=0,2 DO BEGIN
	IF (active(i*i) NE 0) THEN BEGIN
		CASE STRUPCASE(STRMID(res_units,0,1)) OF

		'E' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/mjy) * m2ep
			dbdt(*,i) = db * m2ep
			END

		'W' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/wcm2)
			dbdt(*,i) = db
			END

		'M' : 	BEGIN
			plnck(*,i) = planck(parm(i*i),freq,db, $
					    units=funits,/mjy)
			dbdt(*,i) = db
			END

		ENDCASE
	ENDIF
ENDFOR
		; for temp parms if active parm then calculate and
		; store planck and dbdt


ft = plnck(*,0) + plnck(*,1) * s_fq^parm(3) * parm(2) + $
                  plnck(*,2) * s_fq^parm(6) * parm(5)
		; calculate new fit

last_chi2 = chi2
chi2 = TOTAL(((spec-ft)*sqrt_wt)^2)/ndf
		; calculate new chi_squared

ENDREP UNTIL (ABS((chi2-last_chi2)/last_chi2) LT tol OR n_iter GE 20)
		; repeat if relative change in chi2 gt tolerance


nodust = spec - (plnck(*,1) * s_fq^parm(3) * parm(2) + $
                 plnck(*,2) * s_fq^parm(6) * parm(5))
		; return no dust residual

dust = spec - plnck(*,0)
		; return no cmb residual


i = WHERE(parm * [1,1,1,0,1,1,0] LT 0)
IF (i(0) NE -1 OR n_iter GE 20) THEN fail = 1B ELSE fail = 0B
		; if any parm becomes negative or # inter > 20 then set fail

FOR j=0,n_parm-1 DO parm_sig(act_parm(j)) = SQRT(err_mat(j,j))
		; calculate parameter uncertainties


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


