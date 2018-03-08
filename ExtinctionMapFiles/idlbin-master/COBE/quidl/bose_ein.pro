function bose_ein, Ttemp, mu, nu_or_lambda, derivs, units=units, $
                   mjy=mjy, wcm2=wcm2
;+NAME/ONE LINE DESCRIPTION OF ROUTINE:
;     BOSE_EIN returns a Bose-Einstein spectrum with a chemical potential.
;
;DESCRIPTION:  
;    IDL function to return the spectral radiance of a Bose-Einstein
;    distribution, i.e., a Planck-like curve with a chemical potential
;    term mu.  Output is in units of either MJy/steradian (I_nu)
;    or watts/cm^2/steradian (nu_I_nu).  Inputs are the equivalent
;    blackbody temperature, the chemical potential parameter mu, and 
;    either frequency (in icm or GHz) or wavelength (in microns).  The
;    routine also optionally returns the derivative with respect to 
;    both temperature (in units of MJy/sr/K or W/cm^2/sr/K) and mu.
;    Note that for mu=0 the output spectrum is identical to the Planck
;    function with the same temperature.
;
;CALLING SEQUENCE:  
;     RESULT = BOSE_EIN (temperature, mu, nu_or_lambda  $
;                        [,derivs] [,UNITS=units], $
;                        [,/MJY] [,/WCM2])
;
;ARGUMENTS (I = input, O = output, [] = optional):
;     RESULT        O   flt [arr]  Spectral radiance at each wavelength. 
;                                  Units: W/cm^2/sr/K if /WCM2 specified
;                                         MJy/sr      if /MJY specfied
;     TEMPERATURE   I   flt        Temperature of blackbody, in K.
;     MU            I   flt        B-E chemical potential, unitless
;     NU_OR_LAMBDA  I   flt [arr]  Frequency or wavelength at which to 
;                                  calculate spectrum. Units are as 
;                                  specified with UNITS keyword.
;     DERIVS       [O]  flt [arr,2]Derivative of spectrum with respect to 
;                                  temperature [DERIVS(*,0)] and chemical
;                                  potential [DERIVS(*,1)].
;     UNITS        [I]  str        'Microns', 'icm', or 'GHz' to 
;                                  identify units of NU_OR_LAMBDA. Only 
;                                  first character is required.  If 
;                                  left out, default is 'microns'.
;     /MJY          I   key        Sets output units to MJy/sr
;     /WCM2         I   key        Sets output units to W/cm^2/sr
;
;WARNINGS:
;     1.  One of /MJY or /WCM2 MUST be specified.  
;     2.  Routine gives incorrect results for T < 1 microKelvin and
;            wavelengths shortward of 1.e-10 microns.  (So sue me).
;     3.  We use a dimensionless mu.  It is NOT in units of kT.
;
;EXAMPLE:
;     To produce a 35 K spectrum in MJy/sr at 2, 4, 6, 8, 10 microns
;     with a chemical potential 0.001:
;
;       wavelength = 2. + 2.*findgen(5)
;       temp = 35.
;       mu = 0.001
;       spectrum = BOSE_EIN (temp, mu, wavelength, units='micron', /mjy)
;
;     One could also get back the derivatives into the variable 
;     "derivs" including it in the call:
;
;       spectrum = BOSE_EIN (temp, mu, wavelength, derivs, $
;                            units='m', /mjy)
;
;     Note that dB/dT = derivs(*,0) and dB/dmu= derivs(*,1).
;#
;COMMON BLOCKS:
;     None
;
;PROCEDURE (AND OTHER PROGRAMMING NOTES): 
;     Identifies units using the UNITS keyword, then converts the 
;     supplied independent variable into microns to evaluate the 
;     Planck function.  Uses Rayleigh-Jeans and Wien approximations 
;     for the low- and high-frequency end, respectively.  Reference: 
;     Allen, Astrophysical Quantities, for the Planck formula.
;
;PERTINENT ALGORITHMS, LIBRARY CALLS, ETC.:
;     None
;  
;MODIFICATION HISTORY:
;    Written by Rich Isaacman, General Sciences Corp.  28 Dec 1992 SPR 10397
;
;.TITLE
; Routine BOSE_EIN
;-
;
; Check on input parameters
;
on_error, 2
if n_elements(nu_or_lambda) lt 1 or n_elements(Ttemp) lt 1 or $
     n_elements(mu) lt 1 then $
     message,'CALLING SEQUENCE: spectrum = bose_ein (temp,mu,wavelength,units=<units>)'
if not keyword_set(mjy) and not keyword_set(wcm2) then $
     message, 'Either /MJy or /Wcm2 must be specified!'
if keyword_set(mjy) and keyword_set(wcm2) then $
     message, 'Only one of /MJy or /Wcm2 may be specified!'
;
makederiv = n_params() gt 3         ; see whether derivatives are requested
;
if n_elements(units) lt 1 then begin
   units = 'm'
   message, /continue, "Wavelength units are assumed to be microns."
endif
T = float(ttemp) > 1.e-06             	;force temperature to be real*4
;
; Define some necessary constants
;
c    = 299792.458d0                      ;lightspeed per Physics Today 8/90
hck  = 14387.69d0                        ;h*c/k       "      "     "
thcc = 1.1910439d4                       ;2*h/c^2  Physics Today, Aug.1990 
;
; Convert nu_or_lambda into lambda (in microns) depending on specified units
;
units0 = strupcase (strmid (units,0,1))
case units0 of
   'M':   lambda = nu_or_lambda > 1.e-10                ; microns
   'I':   lambda = 1.e04 / (nu_or_lambda > 1.e-10)      ; icm, avoiding nu=0.
   'G':   lambda = c / (nu_or_lambda > 1.e-10)          ; GHz, avoiding nu=0.
   else:  message, "Units must be specified as 'microns', 'icm', or 'GHz'"
endcase
; 
;  Variable fhz is a scale factor used to go from units of nu_I_nu units to 
;  MJy/sr if the keyword is set.  In that case, its value is
;  frequency in Hz * w/cm^2/Hz ==> MJy
;
if keyword_set(mjy) then fhz = c/lambda * 1.e-15         
if keyword_set(wcm2) then fhz = 1. + fltarr(n_elements(lambda))       
;
;  Introduce dimensionless variable chi, used to check whether we are on 
;  Wien or Rayleigh Jeans tails
;
chi = (hck / lambda / T + mu) > 1.e-10
val = fltarr(n_elements(chi))
if makederiv then derivs = fltarr(n_elements(chi),2)
;
;  Start on Rayleigh Jeans side
;
rj = where (chi lt 0.001)
if rj(0) ne -1 then begin
    val(rj) = thcc / chi(rj) / lambda(rj)^4 / fhz(rj)
    if makederiv then begin
       derivs(rj,0) = val(rj) / T
       derivs(rj,1) = -val(rj) / chi(rj)
    endif
endif
;
;  Now do nonapproximate part
;
exact = where (chi ge 0.001 and chi le 50)
if exact(0) ne -1 then begin
    chi_ex = chi(exact)
    val(exact) = thcc / lambda(exact)^4 / (exp(chi_ex) - 1.) / fhz(exact)
    if makederiv then begin
       derivs(exact,0) = val(exact) * chi_ex / T / (1. - exp(-chi_ex))
       derivs(exact,1) = -val(exact) / (1. - exp(-chi_ex))
    endif
endif
;
;  ...and finally the Wien tail
;
wien = where (chi gt 50.)
if wien(0) ne -1 then begin
    chi_wn = chi(wien)
    val(wien) = thcc / (lambda(wien) > 1.e-6)^4 * exp(-chi_wn) / fhz(wien)
    if makederiv then begin
       derivs(wien,0) = val(wien) * chi_wn / (T > 1.e-6) / (1. - exp(-chi_wn))
       derivs(wien,1) = -val(wien) / (1. - exp(-chi_wn))
    endif
endif
;
;  Redimension to a scalar if only 1 wavelength supplied, then return
;
if n_elements(nu_or_lambda) eq 1 then begin
    val = val(0)
    derivs = reform(derivs)
endif
return, val
;
end
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


