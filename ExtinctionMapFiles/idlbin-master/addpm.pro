;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;+
;
;*NAME: ADDPM
;http://www.astro.washington.edu/docs/idl/cgi-bin/getpro/library43.html?ADDPM
;*PURPOSE: Correctly adds proper motions to input positions.
;
;*CALLING SEQUENCE:
;	ADDPM,ra,dec,pm_ra,pm_dec,v_r,parallax,delta_t,OUTRA,OUTDEC
;*PARAMETERS:
;	INPUT:
;	ra	- (REQ) - (0) - (R,D) - RA in decimel degrees.
;	dec	- (REQ) - (0) - (R,D) - Dec in decimel degrees.
;	pm_ra	- (REQ) - (0) - (R,D) - RA proper motion in arcsec/century.
;	pm_dec	- (REQ) - (0) - (R,D) - DEC proper motion in arcsec/century.
;	v_r	- (REQ) - (0) - (R,D) - Radial velocity in
;	parallax- (REQ) - (0) - (R,D) - Parallax angle in 
;	delta_t	- (REQ) - (0) - (R,D) - Time span in years.
;	OUTPUT:
;   	OUTRA	- (REQ) - (0) - (R,D) - Corrected RA in decimel degrees.
;	OUTDEC	- (REQ) - (0) - (R,D) - Corrected DEC in decimel degrees.
;*EXAMPLES:
;	For Epsilon Ind given in Taff:
;	ADDPM,329.887726,-56.992683,3.93999,-2.5554,-40.4,0.285,50.0,ra,dec
;*RESTRICTIONS:
;	RA and DEC are assumed to be in Degrees.
;	PM_RA and PM_DEC are assumed to be in 
;	V-R assumed to be in
;	PARALLAX assumed to be in 
;	DELTA_T	assumed to be in decimel years.
;*NOTES:
;	This procedure is used with TA_POS.PRO but can be used independantly.
;*PROCEDURE:
;	Taken from Taff, "Computational Spherical Astronomy".  Uses third
;	order polynomial solution (eq. 3.34, page 39).
;*MODIFICATION HISTORY:
;	Ver 1.0	- 12/xx/90 - A. Warnock   - ST Systems Corp.
;	Ver 2.0 - 02/04/91 - J. Blackwell - GSFC - Modified to conform with
;	                                           GHRS DAF standards. 
;       Mar 27 1991      JKF/ACC    - moved to GHRS DAF (IDL Version 2)
;-
;-------------------------------------------------------------------------------
pro addpm,ra,dec,pm_ra,pm_dec,v_r,parallax,delta_t,OUTRA,OUTDEC
;
deg2rad  = 3.14159265335 / 180.0                   ; conversion factor
sec2rad  = deg2rad / 3600.0                        ; conversion factor
r_ra     = ra  * deg2rad                           ; radians
r_dec    = dec * deg2rad                           ; radians
nu       = 1.0227e-4 * parallax * v_r              ; radians/century
t        = delta_t / 100.0                         ; in centuries
mu       = pm_ra  * 100.0 * sec2rad / cos(r_dec)   ; radians/century
mu_p     = pm_dec * 100.0 * sec2rad                ; radians/century
omega_2  = (mu*cos(r_dec))^2 + mu_p^2              ; radians/century
fact1    = nu - mu_p * tan(r_dec)
;
; Compute the coefficients of t for right ascension
;
r1 = mu
r2 = mu * fact1
r3 = mu * ((3.0 * fact1^2) - mu^2) / 3.0
;
; Compute the coefficients of t for declination
;
d1 = mu_p
d2 = ( (mu^2 * sin(r_dec) * cos(r_dec)) + (2 * mu_p * nu)) / 2.0
d3 = (omega_2 * ((4 * nu) - (mu_p * tan(r_dec))) * tan(r_dec)       $
                 - mu_p * (mu^2 + omega_2)                          $
                 + mu_p * (2 * nu - mu_p * tan(r_dec))^2) /4.0      $
                 - mu_p^3 / 12.0
;
; Correct the positions via polynomials
;
outra  = r_ra  + (r1 * t) - (r2 * t^2) + (r3 * t^3)
outdec = r_dec + (d1 * t) - (d2 * t^2) + (d3 * t^3)
;
; Convert results back to decimal degrees
;
outra  = outra / deg2rad
outdec = outdec / deg2rad
;
; Back to the caller
;
return
end
