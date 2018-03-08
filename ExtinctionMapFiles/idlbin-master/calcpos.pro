;+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
;+
;
;*NAME: CALCPOS
;
;*PURPOSE: Updates an input position for precession and proper motion.
;
;*CALLING SEQUENCE:
;	CALCPOS,ra_deg,dec_deg,pm_ra,pm_dec,v_r,parallax,in_equinox,
;	out_equinox,pm_equinox,obs_epoch,pm_epoch,pm_yn,prec_pos,OUTRA,OUTDEC
;*PARAMETERS:
;	INPUT:
;	ra_dec	   - (REQ) - (0) - (R,D) - Input RA, decimal degrees.
;	dec_deg	   - (REQ) - (0) - (R,D) - Input DEC, decimal degrees.
;	pm_ra	   - (REQ) - (0) - (R,D) - Proper motion in RA, arcseconds/yr.
;	pm_dec	   - (REQ) - (0) - (R,D) - Proper motion in DEC, arcseconds/yr.
;	v_r	   - (REQ) - (0) - (R,D) - Radial velocity, km/sec.
;	parallax   - (REQ) - (0) - (R,D) - Parallax, arcseconds.
;	in_equinox - (REQ) - (0) - (R,D) - Starting coordinate equinox
;		                           (e.g., 1950.0 coords).
;	out_equinox- (REQ) - (0) - (R,D) - Ending coordinate equinox
;		                           (e.g., 2000.0 coords).
;	pm_equinox - (REQ) - (0) - (R,D) - Equinox for proper motions
;		                           (e.g., 1950.0 coords).
;	obs_epoch  - (REQ) - (0) - (R,D) - Date of observation (e.g., 1990.9).
;	pm_epoch   - (REQ) - (0) - (R,D) - Date to update proper motions
;		                           from (e.g., 1958.5).
;	pm_yn      - (REQ) - (0) - (S)   - 'Y' if proper motion corrections
;		                           are to be done.
;	prec_pos   - (REQ) - (0) - (S)   - 'Y' if precession is to be performed.
;	OUTPUT:
;	outra	   - (REQ) - (0) - (R,D) - Resulting RA, decimal degrees.
;	outdec	   - (REQ) - (0) - (R,D) - Resulting DEC, decimal degrees.
;*SUBROUTINES CALLED:
;	ZPRECESS (from UIT ASTRO Ver 1.0 User's library PRECESS)
;*EXAMPLES:
;	See CHECKPOS.PRO
;*RESTRICTIONS:
;	Proper motion equinox must be either J2000.0 or the same as
;	the input equinox.
;*PROCEDURE:
;	Precesses position, if necessary.  Proper motion correction is
;	applied either before or after precession, depending on the
;	equinox of the proper motions.
;*MODIFICATION HISTORY:
;	Ver 1.0	- 12/xx/90 - A. Warnock   - ST Systems Corp.
;	Ver 2.0	- 02/11/91 - J. Blackwell - GSFC - Modified to conform with
;	                                    GHRS DAF standards.
;       Mar 28 1991      JKF/ACC    - moved to GHRS DAF (IDL Version 2)
;-
;-------------------------------------------------------------------------------
pro calcpos, ra_deg, dec_deg, pm_ra, pm_dec, v_r, parallax,      $
            in_equinox, out_equinox, pm_equinox,                 $
            obs_epoch, pm_epoch, pm_yn, prec_pos,                $
            outra,outdec
if n_params(0) eq 0 then begin
   print,'  Calling Sequence: CALCPOS, ra_deg,dec_deg,pm_ra,pm_dec,v_r,'
   print,'parallax,in_equinox,out_equinox,pm_equinox,obs_epoch,pm_epoch,'
   print,'pm_yn,prec_pos,OUTRA,OUTDEC'
   retall
endif
;
; Save the input values
;
pm_flag = pm_yn
tmpra   = ra_deg
tmpdec  = dec_deg
outra   = ra_deg
outdec  = dec_deg
in_eq   = in_equinox
;
; See if we need to correct pm's first, and do so
;
if (pm_flag eq 'Y') then begin
   if (pm_equinox eq in_equinox) then begin
   ;
   ; We don't need to precess the input if the equinoxes are the same
   ;
   ;  outra  = outra  + (pm_ra  / 3600.0)*(obs_epoch - pm_epoch)
   ;  outdec = outdec + (pm_dec / 3600.0)*(obs_epoch - pm_epoch)
      delta_t = obs_epoch - pm_epoch
      addpm,tmpra,tmpdec,pm_ra,pm_dec,v_r,parallax,delta_t,outra,outdec
      pm_flag = 'N'
   endif else begin
      if (pm_equinox ne out_equinox) then begin
      ;
      ; If the equinoxes differ, and the equinox for the proper motions
      ; isn't the same as the output equinox, precess the input equinox
      ; to the proper motion equinox so we can add them together.
      ;
      ; Otherwise, we can wait until after the "big" precession to add
      ; in the proper motions
      ;
         zprecess,tmpra,tmpdec,in_eq,pm_equinox
         delta_t = obs_epoch - pm_epoch
         addpm,tmpra,tmpdec,pm_ra,pm_dec,v_r,parallax,delta_t,outra,outdec
         if (prec_pos eq 'Y') then begin
         ;
         ; If we're precessing anyway, just make it from the current
         ; equinox (ie, proper motion equinox)
         ;
            in_eq  = pm_equinox
         endif else begin
         ;
         ; But if no precession is to be done, we'd better put things
         ; back to where they started
         ;
            zprecess,outra,outdec,pm_equinox,in_eq
         endelse
      endif
   endelse
endif
;
; Do precession, if required
;
if ((prec_pos eq 'Y') and (in_eq ne out_equinox)) then begin
   zprecess,outra,outdec,in_eq,out_equinox
endif
;
; See if we needed to correct pm's after precession, and do so
;
if ((pm_flag eq 'Y') and (pm_equinox eq out_equinox)) then begin
   delta_t = obs_epoch - pm_epoch
   addpm,outra,outdec,pm_ra,pm_dec,v_r,parallax,delta_t,tmpra,tmpdec
   outra = tmpra
   outdec = tmpdec
endif
;
; Head back out
;
return
end
