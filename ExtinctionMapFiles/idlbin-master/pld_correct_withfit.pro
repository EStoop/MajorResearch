FUNCTION pld_correct_withfit,pixgrid,flux,sigma_flux,SCALE=scale,REQKEY=reqkey,$
         CORR_UNC=corr_unc,IFIT=ifit,ORDER=order,_EXTRA=extra
;;
;;  corrected_flux = PLD_CORRECT_WITHFIT(pixgrid,flux,sigma_flux)
;;
;;      INPUTS
;;            PIXGRID:  ngrid x ngrid x ndat image cube giving the BCD sub-image of the data
;;                      ngrid can be 3 or 5.  If your data is in the form ngrid^2 x ndat, change it to
;;                      ngrid x ngrid x ndat using REFORM. 
;;                      Each subimage must correspond to the SAME pixels.
;;            FLUX: ndat-element array of aperture fluxes
;;            SIGMA_FLUX: ndat-element array of aperture flux uncertainties
;;        
;;     KEYWORDS
;;            CORR_UNC: ndat-element vector of approximate uncertainties in the corrected data.
;;             SCALE:  Scaling factor to divide into the data before fitting.  This can be used to
;;                     pre-normalize the data to the same scale as another dataset to be fit later.
;;             REQKEY: ndat-element string array giving the reqkeys corresponding to each data point
;;                     If entered, then separate fits and corrections will be carried out on each 
;;                     unique reqkey.  
;;             IFIT:  (nfit < ndat) -element vector giving the indices in the input data to be  
;;                     used in the fit.  All data will be corrected, if possible.  If not set then all data,
;;                     with the exception of NaN's, will be used in the fit.
;;             ORDER:  2 or 1 (default).  Gives the upper order of terms to keep in the pixel
;;                     correlation equation.  
;;


ndat_tot = N_ELEMENTS(flux)
IF N_ELEMENTS(ifit) EQ 0 THEN ifit = WHERE(finite(flux) AND finite(sigma_flux))
IF N_ELEMENTS(REQKEY) EQ ndat_tot THEN BEGIN
  iuniq = UNIQ(reqkey,SORT(reqkey))
  nreq = N_ELEMENTS(iuniq)
  rk = reqkey[iuniq]
  irk = LIST()
  irk_fit = LIST()
  FOR i = 0,nreq-1 DO BEGIN
    irk.ADD,WHERE(reqkey EQ rk[i])
    irk_fit.ADD,WHERE(reqkey[ifit] EQ rk[i])
  ENDFOR
ENDIF ELSE BEGIN
  nreq = 1
  irk = LIST(LINDGEN(ndat_tot))
  irk_fit = LIST(LINDGEN(ndat_tot))
  rk = 0
ENDELSE


IF N_ELEMENTS(sigma_flux) LT ndat_tot THEN sigma_flux = REPLICATE(1d0,ndat_tot)
;coeffs = DBLARR(10,nreq)
fcorrect = MAKE_ARRAY(ndat_tot,VALUE=!VALUES.D_NAN)
corr_unc = MAKE_ARRAY(ndat_tot,VALUE=!VALUES.D_NAN)
FOR j = 0,nreq-1 DO BEGIN
   ON_ERROR, 2
   IF N_ELEMENTS(irk_fit[j]) EQ 0 THEN BEGIN
      PRINT,'*** PLD_CORRECT_WITHFIT:  NO FITTING DATA FOR AOR '+STRNG(rk[j])+'.  VALUES WILL RETURN NaN.****'
      CONTINUE
   ENDIF
   cj = pld_fit(pixgrid[*,*,ifit[irk_fit[j]]],flux[ifit[irk_fit[j]]],sigma_flux[ifit[irk_fit[j]]],sigma_coeffs=sigma_cj,SCALE=scale,/NOCONSTANT,ORDER=order,_EXTRA=extra)
   fcj = pld_correct(pixgrid[*,*,irk[j]],flux[irk[j]],cj,sigma_cj,sigma_flux=sigma_flux[irk[j]],sigma_correct=cuncj,/NOCONSTANT,ORDER=order,_EXTRA=extra)
   fcorrect[irk[j]] = fcj
   corr_unc[irk[j]] = cuncj
ENDFOR

RETURN,fcorrect
END
