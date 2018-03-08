FUNCTION correct_pld,pixgrid,flux,coeffs,sigma_coeffs,sigma_flux=sigma_flux,sigma_correct=sigma_correct,fpredict=fpredict,fpredict_unc=fpredict_unc,$
                     NOCONSTANT=noconstant
  ;; Build the Design matrix
  a = pld_design_matrix(pixgrid,ndat,ncoef,NOCONSTANT=noconstant)
  
  fpredict = DBLARR(ndat)
  fpredict_unc = fpredict
  IF N_ELEMENTS(sigma_coeffs) NE ncoef THEN SIGMA_COEFFS = REPLICATE(1d0,ncoef)
  IF N_ELEMENTS(sigma_flux) NE ndat THEN SIGMA_FLUX = REPLICATE(1d0,ndat)
  FOR j = 0,ndat-1 DO BEGIN
    fpredict[j] = TOTAL(A[*,j]*coeffs,/DOUBLE,/NAN)
    fpredict_unc[j] = SQRT( TOTAL( (A[*,j]*sigma_coeffs)^2,/DOUBLE,/NAN) )
  ENDFOR
  fcorrect = flux/fpredict
  sigma_correct = fcorrect * SQRT( (fpredict_unc/fpredict)^2 + (sigma_flux/flux)^2 )

  RETURN,fcorrect
END
