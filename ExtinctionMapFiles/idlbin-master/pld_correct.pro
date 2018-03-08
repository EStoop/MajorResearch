FUNCTION pld_correct,pixgrid,flux,coeffs,sigma_coeffs,sigma_flux=sigma_flux,sigma_correct=sigma_correct,fpredict=fpredict,$
                     sigma_fpredict=sigma_fpredict,NOCONSTANT=noconstant,ORDER=order
  ;; Build the Design matrix
  a = pld_design_matrix(pixgrid,ndat,ncoef,NOCONSTANT=noconstant,ORDER=order)
  
  fpredict = DBLARR(ndat)
  sigma_fpredict = fpredict
  IF N_ELEMENTS(sigma_coeffs) NE ncoef THEN SIGMA_COEFFS = REPLICATE(1d0,ncoef)
  IF N_ELEMENTS(sigma_flux) NE ndat THEN SIGMA_FLUX = REPLICATE(1d0,ndat)
  FOR j = 0,ndat-1 DO BEGIN
    fpredict[j] = TOTAL(A[*,j]*coeffs,/DOUBLE,/NAN)
    sigma_fpredict[j] = SQRT( TOTAL( (A[*,j]*sigma_coeffs)^2,/DOUBLE,/NAN) )
  ENDFOR
  fcorrect = flux/fpredict
  sigma_correct = fcorrect * SQRT( (sigma_fpredict/fpredict)^2 + (sigma_flux/flux)^2 )
  RETURN,fcorrect
END
