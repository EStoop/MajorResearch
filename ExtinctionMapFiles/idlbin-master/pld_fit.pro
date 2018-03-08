FUNCTION pld_fit,pixgrid,flux,sigma_flux,fpredict=fpredict,fcorrect=fcorrect,sigma_coeffs=sigma_coeffs,SCALE=scale,REQKEY=reqkey,NOCONSTANT=noconstant,ORDER=order

  ;; Build the matrices for linear least squares
  ;; Design matrix

  A = pld_design_matrix(pixgrid,ndat,ncoef,NOCONSTANT=noconstant,ORDER=order)

IF N_ELEMENTS(sigma_flux) NE ndat then sigma_flux = REPLICATE(1.d0,ndat)
sfmatrix = MAKE_ARRAY(ncoef,VALUE=1d0) # sigma_flux
A /= sfmatrix
;;
;; Use Singular Value Decomposition on A
;;
SVDC,A,W,U,V,/DOUBLE
;; Eliminate singular values
ISING = WHERE(W LE 1e-6,nsing)
IF NSING NE 0 THEN w[ising] = 0.d0
;; Solve the system
IF N_ELEMENTS(scale) EQ 0 THEN scale = MEDIAN(flux,/EVEN,/DOUBLE)
coeffs = SVSOL(U,W,V,flux/scale/sigma_flux,/DOUBLE)
;print,'Coeffs = ',coeffs
sigma_coeffs = DBLARR(ncoef)
FOR i = 0,ncoef-1 DO sigma_coeffs[i] = SQRT( TOTAL( (V[*,i]/W[i])^2, /DOUBLE, /NAN ) )
;; Return the predicted flux 
fpredict = DBLARR(ndat)
FOR j = 0,ndat-1 DO fpredict[j] = TOTAL(A[*,j]*coeffs,/DOUBLE,/NAN)
fcorrect = flux/fpredict

RETURN,coeffs
END

