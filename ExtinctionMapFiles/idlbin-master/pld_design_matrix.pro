FUNCTION pld_design_matrix,pixgrid,ndat,ncoef,NOCONSTANT=noconstant,ORDER=order
;; Default is 1st order
   IF N_ELEMENTS(order) EQ 0 THEN order=1
   SZ = size(pixgrid)
   ndat = SZ[0] GE 3 ? SZ[3] : 1
   Ngrid = SZ[1]
;;; Set up the background indices as the 1-pixel boundary of the nxn region
   IF NGRID EQ 5 THEN BEGIN
    ;; 5x5
      bgindices1 = [0,1,2,3,4,0,1,2,3,4,0,0,0,4,4,4]
      bgindices2 = [0,0,0,0,0,4,4,4,4,4,1,2,3,1,2,3]
   ENDIF ELSE BEGIN
    ;; Assume 7x7
      bgindices1 = [0,1,2,3,4,5,6,0,1,2,3,4,5,6,0,0,0,0,0,6,6,6,6,6]
      bgindices2 = [0,0,0,0,0,0,0,6,6,6,6,6,6,6,1,2,3,4,5,1,2,3,4,5]
   ENDELSE
   nbg = 4*ngrid-4
   bgdata = FLTARR(nbg,ndat)
   FOR i = 0,nbg-1 DO bgdata[i,*] = pixgrid[bgindices1[i],bgindices2[i],*]
   bgval = MEAN( bgdata,DIMENSION=1,/NAN,/DOUBLE)
;;; Extract (NGRID-2)x(NGRID-2)x NDAT region, convert to (NGRID^2 - 4*NGRID + 4) x NDAT
   ncoef1 =  ngrid^2-4*ngrid+4   ;; Number of coefficients for a first order fit (no constant)
   pixsrc = REFORM(pixgrid[1:ngrid-2,1:ngrid-2,*],ncoef1,ndat)
;;; Subtract the background value
   FOR j = 0,ndat-1 DO BEGIN
      pixsrc[*,j] -= bgval[j]
      pixsrc[*,j] /= TOTAL(pixsrc[*,j],/DOUBLE,/NAN)  ; normalize by the sum of pixel values
   ENDFOR

   ncoef1final = KEYWORD_SET(NOCONSTANT) ? ncoef1 : ncoef1 + 1   ;;; Account for the constant 
   
   ncoef2 =  ncoef1 + ncoef1*(ncoef1-1)/2.0   ;; Number of (additional) coefficients for a 2nd order fit
   ncoef = ORDER EQ 2 ? ncoef1final + ncoef2 : ncoef1final 
   
;; Design matrix
   A = DBLARR(ncoef,ndat)
   IF ~KEYWORD_SET(NOCONSTANT) THEN A[0,*] = 1d0
   A[ncoef1final-ncoef1:ncoef1final-1,*] = pixsrc
   IF ORDER EQ 2 THEN BEGIN
      A[ncoef1final:ncoef1final+ncoef1-1,*] = pixsrc^2
      k = ncoef1final+ncoef1
      FOR i = 0,ncoef1-2 DO BEGIN
         FOR j = i+1,ncoef1-1 DO BEGIN
            A[k,*] = pixsrc[i,*] * pixsrc[j,*]
            k++
         ENDFOR
      ENDFOR
   ENDIF
   

RETURN,A
END
