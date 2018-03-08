PRO MEDCLIP, Image, med, Sigma, CLIPSIG=clipsig, MAXITER=maxiter, $
    CONVERGE_NUM=converge_num, SUBS=subs
;+
; NAME:
;       MEDCLIP
;
; PURPOSE:
;       Computes an iteratively sigma-clipped median on a data set
; EXPLANATION:
;       Clipping is done about median.
;       Called by SKYADJ_CUBE
;
; CATEGORY:
;       Statistics
;
; CALLING SEQUENCE:
;       MEANCLIP, Data, med, Sigma
;
; INPUT POSITIONAL PARAMETERS:
;       Data:     Input data, any numeric array
;
; OUTPUT POSITIONAL PARAMETERS:
;       med:     N-sigma clipped median.
;       Sigma:    Standard deviation of remaining pixels.
;
; INPUT KEYWORD PARAMETERS:
;       CLIPSIG:  Number of sigma at which to clip.  Default=3
;       MAXITER:  Ceiling on number of clipping iterations.  Default=5
;       CONVERGE_NUM:  If the proportion of rejected pixels is less
;           than this fraction, the iterations stop.  Default=0.02, i.e.,
;           iteration stops if fewer than 2% of pixels excluded.
;
; OUTPUT KEYWORD PARAMETER:
;       SUBS:     Subscript array for pixels finally used.
;
;
; MODIFICATION HISTORY:
;       Based on meanclip written by:   RSH, RITSS, 21 Oct 98
;       20 Jan 99 - Added SUBS, fixed misplaced paren on float call,
;                   improved doc.  RSH
;       Edited by R. Paladini, 14 Aug 2010 to output median instead of mean
;-
IF n_params(0) LT 1 THEN BEGIN
    print, 'CALLING SEQUENCE:  MEANCLIP, Image, med, Sigma'
    print, 'KEYWORD PARAMETERS:  CLIPSIG, MAXITER, CONVERGE_NUM, ' $
        + 'SUBS'
    RETALL
ENDIF

prf = 'MEDCLIP:  '

IF n_elements(maxiter) LT 1 THEN maxiter = 5
;IF n_elements(clipsig) LT 1 THEN clipsig = 3
IF n_elements(clipsig) LT 1 THEN clipsig = 1
IF n_elements(converge_num) LT 1 THEN converge_num = 0.02
clipsig = 1.5

;print,image
subs = where(finite(image),ct)
print, 'finite ',ct
iter=0
REPEAT BEGIN
    skpix = image[subs]
; print,skpix
    iter = iter + 1
    lastct = ct
;   print,' last ct', lastct,ct
    medval = median(skpix)
    sig = stdev(skpix)
    sig1 = clipsig*sig
;    print,'iter med sig sig1', iter, medval, sig, sig1
;   print, abs(skpix-medval) 
    wsm = where(abs(skpix-medval) LT clipsig*sig,ct)
;   print,'clipped ct', ct
;   print,wsm
    IF ct GT 0 THEN subs = subs[wsm]
;  print,'sub',subs
    f1 = float(abs(ct-lastct))/lastct
;  print,' f1=',f1
ENDREP UNTIL (float(abs(ct-lastct))/lastct LT converge_num) $
          OR (iter GT maxiter)


;print,'elems=',n_elements(subs)
;   print,'sub',subs
med = median(image[subs])
sigma = stddev(image[subs],/NaN)
print,'final med sigma',med,sigma


RETURN
END

