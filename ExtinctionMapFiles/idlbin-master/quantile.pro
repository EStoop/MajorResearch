function QUANTILE,data,q,P=p,METHOD=method
;;
;;  Compute the set of quantiles for a dataset DATA.  
;;
;; The output is a vector of q+1 elements giving the p = k/q - ile values (p ranges from 0 to 1 in 1/q intervals), 
;; determined using one of ten methods, set by keyword METHOD, a value from 1 to 10.  Default METHOD=5.  
;; If q is not set, assume we are looking for quartiles (q=4).
;;
;; METHOD = 1,2,5-10 give the correct median for an odd number of data points.  METHOD = 2, 5-10 gives the standard 
;;          approximation of median for an even number of data points (average of the two points flanking the halfway mark).   
;;
;;       METHOD: 1) The inverse of the empirical cumulative distribution function, where quantiles correspond exactly to
;;                  data points.
;;               2) The same as (1), but with averaging at discontinuities.
;;               3) The observation numbered closest to N*p, where N is the number of data points and p is the 
;;                  quantile.
;;               4) Linear interpolation of the empirical distribution function.
;;               5) Piecewise linear function where the knots are the values midway through the steps of the empirical
;;                  distribution function.
;;               6) Linear interpolation of the expectations for the order statistics for the uniform distribution
;;                  on [0,1].
;;               7) Linear interpolation of the modes for the order statistics for the uniform distribution on
;;                  [0,1].
;;               8) Linear interpolation of the approximate medians for order statistics.
;;               9) Approximately unbiased quantile estimates for the expected order statistics if x is normal.
;;              10) The order statistic with the least expected square deviation relative to p.   
;; 
;;
IF N_PARAMS() LT 2 THEN q = 4
N = N_ELEMENTS(data)
IF q gt N THEN BEGIN
  PRINT,'QUANTILE.PRO: more quantiles requested than number of data points.  Returning.'
  RETURN,-1   
ENDIF
IF N LT 5*Q THEN PRINT,'QUANTILE.PRO: small number of data points.  Be warned that the statistics might be highly uncertain.'
IF N_ELEMENTS(METHOD) EQ 0 THEN METHOD=5

;;; Create the quantile array, and the set of quantiles to determine
Qp = MAKE_ARRAY(q+1, TYPE=SIZE(data,/TYPE))
p = DINDGEN(q+1) / q
;;; Sort the dataset.
ds = DATA[SORT(data)]
;; Determine h, the method-dependent set of indices into the sorted dataset, on which we base the quantiles.   
CASE METHOD of
  1: h = N * p + 0.5d0
  2: h = N * p + 0.5d0
  3: h = N * p
  4: h = N * p
  5: h = N * p + 0.5d0
  6: h = (N + 1d0) * p
  7: h = (N - 1d0) * p + 1d0
  8: h = (N + 1d/3d) * p + 1d/3d
  9: h = (N + 1d/4d) * p + 3d/8d
  10: h = (N + 2d) * p - 0.5d
ENDCASE
;; Determine Qp
Qinterp = ds[FLOOR(h) - 1] + (h - FLOOR(h)) * (ds[FLOOR(h)] - ds[FLOOR(h)-1])
n1 = 0
n_n = 0
CASE METHOD of
  1: BEGIN
        Qp = ds[CEIL(h - 0.5d - 1)]
        Qp[0] = ds[0]
     END  
  2: BEGIN
        Qp = 0.5 * (ds[CEIL(h-0.5d) - 1] + ds[FLOOR(h+0.5d) - 1] )
        Qp[0] = ds[0]
        Qp[q] = ds[N-1]
     END
  3: BEGIN
        Qp = ds[ROUND(h)-1]
        i1 = WHERE(p LE 0.5/N,n1)
     END
  4: BEGIN
       Qp = Qinterp
       i1 = WHERE(p LT 1d/N,n1)
       i_n = WHERE(p EQ 1,n_n)
     END
  5: BEGIN
       Qp = Qinterp
       i1 = WHERE(p LT 0.5d/N,n1)
       i_n = WHERE(p GE (N-0.5d)/N,n_n)
     END
  6: BEGIN
       Qp = Qinterp
       i1 = WHERE(p LT 1d/(N+1),n1)
       i_n = WHERE(p GE N/(N+1d),n_n)
     END
  7: BEGIN
       Qp = Qinterp
       i_n = WHERE(p EQ 1,n_n)
     END
  8: BEGIN
       Qp = Qinterp
       i1 = WHERE(p LT (2d/3)/(N+1d/3),n1)
       i_n = WHERE(p GE (N-1d/3)/(N+1d/3),n_n)
     END
  9: BEGIN
       Qp = Qinterp
       i1 = WHERE(p LT (5d/8)/(N+0.25d),n1)
       i_n = WHERE(p GE (N-3d/8)/(N+0.25d),n_n)
     END
 10: BEGIN
       Qp = Qinterp
       i1 = WHERE(p LT 1.5d/(N+2d),n1)
       i_n = WHERE(p GE (N+0.5d)/(N+2d),n_n)
     END
    
ENDCASE
IF n1 NE 0 THEN Qp[i1] = ds[0]
IF n_n NE 0 THEN Qp[i_n] = ds[N-1]

RETURN,Qp
END

