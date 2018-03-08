FUNCTION schechter_rines, X, P

;schechter function (M)defined as

return, (P(0)*0.4*alog(10))*(exp(-10^(0.4*((-21.32)-X))))*(10^(0.4*((-21.32)-X)*((-1.28)+1)))

;P(0) = normalization
; rest are fixed to Rines & geller 2008 for the virgo cluster



END

