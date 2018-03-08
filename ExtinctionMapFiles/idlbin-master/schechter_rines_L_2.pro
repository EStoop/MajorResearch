FUNCTION schechter_rines_L_2, X, P

;schechter function (L)defined as

return, (2.30718E11 / 2.70E10) * ((X/2.70E10)^(-0.278)) * (exp(-X/2.70E10))

;return, (P(0)*0.4*alog(10))*(exp(-10^(0.4*((-21.32)-X))))*(10^(0.4*((-21.32)-X)*((-1.28)+1)))

;P(0) = normalization
; rest are fixed to Rines & geller 2008 for the virgo cluster



END

