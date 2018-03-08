FUNCTION schechter3, X, P

;schechter function (M)defined as
;return, (10^(-0.4*(P(1)+1.)*X))*exp(-10^(0.4*(P(0)-X)))
;return, (-P(2)*0.4*alog(10))*(10^(-0.4*(P(1)+1.)*(X-P(0))))*exp(-10^(0.4*(P(0)-X)))

return, (P(0)*0.4*alog(10))*(exp(-10^(0.4*(19.4-X))))*(10^(0.4*(19.4 -X)*((-0.94)+1)))

;P(0) = mstar
;P(1) = alpha

;schechter function (L)
;return, 

END

