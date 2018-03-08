FUNCTION schechter1, X, P

;schechter function (M)defined as
;return, (10^(-0.4*(P(1)+1.)*X))*exp(-10^(0.4*(P(0)-X)))
;return, (-P(2)*0.4*alog(10))*(10^(-0.4*(P(1)+1.)*(X-P(0))))*exp(-10^(0.4*(P(0)-X)))

;;return, (P(2)*0.4*alog(10))*(exp(-10^(0.4*(P(0)-X))))*(10^(0.4*(P(0)-X)*(P(1)+1)))
return, (P(0)*0.4*alog(10))*(exp(-10^(0.4*((-25.3)-X))))*(10^(0.4*((-25.3)-X)*((-1.03)+1)))

;P(0) = mstar
;P(1) = alpha

;schechter function (L)
;return, 

END

