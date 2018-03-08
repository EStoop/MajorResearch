;******************************************************************************
;This function returns the integral of the differential schecter
;luminosity function, given a lower and upper bound and the values of
;Mstar and alpha.
;******************************************************************************

function int_lumfunc, lower, upper, phistar, Mstar, alpha

MY_QSIMP,'DIFF_LUMFUNC',lower, upper, Mstar, alpha, integral,MAX_ITER=30.0

RETURN, integral * phistar

END
