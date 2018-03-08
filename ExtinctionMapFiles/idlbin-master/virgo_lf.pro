pro virgo_lf

;from Rines & Geller 2008
mag = [9.8,10.3,10.8,11.3,11.8,12.3,12.8,13.3,13.8,14.3,14.8,15.3,15.8,16.3,16.8,17.3]
abs_mag =mag - 31.1
lum = 10^((abs_mag - 4.76)/(-2.5))
print, 'lum' , lum
n = [2,6,10,12,14,18,18,18,13,28,30,40,37,45,55,60]
noise = [1.5,2,4,4,4,4,4,4,5,5,5,5,5,10,10,10]
;they quote
Mstar = -21.32
alpha = -1.28
; I calculated
Lstar = 2.7E10 ;solar luminosities

;keep these fixed, then fit this data with a schechter function to
;determine the normalization

start = [2.E10, -1.28]
result= MPFITFUN('schechter_rines_L',lum,n, noise, start)    ;ICL

;plot, abs_mag, n,/ylog
;oplot, abs_mag,  (result(0)*0.4*alog(10))*(exp(-10^(0.4*((-21.32)-abs_mag))))*(10^(0.4*((-21.32)-abs_mag)*((-1.28)+1)))

plot, lum, n, /ylog, xrange = [2E7,4E10], /xlog
oplot, lum, (result(0) / 2.70E10) * ((lum/2.70E10)^(result(1))) * (exp(-lum/2.70E10)), linestyle = 2

;now what are the units of this thing,
;number of galaxies per area of the survey = 1 Mpc radius = 3.14Mpc^2

; and if I integrate it, what units will I be left with?

;want to integrate over abs_mag from -22 to -11
;use qsimp or int_tabulated

;this doesn't use the schecter function, just the data
;also then doesn't extrapolate
;just use for an idea
;value = int_tabulated(abs_mag, float(n))
;print, value
abs_limit= [-11,-18]
lum_limit = 10^((abs_limit - 4.76)/(-2.5))
print, 'lumlimit', lum_limit
int_value = qsimp('schechter_rines_L_2', lum_limit(0), lum_limit(1))
print, int_value
end
