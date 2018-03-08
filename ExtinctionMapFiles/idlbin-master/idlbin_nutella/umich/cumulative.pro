PRO cumulative
close, /all
I0 = 120.
re = 5.
n = 4.
b = 2.*n - 0.331

outer = 140
rfin = findgen(outer) / 10 + 0.1
L = findgen(outer) / 10
Lfin = (I0 * re^2. * 2.* !PI * n) / (b^(2.*n)) * IGAMMA(2.*n, b*(rfin/re)^(1./n))
Ltot = (I0 * re^2. * 2.* !PI * n) / (b^(2.*n)) * GAMMA(2.*n)
A =  IGAMMA(8., 7.67*(rfin/re)^(0.25))
B = GAMMA(8.)
print, A
print, B
F = 1 - A/B
;print, F
plot, rfin, Lfin
;plot, rfin,L2
plot, rfin, IGAMMA(8.,7.67*rfin^(1/4), itmax = 400, eps = 8E-10)*GAMMA(8.)
END
