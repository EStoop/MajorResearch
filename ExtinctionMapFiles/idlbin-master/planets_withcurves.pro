pro planets_withcurves
ps_start, filename='/Users/jkrick/planets_proposal/2011/planets_withcurves.ps'


;those with secondaries
s_pl_mass = [0.85,0.599,1.057,7.725,10.43,0.847,1.02,2.06,0.59,0.96,0.57,1.72]

s_st_mass = [1.22,0.98,1.29,1.319,1.24,0.84,1.12,1.24,1.01,1.276,0.98,1.32]

plot, s_st_mass, s_pl_mass, psym = 5, xrange = [0.7, 1.4], yrange = [0.1,20], xstyle = 1, ystyle = 1,xtitle = 'Stellar Mass (Solar masses)', ytitle = 'Planetary Mass (Jupiter masses)',/ylog, thick = 3, xthick = 3, ythick = 3, charthick = 3

;those without secondaries

n_pl_mass = [2.2,4.193,2.147,0.46,0.542,0.855,2.712]

n_st_mass= [1.386,1.218,0.916,1.00,1.18,1.022,1.216]

oplot, n_st_mass, n_pl_mass, psym = 6, thick = 3

;those already done
d_pl_mass = [0.714,1.15,0.356,8.74,11.79,0.9,1.41]

d_st_mass = [1.148,0.8,1.3,1.36,1.213,1,1.35]

oplot, d_st_mass, d_pl_mass, psym = 2, thick = 3



;plot up some lines for comparison
a = findgen(10) / 2
b = fltarr(10)
b[*] = 0.55

;oplot, a, b, linestyle = 2



ps_end;, /noprint,/noid

end
