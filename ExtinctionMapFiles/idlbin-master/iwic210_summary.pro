pro iwic210_summary
;all numbers from dark_summary.xlsx

;first check repeatability -3.6, 750
;ps_start, filename = '/Users/jkrick/iwic/iwic_recovery12/repeat.ps'

!P.multi = [0,0,1]
exptime = [0.1,0.4,2,12,100,200]
ch1mean=[4.5683,1.505,0.16984,0.0310308,0.0066448,0.004247]
ch1mean_2 = [4.6263,1.5169,0.172072,0.03139,0.006728,0.0043104]

plot, exptime, ch1mean, psym = 4, xtitle = 'exposure time', ytitle = 'mean of stddev in Mjy/sr', thick = 3, /ylog, /xlog
oplot, exptime, ch1mean_2, psym = 5, thick =3
;oplot, exptime, 1/exptime
;ps_end, /png

end
