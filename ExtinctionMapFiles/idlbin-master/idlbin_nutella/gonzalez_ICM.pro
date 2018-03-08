pro enrichment
;taken from Gonzalez Zaritsky Zabludoff 2007
; and Sivanandam et al 2008 (maybe 2009)

cluster_names = ["A0496", "A1651", "A2811", "A2984", "A3112", "A3693", "A4010", "A4059", "AS84", "AS540"]
Fe_fraction = [29, 20, 20, 188, 36, 28,42,47,69,55]
fracerr=[6,5,8,43,9,12,11,14,17,15]
bcg_icl_fraction = [0,24,33,53,25,26,35,0,33,0]
bcg_icl_err = [0,4,8,5,3,5,8,0,4,0]


icl_fraction = 0.8 * bcg_icl_fraction
icl_err = 0.8 * bcg_icl_err

ps_open, filename='/Users/jkrick/spitzer/virgo_proposal/enrichment.ps',/portrait,/square,/color

ploterror, icl_fraction, Fe_fraction, icl_err, fracerr, psym = 2, thick = 3, xtitle = 'ICL fraction', ytitle = 'Fe fraction', charthick = 3, xthick=3,ythick=3, errthick=3

ps_close, /noprint, /noid
end
