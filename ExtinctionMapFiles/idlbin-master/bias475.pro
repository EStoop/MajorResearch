pro bias475

bias = [450,450,450,450,450,450,450,450,450,450,450,450,450,450,475,475,475,475,475,475,475,475,475,475,475,475,475,475]

meannoise = [4.71025,1.59474,0.21878,0.07866,0.045559,0.007275,0.00419044,4.90694,1.70162,0.203714,0.067596,0.04094,0.01052,0.0071476,4.83539,1.618764,0.2223,0.07824,0.045468,0.0074933,0.004348,5.02886,1.7015448,0.206098,0.0685198,0.042384,0.0108814,0.0075295]

ch = [1,1,1,1,1,1,1,2,2,2,2,2,2,2,1,1,1,1,1,1,1,2,2,2,2,2,2,2]

noisepix = [0,0,0,0,0,0.00185246,0,0,0,0,0,0,0.00109872,0,0,0,0,0,0,0.002531,0,0,0,0,0,0,0.00149657,0]

exptime = [0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200,0.1,0.4,2,6,12,100,200]

ch1_100 = where(ch eq 1 and exptime eq 100)
ch1_12 = where(ch eq 1 and exptime eq 12)
ch1_0p4= where(ch eq 1 and exptime eq 0.4)
ch2_100 = where(ch eq 2 and exptime eq 100)
ch2_12 = where(ch eq 2 and exptime eq 12)
ch2_0p4= where(ch eq 2 and exptime eq 0.4)
!P.multi=[0,1,1]

ps_start, filename= '/Users/jkrick/iwic/iwicDAA/noisypixvsbias475.ps'

plot, bias(ch1_100), noisepix(ch1_100), psym = 2, xrange=[425, 500], xtitle = 'bias', ytitle = 'fraction of pixels 4x zodi noise', title = '29K 100s r36 g368_361 d36', yrange=[0,.004]
oplot, bias(ch1_100), noisepix(ch1_100)
oplot, bias(ch2_100), noisepix(ch2_100), psym = 4
oplot, bias(ch2_100), noisepix(ch2_100)
legend, ['ch1', 'ch2'], psym=[2,4],/left,/top

;put number of CRs  onto plots

oplot, [400,600], [12./(256.*256.), 12./(256.*256.)]
xyouts, 425, 12./(256.*256.), 'number of CR in 12s'
oplot, [400,600], [100./(256.*256.), 100./(256.*256.)]
xyouts, 425, 100./(256.*256.), 'number of CR in 100s'
oplot, [400,600], [200./(256.*256.), 200./(256.*256.)]
xyouts, 425, 200./(256.*256.), 'number of CR in 200s'


ps_end, /png

!P.multi=[0,1,1]
;ps_start, filename= '/Users/jkrick/iwic/iwicDAA/noisevsbias475.ps'
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)

plot, bias(ch1_0p4), meannoise(ch1_0p4)/max(meannoise(ch1_0p4)), psym = 2, xrange=[425, 500], xtitle = 'bias', ytitle = 'fractional mean noise MJy/sr', title = '29K r36 g368_361 d36', yrange=[.95,1.0], charsize=1.5, /nodata
oplot, bias(ch1_0p4), meannoise(ch1_0p4)/max(meannoise(ch1_0p4)), psym = 2, color=greencolor
oplot, bias(ch1_0p4), meannoise(ch1_0p4)/max(meannoise(ch1_0p4)), color=greencolor
oplot, bias(ch2_0p4), meannoise(ch2_0p4)/max(meannoise(ch2_0p4)), psym = 4, color=greencolor
oplot, bias(ch2_0p4), meannoise(ch2_0p4)/max(meannoise(ch2_0p4)), color=greencolor

legend, ['0.4s', '12s','100s'], linestyle = [0,0,0], color = [greencolor, redcolor, bluecolor], /right,/bottom
legend, ['ch1', 'ch2'], psym=[2,4],/left,/top


oplot, bias(ch1_12), meannoise(ch1_12)/max(meannoise(ch1_12)), psym = 2, color=redcolor;, xrange=[425, 500], xtitle = 'bias', ytitle = 'mean noise MJy/sr', title = '29K 12s r36 g368_361 d36', yrange=[.035,.05], charsize=1.5
oplot, bias(ch1_12), meannoise(ch1_12)/max(meannoise(ch1_12)), color=redcolor
oplot, bias(ch2_12), meannoise(ch2_12)/max(meannoise(ch2_12)), psym = 4, color=redcolor
oplot, bias(ch2_12), meannoise(ch2_12)/max(meannoise(ch2_12)), color=redcolor

oplot, bias(ch1_100), meannoise(ch1_100)/max(meannoise(ch1_100)), psym = 2, color=bluecolor;, xrange=[425, 500], xtitle = 'bias', ytitle = 'mean noise MJy/sr', title = '29K 100s r36 g368_361 d36', yrange=[.004,.015], charsize=1.5
oplot, bias(ch1_100), meannoise(ch1_100)/max(meannoise(ch1_100)), color=bluecolor
oplot, bias(ch2_100), meannoise(ch2_100)/max(meannoise(ch2_100)), psym = 4, color=bluecolor
oplot, bias(ch2_100), meannoise(ch2_100)/max(meannoise(ch2_100)), color=bluecolor

;ps_end, /png

end
