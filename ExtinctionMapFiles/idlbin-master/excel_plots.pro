pro excel_plots
!p.multi = [0, 2, 2]

ps_start, filename = '/Users/jkrick/iwic/meanvstemp_r35b75.ps'
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

;12s ch1 750
cernox = [21,22,23.6,23.5,23.5,23.5,25.4,39.7,30,34.8,37.2,33.2,30.7]
diode = [21.6,21.6,21.6,21.1,21.1,21.1,24.6,36.1,28.1,32,34,30.8,28.7]
mean12 = [0.03414,0.03461,0.035809,0.0311878,0.03139,0.0310308,0.0396,0.2603,0.04834,0.087655,0.15095,0.053124,0.046622]

;plot, cernox, mean12, psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'cernox temperature (K)', ytitle ='mean (stddev) in Mjy/sr', title ='Ch1 12s b=750',/ylog
plot, diode, mean12, psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'corrected cernox temperature (K)', ytitle ='mean (stddev) in Mjy/sr', title ='Ch1 12s b=750', /ylog, xrange=[20,40]

;12s ch2 500
diode = [25.4,25.4,25.4,27.2,33.3,28.7,31.8,35,36.5,28,30.1,32,35.1,34]
cernox = [20.8,21.8,26.3,28.8,36.4,30.8,34.6,38.3,40.2,29.9,32.4,34.8,38.4,37.2]
mean12 = [0.04885,0.04651,0.01674,0.0442,0.0474,0.045053,0.046017,0.048545,0.04614,0.0415608,0.042004,0.042928,0.044847,0.044086]


;plot, cernox, mean12, psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'cernox temperature (K)', ytitle ='mean (stddev) in Mjy/sr', title ='Ch2 12s b=500',/ylog
plot, diode, mean12, psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'corrected cernox temperature (K)', ytitle ='mean (stddev) in Mjy/sr', title ='Ch2 12s b=500', /ylog, xrange=[20,40]


;100s
;not enough data 
diode = [21.6,21.6,21.6,21.1,21.1,21.1,30.8,28.7]
cernox = [21,22,23.6,23.5,23.5,23.5,33.2,30.7]
mean100 = [0.070602,0.071144,0.0735162,0.006641,0.006728,0.0066448,0.021193,0.010209]


;0.4s ch1
diode = [21.1,21.1,21.1,24.6,36.1,28.1,32,34,30.8,28.7]
cernox = [23.5,23.5,23.5,25.4,39.7,30,34.8,37.2,33.2,30.7]
mean0p4 = [1.58478,1.5169,0.505,1.6024,2.326,1.74345,2.0204,2.17344,1.75394,1.60866]

;plot, cernox, mean0p4, psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'cernox temperature (K)', ytitle ='mean (stddev) in Mjy/sr', title ='Ch1 0.4s b=750',/ylog, xrange=[20,40]
plot, diode, mean0p4, psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'corrected cernox temperature (K)', ytitle ='mean (stddev) in Mjy/sr', title ='Ch1 0.4s b=750', /ylog, xrange=[20,40]

;0.4s ch2
diode=[27.2,33.3,28.7,31.8,35,36.5,28,30.1,32,35.1,34]
cernox = [28.8,36.4,30.8,34.6,38.3,40.2,29.9,32.4,34.8,38.4,37.2]
mean0p4 = [1.8027,1.9597,1.8212305,1.90297,2.04736,1.9149,1.70037,1.69424,1.76378,1.87157,1.8203]

;plot, cernox, mean0p4, psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'cernox temperature (K)', ytitle ='mean (stddev) in Mjy/sr', title ='Ch2 0.4s b=500',/ylog, xrange=[20,40]
plot, diode, mean0p4, psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'corrected cernox temperature (K)', ytitle ='mean (stddev) in Mjy/sr', title ='Ch2 0.4s b=500', /ylog, xrange=[20,40]



;is there a difference between 30 and 32


diode = [30.8,28.7]
cernox = [33.2,30.7]
mean12 = [0.053124,0.046622]
mean0p4 = [1.75394,1.60866]

;plot, cernox, mean12/max(mean12), psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'cernox temperature (K)', ytitle ='fractional mean (stddev) ', title ='Ch1 ', yrange=[0.8,1.0]
;oplot, cernox, mean0p4/max(mean0p4), psym = 2, thick = 3
;legend, ['12s', '0.4s'], psym=[4,2],/right,/bottom

plot, diode, mean12/max(mean12), psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'corrected cernox temperature (K)', ytitle ='fractional mean (stddev) ', title ='Ch1 ' , yrange=[0.8,1.0]
oplot, diode, mean0p4/max(mean0p4), psym = 2, thick = 3
legend, ['12s', '0.4s'], psym=[4,2],/right,/bottom

diode=[28,30.1,32]
cernox = [29.9,32.4,34.8]
mean0p4=[1.70037,1.69424,1.76378]
mean12=[0.0415608,0.042004,0.042928]


;plot, cernox, mean12/max(mean12), psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'cernox temperature (K)', ytitle ='fractional mean (stddev) ', title ='Ch2 ', yrange=[0.8,1.0]
;oplot, cernox, mean0p4/max(mean0p4), psym = 2, thick = 3
;legend, ['12s', '0.4s'], psym=[4,2],/right,/bottom

plot, diode, mean12/max(mean12), psym = 4, thick =3, charthick = 3, xthick=3,ythick=3,xtitle = 'corrected cernox temperature (K)', ytitle ='fractional mean (stddev) ', title ='Ch2 ' , yrange=[0.8,1.0]
oplot, diode, mean0p4/max(mean0p4), psym = 2, thick = 3
legend, ['12s', '0.4s'], psym=[4,2],/right,/bottom

ps_end, /png

;--------
;what about bias
!p.multi = [0, 2, 1]

ps_start, filename = '/Users/jkrick/iwic/noisypixvsbias.ps'

;ch1 
bias_ch1 = [450,600,750]
mean30_12s=[0.0498516,0.045598,0.046622] ;30.7
mean32_12s = [0.051963,0.049096,0.053124] ;33.2
noisy30_12s = [67,299,991]
noisy32_12s = [96,555,3618]

mean30_100s = [0.00771637,0.007854,0.010209] 
mean32_100s = [0.0077577,0.008179,0.021193]
noisy30_100s = [207,720,5360]
noisy32_100s = [355,2463,8973]

;plot, bias_ch1, mean30_12s, psym = 4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'bias', ytitle ='mean (stddev) in Mjy/sr ', title ='Ch1 ', yrange=[.001,.1],/ylog
;oplot, bias_ch1, mean32_12s, psym = 2,thick=3
;oplot, bias_ch1, mean30_100s, psym=4, thick=3
;oplot, bias_ch1, mean32_100s, psym=2, thick=3
;legend, ['corr. cernox28.7K', 'corr cernox 30.8K'], psym=[4,2],/right,/bottom

plot, bias_ch1, noisy30_12s, psym =4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'bias', ytitle ='number of 10sigma pixels ', title ='Ch1 ', yrange=[10,10000],/ylog
oplot, bias_ch1, noisy32_12s, psym = 2,thick=3
oplot, bias_ch1, noisy30_100s, psym=4, thick=3, color = redcolor
oplot, bias_ch1, noisy32_100s, psym=2, thick=3, color=redcolor
legend, ['corr. cernox 28.7K', 'corr cernox 30.8K'], psym=[4,2],/right,/bottom
xyouts, 600,60,'100s', color=redcolor, charthick=3
xyouts, 600,30,'12s', charthick=3


;ch2
bias_ch2=[500,450,600]
mean_30_12s=[0.042004,0.041059,0.045066] ;32.5
mean_32_12s = [0.042928,0.04229,0.045846] ;34.9
noisy30_12s = [114,76,216]
noisy32_12s=[216,149,377]

mean30_100s = [0.010966,0.0104835,0.0120862]
mean32_100s = [0.010764,0.010387,0.0122286]
noisy30_100s = [225,183,381]
noisy32_100s=[348,273,840]

;plot, bias_ch2, mean30_12s, psym = 4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'bias', ytitle ='mean (stddev) in Mjy/sr ', title ='Ch2 ', yrange=[.01,.1],/ylog, xrange=[400,700]
;oplot, bias_ch2, mean32_12s, psym = 2,thick=3
;oplot, bias_ch2, mean30_100s, psym=4, thick=3
;oplot, bias_ch2, mean32_100s, psym=2, thick=3
;legend, ['corr. cernox 30.1K', 'cernox 32.0K'], psym=[4,2],/right,/center

plot, bias_ch2, noisy30_12s, psym =4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'bias', ytitle ='number of 10sigma pixels ', title ='Ch2 ', xrange=[400,700], yrange=[10,1000],/ylog
oplot, bias_ch2, noisy32_12s, psym = 2,thick=3
oplot, bias_ch2, noisy30_100s, psym=4, thick=3, color = redcolor
oplot, bias_ch2, noisy32_100s, psym=2, thick=3, color = redcolor
legend, ['corr. cernox 30.1K', 'corr cernox 32.0K'], psym=[4,2],/right,/bottom
xyouts, 600,40,'100s', color=redcolor, charthick=3
xyouts, 600,25,'12s', charthick=3




ps_end, /png

;---------------------------------------------------------
;latent as func of temp
!p.multi = [0, 2, 1]
ps_start, filename = '/Users/jkrick/iwic/latentvstemp.ps'

; ch1b = 450
temp1 = [30.8,28.7,31.1,30.8,28.7,31.1,30.8,28.7,31.1,30.8,28.7,31.1,30.8,28.7,31.1,30.8,28.7,31.1,30.8,28.7,31.1]
strength1 = [0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386]
decay1=[0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746]

;plot, temp1, strength1, psym = 4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'corrected cernox temp (deg K)', ytitle ='fractional strength of latent ', title ='Ch1 b = 450 ', yrange=[.0015,.0025]

;plot, temp1, decay1, psym = 4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'corrected cernox temp (deg K)', ytitle ='decay of latent(counts/s)', title ='Ch1 b = 450 ', yrange=[.08,.12]

;ch2 b = 450
temp2= [32,30.2,31.1,32,30.2,31.1,32,30.2,31.1,32,30.2,31.1,32,30.2,31.1,32,30.2,31.1,32,30.2,31.1]
strength2 = [0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386,0.00181787,0.00204855,0.00179386]
decay2 = [0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746,0.099695,0.110081,0.0990746]

plot, temp2, strength2, psym = 4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'corrected cernox temp (deg K)', ytitle ='fractional strength of latent ', title ='Ch2 b = 450 ', xrange=[30,32.5], yrange=[.0015,.0025]

plot, temp2, decay2, psym = 4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'corrected cernox temp (deg K)', ytitle ='decay of latent(counts/s)', title ='Ch2 b = 450 ', xrange=[30,32.5], yrange=[.08,.12]

ps_end, /png

;---------------------------------------------------------
;latent as func of bias
!p.multi = [0, 2,1]
ps_start, filename = '/Users/jkrick/iwic/latentvsbias.ps'

;ch1
bias1 = [450,600,750]
strength_32s = [0.00181787,0.00211641,0.00254605] ;30.8
decay_32s = [0.099695,0.111524,0.1314]

strength_30s = [0.00204855,0.00221165,0.00170674] ;28.7
decay_30s=[0.110081,0.118162,0.0529277]


;plot, bias1, strength_30s, psym =4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'bias', ytitle ='fractional strength of latent ', title ='Ch1 ',yrange=[0.001,0.003]
;oplot, bias1, strength_32s, psym = 2,thick=3
;legend, ['corr. cernox 28.7K', 'corr cernox 30.8K'], psym=[4,2],/right,/bottom
;plot, bias1, decay_30s, psym =4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'bias', ytitle ='decay of latent(counts/s)', title ='Ch1 ', yrange=[.02,.14]
;oplot, bias1, decay_32s, psym = 2,thick=3
;legend, ['corr. cernox 28.7K', 'corr cernox 30.8K'], psym=[4,2],/right,/bottom


;ch2
bias2 = [500,450,600]
strength_30s = [0.00221165, 0.00204855, 0.00170674] ;30.2
decay_30s = [0.118162,0.110081,0.0529277]

strength_32s=[0.00211641,0.00181787,0.00254605] ;32
decay_32s =[0.111524,0.099695,0.1314]



plot, bias2, strength_30s, psym =4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'bias', ytitle ='fractional strength of latent ', title ='Ch2 ', xrange=[400,650],yrange=[0.001,0.003]
oplot, bias2, strength_32s, psym = 2,thick=3
legend, ['corr. cernox 30.1K', 'corr cernox 32.0K'], psym=[4,2],/right,/bottom
plot, bias2, decay_30s, psym =4, thick = 3, charthick = 3, xthick=3,ythick=3,xtitle = 'bias', ytitle ='decay of latent(counts/s)', title ='Ch2 ', xrange=[400,650], yrange=[.02,.14]
oplot, bias2, decay_32s, psym = 2,thick=3
;legend, ['corr. cernox 30.1K', 'corr cernox 32.0K'], psym=[4,2],/right,/bottom


ps_end, /png

end
