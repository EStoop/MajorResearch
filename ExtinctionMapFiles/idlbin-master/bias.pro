pro bias

!p.multi = [0,1,3]

;ps_start, filename = '/Users/jkrick/iwic/allbiases.ps'
ps_open, filename= '/Users/jkrick/iwic/allbiases.ps',/portrait,/square,/color
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

;bias 450, ch1, corrected cernox = 32.1 exptime = 12 & 100,vgg1 = -3.75
;vdduc = -3.5
mean12 = [0.048765,0.0497592]
noisy12 = [114,106]
latentstrength = [0.00264021,0.00266504]
mean100=[0.00754524,0.007676]
noisy100 = [431,395]
vreset = [-3.35, -3.4]
!P.charsize = 2
plot, vreset, mean100, psym = 4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vreset', ytitle = 'mean(stddev) Mjy/sr', title = 'ch1, b450, corr. cernox=32.1, vgg1=-3.75, vdduc=-3.5', yrange=[.001,.1], /ylog
oplot, vreset, mean12, psym =2
legend, ['100s', '12s'], psym=[4,2],/center,/bottom

plot, vreset, noisy100, psym=4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vreset', ytitle = 'number of 10sigma pixels', title = 'ch1, b450, corr. cernox=32.1, vgg1=-3.75, vdduc=-3.5', yrange=[1,1000], /ylog
oplot, vreset, noisy12, psym = 2

plot, vreset, latentstrength, psym=4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vreset', ytitle = 'fractional latent strength', title = 'ch1, b450, corr. cernox=32.1, vgg1=-3.75, vdduc=-3.5', yrange=[.001, .004]

;bias 450, ch1, corrected cernox = 32.1, exptime = 12 & 100, vdduc = -3.5,
;vreset = [-3.5]
vgg1 = [-3.55,-4.0]
mean12 = [0.040015347,0.05897]
noisy12 =[150,68]
latentstrength = [0.00297933,0.00186678]
mean100 =[0.006395,0.008079]
noisy100 = [530,324]

plot, vgg1, mean100, psym = 4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vgg1', ytitle = 'mean(stddev) Mjy/sr', title = 'ch1, b450, corr. cernox=32.1,vreset=-3.5, vdduc=-3.5', yrange=[.001,.1], /ylog, xrange=[-4.2,-3.4]
oplot, vgg1, mean12, psym =2
legend, ['100s', '12s'], psym=[4,2],/center,/bottom

plot, vgg1, noisy100, psym=4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vgg1', ytitle = 'number of 10sigma pixels', title = 'ch1, b450, corr. cernox=32.1, vreset=-3.5, vdduc=-3.5', yrange=[1,1000], /ylog,xrange=[-4.2,-3.4]
oplot, vgg1, noisy12, psym = 2

plot, vgg1, latentstrength, psym=4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vgg1', ytitle = 'fractional latent strength', title = 'ch1, b450, corr. cernox=32.1,vreset=-3.5, vdduc=-3.5', yrange=[.001, .004],xrange=[-4.2,-3.4]

;-------
;bias 450, ch2, corrected cernox = 32.1 exptime = 12 & 100,vgg1 = -3.7
;vdduc = -3.5
vreset = [-3.5,-3.35,-3.4]
mean12 = [0.04229,0.0421009,0.04233826]
noisy12 =[149,156,154]
latentstrength = [0.00181787,0.00264021,0.00266504]
mean100 =[0.010387,0.010485,0.010494]
noisy100 = [273,284,269]

plot, vreset, mean100, psym = 4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vreset', ytitle = 'mean(stddev) Mjy/sr', title = 'ch2, b450, corr. cernox=32.1, vgg1=-3.75, vdduc=-3.5', yrange=[.001,.1], /ylog, xrange=[-3.6,-3.2]
oplot, vreset, mean12, psym =2
legend, ['100s', '12s'], psym=[4,2],/right,/bottom

plot, vreset, noisy100, psym=4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vreset', ytitle = 'number of 10sigma pixels', title = 'ch2, b450, corr. cernox=32.1, vgg1=-3.75, vdduc=-3.5', yrange=[100,1000], /ylog, xrange=[-3.6,-3.2]
oplot, vreset, noisy12, psym = 2

plot, vreset, latentstrength, psym=4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vreset', ytitle = 'fractional latent strength', title = 'ch2, b450, corr. cernox=32.1, vgg1=-3.75, vdduc=-3.5', yrange=[.001, .003], xrange=[-3.6,-3.2]

;bias 450, ch2, corrected cernox = 32.1, exptime = 12 & 100, vdduc = -3.5,
;vreset = [-3.5]
vgg1 = [-3.7,-3.5,-3.95]
mean12 =[0.04229,0.040578,0.04597]
noisy12 = [149,144,155]
latentstrength = [0.00181787,0.00297933,0.00186678]
mean100 =[0.010387,0.010399,0.01075]
noisy100 = [273,265,273]

plot, vgg1, mean100, psym = 4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vgg1', ytitle = 'mean(stddev) Mjy/sr', title = 'ch2, b450, corr. cernox=32.1,vreset=-3.5, vdduc=-3.5', yrange=[.001,.1], /ylog, xrange=[-4.2,-3.4]
oplot, vgg1, mean12, psym =2
legend, ['100s', '12s'], psym=[4,2],/left,/bottom

plot, vgg1, noisy100, psym=4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vgg1', ytitle = 'number of 10sigma pixels', title = 'ch2, b450, corr. cernox=32.1, vreset=-3.5, vdduc=-3.5', yrange=[10,1000], /ylog,xrange=[-4.2,-3.4]
oplot, vgg1, noisy12, psym = 2

plot, vgg1, latentstrength, psym=4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vgg1', ytitle = 'fractional latent strength', title = 'ch2, b450, corr. cernox=32.1,vreset=-3.5, vdduc=-3.5', yrange=[.001, .003],xrange=[-4.2,-3.4]

;bias 450, ch2, corrected cernox = 32.1, exptime = 12 & 100, 
;vreset = [-3.5], vgg1 = -3.7

vdduc = [-3.5,-3.1]
mean12 = [0.04229,0.041539]
mean100 = [0.010387,0.0103689]
noisy12 = [149,187]
noisy100 = [273,325]
latentstrength = [0.00181787,0.00235575]


plot, vdduc, mean100, psym = 4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vdduc', ytitle = 'mean(stddev) Mjy/sr', title = 'ch2, b450, corr. cernox=32.1,vreset=-3.5, vgg1=-3.7', yrange=[.001,.1], /ylog, xrange=[-3.6,-3.0]
oplot, vdduc, mean12, psym =2
legend, ['100s', '12s'], psym=[4,2],/center,/bottom

plot, vdduc, noisy100, psym=4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vdduc', ytitle = 'number of 10sigma pixels', title = 'ch2, b450, corr. cernox=32.1, vreset=-3.5, vgg1=-3.7', yrange=[10,1000], /ylog,xrange=[-3.6,-3.0]
oplot, vdduc, noisy12, psym = 2

plot, vdduc, latentstrength, psym=4, thick=3, charthick=3, xthick=3, ythick=3, xtitle = 'vdduc', ytitle = 'fractional latent strength', title = 'ch2, b450, corr. cernox=32.1,vreset=-3.5, vgg1=-3.7', yrange=[.001, .003],xrange=[-3.6,-3.0]

ps_close, /noprint,/noid
;ps_end

end
