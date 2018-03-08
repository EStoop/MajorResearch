pro bias1
ps_start, filename='/Users/jkrick/iwic/noisevsbias1.ps';,/portrait,/square,/color

!P.multi=[0,3,2]


;---------------------
;ch1 0.4s 30k
bias_ch1_0p4s_30k = [400,400,450,450,450,450,450,450,550,750,750]
vreset_ch1_0p4s_30k =[-3.5,-3.5,-3.5,-3.5,-3.5,-3.4,-3.5,-3.6,-3.5,-3.5,-3.6]
vdduc_ch1_0p4s_30k = [-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.1]
vgg1_ch1_0p4s_30k = [-3.75,-3.65,-3.55,-3.65,-3.65,-3.65,-3.75,-3.65,-3.75,-4,-3.75]
mean_ch1_0p4s_30k = [1.8197,1.58985,1.4762,1.5701,1.55598,1.53718,1.6533,1.5657,1.6491,1.8851,1.61431]
noisy_ch1_0p4s_30k = [1,1,4,2,2,3,4,2,6,25,55]
latent_ch1_0p4s_30k = [0.00205113,0.00180663,0.00260332,0.00222764,0.00224228,0.00235581,0.00161768,0.00184746,0.00172491,0.00136504,0.00157485]

resetgood_ch1_0p4s_30k = where(bias_ch1_0p4s_30k eq 450 and vdduc_ch1_0p4s_30k eq -3.5 and vgg1_ch1_0p4s_30k eq -3.65)
dducgood_ch1_0p4s_30k = where(bias_ch1_0p4s_30k eq 450 and vgg1_ch1_0p4s_30k eq -3.65 and vreset_ch1_0p4s_30k eq -3.5)
vggood_ch1_0p4s_30k = where(bias_ch1_0p4s_30k eq 450 and vdduc_ch1_0p4s_30k eq -3.5 and vreset_ch1_0p4s_30k eq -3.5)

;ch1 0.4s 32k

bias_ch1_0p4s_32k =[450,450,450,450,450,450,450,750,750,750,750]
vreset_ch1_0p4s_32k =[-3.35,-3.4,-3.4,-3.5,-3.5,-3.5,-3.6,-3.5,-3.5,-3.6,-3.6]
vdduc_ch1_0p4s_32k =[-3.5,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5,-3.1,-3.5]
vgg1_ch1_0p4s_32k=[-3.75,-3.75,-3.7,-3.55,-3.75,-4,-3.75,-3.55,-4,-3.75,-3.75]
mean_ch1_0p4s_32k=[1.84901,1.86707,1.75783,1.57114,1.87693,2.49029,1.89255,1.58897,2.07377,1.72553,1.7564]
noisy_ch1_0p4s_32k=[3,3,3,5,3,4,3,114,56,110,85]
latent_ch1_0p4s_32k =[0.00264021,0.00266504,0.00257865,0.00297933,0.00235575,0.00186678,0.0018307,0.00275751,0.00243414,0.00237035,0.00236514]

resetgood_ch1_0p4s_32k = where(bias_ch1_0p4s_32k eq 450 and vdduc_ch1_0p4s_32k eq -3.5 and vgg1_ch1_0p4s_32k eq -3.75)
dducgood_ch1_0p4s_32k = where(bias_ch1_0p4s_32k eq 450 and vgg1_ch1_0p4s_32k eq -3.75 and vreset_ch1_0p4s_32k eq -3.5)
vggood_ch1_0p4s_32k = where(bias_ch1_0p4s_32k eq 450 and vdduc_ch1_0p4s_32k eq -3.5 and vreset_ch1_0p4s_32k eq -3.5)



;---------------------

;ch1 12s 30k
bias_ch1_12s_30k = [750,450,750,400,550,450,400,450,450,450,450]
vreset_ch1_12s_30k =[-3.6,-3.6,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.4]
vdduc_ch1_12s_30k = [-3.1,-3.5,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5,-3.5]
vgg1_ch1_12s_30k = [-3.75,-3.65,-4,-3.75,-3.75,-3.75,-3.65,-3.65,-3.65,-3.55,-3.65]
mean_ch1_12s_30k = [0.046956,0.044566,0.054018,0.05078,0.04654,0.047469,0.045485,0.043763,0.04359,0.038882,0.042194]
noisy_ch1_12s_30k = [1510,59,600,15,131,85,37,53,64,31,61]
latent_ch1_12s_30k = [0.00157485,0.00184746,0.00136504,0.00205113,0.00172491,0.00161768,0.00180663,0.00222764,0.00224228,0.00260332,0.00235581]

resetgood_ch1_12s_30k = where(bias_ch1_12s_30k eq 450 and vdduc_ch1_12s_30k eq -3.5 and vgg1_ch1_12s_30k eq -3.65)
dducgood_ch1_12s_30k = where(bias_ch1_12s_30k eq 450 and vgg1_ch1_12s_30k eq -3.65 and vreset_ch1_12s_30k eq -3.5)
vggood_ch1_12s_30k = where(bias_ch1_12s_30k eq 450 and vdduc_ch1_12s_30k eq -3.5 and vreset_ch1_12s_30k eq -3.5)

;ch1 12s 32K
bias_ch1_12s_32k =[450,750,750,450,750,450,450,750,450,450,450]
vreset_ch1_12s_32k =[-3.6,-3.6,-3.6,-3.5,-3.5,-3.5,-3.5,-3.5,-3.4,-3.4,-3.35]
vdduc_ch1_12s_32k =[-3.5,-3.5,-3.1,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5,-3.5]
vgg1_ch1_12s_32k=[-3.75,-3.75,-3.75,-4,-4,-3.75,-3.55,-3.55,-3.75,-3.7,-3.75]
mean_ch1_12s_32k=[0.051127,0.05404,0.05775,0.05897,0.05831,0.05039,0.040015347,0.050633,0.0497592,0.046519,0.048765]
noisy_ch1_12s_32k=[87,2598,3745,51,1317,116,107,4545,94,98,100]
latent_ch1_12s_32k =[0.0018307,0.00236514,0.00237035,0.00186678,0.00243414,0.00235575,0.00297933,0.00275751,0.00266504,0.00257865,0.00264021]

resetgood_ch1_12s_32k = where(bias_ch1_12s_32k eq 450 and vdduc_ch1_12s_32k eq -3.5 and vgg1_ch1_12s_32k eq -3.75)
dducgood_ch1_12s_32k = where(bias_ch1_12s_32k eq 450 and vgg1_ch1_12s_32k eq -3.75 and vreset_ch1_12s_32k eq -3.5)
vggood_ch1_12s_32k = where(bias_ch1_12s_32k eq 450 and vdduc_ch1_12s_32k eq -3.5 and vreset_ch1_12s_32k eq -3.5)

;---------------------

;ch1 100s 30k
bias_ch1_100s_30k = [750,450,750,400,550,450,400,450,450,450,450]
vreset_ch1_100s_30k = [-3.6,-3.6,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.4]
vdduc_ch1_100s_30k = [-3.1,-3.5,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5,-3.5]
vgg1_ch1_100s_30k = [-3.75,-3.65,-4,-3.75,-3.75,-3.75,-3.65,-3.65,-3.65,-3.55,-3.65]
mean_ch1_100s_30k = [0.01342,0.006956,0.009518,0.007386,0.00753,0.00743,0.00698,0.006961,0.00694,0.0063999,0.006775]
noisy_ch1_100s_30k = [5522,207,2601,58,434,271,101,184,208,194,212]
latent_ch1_100s_30k = [0.00157485,0.00184746,0.00136504,0.00205113,0.00172491,0.00161768,0.00180663,0.00222764,0.00224228,0.00260332,0.00235581]
noise4x_ch1_100s_30k=[0.118043,0.00276443,0.051382,0.00132706,0.00782853,0.00333561,0.0014235,0.00265115,0.00272027,0.00321454,0.00298852]

resetgood_ch1_100s_30k = where(bias_ch1_100s_30k eq 450 and vdduc_ch1_100s_30k eq -3.5 and vgg1_ch1_100s_30k eq -3.65)
dducgood_ch1_100s_30k = where(bias_ch1_100s_30k eq 450 and vgg1_ch1_100s_30k eq -3.65 and vreset_ch1_100s_30k eq -3.5)
vggood_ch1_100s_30k = where(bias_ch1_100s_30k eq 450 and vdduc_ch1_100s_30k eq -3.5 and vreset_ch1_100s_30k eq -3.5)

;ch1 100s 32K
bias_ch1_100s_32k = [450,750,750,450,750,450,450,750,450,450,450]
vreset_ch1_100s_32k = [-3.6,-3.6,-3.6,-3.5,-3.5,-3.5,-3.5,-3.5,-3.4,-3.4,-3.35]
vdduc_ch1_100s_32k =[-3.5,-3.5,-3.1,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5,-3.5]
vgg1_ch1_100s_32k = [-3.75,-3.75,-3.75,-4,-4,-3.75,-3.55,-3.55,-3.75,-3.7,-3.75]
mean_ch1_100s_32k= [0.00765,0.02077,0.033288,0.008079,0.01219,0.007679,0.006395,0.028852,0.007676,0.007244,0.00754524]
noisy_ch1_100s_32k=[281,7532,8229,238,5699,399,368,8349,325,352,350]
latent_ch1_100s_32k = [0.0018307,0.00236514,0.00237035,0.00186678,0.00243414,0.00235575,0.00297933,0.00275751,0.00266504,0.00257865,0.00264021]
noise4x_ch1_100s_32k=[0.0050177,0.202771,0.253307,0.00491885,0.138904,0.00730718,0.00717723,0.258647,0.0054736,0.00603554,0.00578842]

resetgood_ch1_100s_32k = where(bias_ch1_100s_32k eq 450 and vdduc_ch1_100s_32k eq -3.5 and vgg1_ch1_100s_32k eq -3.75)
dducgood_ch1_100s_32k = where(bias_ch1_100s_32k eq 450 and vgg1_ch1_100s_32k eq -3.75 and vreset_ch1_100s_32k eq -3.5)
vggood_ch1_100s_32k = where(bias_ch1_100s_32k eq 450 and vdduc_ch1_100s_32k eq -3.5 and vreset_ch1_100s_32k eq -3.5)

;------------------------------------------------------------------------------------------------------
;0.4s plots

plot, vreset_ch1_0p4s_30k(resetgood_ch1_0p4s_30k), mean_ch1_0p4s_30k(resetgood_ch1_0p4s_30k), psym = 4, xrange=[-3.7,-3.3], yrange=[1.0,3.0], xtitle ='Vreset', ytitle ='mean(stddev) MJy/sr', title = '0p4s Ch1 b450 g365(30k) g375(32k) d35',charsize=1.5
oplot, vreset_ch1_0p4s_30k(resetgood_ch1_0p4s_30k), mean_ch1_0p4s_30k(resetgood_ch1_0p4s_30k)
oplot,  vreset_ch1_0p4s_32k(resetgood_ch1_0p4s_32k), mean_ch1_0p4s_32k(resetgood_ch1_0p4s_32k)
oplot,  vreset_ch1_0p4s_32k(resetgood_ch1_0p4s_32k), mean_ch1_0p4s_32k(resetgood_ch1_0p4s_32k),psym = 2

plot, vreset_ch1_0p4s_30k(resetgood_ch1_0p4s_30k), noisy_ch1_0p4s_30k(resetgood_ch1_0p4s_30k), psym = 4, xrange=[-3.7,-3.3], yrange=[0,10], xtitle ='Vreset', ytitle ='10sigma noisypixels', title = '0p4s Ch1 b450 g365(30k) g375(32k) d35', charsize=1.5
oplot, vreset_ch1_0p4s_30k(resetgood_ch1_0p4s_30k), noisy_ch1_0p4s_30k(resetgood_ch1_0p4s_30k)
oplot,  vreset_ch1_0p4s_32k(resetgood_ch1_0p4s_32k), noisy_ch1_0p4s_32k(resetgood_ch1_0p4s_32k)
oplot,  vreset_ch1_0p4s_32k(resetgood_ch1_0p4s_32k), noisy_ch1_0p4s_32k(resetgood_ch1_0p4s_32k),psym = 2

plot, vreset_ch1_0p4s_30k(resetgood_ch1_0p4s_30k), latent_ch1_0p4s_30k(resetgood_ch1_0p4s_30k), psym = 4, xrange=[-3.7,-3.3], xtitle ='Vreset', ytitle ='latent strength', title = '0p4s Ch1 b450 g365(30k) g375(32k) d35', yrange=[.001,.004],charsize=1.5
oplot, vreset_ch1_0p4s_30k(resetgood_ch1_0p4s_30k), latent_ch1_0p4s_30k(resetgood_ch1_0p4s_30k)
oplot,  vreset_ch1_0p4s_32k(resetgood_ch1_0p4s_32k), latent_ch1_0p4s_32k(resetgood_ch1_0p4s_32k)
oplot,  vreset_ch1_0p4s_32k(resetgood_ch1_0p4s_32k), latent_ch1_0p4s_32k(resetgood_ch1_0p4s_32k),psym = 2


;plot, vdduc_ch1_0p4s_32k(dducgood_ch1_0p4s_32k), mean_ch1_0p4s_32k(dducgood_ch1_0p4s_32k), psym = 2, xrange=[-3.6,-3.0],  yrange=[1.5,2.0], xtitle ='Vdduc', ytitle ='mean(stddev) MJy/sr', title = '0p4s Ch1 b450 g37 r35', charsize=1.5
;oplot, vdduc_ch1_0p4s_32k(dducgood_ch1_0p4s_32k), mean_ch1_0p4s_32k(dducgood_ch1_0p4s_32k)
;legend, ['30k', '32k'], psym=[4,2],/center,/bottom

;plot, vdduc_ch1_0p4s_32k(dducgood_ch1_0p4s_32k), noisy_ch1_0p4s_32k(dducgood_ch1_0p4s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[0,20],xtitle ='Vdduc', ytitle ='10sigma noisypixels', title = '0p4s Ch1 b450 g37 r35', charsize=1.5
;oplot, vdduc_ch1_0p4s_32k(dducgood_ch1_0p4s_32k), noisy_ch1_0p4s_32k(dducgood_ch1_0p4s_32k)
;legend, ['30k', '32k'], psym=[4,2],/center,/bottom

;plot, vdduc_ch1_0p4s_32k(dducgood_ch1_0p4s_32k), latent_ch1_0p4s_32k(dducgood_ch1_0p4s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[.001,.004], xtitle ='Vdduc', ytitle ='latent strength', title = '0p4s Ch1 b450 g37 r35', charsize=1.5
;oplot, vdduc_ch1_0p4s_32k(dducgood_ch1_0p4s_32k), latent_ch1_0p4s_32k(dducgood_ch1_0p4s_32k)
;legend, ['30k', '32k'], psym=[4,2],/center,/bottom



plot, vgg1_ch1_0p4s_32k(vggood_ch1_0p4s_32k), mean_ch1_0p4s_32k(vggood_ch1_0p4s_32k), psym = 2, xtitle ='Vgg1', ytitle ='mean(stddev) MJy/sr', title = '0p4s Ch1 b450 r35 d35', xrange=[-4.1,-3.4],   yrange=[1.0,3.0], charsize=1.5
oplot, vgg1_ch1_0p4s_32k(vggood_ch1_0p4s_32k), mean_ch1_0p4s_32k(vggood_ch1_0p4s_32k)
oplot,  vgg1_ch1_0p4s_30k(vggood_ch1_0p4s_30k), mean_ch1_0p4s_30k(vggood_ch1_0p4s_30k), psym = 4
oplot,  vgg1_ch1_0p4s_30k(vggood_ch1_0p4s_30k), mean_ch1_0p4s_30k(vggood_ch1_0p4s_30k)
legend, ['30k', '32k'], psym=[4,2],/center,/top

plot, vgg1_ch1_0p4s_32k(vggood_ch1_0p4s_32k), noisy_ch1_0p4s_32k(vggood_ch1_0p4s_32k), psym = 2, xtitle ='Vgg1', ytitle ='10sigma noisypixels', title = '0p4s Ch1 b450 r35 d35', xrange=[-4.1,-3.4],yrange=[0,10],charsize=1.5
oplot, vgg1_ch1_0p4s_32k(vggood_ch1_0p4s_32k), noisy_ch1_0p4s_32k(vggood_ch1_0p4s_32k)
oplot,  vgg1_ch1_0p4s_30k(vggood_ch1_0p4s_30k), noisy_ch1_0p4s_30k(vggood_ch1_0p4s_30k), psym = 4
oplot,  vgg1_ch1_0p4s_30k(vggood_ch1_0p4s_30k), noisy_ch1_0p4s_30k(vggood_ch1_0p4s_30k)
legend, ['30k', '32k'], psym=[4,2],/center,/top

plot, vgg1_ch1_0p4s_32k(vggood_ch1_0p4s_32k), latent_ch1_0p4s_32k(vggood_ch1_0p4s_32k), psym = 2, xtitle ='Vgg1', ytitle ='latent strength', title = '0p4s Ch1 b450 r35 d35', xrange=[-4.1,-3.4], yrange=[.001,.004], charsize=1.5
oplot, vgg1_ch1_0p4s_32k(vggood_ch1_0p4s_32k), latent_ch1_0p4s_32k(vggood_ch1_0p4s_32k)
oplot,  vgg1_ch1_0p4s_30k(vggood_ch1_0p4s_30k), latent_ch1_0p4s_30k(vggood_ch1_0p4s_30k), psym = 4
oplot,  vgg1_ch1_0p4s_30k(vggood_ch1_0p4s_30k), latent_ch1_0p4s_30k(vggood_ch1_0p4s_30k)
legend, ['30k', '32k'], psym=[4,2],/center,/top


;------------------------------------------------------------------------------------------------------
;12s plots

plot, vreset_ch1_12s_30k(resetgood_ch1_12s_30k), mean_ch1_12s_30k(resetgood_ch1_12s_30k), psym = 4, xrange=[-3.7,-3.3], yrange=[0.03,0.08], xtitle ='Vreset', ytitle ='mean(stddev) MJy/sr', title = '12s Ch1 b450 g365(30k) g375(32k) d35',charsize=1.5
oplot, vreset_ch1_12s_30k(resetgood_ch1_12s_30k), mean_ch1_12s_30k(resetgood_ch1_12s_30k)
oplot,  vreset_ch1_12s_32k(resetgood_ch1_12s_32k), mean_ch1_12s_32k(resetgood_ch1_12s_32k)
oplot,  vreset_ch1_12s_32k(resetgood_ch1_12s_32k), mean_ch1_12s_32k(resetgood_ch1_12s_32k),psym = 2

plot, vreset_ch1_12s_30k(resetgood_ch1_12s_30k), noisy_ch1_12s_30k(resetgood_ch1_12s_30k), psym = 4, xrange=[-3.7,-3.3], xtitle ='Vreset', ytitle ='10sigma noisypixels', title = '12s Ch1 b450 g365(30k) g375(32k) d35', charsize=1.5
oplot, vreset_ch1_12s_30k(resetgood_ch1_12s_30k), noisy_ch1_12s_30k(resetgood_ch1_12s_30k)
oplot,  vreset_ch1_12s_32k(resetgood_ch1_12s_32k), noisy_ch1_12s_32k(resetgood_ch1_12s_32k)
oplot,  vreset_ch1_12s_32k(resetgood_ch1_12s_32k), noisy_ch1_12s_32k(resetgood_ch1_12s_32k),psym = 2

plot, vreset_ch1_12s_30k(resetgood_ch1_12s_30k), latent_ch1_12s_30k(resetgood_ch1_12s_30k), psym = 4, xrange=[-3.7,-3.3], xtitle ='Vreset', ytitle ='latent strength', title = '12s Ch1 b450 g365(30k) g375(32k) d35', yrange=[.001,.004],charsize=1.5
oplot, vreset_ch1_12s_30k(resetgood_ch1_12s_30k), latent_ch1_12s_30k(resetgood_ch1_12s_30k)
oplot,  vreset_ch1_12s_32k(resetgood_ch1_12s_32k), latent_ch1_12s_32k(resetgood_ch1_12s_32k)
oplot,  vreset_ch1_12s_32k(resetgood_ch1_12s_32k), latent_ch1_12s_32k(resetgood_ch1_12s_32k),psym = 2


;plot, vdduc_ch1_12s_32k(dducgood_ch1_12s_32k), mean_ch1_12s_32k(dducgood_ch1_12s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[.03,.05], xtitle ='Vdduc', ytitle ='mean(stddev) MJy/sr', title = '12s Ch1 b450 g37 r35', charsize=1.5
;oplot, vdduc_ch1_12s_32k(dducgood_ch1_12s_32k), mean_ch1_12s_32k(dducgood_ch1_12s_32k)
;legend, ['30k', '32k'], psym=[4,2],/center,/bottom

;plot, vdduc_ch1_12s_32k(dducgood_ch1_12s_32k), noisy_ch1_12s_32k(dducgood_ch1_12s_32k), psym = 2, xrange=[-3.6,-3.0], xtitle ='Vdduc', ytitle ='10sigma noisypixels', title = '12s Ch1 b450 g37 r35', charsize=1.5
;oplot, vdduc_ch1_12s_32k(dducgood_ch1_12s_32k), noisy_ch1_12s_32k(dducgood_ch1_12s_32k)
;legend, ['30k', '32k'], psym=[4,2],/center,/bottom

;plot, vdduc_ch1_12s_32k(dducgood_ch1_12s_32k), latent_ch1_12s_32k(dducgood_ch1_12s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[.001,.004], xtitle ='Vdduc', ytitle ='latent strength', title = '12s Ch1 b450 g37 r35', charsize=1.5
;oplot, vdduc_ch1_12s_32k(dducgood_ch1_12s_32k), latent_ch1_12s_32k(dducgood_ch1_12s_32k)
;legend, ['30k', '32k'], psym=[4,2],/center,/bottom



plot, vgg1_ch1_12s_32k(vggood_ch1_12s_32k), mean_ch1_12s_32k(vggood_ch1_12s_32k), psym = 2, xtitle ='Vgg1', ytitle ='mean(stddev) MJy/sr', title = '12s Ch1 b450 r35 d35', xrange=[-4.1,-3.4],  yrange=[0.04,0.08],charsize=1.5
oplot, vgg1_ch1_12s_32k(vggood_ch1_12s_32k), mean_ch1_12s_32k(vggood_ch1_12s_32k)
oplot,  vgg1_ch1_12s_30k(vggood_ch1_12s_30k), mean_ch1_12s_30k(vggood_ch1_12s_30k), psym = 4
oplot,  vgg1_ch1_12s_30k(vggood_ch1_12s_30k), mean_ch1_12s_30k(vggood_ch1_12s_30k)
legend, ['30k', '32k'], psym=[4,2],/center,/top

plot, vgg1_ch1_12s_32k(vggood_ch1_12s_32k), noisy_ch1_12s_32k(vggood_ch1_12s_32k), psym = 2, xtitle ='Vgg1', ytitle ='10sigma noisypixels', title = '12s Ch1 b450 r35 d35', xrange=[-4.1,-3.4], charsize=1.5
oplot, vgg1_ch1_12s_32k(vggood_ch1_12s_32k), noisy_ch1_12s_32k(vggood_ch1_12s_32k)
oplot,  vgg1_ch1_12s_30k(vggood_ch1_12s_30k), noisy_ch1_12s_30k(vggood_ch1_12s_30k), psym = 4
oplot,  vgg1_ch1_12s_30k(vggood_ch1_12s_30k), noisy_ch1_12s_30k(vggood_ch1_12s_30k)
legend, ['30k', '32k'], psym=[4,2],/center,/top

plot, vgg1_ch1_12s_32k(vggood_ch1_12s_32k), latent_ch1_12s_32k(vggood_ch1_12s_32k), psym = 2, xtitle ='Vgg1', ytitle ='latent strength', title = '12s Ch1 b450 r35 d35', xrange=[-4.1,-3.4], yrange=[.001,.004], charsize=1.5
oplot, vgg1_ch1_12s_32k(vggood_ch1_12s_32k), latent_ch1_12s_32k(vggood_ch1_12s_32k)
oplot,  vgg1_ch1_12s_30k(vggood_ch1_12s_30k), latent_ch1_12s_30k(vggood_ch1_12s_30k), psym = 4
oplot,  vgg1_ch1_12s_30k(vggood_ch1_12s_30k), latent_ch1_12s_30k(vggood_ch1_12s_30k)
legend, ['30k', '32k'], psym=[4,2],/center,/top

;------------------------------------------------------------------------------------------------------
;100s plots

plot, vreset_ch1_100s_30k(resetgood_ch1_100s_30k), mean_ch1_100s_30k(resetgood_ch1_100s_30k), psym = 4, xrange=[-3.7,-3.3],  yrange=[.005,.009],xtitle ='Vreset', ytitle ='mean(stddev) MJy/sr', title = '100s Ch1 b450 g365(30k) g375(32k) d35',  charsize=1.5
oplot, vreset_ch1_100s_30k(resetgood_ch1_100s_30k), mean_ch1_100s_30k(resetgood_ch1_100s_30k)
oplot,  vreset_ch1_100s_32k(resetgood_ch1_100s_32k), mean_ch1_100s_32k(resetgood_ch1_100s_32k)
oplot,  vreset_ch1_100s_32k(resetgood_ch1_100s_32k), mean_ch1_100s_32k(resetgood_ch1_100s_32k),psym = 2

plot, vreset_ch1_100s_30k(resetgood_ch1_100s_30k), noise4x_ch1_100s_30k(resetgood_ch1_100s_30k), psym = 4, xrange=[-3.7,-3.3], yrange=[0,.008],xtitle ='Vreset', ytitle ='fraction noisypixels', title = '100s Ch1 b450 g365(30k) g375(32k) d35', charsize=1.5
oplot, vreset_ch1_100s_30k(resetgood_ch1_100s_30k), noise4x_ch1_100s_30k(resetgood_ch1_100s_30k)
oplot,  vreset_ch1_100s_32k(resetgood_ch1_100s_32k), noise4x_ch1_100s_32k(resetgood_ch1_100s_32k)
oplot,  vreset_ch1_100s_32k(resetgood_ch1_100s_32k), noise4x_ch1_100s_32k(resetgood_ch1_100s_32k),psym = 2

plot, vreset_ch1_100s_30k(resetgood_ch1_100s_30k), latent_ch1_100s_30k(resetgood_ch1_100s_30k), psym = 4, xrange=[-3.7,-3.3], xtitle ='Vreset', ytitle ='latent strength', title = '100s Ch1 b450 g365(30k) g37(32k) d35', yrange=[.001,.004],charsize=1.5
oplot, vreset_ch1_100s_30k(resetgood_ch1_100s_30k), latent_ch1_100s_30k(resetgood_ch1_100s_30k)
oplot,  vreset_ch1_100s_32k(resetgood_ch1_100s_32k), latent_ch1_100s_32k(resetgood_ch1_100s_32k)
oplot,  vreset_ch1_100s_32k(resetgood_ch1_100s_32k), latent_ch1_100s_32k(resetgood_ch1_100s_32k),psym = 2


;plot, vdduc_ch1_100s_32k(dducgood_ch1_100s_32k), mean_ch1_100s_32k(dducgood_ch1_100s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[.008,.012], xtitle ='Vdduc', ytitle ='mean(stddev) MJy/sr', title = '100s Ch1 b450 g37 r35', charsize=1.5
;oplot, vdduc_ch1_100s_32k(dducgood_ch1_100s_32k), mean_ch1_100s_32k(dducgood_ch1_100s_32k)
;legend, ['30k', '32k'], psym=[4,2],/center,/bottom

;plot, vdduc_ch1_100s_32k(dducgood_ch1_100s_32k), noisy_ch1_100s_32k(dducgood_ch1_100s_32k), psym = 2, xrange=[-3.6,-3.0], xtitle ='Vdduc', ytitle ='fraction noisypixels', title = '100s Ch1 b450 g37 r35', charsize=1.5
;oplot, vdduc_ch1_100s_32k(dducgood_ch1_100s_32k), noisy_ch1_100s_32k(dducgood_ch1_100s_32k)
;legend, ['30k', '32k'], psym=[4,2],/center,/bottom

;plot, vdduc_ch1_100s_32k(dducgood_ch1_100s_32k), latent_ch1_100s_32k(dducgood_ch1_100s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[.001,.004], xtitle ='Vdduc', ytitle ='latent strength', title = '100s Ch1 b450 g37 r35', charsize=1.5
;oplot, vdduc_ch1_100s_32k(dducgood_ch1_100s_32k), latent_ch1_100s_32k(dducgood_ch1_100s_32k)
;legend, ['30k', '32k'], psym=[4,2],/center,/bottom



plot, vgg1_ch1_100s_32k(vggood_ch1_100s_32k), mean_ch1_100s_32k(vggood_ch1_100s_32k), psym = 2, xtitle ='Vgg1', ytitle ='mean(stddev) MJy/sr', title = '100s Ch1 b450 r35 d35', xrange=[-4.1,-3.4], yrange=[.005,.009],charsize=1.5
oplot, vgg1_ch1_100s_32k(vggood_ch1_100s_32k), mean_ch1_100s_32k(vggood_ch1_100s_32k)
oplot,  vgg1_ch1_100s_30k(vggood_ch1_100s_30k), mean_ch1_100s_30k(vggood_ch1_100s_30k), psym = 4
oplot,  vgg1_ch1_100s_30k(vggood_ch1_100s_30k), mean_ch1_100s_30k(vggood_ch1_100s_30k)
legend, ['30k', '32k'], psym=[4,2],/center,/top

plot, vgg1_ch1_100s_32k(vggood_ch1_100s_32k), noise4x_ch1_100s_32k(vggood_ch1_100s_32k), psym = 2, xtitle ='Vgg1', ytitle ='fraction noisypixels', title = '100s Ch1 b450 r35 d35', xrange=[-4.1,-3.4], charsize=1.5
oplot, vgg1_ch1_100s_32k(vggood_ch1_100s_32k), noise4x_ch1_100s_32k(vggood_ch1_100s_32k)
oplot,  vgg1_ch1_100s_30k(vggood_ch1_100s_30k), noise4x_ch1_100s_30k(vggood_ch1_100s_30k), psym = 4
oplot,  vgg1_ch1_100s_30k(vggood_ch1_100s_30k), noise4x_ch1_100s_30k(vggood_ch1_100s_30k)
legend, ['30k', '32k'], psym=[4,2],/center,/bottom

plot, vgg1_ch1_100s_32k(vggood_ch1_100s_32k), latent_ch1_100s_32k(vggood_ch1_100s_32k), psym = 2, xtitle ='Vgg1', ytitle ='latent strength', title = '100s Ch1 b450 r35 d35', xrange=[-4.1,-3.4], yrange=[.001,.004], charsize=1.5
oplot, vgg1_ch1_100s_32k(vggood_ch1_100s_32k), latent_ch1_100s_32k(vggood_ch1_100s_32k)
oplot,  vgg1_ch1_100s_30k(vggood_ch1_100s_30k), latent_ch1_100s_30k(vggood_ch1_100s_30k), psym = 4
oplot,  vgg1_ch1_100s_30k(vggood_ch1_100s_30k), latent_ch1_100s_30k(vggood_ch1_100s_30k)
legend, ['30k', '32k'], psym=[4,2],/center,/top

ps_end;, /noprint,/noid

end
