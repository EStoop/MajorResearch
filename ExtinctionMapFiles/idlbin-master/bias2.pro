pro bias2
ps_start, filename='/Users/jkrick/iwic/noisevsbias2.ps';,/landscape;,/square,/color

!P.multi=[0,3,3]

;need to add noisy pixel 4x values, only have those for 100s
                                ;probably should add these by hand to
                                ;this sheet instead of trying to
                                ;resort, is only a few values.
;need to change y range on noisy pixel values for 100s



;---------------------
;ch2 0.4s 30k
bias_ch2_0p4s_30k = [450,450,450,450,450,450,500,550,550,600,650,650,450,600]
vreset_ch2_0p4s_30k =[-3.5,-3.5,-3.5,-3.5,-3.4,-3.5,-3.5,-3.5,-3.6,-3.6,-3.5,-3.5,-3.5,-3.5]
vdduc_ch2_0p4s_30k = [-3.5,-3.5,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5]
vgg1_ch2_0p4s_30k = [-3.5,-3.6,-3.95,-3.6,-3.6,-3.7,-3.7,-3.7,-3.6,-3.7,-3.7,-3.6,-3.7,-3.7]
mean_ch2_0p4s_30k = [1.6413906,1.6605,1.6848,1.67337,1.69301,1.64996,1.69424,1.7673,1.7957,1.77615,1.8508,1.82627,1.67454,1.80202]
noisy_ch2_0p4s_30k = [6,5,5,6,15,5,8,23,13,22,34,33,6,19]
latent_ch2_0p4s_30k = [0.00260332,0.00222764,0.00136504,0.00224228,0.00235581,0.00161768,0.00221165,0.00172491,0.00184746,0.00157485,0.00205113,0.00180663,0.00204855,0.00170674]

resetgood_ch2_0p4s_30k = where(bias_ch2_0p4s_30k eq 450 and vdduc_ch2_0p4s_30k eq -3.5 and vgg1_ch2_0p4s_30k eq -3.6)
dducgood_ch2_0p4s_30k = where(bias_ch2_0p4s_30k eq 450 and vgg1_ch2_0p4s_30k eq -3.7 and vreset_ch2_0p4s_30k eq -3.5)
vggood_ch2_0p4s_30k = where(bias_ch2_0p4s_30k eq 450 and vdduc_ch2_0p4s_30k eq -3.5 and vreset_ch2_0p4s_30k eq -3.5)

;ch2 0.4s 32k

bias_ch2_0p4s_32k =[450,450,450,450,450,600,600,600,600,600,450,500,600]
vreset_ch2_0p4s_32k =[-3.35,-3.4,-3.5,-3.5,-3.5,-3.4,-3.5,-3.5,-3.6,-3.6,-3.5,-3.5,-3.5]
vdduc_ch2_0p4s_32k =[-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5]
vgg1_ch2_0p4s_32k=[-3.7,-3.7,-3.5,-3.7,-3.95,-3.7,-3.5,-3.95,-3.7,-3.7,-3.7,-3.7,-3.7]
mean_ch2_0p4s_32k=[1.74157,1.73517,1.7322433,1.71482,1.779306,1.86604,1.855409,1.84879,1.8504,1.8611,1.7247,1.76378,1.85207]
noisy_ch2_0p4s_32k=[13,9,7,8,6,28,25,22,29,25,8,16,31]
latent_ch2_0p4s_32k =[0.00264021,0.00266504,0.00297933,0.00235575,0.00186678,0.00257865,0.00275751,0.00243414,0.00237035,0.00236514,0.00181787,0.00211641,0.00254605]

resetgood_ch2_0p4s_32k = where(bias_ch2_0p4s_32k eq 450 and vdduc_ch2_0p4s_32k eq -3.5 and vgg1_ch2_0p4s_32k eq -3.7)
dducgood_ch2_0p4s_32k = where(bias_ch2_0p4s_32k eq 450 and vgg1_ch2_0p4s_32k eq -3.7 and vreset_ch2_0p4s_32k eq -3.5)
vggood_ch2_0p4s_32k = where(bias_ch2_0p4s_32k eq 450 and vdduc_ch2_0p4s_32k eq -3.5 and vreset_ch2_0p4s_32k eq -3.5)



;---------------------

;ch2 12s 30k
bias_ch2_12s_30k = [450,450,450,450,450,500,550,550,600,650,650,450]
vreset_ch2_12s_30k =[-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.6,-3.6,-3.5,-3.5,-3.4]
vdduc_ch2_12s_30k = [-3.5,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5]
vgg1_ch2_12s_30k = [-3.5,-3.6,-3.95,-3.6,-3.7,-3.7,-3.7,-3.6,-3.7,-3.7,-3.6,-3.6]
mean_ch2_12s_30k = [0.039556,0.0406,0.04531,0.041156,0.041406,0.042004,0.04335,0.0438719,0.044964,0.04562,0.045195,0.041526]
noisy_ch2_12s_30k = [40,41,40,38,56,61,98,94,148,177,176,43]
latent_ch2_12s_30k = [0.00260332,0.00222764,0.00136504,0.00224228,0.00161768,0.00221165,0.00172491,0.00184746,0.00157485,0.00205113,0.00180663,0.00235581]

resetgood_ch2_12s_30k = where(bias_ch2_12s_30k eq 450 and vdduc_ch2_12s_30k eq -3.5 and vgg1_ch2_12s_30k eq -3.6)
dducgood_ch2_12s_30k = where(bias_ch2_12s_30k eq 450 and vgg1_ch2_12s_30k eq -3.7 and vreset_ch2_12s_30k eq -3.5)
vggood_ch2_12s_30k = where(bias_ch2_12s_30k eq 450 and vdduc_ch2_12s_30k eq -3.5 and vreset_ch2_12s_30k eq -3.5)

;ch2 12s 32K
bias_ch2_12s_32k =[450,500,600,450,450,450,450,450,450,600,600,600,600,600]
vreset_ch2_12s_32k =[-3.5,-3.5,-3.5,-3.35,-3.4,-3.5,-3.5,-3.5,-3.6,-3.4,-3.5,-3.5,-3.6,-3.6]
vdduc_ch2_12s_32k =[-3.5,-3.5,-3.5,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5,-3.5,-3.1,-3.5]
vgg1_ch2_12s_32k=[-3.7,-3.7,-3.7,-3.7,-3.7,-3.5,-3.7,-3.95,-3.7,-3.7,-3.5,-3.95,-3.7,-3.7]
mean_ch2_12s_32k=[0.04229,0.042928,0.045846,0.0421009,0.04233826,0.040578,0.041539,0.04597,0.043886,0.045344,0.04397,0.047223,0.04601,0.04593]
noisy_ch2_12s_32k=[87,114,220,85,84,81,103,82,558,232,219,209,282,207]
latent_ch2_12s_32k =[0.00181787,0.00211641,0.00254605,0.00264021,0.00266504,0.00297933,0.00235575,0.00186678,0.0018307,0.00257865,0.00275751,0.00243414,0.00237035,0.00236514]

resetgood_ch2_12s_32k = where(bias_ch2_12s_32k eq 450 and vdduc_ch2_12s_32k eq -3.5 and vgg1_ch2_12s_32k eq -3.7)
dducgood_ch2_12s_32k = where(bias_ch2_12s_32k eq 450 and vgg1_ch2_12s_32k eq -3.7 and vreset_ch2_12s_32k eq -3.5)
vggood_ch2_12s_32k = where(bias_ch2_12s_32k eq 450 and vdduc_ch2_12s_32k eq -3.5 and vreset_ch2_12s_32k eq -3.5)

;---------------------

;ch2 100s 30k
bias_ch2_100s_30k = [450,450,450,450,450,450,500,550,550,600,650,650]
vreset_ch2_100s_30k = [-3.5,-3.5,-3.5,-3.5,-3.4,-3.5,-3.5,-3.5,-3.6,-3.6,-3.5,-3.5]
vdduc_ch2_100s_30k = [-3.5,-3.5,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5]
vgg1_ch2_100s_30k = [-3.5,-3.6,-3.95,-3.6,-3.6,-3.7,-3.7,-3.7,-3.6,-3.7,-3.7,-3.6]
mean_ch2_100s_30k = [0.010124,0.010358,0.010758,0.01025,0.01046,0.01034,0.010966,0.011533,0.0114,0.012037,0.0125505,0.012571]
noisy_ch2_100s_30k = [97,100,95,98,101,122,125,183,180,329,411,375]
latent_ch2_100s_30k = [0.00260332,0.00222764,0.00136504,0.00224228,0.00235581,0.00161768,0.00221165,0.00172491,0.00184746,0.00157485,0.00205113,0.00180663]
noise4x_ch2_100s_30k=[0.00159075,0.00166968,0.00156472,0.00365443,0.00171726,0.00224172,0,0.00372114,0.0015977,0.00568311,0.00729459,0.00697771]

resetgood_ch2_100s_30k = where(bias_ch2_100s_30k eq 450 and vdduc_ch2_100s_30k eq -3.5 and vgg1_ch2_100s_30k eq -3.6)
dducgood_ch2_100s_30k = where(bias_ch2_100s_30k eq 450 and vgg1_ch2_100s_30k eq -3.7 and vreset_ch2_100s_30k eq -3.5)
vggood_ch2_100s_30k = where(bias_ch2_100s_30k eq 450 and vdduc_ch2_100s_30k eq -3.5 and vreset_ch2_100s_30k eq -3.5)

;ch2 100s 32K
bias_ch2_100s_32k = [450,450,450,450,450,450,600,600,600,600,600,450,500,600]
vreset_ch2_100s_32k = [-3.35,-3.4,-3.5,-3.5,-3.5,-3.6,-3.4,-3.5,-3.5,-3.6,-3.6,-3.5,-3.5,-3.5]
vdduc_ch2_100s_32k =[-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5,-3.5,-3.1,-3.5,-3.5,-3.5,-3.5]
vgg1_ch2_100s_32k = [-3.7,-3.7,-3.5,-3.7,-3.95,-3.7,-3.7,-3.5,-3.95,-3.7,-3.7,-3.7,-3.7,-3.7]
mean_ch2_100s_32k= [0.010485,0.010494,0.010399,0.0103689,0.01075,0.01057,0.012354,0.011984,0.01224,0.012248,0.012007,0.010387,0.010764,0.0122286]
noisy_ch2_100s_32k=[157,149,149,182,150,279,544,506,456,767,463,162,196,525]
latent_ch2_100s_32k = [0.00264021,0.00266504,0.00297933,0.00235575,0.00186678,0.0018307,0.00257865,0.00275751,0.00243414,0.00237035,0.00236514,0.00181787,0.00211641,0.00254605]
noise4x_ch2_100s_32k=[0.0034709,0.00323427,0.0031826,0.00399558,0.00323423,0.00331001,0.00843324,0.00829871,0.00760383,0.0105887,0.00802234,0,0,0]

resetgood_ch2_100s_32k = where(bias_ch2_100s_32k eq 450 and vdduc_ch2_100s_32k eq -3.5 and vgg1_ch2_100s_32k eq -3.7)
dducgood_ch2_100s_32k = where(bias_ch2_100s_32k eq 450 and vgg1_ch2_100s_32k eq -3.7 and vreset_ch2_100s_32k eq -3.5)
vggood_ch2_100s_32k = where(bias_ch2_100s_32k eq 450 and vdduc_ch2_100s_32k eq -3.5 and vreset_ch2_100s_32k eq -3.5)

;------------------------------------------------------------------------------------------------------
;0.4s plots

plot, vreset_ch2_0p4s_30k(resetgood_ch2_0p4s_30k), mean_ch2_0p4s_30k(resetgood_ch2_0p4s_30k), psym = 4, xrange=[-3.7,-3.3], yrange=[1.5,2.0], xtitle ='Vreset', ytitle ='mean(stddev) MJy/sr', title = '0p4s Ch2 b450 g36(30k) g37(32k) d35',charsize=1.5
oplot, vreset_ch2_0p4s_30k(resetgood_ch2_0p4s_30k), mean_ch2_0p4s_30k(resetgood_ch2_0p4s_30k)
oplot,  vreset_ch2_0p4s_32k(resetgood_ch2_0p4s_32k), mean_ch2_0p4s_32k(resetgood_ch2_0p4s_32k)
oplot,  vreset_ch2_0p4s_32k(resetgood_ch2_0p4s_32k), mean_ch2_0p4s_32k(resetgood_ch2_0p4s_32k),psym = 2

plot, vreset_ch2_0p4s_30k(resetgood_ch2_0p4s_30k), noisy_ch2_0p4s_30k(resetgood_ch2_0p4s_30k), psym = 4, xrange=[-3.7,-3.3], yrange=[0,20],xtitle ='Vreset', ytitle ='10sigma noisypixels', title = '0p4s Ch2 b450 g36(30k) g37(32k) d35', charsize=1.5
oplot, vreset_ch2_0p4s_30k(resetgood_ch2_0p4s_30k), noisy_ch2_0p4s_30k(resetgood_ch2_0p4s_30k)
oplot,  vreset_ch2_0p4s_32k(resetgood_ch2_0p4s_32k), noisy_ch2_0p4s_32k(resetgood_ch2_0p4s_32k)
oplot,  vreset_ch2_0p4s_32k(resetgood_ch2_0p4s_32k), noisy_ch2_0p4s_32k(resetgood_ch2_0p4s_32k),psym = 2

plot, vreset_ch2_0p4s_30k(resetgood_ch2_0p4s_30k), latent_ch2_0p4s_30k(resetgood_ch2_0p4s_30k), psym = 4, xrange=[-3.7,-3.3], xtitle ='Vreset', ytitle ='latent strength', title = '0p4s Ch2 b450 g36(30k) g37(32k) d35', yrange=[.001,.004],charsize=1.5
oplot, vreset_ch2_0p4s_30k(resetgood_ch2_0p4s_30k), latent_ch2_0p4s_30k(resetgood_ch2_0p4s_30k)
oplot,  vreset_ch2_0p4s_32k(resetgood_ch2_0p4s_32k), latent_ch2_0p4s_32k(resetgood_ch2_0p4s_32k)
oplot,  vreset_ch2_0p4s_32k(resetgood_ch2_0p4s_32k), latent_ch2_0p4s_32k(resetgood_ch2_0p4s_32k),psym = 2


plot, vdduc_ch2_0p4s_32k(dducgood_ch2_0p4s_32k), mean_ch2_0p4s_32k(dducgood_ch2_0p4s_32k), psym = 2, xrange=[-3.6,-3.0],  yrange=[1.5,2.0], xtitle ='Vdduc', ytitle ='mean(stddev) MJy/sr', title = '0p4s Ch2 b450 g37 r35', charsize=1.5
oplot, vdduc_ch2_0p4s_32k(dducgood_ch2_0p4s_32k), mean_ch2_0p4s_32k(dducgood_ch2_0p4s_32k)
legend, ['30k', '32k'], psym=[4,2],/center,/bottom

plot, vdduc_ch2_0p4s_32k(dducgood_ch2_0p4s_32k), noisy_ch2_0p4s_32k(dducgood_ch2_0p4s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[0,20],xtitle ='Vdduc', ytitle ='10sigma noisypixels', title = '0p4s Ch2 b450 g37 r35', charsize=1.5
oplot, vdduc_ch2_0p4s_32k(dducgood_ch2_0p4s_32k), noisy_ch2_0p4s_32k(dducgood_ch2_0p4s_32k)
legend, ['30k', '32k'], psym=[4,2],/center,/bottom

plot, vdduc_ch2_0p4s_32k(dducgood_ch2_0p4s_32k), latent_ch2_0p4s_32k(dducgood_ch2_0p4s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[.001,.004], xtitle ='Vdduc', ytitle ='latent strength', title = '0p4s Ch2 b450 g37 r35', charsize=1.5
oplot, vdduc_ch2_0p4s_32k(dducgood_ch2_0p4s_32k), latent_ch2_0p4s_32k(dducgood_ch2_0p4s_32k)
legend, ['30k', '32k'], psym=[4,2],/center,/bottom



plot, vgg1_ch2_0p4s_32k(vggood_ch2_0p4s_32k), mean_ch2_0p4s_32k(vggood_ch2_0p4s_32k), psym = 2, xtitle ='Vgg1', ytitle ='mean(stddev) MJy/sr', title = '0p4s Ch2 b450 r35 d35', xrange=[-4.0,-3.4],  charsize=1.5
oplot, vgg1_ch2_0p4s_32k(vggood_ch2_0p4s_32k), mean_ch2_0p4s_32k(vggood_ch2_0p4s_32k)
oplot,  vgg1_ch2_0p4s_30k(vggood_ch2_0p4s_30k), mean_ch2_0p4s_30k(vggood_ch2_0p4s_30k), psym = 4
oplot,  vgg1_ch2_0p4s_30k(vggood_ch2_0p4s_30k), mean_ch2_0p4s_30k(vggood_ch2_0p4s_30k)

plot, vgg1_ch2_0p4s_32k(vggood_ch2_0p4s_32k), noisy_ch2_0p4s_32k(vggood_ch2_0p4s_32k), psym = 2, xtitle ='Vgg1', ytitle ='10sigma noisypixels', title = '0p4s Ch2 b450 r35 d35', xrange=[-4.0,-3.4], yrange=[0,20],charsize=1.5
oplot, vgg1_ch2_0p4s_32k(vggood_ch2_0p4s_32k), noisy_ch2_0p4s_32k(vggood_ch2_0p4s_32k)
oplot,  vgg1_ch2_0p4s_30k(vggood_ch2_0p4s_30k), noisy_ch2_0p4s_30k(vggood_ch2_0p4s_30k), psym = 4
oplot,  vgg1_ch2_0p4s_30k(vggood_ch2_0p4s_30k), noisy_ch2_0p4s_30k(vggood_ch2_0p4s_30k)

plot, vgg1_ch2_0p4s_32k(vggood_ch2_0p4s_32k), latent_ch2_0p4s_32k(vggood_ch2_0p4s_32k), psym = 2, xtitle ='Vgg1', ytitle ='latent strength', title = '0p4s Ch2 b450 r35 d35', xrange=[-4.0,-3.4], yrange=[.001,.004], charsize=1.5
oplot, vgg1_ch2_0p4s_32k(vggood_ch2_0p4s_32k), latent_ch2_0p4s_32k(vggood_ch2_0p4s_32k)
oplot,  vgg1_ch2_0p4s_30k(vggood_ch2_0p4s_30k), latent_ch2_0p4s_30k(vggood_ch2_0p4s_30k), psym = 4
oplot,  vgg1_ch2_0p4s_30k(vggood_ch2_0p4s_30k), latent_ch2_0p4s_30k(vggood_ch2_0p4s_30k)


;------------------------------------------------------------------------------------------------------
;12s plots

plot, vreset_ch2_12s_30k(resetgood_ch2_12s_30k), mean_ch2_12s_30k(resetgood_ch2_12s_30k), psym = 4, xrange=[-3.7,-3.3], xtitle ='Vreset', ytitle ='mean(stddev) MJy/sr', title = '12s Ch2 b450 g36(30k) g37(32k) d35', yrange=[.03,.05], charsize=1.5
oplot, vreset_ch2_12s_30k(resetgood_ch2_12s_30k), mean_ch2_12s_30k(resetgood_ch2_12s_30k)
oplot,  vreset_ch2_12s_32k(resetgood_ch2_12s_32k), mean_ch2_12s_32k(resetgood_ch2_12s_32k)
oplot,  vreset_ch2_12s_32k(resetgood_ch2_12s_32k), mean_ch2_12s_32k(resetgood_ch2_12s_32k),psym = 2

plot, vreset_ch2_12s_30k(resetgood_ch2_12s_30k), noisy_ch2_12s_30k(resetgood_ch2_12s_30k), psym = 4, xrange=[-3.7,-3.3], xtitle ='Vreset', ytitle ='10sigma noisypixels', title = '12s Ch2 b450 g36(30k) g37(32k) d35', charsize=1.5
oplot, vreset_ch2_12s_30k(resetgood_ch2_12s_30k), noisy_ch2_12s_30k(resetgood_ch2_12s_30k)
oplot,  vreset_ch2_12s_32k(resetgood_ch2_12s_32k), noisy_ch2_12s_32k(resetgood_ch2_12s_32k)
oplot,  vreset_ch2_12s_32k(resetgood_ch2_12s_32k), noisy_ch2_12s_32k(resetgood_ch2_12s_32k),psym = 2

plot, vreset_ch2_12s_30k(resetgood_ch2_12s_30k), latent_ch2_12s_30k(resetgood_ch2_12s_30k), psym = 4, xrange=[-3.7,-3.3], xtitle ='Vreset', ytitle ='latent strength', title = '12s Ch2 b450 g36(30k) g37(32k) d35', yrange=[.001,.004],charsize=1.5
oplot, vreset_ch2_12s_30k(resetgood_ch2_12s_30k), latent_ch2_12s_30k(resetgood_ch2_12s_30k)
oplot,  vreset_ch2_12s_32k(resetgood_ch2_12s_32k), latent_ch2_12s_32k(resetgood_ch2_12s_32k)
oplot,  vreset_ch2_12s_32k(resetgood_ch2_12s_32k), latent_ch2_12s_32k(resetgood_ch2_12s_32k),psym = 2


plot, vdduc_ch2_12s_32k(dducgood_ch2_12s_32k), mean_ch2_12s_32k(dducgood_ch2_12s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[.03,.05], xtitle ='Vdduc', ytitle ='mean(stddev) MJy/sr', title = '12s Ch2 b450 g37 r35', charsize=1.5
oplot, vdduc_ch2_12s_32k(dducgood_ch2_12s_32k), mean_ch2_12s_32k(dducgood_ch2_12s_32k)
legend, ['30k', '32k'], psym=[4,2],/center,/bottom

plot, vdduc_ch2_12s_32k(dducgood_ch2_12s_32k), noisy_ch2_12s_32k(dducgood_ch2_12s_32k), psym = 2, xrange=[-3.6,-3.0], xtitle ='Vdduc', ytitle ='10sigma noisypixels', title = '12s Ch2 b450 g37 r35', charsize=1.5
oplot, vdduc_ch2_12s_32k(dducgood_ch2_12s_32k), noisy_ch2_12s_32k(dducgood_ch2_12s_32k)
legend, ['30k', '32k'], psym=[4,2],/center,/bottom

plot, vdduc_ch2_12s_32k(dducgood_ch2_12s_32k), latent_ch2_12s_32k(dducgood_ch2_12s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[.001,.004], xtitle ='Vdduc', ytitle ='latent strength', title = '12s Ch2 b450 g37 r35', charsize=1.5
oplot, vdduc_ch2_12s_32k(dducgood_ch2_12s_32k), latent_ch2_12s_32k(dducgood_ch2_12s_32k)
legend, ['30k', '32k'], psym=[4,2],/center,/bottom



plot, vgg1_ch2_12s_32k(vggood_ch2_12s_32k), mean_ch2_12s_32k(vggood_ch2_12s_32k), psym = 2, xtitle ='Vgg1', ytitle ='mean(stddev) MJy/sr', title = '12s Ch2 b450 r35 d35', xrange=[-4.0,-3.4], yrange=[.03,.05], charsize=1.5
oplot, vgg1_ch2_12s_32k(vggood_ch2_12s_32k), mean_ch2_12s_32k(vggood_ch2_12s_32k)
oplot,  vgg1_ch2_12s_30k(vggood_ch2_12s_30k), mean_ch2_12s_30k(vggood_ch2_12s_30k), psym = 4
oplot,  vgg1_ch2_12s_30k(vggood_ch2_12s_30k), mean_ch2_12s_30k(vggood_ch2_12s_30k)

plot, vgg1_ch2_12s_32k(vggood_ch2_12s_32k), noisy_ch2_12s_32k(vggood_ch2_12s_32k), psym = 2, xtitle ='Vgg1', ytitle ='10sigma noisypixels', title = '12s Ch2 b450 r35 d35', xrange=[-4.0,-3.4], charsize=1.5
oplot, vgg1_ch2_12s_32k(vggood_ch2_12s_32k), noisy_ch2_12s_32k(vggood_ch2_12s_32k)
oplot,  vgg1_ch2_12s_30k(vggood_ch2_12s_30k), noisy_ch2_12s_30k(vggood_ch2_12s_30k), psym = 4
oplot,  vgg1_ch2_12s_30k(vggood_ch2_12s_30k), noisy_ch2_12s_30k(vggood_ch2_12s_30k)

plot, vgg1_ch2_12s_32k(vggood_ch2_12s_32k), latent_ch2_12s_32k(vggood_ch2_12s_32k), psym = 2, xtitle ='Vgg1', ytitle ='latent strength', title = '12s Ch2 b450 r35 d35', xrange=[-4.0,-3.4], yrange=[.001,.004], charsize=1.5
oplot, vgg1_ch2_12s_32k(vggood_ch2_12s_32k), latent_ch2_12s_32k(vggood_ch2_12s_32k)
oplot,  vgg1_ch2_12s_30k(vggood_ch2_12s_30k), latent_ch2_12s_30k(vggood_ch2_12s_30k), psym = 4
oplot,  vgg1_ch2_12s_30k(vggood_ch2_12s_30k), latent_ch2_12s_30k(vggood_ch2_12s_30k)
;------------------------------------------------------------------------------------------------------
;100s plots

plot, vreset_ch2_100s_30k(resetgood_ch2_100s_30k), mean_ch2_100s_30k(resetgood_ch2_100s_30k), psym = 4, xrange=[-3.7,-3.3], yrange=[.008,.012],xtitle ='Vreset', ytitle ='mean(stddev) MJy/sr', title = '100s Ch2 b450 g36(30k) g37(32k) d35',  charsize=1.5
oplot, vreset_ch2_100s_30k(resetgood_ch2_100s_30k), mean_ch2_100s_30k(resetgood_ch2_100s_30k)
oplot,  vreset_ch2_100s_32k(resetgood_ch2_100s_32k), mean_ch2_100s_32k(resetgood_ch2_100s_32k)
oplot,  vreset_ch2_100s_32k(resetgood_ch2_100s_32k), mean_ch2_100s_32k(resetgood_ch2_100s_32k),psym = 2

plot, vreset_ch2_100s_30k(resetgood_ch2_100s_30k), noise4x_ch2_100s_30k(resetgood_ch2_100s_30k), psym = 4, xrange=[-3.7,-3.3], xtitle ='Vreset', ytitle ='fraction noisypixels', title = '100s Ch2 b450 g36(30k) g37(32k) d35', charsize=1.5
;oplot, vreset_ch2_100s_30k(resetgood_ch2_100s_30k), noise4x_ch2_100s_30k(resetgood_ch2_100s_30k)
;oplot,  vreset_ch2_100s_32k(resetgood_ch2_100s_32k), noise4x_ch2_100s_32k(resetgood_ch2_100s_32k)
oplot,  vreset_ch2_100s_32k(resetgood_ch2_100s_32k), noise4x_ch2_100s_32k(resetgood_ch2_100s_32k),psym = 2

plot, vreset_ch2_100s_30k(resetgood_ch2_100s_30k), latent_ch2_100s_30k(resetgood_ch2_100s_30k), psym = 4, xrange=[-3.7,-3.3], xtitle ='Vreset', ytitle ='latent strength', title = '100s Ch2 b450 g36(30k) g37(32k) d35', yrange=[.001,.004],charsize=1.5
oplot, vreset_ch2_100s_30k(resetgood_ch2_100s_30k), latent_ch2_100s_30k(resetgood_ch2_100s_30k)
oplot,  vreset_ch2_100s_32k(resetgood_ch2_100s_32k), latent_ch2_100s_32k(resetgood_ch2_100s_32k)
oplot,  vreset_ch2_100s_32k(resetgood_ch2_100s_32k), latent_ch2_100s_32k(resetgood_ch2_100s_32k),psym = 2


plot, vdduc_ch2_100s_32k(dducgood_ch2_100s_32k), mean_ch2_100s_32k(dducgood_ch2_100s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[.008,.012], xtitle ='Vdduc', ytitle ='mean(stddev) MJy/sr', title = '100s Ch2 b450 g37 r35', charsize=1.5
oplot, vdduc_ch2_100s_32k(dducgood_ch2_100s_32k), mean_ch2_100s_32k(dducgood_ch2_100s_32k)
legend, ['30k', '32k'], psym=[4,2],/center,/bottom

plot, vdduc_ch2_100s_32k(dducgood_ch2_100s_32k), noise4x_ch2_100s_32k(dducgood_ch2_100s_32k), psym = 2, xrange=[-3.6,-3.0], xtitle ='Vdduc', ytitle ='fraction noisypixels', title = '100s Ch2 b450 g37 r35', charsize=1.5
;oplot, vdduc_ch2_100s_32k(dducgood_ch2_100s_32k), noise4x_ch2_100s_32k(dducgood_ch2_100s_32k)
legend, ['30k', '32k'], psym=[4,2],/center,/bottom

plot, vdduc_ch2_100s_32k(dducgood_ch2_100s_32k), latent_ch2_100s_32k(dducgood_ch2_100s_32k), psym = 2, xrange=[-3.6,-3.0], yrange=[.001,.004], xtitle ='Vdduc', ytitle ='latent strength', title = '100s Ch2 b450 g37 r35', charsize=1.5
oplot, vdduc_ch2_100s_32k(dducgood_ch2_100s_32k), latent_ch2_100s_32k(dducgood_ch2_100s_32k)
legend, ['30k', '32k'], psym=[4,2],/center,/bottom



plot, vgg1_ch2_100s_32k(vggood_ch2_100s_32k), mean_ch2_100s_32k(vggood_ch2_100s_32k), psym = 2, xtitle ='Vgg1', ytitle ='mean(stddev) MJy/sr', title = '100s Ch2 b450 r35 d35', xrange=[-4.0,-3.4], yrange=[.008,.012], charsize=1.5
oplot, vgg1_ch2_100s_32k(vggood_ch2_100s_32k), mean_ch2_100s_32k(vggood_ch2_100s_32k)
oplot,  vgg1_ch2_100s_30k(vggood_ch2_100s_30k), mean_ch2_100s_30k(vggood_ch2_100s_30k), psym = 4
oplot,  vgg1_ch2_100s_30k(vggood_ch2_100s_30k), mean_ch2_100s_30k(vggood_ch2_100s_30k)

plot, vgg1_ch2_100s_32k(vggood_ch2_100s_32k), noise4x_ch2_100s_32k(vggood_ch2_100s_32k), psym = 2, xtitle ='Vgg1', ytitle ='fraction noisypixels', title = '100s Ch2 b450 r35 d35', xrange=[-4.0,-3.4], charsize=1.5
;oplot, vgg1_ch2_100s_32k(vggood_ch2_100s_32k), noise4x_ch2_100s_32k(vggood_ch2_100s_32k)
oplot,  vgg1_ch2_100s_30k(vggood_ch2_100s_30k), noise4x_ch2_100s_30k(vggood_ch2_100s_30k), psym = 4
;oplot,  vgg1_ch2_100s_30k(vggood_ch2_100s_30k), noise4x_ch2_100s_30k(vggood_ch2_100s_30k)

plot, vgg1_ch2_100s_32k(vggood_ch2_100s_32k), latent_ch2_100s_32k(vggood_ch2_100s_32k), psym = 2, xtitle ='Vgg1', ytitle ='latent strength', title = '100s Ch2 b450 r35 d35', xrange=[-4.0,-3.4], yrange=[.001,.004], charsize=1.5
oplot, vgg1_ch2_100s_32k(vggood_ch2_100s_32k), latent_ch2_100s_32k(vggood_ch2_100s_32k)
oplot,  vgg1_ch2_100s_30k(vggood_ch2_100s_30k), latent_ch2_100s_30k(vggood_ch2_100s_30k), psym = 4
oplot,  vgg1_ch2_100s_30k(vggood_ch2_100s_30k), latent_ch2_100s_30k(vggood_ch2_100s_30k)
ps_end

end
