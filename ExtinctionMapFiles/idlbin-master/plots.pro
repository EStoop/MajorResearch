pro plots_now

ps_start, filename='/Users/jkrick/iwic/noisypixvsbias_temp.ps'
!P.multi=[0,2,1]
;ch1
bias30k=[750,450,750,400,550,450,400,450,450,450,450]
noise4x30k=[0.118043,0.00276443,0.051382,0.00132706,0.00782853,0.00333561,0.0014235,0.00265115,0.00272027,0.00321454,0.00298852]


bias32k=[450,750,750,450,750,450,450,750,450,450,450]
noise4x32k=[0.0050177,0.202771,0.253307,0.00491885,0.138904,0.00730718,0.00717723,0.258647,0.0054736,0.00603554,0.00578842]


plot, bias30k, noise4x30k, psym = 4, yrange=[.001,1], /ylog, xrange=[350,800], xtitle ='Bias', ytitle ='fraction of pixels greater than 4x zodiacal light', title = 'Ch1 100s'
oplot, bias32k,noise4x32k, psym = 2
legend, ['30k', '32k'], psym=[4,2],/top,/left



;ch2
bias30k=[450,450,450,450,450,450,500,550,550,600,650,650]
noise4x30k=[0.00159075,0.00166968,0.00156472,0.00365443,0.00171726,0.00224172,0,0.00372114,0.0015977,0.00568311,0.00729459,0.00697771]

bias32k=[450,450,450,450,450,450,600,600,600,600,600]
noise4x32k = [0.0034709,0.00323427,0.0031826,0.00399558,0.00323423,0.00331001,0.00843324,0.00829871,0.00760383,0.0105887,0.00802234]

plot, bias30k, noise4x30k, psym = 4, yrange=[.001,1], /ylog, xrange=[350,800], xtitle ='Bias', ytitle ='fraction of pixels greater than 4x zodiacal light', title = 'Ch2 100s'
oplot, bias32k,noise4x32k, psym = 2
legend, ['30k', '32k'], psym=[4,2],/top,/left

ps_end, /png

end
