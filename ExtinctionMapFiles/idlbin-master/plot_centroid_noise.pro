pro plot_centroid_noise

;measure real sigma clouds

;ch2sub
ch2sub = [0.049,0.019,0.031,0.029,0.025,0.034,0.040,0.026,0.025,0.031,0.026,0.028,0.022,0.030,0.025,0.085,0.031,0.039,0.019,0.049,0.055,0.028,0.063,0.035,0.032,0.077,0.023,0.030, 0.129  , 0.021  , 0.021  , 0.029  , 0.214  , 0.082  , 0.060  , 0.031  , 0.120  , 0.025  ]

;ch2full
ch2full = [0.054 ,0.052 ,0.057 ,0.052 ,0.057 ,0.060 ,0.047 ,0.055 , 0.071, 0.060,0.058 ,0.061 ,0.053 ,0.062 ,0.054 ,0.056 ,0.054 , 0.054, 0.057, 0.062]

;ch1sub
ch1sub = [0.098 ,0.088 ,0.094 ,0.090 ,0.094 ,0.095 ,0.092 ,0.092 ,0.102 ,0.096 ]


 electrons = [56900., 64066., 212558.]
 b = [0.092, 0.06, 0.03]
b = [mean(ch1sub), mean(ch2full), mean(ch2sub)] 
c = plot(electrons, b, 'r4D-', xtitle = 'Number of Source Electrons', ytitle = 'Cloud Sigma', name = 'data', /ylog,yrange = [0.01, 0.1])
 c1 = plot(electrons, 0.5E4*(1/electrons),'b4o-', /overplot, name = 'prediction')
 
l = legend(target = [c, c1], position = [0.25 ,0.35], /normal)


end

