pro plot_barrel_shift

;make final plots
frametime = ['s0.02','s0.1','f0.4','s0.4','h0.6','h0.6','h1.2','f2','h6','f12','h12','f30','h30']
ch1_pipeline_shift1 = [0.7,1.1,0.1,1.3,0.2,-0.01,0.06,1.5,1.3,0.3,0.04,0.8,-0.2]
ch1_pipeline_shift2 = [3.2,5.3,2.8,5.6,2.5,2.8,1.1,2,2.5,1.7,2.3,2.1,1.7]
ch2_pipeline_shift1 = [1.4,1.9,1.3,2.5,-0.5,-0.7,2,2.4,2.1,-1.8,-1.8,1.1,1.1]
ch2_pipeline_shift2 = [5.4,8.7,5.6,11.6,3.7,3.5,4.7,9.7,11.6,5,4.8,3.2,2.9]

frametime = ['s0.02','s0.1','f0.4','s0.4','h0.6','h0.6','h1.2','f2','h6','f12','h12','f30','h30']
ch2_predict_shift1 = [0.3,0.6,0.6,1.1,0.5,0.5,0.9,2.1,1.7,1.0,1.0,1.3,1.3]
ch2_predict_shift2 = [1.5, 3, 2.9, 5.4, 2.6, 2.6, 4.2, 9.9, 8.3, 4.9, 4.9, 0, 0]
ch1_predict_shift1 =[0.3,0.7,0.3,1,0.3,0.3,0.5,1,1.1,0.7,0.7,0.5,0.5]
ch1_predict_shift2 = [1.7,3.2,1.5,5.2,1.5,1.5,2.5,4.7,5.6,3.6,3.6,2.6,2.6]
;ch1_predict_shift1 = [0.3,0.6,0.3,1.1,0.3,0.3,0.6,1,1.5,1.2,1.2,1.3,1.3]
;ch1_predict_shift2 = [1.5, 3, 1.5, 5.4, 1.5, 1.5, 3, 5, 0, 6, 6, 0, 0]
;ch2_predict_shift1 =[0.3,0.7,0.3,1,0.3,0.3,0.5,1,1.1,0.7,0.7,0.5,0.5]
;ch2_predict_shift2 = [1.7,3.2,1.5,5.2,1.5,1.5,2.5,4.7,5.6,3.6,3.6,2.6,2.6]


x = findgen(n_elements(frametime))

;raw
frametime =                      ['s0.02','s0.1','f0.4','s0.4','h0.6', 'h0.6',  'h1.2','f2','h6','f12','h12','f30','h30']
ch1_pipeline_shift1_raw = [0.5,     0.5,  -0.1,    0.6,-0.1,  -0.007,0.02,1.5,1.0,  0.3, -0.3, 0.1, -0.3]
ch1_pipeline_shift2_raw = [2.2,     2.8,    1.8,    3.2,   2.6 ,  1.58,  -.39,  2.1,1.6,  1.6, 1.9 ,1.6, 1.4]
ch2_pipeline_shift1_raw = [0.7,     1.0,    0.5,    1.5,   0.6 ,  0.66,    .77,  1.6,1.5,  0.9, 0.7 ,0.4, 0.3]
ch2_pipeline_shift2_raw = [3.6,     4.9,    3.0,    7.0,   3.5 ,  2.59,   4.45, 7.2,7.0,  4.0, 4.3 ,2.6, 3.1]
;ch1_pipeline_shift1_raw = [0.5,     0.5,  -0.1,    0.6,-0.1,  -0.007,0.02,1.5,1.0,  0.3, -0.3, 0.1, -0.3]
;ch1_pipeline_shift2_raw = [2.2,     2.8,    1.8,    3.2,   2.6 ,  1.58,  -.39,  2.1,1.6,  1.6, 1.9 ,1.6, 1.4]
;ch2_pipeline_shift1_raw = [0.7,     1.0,    0.5,    1.5,   0.6 ,  0.66,    .77,  1.6,1.5,  0.9, 0.7 ,0.4, 0.3]
;ch2_pipeline_shift2_raw = [3.6,     4.9,    3.0,    7.0,   3.5 ,  2.59,   4.45, 7.2,7.0,  4.0, 4.3 ,2.6, 3.1]


ch1_shift1_obs_y = [0.9, 0.3,0.6]
ch1_shift1_obs_x = [2,8,9]
ch1_shift2_obs_y = [1.6]
ch1_shift2_obs_x = [2]
ch2_shift1_obs_y = [0.5, 6.3,2.3]
ch2_shift1_obs_x = [2,8,9]
ch2_shift2_obs_y = [2.7]
ch2_shift2_obs_x = [2]

frametime_ch1 = ['s0.02','s0.1','f0.4','s0.4','h0.6','h0.6','h1.2','f2','h6','f12','h12']
x_ch1 = findgen(n_elements(frametime_ch1))
y_ch1 = [ch1_pipeline_shift2(0), ch1_pipeline_shift1(1), ch1_pipeline_shift2(2), ch1_pipeline_shift1(3), ch1_pipeline_shift2(4), ch1_pipeline_shift2(5),ch1_pipeline_shift1(6), ch1_pipeline_shift1(7) , ch1_pipeline_shift1(8)  , ch1_pipeline_shift1(9)  , ch1_pipeline_shift1(10) ]
frametime_ch2 = ['s0.02','s0.1','f0.4','s0.4','h0.6','h0.6','h1.2','f2','h6','f12','h12','f30','h30']
x_ch2 = findgen(n_elements(frametime_ch2))
y_ch2 = [ch2_pipeline_shift2(0), ch2_pipeline_shift1(1), ch2_pipeline_shift2(2), ch2_pipeline_shift1(3), ch2_pipeline_shift2(4), ch2_pipeline_shift2(5),ch2_pipeline_shift1(6), ch2_pipeline_shift1(7) , ch2_pipeline_shift1(8)  , ch2_pipeline_shift1(9)  , ch2_pipeline_shift1(10) , ch2_pipeline_shift1(11) , ch2_pipeline_shift1(12) ]

print, 'y', y_ch1
print, 'y', y_ch2

e1 = plot(x, ch1_pipeline_shift1, '2s', color = 'orange', yrange = [-5, 15], sym_filled = 1, ytitle = 'Percentage change from nominal barrel shift', name = 'shift 1 bcd', title = 'ch1', xtickname = frametime, xminor = 0)
e2 = plot(x, ch1_pipeline_shift2, '2s', color = 'red', sym_filled = 1,/overplot, name = 'shift 2 bcd')
e3 = plot(x, ch1_predict_shift1, '2', color = 'orange', /overplot, name = 'shift 1 predict')
e4 = plot(x, ch1_predict_shift2, '2', color = 'red', /overplot, name = 'shift 2 predict')
t1 = plot(ch1_shift1_obs_x,ch1_shift1_obs_y, '2 ', symbol = 'star', name= 'shift 1 observed', color = 'orange',/overplot, sym_filled = 1)
t2 = plot(ch1_shift2_obs_x,ch1_shift2_obs_y, '2 ', symbol = 'star', name= 'shift 1 observed', color = 'red',/overplot, sym_filled = 1)
e7 = plot(x, ch1_pipeline_shift1_raw, '2s', color = 'grey', /overplot, sym_filled = 1, name = 'shift 1 raw')
e8 = plot(x, ch1_pipeline_shift2_raw, '2s', color = 'black', /overplot, sym_filled = 1, name = 'shift 2 raw')
l = legend(target = [e3,e4,t1,t2], position = [8,14],/data)
l2 = legend(target = [e1,e2,e7, e8], position = [4, 14], /data)
t = plot(x_ch1, y_ch1, '2o', color = 'grey', sym_size=3, /overplot)



e5 = plot(x, ch2_pipeline_shift1, '2D', color = 'orange', sym_filled = 1,name = 'shift 1 bcd', yrange = [-5,15],ytitle = 'Percentage change from nominal barrel shift',title = 'ch2' , xtickname = frametime, xminor = 0)
e6 = plot(x, ch2_pipeline_shift2, '2D', color = 'red', sym_filled = 1,/overplot, name = 'shift 2 bcd')
e7 = plot(x, ch2_predict_shift1, '2', color = 'orange', /overplot, name = 'shift 1 predict')
e8 = plot(x, ch2_predict_shift2, '2', color = 'red', /overplot, name = 'shift 2 predict')
t3 = plot(ch2_shift1_obs_x,ch2_shift1_obs_y, '2 ', symbol = 'star', name= 'shift 2 observed', color = 'orange',/overplot, sym_filled = 1)
t4 = plot(ch2_shift2_obs_x,ch2_shift2_obs_y, '2 ', symbol = 'star', name= 'shift 2 observed', color = 'red',/overplot, sym_filled = 1)
;add the raw data points
e17 = plot(x, ch2_pipeline_shift1_raw, '2D', color = 'grey', /overplot, sym_filled = 1, name = 'shift 1 raw')
e18 = plot(x, ch2_pipeline_shift2_raw, '2D', color = 'black', /overplot, sym_filled = 1, name = 'shift 2 raw')

l = legend(target = [e7,e8,t3,t4], position = [8,14],/data)
l = legend(target = [e5,e6,e17,e18], position = [3.5,14],/data)

;now track the recommended usage:

t = plot(x_ch2, y_ch2, '2o', color = 'grey', sym_size=3, /overplot)


end
