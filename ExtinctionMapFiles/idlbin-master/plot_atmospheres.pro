pro plot_atmospheres

  deltaz = [-1.48, -0.49, -1.97, 0.19, 1.25, 1.65, 1.37, -0.56]
  deltazerr = [0.71, 0.36, 1.32, 0.93, 0.77, 1.47, 0.79, 0.50]
  name = ['WASP-17b', 'HD209458b', 'WASP-19b', 'HAT-P-1b', 'WASP-31b', 'WASP-12b', 'HAT-P-12b', 'HD189733b']
  color_name = ['pink', 'plum', 'pink', 'plum', 'purple', 'purple', 'purple', 'plum']
  H2O = [94, 32, 105, 68, 31, 38, 17, 53.6]
  H2Oerr = [29, 5, 20, 19, 12, 34, 23, 9.6]

  hazy_x = [ -2.6516, -2.4295, -2.1935, -1.6798, -1.3731, -0.9427, -0.3596, 0.4873, 0.1402, 0.9871, 1.7923, 2.7502, 3.4861]
  hazy_y = [99.9475,  92.5984,  86.0187,  71.6880,  64.3389,  57.3573,  49.2733,  38.9846,  43.0266,  34.9426,  29.4308,  22.8166,  17.6722]

  h_x = [-2.4295,-1.3731,0.4873, 3.4861]
  h_y = [92.5984,64.3389,38.9846, 17.67]
  h_name = ['x10', 'x100', 'x1000']
  h_u =[94,66,41, 19]
  h_l =[92,62,36, 14]

  hazy_y = hazy_y(sort(hazy_x))
  hazy_x = hazy_x(sort(hazy_x))

  solar_x = [-2.6642, -2.7475, -2.5809, -1.9574, -1.3466]
  solar_y = [99.9490,	87.8560,	68.7483,	45.9662,	23.5515]
  solar_u = [105,	93,	73,	48,	23.5515]
  solar_l = [92,	83,	65,	44.,	23.5515]
  solar_name = ['x0.1', 'x0.01','x0.001']
  
  cloudy_x = [-2.6516, -1.3882, -0.5413, -0.1804, 1.0657e-4]
  cloudy_y = [	100,	56,	25,	8,	0]
  cloudy_u = [	101,	58,	28,	12, 5];cloudy_y + 3.
  cloudy_l = [	98.9475,	54.2205,	21.6194,	3.7165,	-5.3675]
  cloudy_t = ['x1', 'x10', 'x100']

  ;;plot data points
  p1 = errorplot(deltaz, H2O, deltazerr, H2Oerr, '1o', sym_filled = 1, xrange = [-3.5, 3.5], $
                 yrange = [0, 140], xtitle = 'Difference between HST & Spitzer Planetary Radii', $
                 ytitle = '$1.4\mum H_{2}O$ amplitude (%)', xstyle = 1)
  p1.xrange = [-3.5, 3.5]

  
  ;;plot model lines
  p1 = plot(cloudy_x, cloudy_y ,thick = 2,  color = 'gray', overplot = p1)
  ;;p1 = plot(cloudy_x[0:3], cloudy_y[0:3], '1tu', overplot = p1, thick = 2,  sym_filled = 1, color = 'gray')
  g1 = polygon([cloudy_x, reverse(cloudy_x)], [cloudy_u, reverse(cloudy_l)], /data, /fill_background, fill_color = 'gray', overplot = p1)
  a1 = arrow([-0.2, 0], [10,0], color = 'gray', head_size = 2.5, /data, overplot = p1)
  t1 = text(-1.7, 11, 'Increasing', /data, color = 'gray', overplot = p1, font_style = 'bold')
  t1 = text(-1.7, 6, 'Clouds', /data, color = 'gray', overplot = p1, font_style = 'bold')
  
  p1 = plot(hazy_x, hazy_y, color = 'dark orchid', overplot = p1, thick = 2)
;;  p1 = plot(h_x, h_y, '1tu', sym_filled =1, color = 'dark orchid', overplot = p1)
;;  t1 = text(h_x, h_y, h_name, color = 'dark orchid',/data, overplot = p1)
  g1 = polygon([h_x, reverse(h_x)], [h_u, reverse(h_l)], /data, /fill_background, fill_color = 'dark orchid', overplot = p1)
  a1 = arrow([3.0, 3.6], [21, 16], color = 'dark orchid', head_size = 2.5, /data, overplot = p1)
  t1= text(2.4, 10, 'Increasing', color = 'dark orchid', /data, font_style = 'bold', overplot = p1)
  t1= text(2.4, 05, 'Haze', color = 'dark orchid', /data, font_style = 'bold', overplot = p1)
  
  p1 = plot(solar_x, solar_y, color = 'red', overplot = p1, thick = 2)
;;  p1 = plot(solar_x[0:3], solar_y[0:3],'1tu', sym_filled = 1, color = 'red', overplot = p1, thick = 2)
;;  t1 = text(solar_x[1:3], solar_y[1:3], solar_name, color = 'red',/data, overplot = p1)
  g1 = polygon([solar_x, reverse(solar_x)], [solar_u, reverse(solar_l)], /data, /fill_background, fill_color = 'red', overplot = p1)
  a1 = arrow([-1.4, -1.3], [25, 20], color = 'red', head_size = 2.0, /data, overplot = p1)
  t1= text(-3.1, 31, 'Decreasing',color = 'red', /data, font_style = 'bold',overplot = p1)
  t1= text(-3.1, 26, 'Solar abund.',color = 'red', /data, font_style = 'bold', overplot = p1)

;;  t1 = text(-3.3, 120, 'Clear', color = 'black', /data, font_style = 'bold', overplot = p1)

  ;;re-plot data points
  
  p1 = errorplot(deltaz, H2O, deltazerr, H2Oerr, '2o',  VERT_COLORS =[[250,159,181], [221,52,151], [250,159,181], [221,52,151], [122,1,119],[122,1,119],[122,1,119], [221,52,151]], sym_filled = 1,overplot = p1)
  t1 = text(deltaz+0.1, H2O+1., name, /data,  font_style = 'bold',overplot = p1)

;  pink = [250,159,181]
;  plum = [221,52,151]
;  purp = [122,1,119]
  
end
