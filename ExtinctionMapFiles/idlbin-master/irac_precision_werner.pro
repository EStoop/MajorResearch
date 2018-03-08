pro irac_precision_werner
;make a plot of the noise levels in an exoplanet staring mode observation
;some things are hard coded in so be careful to expand this to other observations

;maybe they should all be bands to show a range of eg. eclipse depths.
  xmax = 1E3
  bin_scale = findgen(xmax) + 1
  root_n = sqrt(bin_scale)
 
;XXXhow many electrons are in a source with magnitude X?
 ;start with mJy
  ch = '2'
  if ch eq '2' then begin
     gain = 3.71
     pixel_scale = 1.22
     flux_conv = .1469
  endif

  ;try HD189733  ; K = 5.5, K0 star
  source_mjy = 1094. ;mJy  from Star-pet
  exptime = 0.1
     

                                
  starname = ['K=11', 'K=12', 'K=8.6']  ; assuming M5 stars
  source_mjy =[ 7.67, 3.05, 68]  ; ch2 only
  exptime = 30.
  colorarr = ['red', 'blue', 'cyan']
  for s = 0, n_elements(source_mjy) - 1 do begin
     
     source_electrons = mjy_to_electron( source_mjy(s), pixel_scale, gain, exptime, flux_conv)
     
     sigma_poisson = sqrt(source_electrons)
     print, 'poisson noise', sigma_poisson
     

;ok, but want to make the y-axis in percentage of the source flux
     y =  (sigma_poisson / root_n) / source_electrons
     p = plot( bin_scale, y, thick = 3, xtitle = 'Binning Scale (N frames)', ytitle = 'Noise', axis_style = 1,/xlog,/ylog, xrange = [1, xmax], margin = 0.2, overplot = p, color = colorarr(s), name = starname(s), yrange = [1E-5, 1E-2] )
;  t = text(3., 200., 'Source Poisson', color = 'black', /current,/data)
     xaxis = axis('x', location = [0,max(p.yrange)], coord_transform = [0, exptime/60.],target = p, textpos = 1, title = 'Binning Scale (Min.)')
     xaxis = axis('y', location = [max(p.xrange),0], target = p, tickdir = 0, textpos = 0, showtext = 0)


;add a turnover at a certain binning scale
;XXX what binning scale?
;XXX what strength turnover?



  endfor
;highlight predicted transit depths
     y = [5000E-6, 5000E-6]
     y2 = [1000E-6, 1000E-6]
     x = [1,  1000]
     ;ph = plot(x, y,  thick = 2, overplot = p)
     ph = plot(x, y2,  thick = 2, overplot = p)
     ;t = text(10, 5050E-6, '0.5% depth',/data) 
     t = text(10, 1000E-6, '1000ppm',/data) 

;make a legend
     t = text(100,.0007, 'K=12', color ='blue', /data)
     t = text(50,.0003, 'K=11', color ='red', /data)
     t = text(20,.0002, 'K=8.6', color ='cyan', /data)

end
