pro irac_precision_kepler78
;make a plot of the noise levels in an exoplanet staring mode observation
;some things are hard coded in so be careful to expand this to other observations

;maybe they should all be bands to show a range of eg. eclipse depths.
xmax = 1E4
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
  exptime = 2.
     

                                ;Kepler-78  ; K = 9.5 Gstar
  starname = ['Kepler-78 ch2', 'Kepler-78 ch1']
  source_mjy =[ 27, 42.6]
  exptime = 2.
  colorarr = ['red', 'blue']
  for s = 0, n_elements(source_mjy) - 1 do begin
     
     source_electrons = mjy_to_electron( source_mjy(s), pixel_scale, gain, exptime, flux_conv)
     
     sigma_poisson = sqrt(source_electrons)
     print, 'poisson noise', sigma_poisson
     

;ok, but want to make the y-axis in percentage of the source flux
     y =  ((sigma_poisson / root_n) / source_electrons) * 1E6
     p = plot( bin_scale, y, thick = 3, xtitle = 'Binning Scale (N frames)', ytitle = 'PPM', axis_style = 1,/xlog,/ylog, xrange = [1, xmax], margin = 0.2, overplot = p, color = colorarr(s), name = starname(s) )
;  t = text(3., 200., 'Source Poisson', color = 'black', /current,/data)
     xaxis = axis('x', location = [0,max(p.yrange)], coord_transform = [0, exptime/60.],target = p, textpos = 1, title = 'Binning Scale (Min.)')
     xaxis = axis('y', location = [max(p.xrange),0], target = p, tickdir = 0, textpos = 0, showtext = 0)


;add a turnover at a certain binning scale
;XXX what binning scale?
;XXX what strength turnover?


;make a legend
     t = text(10, (2E-3 + s*(1E-3)) * 1E6, starname(s), color = colorarr(s), /data)

;highlight eclipse duration
     y = [1E-5, 1E-4, 1E-3, 1E-2] * 1E6
     x = [1200, 1200, 1200, 1200]
     ph = plot(x, y, linestyle = 2, thick = 2, overplot = p)

;highlight predicted transit depth
     y = [220E-6, 220E-6, 220E-6, 220E-6] * 1E6
     x = [1, 10, 100, 10000]
     ph = plot(x, y,  thick = 2, overplot = p)
     t = text(1.1, 250, 'predicted transit depth',/data) 

     t = text(1700, 250, 'transit duration',/data, orientation = 90) 


  endfor


end
