pro ch2_pixelslew
;ps_open, filename='/Users/jkrick/nutella/spitzer/planets_proposal/centering.eps',/color,/encapsulated
ps_start, filename='/Users/jkrick/nutella/spitzer/planets_proposal/centering.eps'

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)

pixel_dist = [-1.0,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0.0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
xhist = [3.0,5.0,3.0,10.0,9.0,15.0,28.0,42.0,55.0,85.0,136.0,144.0,143.0,135.0,108.0,68.0,65.0,42.0,28.0,16.0,14.0]
yhist = [12.0,17.0,44.0,40.0,53.0,64.0,85.0,81.0,129.0,136.0,136.0,121.0,78.0,64.0,37.0,18.0,16.0,11.0,6.0,4.0,2.0]

plot, pixel_dist, yhist, xtitle = 'Distance from Desired Position in Pixels', ytitle = 'Number', /nodata, thick = 3, charthick = 3, xthick = 3, ythick = 3
oplot, pixel_dist+.108, yhist, color = redcolor, thick = 3
oplot,  pixel_dist - .18, xhist, color = bluecolor, thick = 3

;y
start = [0.1,0.2, 5000.]
noise = fltarr(n_elements(yhist))
noise[*] = 1                                                    ;equally weight the values
resulty= MPFITFUN('mygauss',pixel_dist,yhist, noise, start)    ;./quiet    

oplot, pixel_dist +.108, (resulty(2))/sqrt(2.*!Pi) * exp(-0.5*((pixel_dist - (resulty(0)))/(resulty(1)))^2.), $
          linestyle =  2 , color = redcolor, thick  = 3

;x
resultx= MPFITFUN('mygauss',pixel_dist,xhist, noise, start)    ;./quiet    


oplot, pixel_dist - 0.18, (resultx(2))/sqrt(2.*!Pi) * exp(-0.5*((pixel_dist - (resultx(0)))/(resultx(1)))^2.), $
          linestyle =  2 , color = bluecolor, thick  = 3
ps_end
;ps_close, /noprint,/noid

end
