pro eclipse_test, planetname
;restore,'/Users/jkrick/irac_warm/hd209458/hd209458_plot_Ch4.sav'
;restore, '/Users/jkrick/irac_warm/pcrs_planets/hat22/hat22_plot_ch2.sav'
restore, '/Users/jkrick/external/irac_warm/XO3/XO3_plot_ch2.sav'

x = (bin_timearr - bin_timearr(0))/60./60.
y = bin_corrfluxp / median(bin_corrfluxp)            
yerr = bin_corrfluxerrp 

t1 = 3.; hrs start of ingress
dt = .5 ; hrs duration of ingress (or egress)
;t2 = t1+dt ; hrs end of ingress
;t4 = 3.36 + t1 ; hrs, end of egress
t3 =5.- dt ; hrs start of egress

  pa0 = [median(y), 0.001,t1,t3,dt] ;f0, d, t1, t3, dt

  afargs = {FLUX:y, X:x, ERR:yerr}
  pa = mpfit('eclipse_trapezoid', pa0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg);, parinfo = parinfo) ;,/quiet)
  print, 'status', status
  print, 'errmsg', errmsg
  achi = achi / adof
  print, 'reduced chi squared', achi
 
  model_fit =eclipse_trapezoid(pa, flux = y, x =x, err=yerr)
 ; print, 'model)fit', model_fit[0:100]

;normalize more exactly

  

 ; a = plot(x, y, '1.', color = 'black', xtitle = 'Time (hrs)', ytitle = 'Flux (mJy)', title = 'Ch2')
  a = errorplot(x, y/ pa(0),yerr, '1s', sym_size = 0.3, sym_filled = 1,$
                   color = 'black', xtitle = 'Time', ytitle = 'Normalized Flux', title = planetname, $
                   name = 'flux', yrange =[0.997, 1.003], xrange = [0,7], axis_style = 1,  $
                   xstyle = 1, errorbar_capsize = 0.025)

  flat = where(x lt pa(2), flatcount)
  x1 = x[flat]
  y1 = fltarr(flatcount) + 1.0;  pa(0)
  print, 'x,y', x1[15], y1[15]
  f1 = plot(x1,   y1  , '2b',/overplot)

  ingress = where(x ge pa(2) and x lt (pa(2) + pa(4)))
;  x2 = x[ingress]
;  y2 = y[ingress]- pa[1] * (1.0 - ((pa[2] + pa[4]) - x[ingress]) / pa[4])
;  y2 = (-pa[1]/pa[4])*x2 + 1.0
  x2 = fltarr(2)
  x2(0) = pa(2)
  x2(1) = pa(2) + pa(4)
  y2 = fltarr(2)
  y2(0) = 1.0
  y2(1) = 1.0 - pa[1]

  print, 'x2,y2', x2, y2
  f2 = plot(x2, y2, '2b',/overplot)
 
; In eclipse
   eclipse = where(x ge (pa[2] + pa[4]) and x lt pa[3], count)
   y3 = fltarr(count) + 1.0 - pa[1]; pa[0] - pa[1]
   x3 = x[eclipse]
   f3 = plot(x3, y3, '2b',/overplot)

; Leaving eclipse
   egress = where(x ge pa[3] and x lt (pa[3] + pa[4]), count)
;   y4 = y[egress] + pa[1] * ((x[egress] - pa[3]) / pa[4] - 1.0)
;   x4 = x[egress]
   x4 = fltarr(2)
   x4(0) = pa(3)
   x4(1) = pa(3) + pa(4)
   y4 = fltarr(2)
   y4(0) = 1.0 - pa[1]
   y4(1) = 1.0 
   
   f4 = plot(x4, y4, '2b',/overplot)

   flat = where(x ge (pa[3] + pa[4]), flatcount)
   x5 = x[flat]
   y5 = fltarr(flatcount) + 1.0; pa(0)
   print, 'x,y', x1[15], y1[15]
   f5 = plot(x5,   y5  , '2b',/overplot)
   
end


function eclipse_trapezoid, p, FLUX=flux, X=x, ERR=err
; p is parameters, f0, d, t1, t3, dt
; t is time
; Start of ingress is t1, end of ingress is t2, start of egress is t3, end of egress is t4
; assume ingress and egress have same durations, dt 
   n = n_elements(x)
   model = fltarr(n) + p[0]

   t2 = p[2] + p[4]
   t4 = p[3] + p[4]

; Out of eclipse
    ptr = where(x le p[2], count)
    if (count gt 0) then model[ptr] = p[0]

    ptr = where(x gt t4, count)
    if (count gt 0) then model[ptr] = p[0]

; Beginning eclipse
   ptr = where(x gt p[2] and x le t2, count)
   if (count gt 0) then $
       model[ptr] = model[ptr] - p[1] * (1.0 - (t2 - x[ptr]) / p[4])

; In eclipse
   ptr = where(x gt t2 and x le p[3], count)
   if (count gt 0) then model[ptr] = model[ptr] - p[1]

; Leaving eclipse
   ptr = where(x gt p[3] and x le t4, count)
   if (count gt 0) then $
       model[ptr] = model[ptr] + p[1] * ((x[ptr] - p[3]) / p[4] - 1.0)

 model = (flux - model) / err

return, model
end
