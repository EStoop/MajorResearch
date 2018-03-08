pro hd158460_stare_fit
  
  restore,  '/Users/jkrick/irac_warm/snapshots/HD158460_stare_v2.sav'
    
for var = 0.0006, 0.003, 0.0002 do begin
   print, 'working on ', var

  for a = 0, 1 do begin
     print, 'working on ', aorname(a)
     if a eq 0 then begin
        xarr = xarr_0
        yarr = yarr_0
        fluxarr = fluxarr_0
        timearr = timearr_0
        fluxerrarr = fluxerrarr_0
     endif
     
     if a eq 1 then begin
        xarr = xarr_1
        yarr = yarr_1
        fluxarr = fluxarr_1
        timearr = timearr_1
        fluxerrarr = fluxerrarr_1
     endif
     
;get rid of outliers in the pointing
     if a eq 1 then good = where(xarr gt 14.94 and xarr lt 15.12 and yarr gt 14.9 and yarr lt 15.2 )
     if a eq 0 then good = where(xarr gt 14.84 and xarr lt 15.15 and yarr gt 14.55 and yarr lt 14.9 )
     xarr = xarr[good]
     yarr = yarr[good]
     timearr = timearr[good]
     fluxerrarr = fluxerrarr[good]
     fluxarr = fluxarr[good]
     
;get rid of first 30 min of observations
     if a eq 1 then begin
        good = where(timearr gt 0.5)
        xarr = xarr[good]
        yarr = yarr[good]
        timearr = timearr[good]
        fluxerrarr = fluxerrarr[good]
        fluxarr = fluxarr[good]
     endif
     
     x = timearr / max(timearr)  ; fix to 1 phase
     y = fluxarr / fluxarr[1]    ;normalize
     yerr = fluxerrarr / fluxarr[1] ;normalize
     
;the sin curve
     amplitude =  var
     period =  1/(2*!Pi)        ; 1 phase
     phase = 0.
     y_sin = amplitude*sin(x/period + phase ) 
     noise = randomn(seed, n_elements(x), /Double)
     yerr_sin= (noise- mean(noise))/stddev(noise)*0.0002 + 0.0002
     
     y_sum = y + y_sin
     yerr_sum = sqrt( yerr^2 + yerr_sin^2)
     
;the plots!
 ;    p1 = errorplot(x, y, yerr, '+2', yrange = [0.96,1.02], title = aorname(a), xtitle = 'Phase', ytitle = 'Normalized Flux ') ;, errorbar_capsize = .001)
 ;    p3 = errorplot(x, y_sum,yerr_sum, '+b2',/overplot)
     

;test input values
     ;print, 'mean y', mean(y)
     ;print, 'mean yerr', mean(yerr)
     
;can I recover the sin curve? Maybe?
;now need to make a hybrid function
; Initial guess	
     pa0 = [median(y_sum), 1.,1.,1.,1.,1.,amplitude,period,phase]
     amx = median(xarr)         ;xa[agptr])
     amy = median(yarr)         ; ya[agptr])	
     adx = xarr - amx           ; xa[agptr] - amx
     ady = yarr - amy           ; ya[agptr] - amy
     
     parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(pa0))
                                ;limit the range of the sin cuve phase
     parinfo[8].fixed = 1 
     parinfo[7].fixed = 1          
     parinfo[6].limited[0] = 1
     parinfo[6].limits[0] = 0.0

     
     afargs = {FLUX:y_sum, DX:adx, DY:ady, T:x, ERR:yerr_sum}
     pa = mpfit('fpa1_xfunc3', pa0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg, parinfo = parinfo,/quiet)
     print, 'status', status
     print, errmsg
     achi = achi / adof
     print, 'reduced chi squared', achi
     print, 'input amplitude, period, phase', amplitude, period*2.*!pi, phase/(2.*!pi)
     print, 'fitted amplitude', pa[6], spa[6]
     print, 'fitted period (hours)', pa[7]*2.*!pi, spa[7]*2.*!pi
     print, 'fitted phase [0,1]', pa[8]/(2.*!pi), spa[8]/(2.*!pi)
     
     plus = amplitude + 2*spa[6]
     minus = amplitude - 2*spa[6]
     if pa[6] lt plus and pa[6] gt minus then metric = 'within 2sigma' else metric = 'not 2sigma'
     nsigma = (abs(amplitude - pa[6])) / spa[6]
     
     ;print, 'number of sigma away from input   ', nsigma
     print, 'percent from input ',100.* (abs( amplitude - pa[6]) / amplitude)
     print, '------------------------------------------------'
     ;t = text(0.2, 0.965, 'chi sq' + string(achi), /data, /current)
     ;t = text(0.2, 0.962, 'ampl.' + string(pa[6]) + string(spa[6]), /data, /current)
     
     
  ;look at correlations amongs parameters
  ;PCOR = COV * 0
  ;n = n_elements(spa)
  ;FOR i = 0, n-1 DO FOR j = 0, n-1 DO PCOR[i,j] = COV[i,j]/sqrt(COV[i,i]*COV[j,j])
  ;print, 'pcor', pcor
  
;     f1 = plot(x,   pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + $
 ;                            pa[5] * ady * ady ) + pa[6]*sin(x/pa[7] + pa[8]) , '6r1.',/overplot)
     
  ;plot residuals
     sub = y_sum - (pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + $
                             pa[5] * ady * ady ) + pa[6]*sin(x/pa[7] + pa[8]))
;     f2 = plot(x, sub+0.975, '6g1.',/overplot)
     
     
;fit residuals to know what the errors look like.
;     plothist, sub, xhist, yhist,bin=0.0010, /noprint,/noplot
;     r0 = barplot(xhist, yhist, fill_color = 'green', xrange = [-0.04, 0.04], xtitle = 'Residuals', ytitle = 'Number',title = aorname(a))
     
;     noise = fltarr(n_elements(yhist)) + 1.0 ;  
;     by = where(yhist lt 100.)
;     noise(by) = 100.
;     ;print, 'noise', noise
;     start = [0.0001, 0.01, 10000.]
;     result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) 
;     r1 = plot( xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.),  '3',/overplot)
     
;     t_sig = text(0.02, 2000, string(result(1)),'3',/data)

  endfor  ; for each aor
endfor ; for each variable

end

function fpa1_xfunc2, p, FLUX=flux, DX=dx, DY=dy, T=t, ERR=err
; p is parameters, f0, a1..a5
; t is time
    n = n_elements(dx)
    model = dblarr(n) + p[0]

    scale = p[1] * dx + p[2] * dy + p[3] * dx * dy + p[4] * dx * dx + $
            p[5] * dy * dy 
    dt = t - min(t)
    if (n_elements(p) ge 7) then scale2 = 1.0 + p[6] * dt else scale2 = 1.0
    if (n_elements(p) eq 8) then scale2 = scale2 + p[7] * dt * dt
    
    model = model * (1. + scale) * scale2
    model = (flux - model) / err

return, model
end

function fpa1_xfunc3, p, FLUX=flux, DX=dx, DY=dy, T=t, ERR=err
; p is parameters, f0, a1..a5
; t is time
    n = n_elements(dx)
    model = dblarr(n) + p[0]

    scale = p[1] * dx + p[2] * dy + p[3] * dx * dy + p[4] * dx * dx + $
            p[5] * dy * dy 
    dt = t - min(t)
    
    model = model * (1. + scale) 
    
    ;now add in the sin curve (XXX)
    model = model + p[6]*sin(t/p[7] + p[8]) 
 

   model = (flux - model) / err

return, model
end
