pro plot_hd209458_v2,ch, binning = binning, sine=sine
period = 3.524749

 if ch eq 2 then begin
    restore,  '/Users/jkrick/irac_warm/hd209458/hd209458_ch2.sav'
    yyrange = [14,15]
    xyrange =  [14.5, 15.0]
    fyrange = [0.47,0.50]
    add = 0.48
    titlename = 'Ch2'
    xrange = xyrange
    yrange = yyrange
    tt = 48.5
    timestart = 1800.
 endif
 if ch eq 1 then begin
    restore,  '/Users/jkrick/irac_warm/hd209458/hd209458_ch1_new.sav'
    yyrange = [14.0,15]
    xyrange =  [14.2, 15.2]
    fyrange = [0.7,0.83]
    add = 0.73
    titlename = 'Ch1'
    xrange = xyrange
    yrange = yyrange
    tt = 52.
    timestart = 0.
 endif
 if ch eq 4 then begin
    restore,  '/Users/jkrick/irac_warm/hd209458/hd209458_ch4.sav'
    yyrange = [13.8,14.8]
    xyrange =  [14.5, 15.5]
    fyrange = [0.16,0.17]
    add = 0.163
    ;make these up so the code doesn't break
    sortcorrfluxarr = sortfluxarr
    sortcorrfluxerrarr = sortfluxerrarr
    titlename = 'Ch4'
    xrange = [14.5, 14.8]
    yrange = [13.8, 14.8]
    tt = 6.0
    timestart = 3000.
 endif

 
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;get rid of outliers
 print, 'start sortfluxarr', n_elements(sortfluxarr)
 ;chop the first 30 min.
 good = where(sorttime - sorttime(0) gt timestart)
 sortxarr = sortxarr[good]
 sortyarr = sortyarr[good]
 sorttime = sorttime[good]
 sortfluxerr = sortfluxerrarr[good]
 sortfluxarr = sortfluxarr[good]
 sortcorrfluxarr = sortcorrfluxarr[good]
 sortcorrfluxerrarr = sortcorrfluxerrarr[good]
 sortbmjd = sortbmjd[good]
 sortutcs = sortutcs[good]
 print, 'stmiddleart sortfluxarr', n_elements(sortfluxarr)

 ; the position outliers
  good = where(sortxarr lt mean(sortxarr) +3.0*stddev(sortxarr) and sortxarr gt mean(sortxarr) -3.0*stddev(sortxarr) and sortyarr lt mean(sortyarr) +3.0*stddev(sortyarr) and sortyarr gt mean(sortyarr) - 3.0*stddev(sortyarr),ngood_pmap, complement=bad) 
  print, 'bad ',n_elements(bad), n_elements(good)
  xarr = sortxarr[good]
  yarr = sortyarr[good]
  timearr = sorttime[good]
  fluxerr = sortfluxerrarr[good]
  flux = sortfluxarr[good]
  corrflux = sortcorrfluxarr[good]
  corrfluxerr = sortcorrfluxerrarr[good]
  bmjd = sortbmjd[good]
  utcs = sortutcs[good]
 print, 'middle flux', n_elements(flux)

  print, 'max x, y', max(xarr), min(xarr), max(yarr), min(yarr)

  ;try getting rid of flux outliers.
  ;do some running mean with clipping
  start = 0
  print, 'nflux', n_elements(flux)
  for ni = 100, n_elements(flux) -1,100 do begin
     meanclip,flux[start:ni], m, s, subs = subs;,/verbose
    ; print, 'good', subs+start
     ;now keep a list of the good ones
     if ni eq 100 then good_ni = subs+start else good_ni = [good_ni, subs+start]
     start = ni + 1
  endfor
  print, 'n good_ni', n_elements(good_ni)

  ;see if it worked
  xarr = xarr[good_ni]
  yarr = yarr[good_ni]
  timearr = timearr[good_ni]
  fluxerr = fluxerr[good_ni]
  flux = flux[good_ni]
  corrflux = corrflux[good_ni]
  corrfluxerr = corrfluxerr[good_ni]
  bmjd = bmjd[good_ni]
  utcs = utcs[good_ni]
  print, 'nflux', n_elements(flux)

;-------------------------------------------------------------------------------------
  ;now find the factors for binning
  ;nframes = find_factors(n_elements(flux), 60)
  ;nfits = n_elements(flux) / nframes
  
  ;instead just pick my own binning level so that it is consistent between channels
  nframes = 63L
  nfits = long(n_elements(flux)) / nframes
  print, 'nframes, nfits', nframes, nfits
  help, nfits

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;the plotting

;put time in hours
  timearr = (timearr - timearr(0))/60./60.
  
  ;now make the plots
  ;-----------------
  ;y vs. time
  ;st = plot(timearr, yarr,'1o',xtitle = 'Time(hours)',ytitle = 'Y pix',   color = 'black', sym_filled = 1, sym_size = 0.4, yrange = yyrange, title = titlename) 
  ;-----------------
  ;x vs time
  ;st2 = plot(timearr, xarr,'1o', xtitle = 'Time(hours)',ytitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4, yrange =xyrange, title = titlename) ;
  ;-----------------
  ;x vs y
 ; st2 = plot(xarr, yarr,'1o', ytitle = 'Y pix',xtitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4, title = titlename, yrange = yrange, xrange = xrange) ;
   ;-----------------
;-------------------------------------------------------------------------------------
;can I find the breaks in time?
  count = 0
  print, 'starting breaks', timearr(10) - timearr(9)
  breakpt = fltarr(4)
  for ti = 0, n_elements(timearr) - 5 do begin
     if timearr(ti+1) - timearr(ti) gt 1. then begin
        print, 'found a break', timearr(ti)
        breakpt[count] = ti
        count = count + 1
     endif
  endfor
  breakpt = breakpt[0:count-1]
  nsegments = count - 1

;test this
  print,'testing timearr',  timearr[breakpt[0]], timearr[breakpt[0]+1], timearr[breakpt[1]], timearr[breakpt[1] + 1], timearr[breakpt[2]], timearr[breakpt[2]+1]

  ;now make separate arrays
 xarr_1 = xarr[0:breakpt[0]]
 xarr_2 = xarr[breakpt[0]:breakpt[1]]
 xarr_3 = xarr[breakpt[1]:breakpt[2]]
 xarr_4 = xarr[breakpt[2]:breakpt[3]]
 xarr_5 = xarr[breakpt[3]:*]

 yarr_1 = yarr[0:breakpt[0]]
 yarr_2 = yarr[breakpt[0]:breakpt[1]]
 yarr_3 = yarr[breakpt[1]:breakpt[2]]
 yarr_4 = yarr[breakpt[2]:breakpt[3]]
 yarr_5 = yarr[breakpt[3]:*]

 timearr_1 = timearr[0:breakpt[0]]
 timearr_2 = timearr[breakpt[0]:breakpt[1]]
 timearr_3 = timearr[breakpt[1]:breakpt[2]]
 timearr_4 = timearr[breakpt[2]:breakpt[3]]
 timearr_5 = timearr[breakpt[3]:*]

 flux_1 = flux[0:breakpt[0]]
 flux_2 = flux[breakpt[0]:breakpt[1]]
 flux_3 = flux[breakpt[1]:breakpt[2]]
 flux_4 = flux[breakpt[2]:breakpt[3]]
 flux_5 = flux[breakpt[3]:*]

 fluxerr_1 = fluxerr[0:breakpt[0]]
 fluxerr_2 = fluxerr[breakpt[0]:breakpt[1]]
 fluxerr_3 = fluxerr[breakpt[1]:breakpt[2]]
 fluxerr_4 = fluxerr[breakpt[2]:breakpt[3]]
 fluxerr_5 = fluxerr[breakpt[3]:*]

 bmjd_1 = bmjd[0:breakpt[0]]
 bmjd_2 = bmjd[breakpt[0]:breakpt[1]]
 bmjd_3 = bmjd[breakpt[1]:breakpt[2]]
 bmjd_4 = bmjd[breakpt[2]:breakpt[3]]
 bmjd_5 = bmjd[breakpt[3]:*]

 utcs_1 = utcs[0:breakpt[0]]
 utcs_2 = utcs[breakpt[0]:breakpt[1]]
 utcs_3 = utcs[breakpt[1]:breakpt[2]]
 utcs_4 = utcs[breakpt[2]:breakpt[3]]
 utcs_5 = utcs[breakpt[3]:*]

;--------------------------------------------------------------------------
;--------------------------------------------------------------------------
;now work on the data fitting with a polynomial for comparison:
;can I get rid of pixel phase effect?

  x = timearr                   ;/ max(timearr)       
  y = flux                  ;/ fluxarr[1]          ;normalize
  yerr = fluxerr            ;/ fluxarr[1]    ;normalize
  amx = median(xarr)        ;xa[agptr])
  amy = median(yarr)        ; ya[agptr])	
  adx = xarr - amx          ; xa[agptr] - amx
  ady = yarr - amy          ; ya[agptr] - amy 

  if ch eq 1 then begin
     x = timearr_3
     y = flux_3
     yerr = fluxerr_3
     amx = median(xarr_3)         
     amy = median(yarr_3)        
     adx = xarr_3 - amx           
     ady = yarr_3 - amy    
     nfits = n_elements(x) / nframes
     xarr = xarr_3
     yarr = yarr_3
     bmjd = bmjd_3
     utcs= utcs_3
  endif


  print, 'n x, y, yerr', n_elements(x), n_elements(y), n_elements(yerr), n_elements(flux), nfits
;the plots!
 ; p1 = errorplot(x, y, yerr, '1-s',  xtitle = 'Time(hours)', ytitle = ' Flux', title = 'Data Fit', errorbar_capsize = .05 , sym_size = 0.4, sym_filled = 1.,yrange = [1.4, 1.5]) ;, xrange = [2,4]
  ;p1 = plot(x, y,  '1s',  xtitle = 'Time(hours)', ytitle = ' Flux',  sym_size = 0.4, sym_filled = 1.,yrange = fyrange, title = titlename);, xrange = [2,4]

; Initial guesses	
  if keyword_set(sine) then begin  ;use a fitting function with a sine curve or not
     known_period = 84.593976; hours  304538.3 ;seconds
     pa0 = [median(y), 1.,1.,1.,1.,1.,0.001,known_period*(2*!Pi),0.9]
     func = 'fpa1_xfunc2'
  endif else begin
     pa0 = [median(y), 1.,1.,1.,1.,1.]
     func = 'fpa1_xfunc3'

  endelse

  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(pa0))
   ;limit the range of the sin cuve phase
  if keyword_set(sine) then begin  ;use a fitting function with a sine curve
     parinfo[7].fixed = 1          
     parinfo[6].limited[0] = 1
     parinfo[6].limits[0] = 0.0
  endif


  afargs = {FLUX:y, DX:adx, DY:ady, T:x, ERR:yerr}
  pa = mpfit('fpa1_xfunc2', pa0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg, parinfo = parinfo) ;,/quiet)
  print, 'status', status
  print, errmsg
  achi = achi / adof
  print, 'reduced chi squared', achi
  if keyword_set(sine) then begin  ;use a fitting function with a sine curve or not
     model_fit = pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + pa[5] * ady * ady)+ pa[6]*sin(x/pa[7] + pa[8]) 
  endif else begin
     model_fit = pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + pa[5] * ady * ady)
  endelse

  
  ;plot residuals
  sub = y - model_fit
 ; f1 = plot(x,   model  , '6b1.',/overplot)
  ;f2 = plot(x, sub+add, '1sr',/overplot , sym_size = 0.4, sym_filled = 1.)

  print, 'nsub', n_elements(sub)
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
 ;binning after fitting
  if keyword_set(binning) then begin
     
;try binning the fluxes by subarray image
     bin_flux = dblarr(nfits)
     bin_fluxerr = dblarr(nfits)
     ;bin_corrflux= dblarr(nfits)
     ;bin_corrfluxerr= dblarr(nfits)
     bin_time= dblarr(nfits)
     bin_xcen= dblarr(nfits)
     bin_ycen= dblarr(nfits)
     bin_bmjd = dblarr(nfits)
     bin_utcs = dblarr(nfits)
     bin_sub = dblarr(nfits)
     bin_model = dblarr(nfits)

     print, 'nfits', nfits
     help, nfits
     for si = 0L, long(nfits) - 1L do begin
      ;  print, 'working on si', si, n_elements(y), si*nframes, si*nframes + (nframes - 1)
        idata = y[si*nframes:si*nframes + (nframes - 1)]
        idataerr = yerr[si*nframes:si*nframes + (nframes - 1)]
        bin_flux[si] = mean(idata,/nan)
        bin_fluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        ;bin_corrflux[si] = mean(corrflux[si*nframes:si*nframes + (nframes - 1)],/nan)
        ;bin_corrfluxerr[si] = sqrt(total((corrfluxerr[si*nframes:si*nframes + (nframes - 1)])^2))/ (n_elements(corrfluxerr[si*nframes:si*nframes + (nframes - 1)]))
        bin_time[si] = mean(x[si*nframes:si*nframes + (nframes - 1)])
        bin_xcen[si] = mean(xarr[si*nframes:si*nframes + (nframes - 1)])
        bin_ycen[si] = mean(yarr[si*nframes:si*nframes + (nframes - 1)])
        bin_bmjd[si] = mean(bmjd[si*nframes:si*nframes + (nframes - 1)])
        bin_utcs[si] = mean(utcs[si*nframes:si*nframes + (nframes - 1)])
        bin_sub[si] = mean(sub[si*nframes:si*nframes + (nframes - 1)])
        bin_model[si] =  mean(model_fit[si*nframes:si*nframes + (nframes - 1)])
     endfor                     ;for each fits image
     
  
;-------------------------------------------------------------------------------------
 ;plot the binned stuff

;set time = 0 in the middle of the transit
     bin_time = bin_time -tt

;set time in orbital phase from 0 at transit to 1 at next transit
     phase_time = bin_time / (period*24.)

  ;y vs. time
     st = plot(bin_time, bin_ycen,'1o',xtitle = 'Time from transit (hours)',ytitle = 'Y pix',   color = 'black', sym_filled = 1, sym_size = 0.4, yrange = yyrange, title = titlename) ;, axis_style = 1) 
;tjunk = plot(bin_time,bin_ycen,/current, xrange = st.xrange, yrange = yyrange, nodata = 1)
;topx = axis('X', location = [0, 14.8],  target = tjunk, title = 'Time from transit (hours)')

  ;-----------------
  ;x vs time
     st1 = plot(bin_time, bin_xcen,'1o', xtitle = 'Time from transit (hours)',ytitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4, yrange =xyrange, title = titlename) ;
  ;-----------------
  ;x vs y
     st2 = plot(bin_xcen, bin_ycen,'1o', ytitle = 'Y pix',xtitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4, title = titlename, yrange = yrange, xrange = xrange) ;
   ;-----------------

     x = phase_time             ;/ max(timearr)        ; fix to 1 phase
     y = bin_flux               ;/ fluxarr[1]          ;normalize
     yerr = bin_fluxerr         ;/ fluxarr[1]    ;normalize
     print, 'n x, y, yerr', n_elements(x), n_elements(y), n_elements(yerr)
     p1 = plot(x, y,  '1s',  xtitle = 'Orbital Phase', ytitle = ' Flux',  sym_size = 0.2, sym_filled = 1.,yrange = fyrange, title = titlename, xrange = [-0.7,0.7]) ;, xrange = [2,4]
     f1 = plot(x,bin_model,  '1sb',/overplot , sym_size = 0.2, sym_filled = 1.)
     f2 = plot(x, bin_sub+add, '1sr',/overplot , sym_size = 0.2, sym_filled = 1.)

     p1.save, "/Users/jkrick/irac_warm/hd209458/fluxvsphase_"+titlename +".png"
     st.save, "/Users/jkrick/irac_warm/hd209458/yvstime_"+titlename +".png"
     st1.save,  "/Users/jkrick/irac_warm/hd209458/xvstime_"+titlename +".png"
     st2.save,  "/Users/jkrick/irac_warm/hd209458/xvsy_"+titlename +".png"
     
  endif                         ;binning

;save, /all, filename = '/Users/jkrick/irac_warm/hd209458/hd209458_plot_'+titlename+'.sav'

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
    
    ;now add in the sin curve 
    model = model + p[6]*sin(t/p[7] + p[8]) 
 

   model = (flux - model) / err

return, model
end

function funcwithtransit, p, FLUX=flux, DX=dx, DY=dy, T=t, ERR=err
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
 
    ;need to add in parameters for the transit and eclipse



   model = (flux - model) / err

return, model
end
