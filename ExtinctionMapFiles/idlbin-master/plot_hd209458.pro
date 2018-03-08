pro plot_hd209458,ch, binning = binning
 if ch eq 2 then begin
    restore,  '/Users/jkrick/irac_warm/hd209458/hd209458_ch2.sav'
    nfits = 12429.
    nframes = 63.
    yyrange = [14,15]
    xyrange =  [14.5, 15.0]
    fyrange = [0.47,0.50]
    add = 0.48
    titlename = 'Ch2'
    xrange = xyrange
    yrange = yyrange
endif
 if ch eq 1 then begin
    restore,  '/Users/jkrick/irac_warm/hd209458/hd209458_ch1.sav'
    nfits = 40633.
    nframes = 63.
    yyrange = [14.0,15]
    xyrange =  [14.2, 15.2]
    fyrange = [0.7,0.83]
    add = 0.73
    titlename = 'Ch1'
    xrange = xyrange
    yrange = yyrange
 endif
 if ch eq 4 then begin
    restore,  '/Users/jkrick/irac_warm/hd209458/hd209458_ch4.sav'
    nfits = 6870.
    nframes = 63.
    yyrange = [13.8,14.8]
    xyrange =  [14.5, 15.5]
    fyrange = [0.16,0.17]
    add = 0.162
    ;make these up so the code doesn't break
    sortcorrfluxarr = sortfluxarr
    sortcorrfluxerrarr = sortfluxerrarr
    titlename = 'Ch4'
    xrange = [14.56, 14.72]
    yrange = [13.8, 14.8]
 endif

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
 ;binning
  if keyword_set(binning) then begin
     
;try binning the fluxes by subarray image
     bin_flux = dblarr(nfits)
     bin_fluxerr = dblarr(nfits)
     bin_corrflux= dblarr(nfits)
     bin_corrfluxerr= dblarr(nfits)
     bin_time= dblarr(nfits)
     bin_xcen= dblarr(nfits)
     bin_ycen= dblarr(nfits)
     for si = 0., nfits - 1 do begin
                                ;print, 'working on si', si
        idata = sortfluxarr[si*nframes:si*nframes + (nframes - 1)]
        idataerr = sortfluxerrarr[si*nframes:si*nframes + (nframes - 1)]
        bin_flux[si] = mean(idata,/nan)
        bin_fluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        bin_corrflux[si] = mean(sortcorrfluxarr[si*nframes:si*nframes + (nframes - 1)],/nan)
        bin_corrfluxerr[si] = sqrt(total((sortcorrfluxerrarr[si*nframes:si*nframes + (nframes - 1)])^2))/ (n_elements(sortcorrfluxerrarr[si*nframes:si*nframes + (nframes - 1)]))
        bin_time[si] = mean(sorttime[si*nframes:si*nframes + (nframes - 1)])
        bin_xcen[si] = mean(sortxarr[si*nframes:si*nframes + (nframes - 1)])
        bin_ycen[si] = mean(sortyarr[si*nframes:si*nframes + (nframes - 1)])
        
           
     endfor                     ;for each fits image
;rename to fit with the nonbinning nomenclature
     sorttime = bin_time
     sortxarr = bin_xcen
     sortyarr = bin_ycen
     sortfluxarr = bin_flux
     sortfluxerrarr = bin_fluxerr
     sortcorrfluxarr = bin_corrflux
     sortcorrfluxerrarr = bin_corrfluxerr
     
  endif                         ;binning
  
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;get rid of outliers
  print, 'nflux', n_elements(sortfluxarr)

  good = where(sortxarr lt mean(sortxarr) + 2.5*stddev(sortxarr) and sortxarr gt mean(sortxarr) -2.5*stddev(sortxarr) and sortyarr lt mean(sortyarr) +3.0*stddev(sortyarr) and sortyarr gt mean(sortyarr) - 3.0*stddev(sortyarr),ngood_pmap, complement=bad) 
  print, 'bad ',n_elements(bad)
  xarr = sortxarr[good]
  yarr = sortyarr[good]
  timearr = sorttime[good]
  fluxerr = sortfluxerrarr[good]
  flux = sortfluxarr[good]
  corrflux = sortcorrfluxarr[good]
  corrfluxerr = sortcorrfluxerrarr[good]
 
  print, 'max x, y', max(xarr), min(xarr), max(yarr), min(yarr)

  ;try getting rid of flux outliers.
  ;do some running mean with clipping
  start = 0
  print, 'nflux', n_elements(flux)
  for ni = 100, n_elements(flux) -1,100 do begin
     meanclip,flux[start:ni], m, s, subs = subs;,/verbose
     ;print, 'good', subs+start
     ;now keep a list of the good ones
     if ni eq 100 then good_ni = subs+start else good_ni = [good_ni, subs+start]
     start = ni
  endfor

  ;see if it worked
  xarr = xarr[good_ni]
  yarr = yarr[good_ni]
  timearr = timearr[good_ni]
  fluxerr = fluxerr[good_ni]
  flux = flux[good_ni]
  corrflux = corrflux[good_ni]
  corrfluxerr = corrfluxerr[good_ni]

  print, 'nflux', n_elements(flux)

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;the plotting
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;sorttime = timearr[a]
;sortxarr = xarr[a]
;sortyarr = yarr[a]
;sortfluxarr = fluxarr[a]
;sortfluxerrarr = fluxerrarr[a]
;sortcorrfluxarr = corrfluxarr[a]
;sortcorrfluxerrarr = corrfluxerrarr[a]

;put time in hours
  timearr = (timearr - timearr(0))/60./60.
  
  ;now make the plots
  ;-----------------
  ;y vs. time
  st = plot(timearr, yarr,'1o',xtitle = 'Time(hours)',ytitle = 'Y pix',   color = 'black', sym_filled = 1, sym_size = 0.4, yrange = yyrange, title = titlename) 
  ;-----------------
  ;x vs time
  st2 = plot(timearr, xarr,'1o', xtitle = 'Time(hours)',ytitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4, yrange =xyrange, title = titlename) ;
  ;-----------------
  ;x vs y
  st2 = plot(xarr, yarr,'1o', ytitle = 'Y pix',xtitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4, title = titlename, yrange = yrange, xrange = xrange) ;
   ;-----------------
  ;flux vs time
   ;st3 = plot(timearr, flux,'1-', xtitle = 'Time(hours)',ytitle = 'Flux (Jy)',  title = 'Pmap', color = 'light gray', yrange = [0.4,0.6] ); , xrange = [2, 4]) ;
  ;st3 = plot(timearr, flux,'1o',  color = 'black' , sym_filled = 1, sym_size = 0.3 , xtitle = 'Time(hours)',ytitle = 'Flux (Jy)',yrange = fyrange) ; , xrange = [2, 4])                  ;
 ;     st4 = plot(timearr, corrflux +0.01, '1-',  color = 'light gray',/overplot) ; 
 ;     st4 = plot(timearr, corrflux +0.01, '1o',  color = 'red', sym_filled = 1, sym_size = 0.3 ,/overplot) ; 
 ;-----------------
 ;Normalized flux vs time
 ; st6 = plot(timearr, corrflux/corrflux(0), '1-o',  color = 'red', sym_filled = 1, sym_size = 0.4,NAME = 'HD7924', yrange = [0.998, 1.0023], xrange = [0,8],xtitle = 'Time (hours)', ytitle = 'Normalized Flux ', title = 'PMAP')
;  st6 = errorplot(timearr, corrflux/median(corrflux), corrfluxerr/median(corrflux), '1-o', color = 'red', sym_filled = 1, sym_size = 0.4,NAME = 'pmap', xtitle = 'Time (hours)', ytitle = 'Normalized Flux ', title = 'HD7924')
;this is the uncorrected flux
;  p1 = errorplot(timearr, flux/median(flux), fluxerr/median(flux), '1-o', color = 'black', sym_size = 0.4, sym_filled = 1., NAME = 'uncorrected', ytitle = 'Normalized Flux ', xtitle = 'Time (hours)', title = 'HD209458') ;, xrange = [2,4]

    
;--------------------------------------------------------------------------
;--------------------------------------------------------------------------
;now work on the data fitting with a polynomial for comparison:
;can I get rid of pixel phase effect?

  x = timearr                   ;/ max(timearr)        ; fix to 1 phase
  y = flux                   ;/ fluxarr[1]          ;normalize
  yerr = fluxerr            ;/ fluxarr[1]    ;normalize
  print, 'n x, y, yerr', n_elements(x), n_elements(y), n_elements(yerr)
;the plots!
 ; p1 = errorplot(x, y, yerr, '1-s',  xtitle = 'Time(hours)', ytitle = ' Flux', title = 'Data Fit', errorbar_capsize = .05 , sym_size = 0.4, sym_filled = 1.,yrange = [1.4, 1.5]) ;, xrange = [2,4]
  p1 = plot(x, y,  '1s',  xtitle = 'Time(hours)', ytitle = ' Flux',  sym_size = 0.4, sym_filled = 1.,yrange = fyrange, title = titlename);, xrange = [2,4]

; Initial guess	
;  pa0 = [median(y), 1.,1.,1.,1.,1.,0.006,20.618,0.5]
  pa0 = [median(y), 1.,1.,1.,1.,1.]
  amx = median(xarr)            ;xa[agptr])
  amy = median(yarr)            ; ya[agptr])	
  adx = xarr - amx              ; xa[agptr] - amx
  ady = yarr - amy              ; ya[agptr] - amy 
  
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(pa0))
                                ;limit the range of the sin cuve phase
  ;parinfo[8].fixed = 1 
  ;parinfo[7].fixed = 1          
  ;parinfo[6].limited[0] = 1
  ;parinfo[6].limits[0] = 0.0
  

  afargs = {FLUX:y, DX:adx, DY:ady, T:x, ERR:yerr}
  pa = mpfit('fpa1_xfunc2', pa0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg, parinfo = parinfo) ;,/quiet)
  print, 'status', status
  print, errmsg
  achi = achi / adof
  print, 'reduced chi squared', achi
  
  
  f1 = plot(x,   pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + $
                          pa[5] * ady * ady )  , '6b1.',/overplot)
  
  ;plot residuals
  sub = y - (pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + $
                      pa[5] * ady * ady ) )
  print, 'mean sub', mean(sub) 

    ;sub = y - (pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + $
  ;                           pa[5] * ady * ady ) + pa[6]*sin(x/pa[7] + pa[8]))
  ;f2 = errorplot(x, sub+1.427, yerr, '1-sr',/overplot ,errorbar_capsize = .05 , sym_size = 0.4, sym_filled = 1.)
  f2 = plot(x, sub+add, '1sr',/overplot , sym_size = 0.4, sym_filled = 1.)
  
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

