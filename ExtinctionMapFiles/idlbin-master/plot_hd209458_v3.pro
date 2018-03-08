pro plot_hd209458_v3,ch, binning = binning, sine=sine
period = 3.524749

 if ch eq 2 then begin
    restore,  '/Users/jkrick/irac_warm/hd209458/hd209458_ch2.sav'
    yyrange = [14,15]
    xyrange =  [14.5, 15.0]
    fyrange = [0.475,0.50]
    add = 0.48
    titlename = 'Ch2'
    xrange = xyrange
    yrange = yyrange
    tt = 48.5
    timestart = 1800.
 endif
 if ch eq 1 then begin
    restore,  '/Users/jkrick/irac_warm/hd209458/hd209458_ch1.sav'
    yyrange = [14.0,15]
    xyrange =  [14.2, 15.2]
    fyrange = [0.72,0.83]
    add = 0.74
    titlename = 'Ch1'
    xrange = xyrange
    yrange = yyrange
    tt = 52.
    timestart = 0; 13000.
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
    timestart = 16000.
 endif

 
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;get rid of outliers
 print, 'start fluxarr', n_elements(sortfluxarr)
 ;chop the first 30 min.
; good = where(sorttime - sorttime(0) gt timestart)
; sortxarr = sortxarr[good]
; sortyarr = sortyarr[good]
; sorttime = sorttime[good]
; sortfluxerr = sortfluxerrarr[good]
; sortfluxarr = sortfluxarr[good]
; sortcorrfluxarr = sortcorrfluxarr[good]
; sortcorrfluxerrarr = sortcorrfluxerrarr[good]
; sortbmjd = sortbmjd[good]
; sortutcs = sortutcs[good]
; print, 'stmiddleart sortfluxarr', n_elements(sortfluxarr)

 ; the position outliers
  good = where(sortxarr lt mean(sortxarr) +3.0*stddev(sortxarr) and sortxarr gt mean(sortxarr) -3.0*stddev(sortxarr) and sortyarr lt mean(sortyarr) +3.0*stddev(sortyarr) and sortyarr gt mean(sortyarr) - 3.0*stddev(sortyarr),ngood_pmap, complement=bad) 
  print, 'bad position',n_elements(bad), n_elements(good)
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

  ;print, 'max x, y', max(xarr), min(xarr), max(yarr), min(yarr)

  ;try getting rid of flux outliers.
  ;do some running mean with clipping
 ;also cut the transit
  start = 0
  print, 'nflux', n_elements(flux)
  for ni = 100, n_elements(flux) -1,100 do begin
     meanclip,flux[start:ni], m, s, subs = subs;,/verbose
    ; print, 'good', subs+start
     ;now keep a list of the good ones
     if ni eq 100 then good_ni = subs+start else good_ni = [good_ni, subs+start]
     start = ni + 1
  endfor
  print, 'n good fluxclip', n_elements(good_ni)

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


;put time in hours
  timearr = (timearr - timearr(0))/60./60.
    print, timearr[0:100], format = '(F0)'


;try chopping out the transit
  ;transit = where(flux lt 0.
;can I find the edges of the transit?
 

;-------------------------------------------------------------------------------------
  ;now find the factors for binning
  ;nframes = find_factors(n_elements(flux), 60)
  ;nfits = n_elements(flux) / nframes
  
  ;instead just pick my own binning level so that it is consistent between channels
  bin_level = 120L
  nframes = bin_level
  nfits = long(n_elements(flux)) / nframes
  print, 'nframes, nfits', nframes, nfits
  help, nfits

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
;the plotting



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
print, 'ntime, flux', n_elements(timearr), n_elements(flux), n_elements(xarr)

  j2 = plot(timearr, flux)
;-------------------------------------------------------------------------------------


;  if ch eq 2 then begin
;     ;add the transit as a break in the light curve
;     starttransit = max(where(timearr lt 46.))
;     stoptransit = min(where(timearr gt 50.1))
;     ;flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
;     xarr = [xarr[0:starttransit],xarr[stoptransit:n_elements(timearr)-1]]
;     yarr = [yarr[0:starttransit],yarr[stoptransit:n_elements(timearr)-1]]
;     fluxerr = [fluxerr[0:starttransit],fluxerr[stoptransit:n_elements(timearr)-1]]
;     flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
;     corrflux = [corrflux[0:starttransit],corrflux[stoptransit:n_elements(timearr)-1]]
;     corrfluxerr = [corrfluxerr[0:starttransit],corrfluxerr[stoptransit:n_elements(timearr)-1]]
;     bmjd = [bmjd[0:starttransit],bmjd[stoptransit:n_elements(timearr)-1]]
;     utcs = [utcs[0:starttransit],utcs[stoptransit:n_elements(timearr)-1]]
;     timearr = [timearr[0:starttransit],timearr[stoptransit:n_elements(timearr)-1]]
;
;     starteclipse= max(where(timearr lt 89.))
;     stopeclipse = min(where(timearr gt 93.1))
;     xarr = [xarr[0:starteclipse],xarr[stopeclipse:n_elements(timearr)-1]]
;     yarr = [yarr[0:starteclipse],yarr[stopeclipse:n_elements(timearr)-1]]
;     fluxerr = [fluxerr[0:starteclipse],fluxerr[stopeclipse:n_elements(timearr)-1]]
;     flux = [flux[0:starteclipse],flux[stopeclipse:n_elements(timearr)-1]]
;     corrflux = [corrflux[0:starteclipse],corrflux[stopeclipse:n_elements(timearr)-1]]
;     corrfluxerr = [corrfluxerr[0:starteclipse],corrfluxerr[stopeclipse:n_elements(timearr)-1]]
;     bmjd = [bmjd[0:starteclipse],bmjd[stopeclipse:n_elements(timearr)-1]]
;     utcs = [utcs[0:starteclipse],utcs[stopeclipse:n_elements(timearr)-1]]
;     timearr = [timearr[0:starteclipse],timearr[stopeclipse:n_elements(timearr)-1]]
;
;     starteclipse= max(where(timearr lt 4.4))
;     stopeclipse = min(where(timearr gt 8.5))
;     xarr = [xarr[0:starteclipse],xarr[stopeclipse:n_elements(timearr)-1]]
;     yarr = [yarr[0:starteclipse],yarr[stopeclipse:n_elements(timearr)-1]]
;     fluxerr = [fluxerr[0:starteclipse],fluxerr[stopeclipse:n_elements(timearr)-1]]
;     flux = [flux[0:starteclipse],flux[stopeclipse:n_elements(timearr)-1]]
;     corrflux = [corrflux[0:starteclipse],corrflux[stopeclipse:n_elements(timearr)-1]]
;     corrfluxerr = [corrfluxerr[0:starteclipse],corrfluxerr[stopeclipse:n_elements(timearr)-1]]
;     bmjd = [bmjd[0:starteclipse],bmjd[stopeclipse:n_elements(timearr)-1]]
;     utcs = [utcs[0:starteclipse],utcs[stopeclipse:n_elements(timearr)-1]]
;     timearr = [timearr[0:starteclipse],timearr[stopeclipse:n_elements(timearr)-1]]
;  endif
;  
; if ch eq 4 then begin
;     ;add the transit as a break in the light curve
;     starttransit = max(where(timearr lt 4.5))
;     stoptransit = min(where(timearr gt 8.0))
;     ;flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
;     xarr = [xarr[0:starttransit],xarr[stoptransit:n_elements(timearr)-1]]
;     yarr = [yarr[0:starttransit],yarr[stoptransit:n_elements(timearr)-1]]
;     fluxerr = [fluxerr[0:starttransit],fluxerr[stoptransit:n_elements(timearr)-1]]
;     flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
;     corrflux = [corrflux[0:starttransit],corrflux[stoptransit:n_elements(timearr)-1]]
;     corrfluxerr = [corrfluxerr[0:starttransit],corrfluxerr[stoptransit:n_elements(timearr)-1]]
;     bmjd = [bmjd[0:starttransit],bmjd[stoptransit:n_elements(timearr)-1]]
;     utcs = [utcs[0:starttransit],utcs[stoptransit:n_elements(timearr)-1]]
;     timearr = [timearr[0:starttransit],timearr[stoptransit:n_elements(timearr)-1]]
;
;     starteclipse= max(where(timearr lt 46.8))
;     stopeclipse = min(where(timearr gt 50.3))
;     xarr = [xarr[0:starteclipse],xarr[stopeclipse:n_elements(timearr)-1]]
;     yarr = [yarr[0:starteclipse],yarr[stopeclipse:n_elements(timearr)-1]]
;     fluxerr = [fluxerr[0:starteclipse],fluxerr[stopeclipse:n_elements(timearr)-1]]
;     flux = [flux[0:starteclipse],flux[stopeclipse:n_elements(timearr)-1]]
;     corrflux = [corrflux[0:starteclipse],corrflux[stopeclipse:n_elements(timearr)-1]]
;     corrfluxerr = [corrfluxerr[0:starteclipse],corrfluxerr[stopeclipse:n_elements(timearr)-1]]
;     bmjd = [bmjd[0:starteclipse],bmjd[stopeclipse:n_elements(timearr)-1]]
;     utcs = [utcs[0:starteclipse],utcs[stopeclipse:n_elements(timearr)-1]]
;     timearr = [timearr[0:starteclipse],timearr[stopeclipse:n_elements(timearr)-1]]
;  endif
;
; if ch eq 1 then begin
;     ;add the transit as a break in the light curve
;     starttransit = max(where(timearr lt 50.5))
;     stoptransit = min(where(timearr gt 54.0))
;     ;flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
;     xarr = [xarr[0:starttransit],xarr[stoptransit:n_elements(timearr)-1]]
;     yarr = [yarr[0:starttransit],yarr[stoptransit:n_elements(timearr)-1]]
;     fluxerr = [fluxerr[0:starttransit],fluxerr[stoptransit:n_elements(timearr)-1]]
;     flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
;     corrflux = [corrflux[0:starttransit],corrflux[stoptransit:n_elements(timearr)-1]]
;     corrfluxerr = [corrfluxerr[0:starttransit],corrfluxerr[stoptransit:n_elements(timearr)-1]]
;     bmjd = [bmjd[0:starttransit],bmjd[stoptransit:n_elements(timearr)-1]]
;     utcs = [utcs[0:starttransit],utcs[stoptransit:n_elements(timearr)-1]]
;     timearr = [timearr[0:starttransit],timearr[stoptransit:n_elements(timearr)-1]]
;
;     starteclipse= max(where(timearr lt 92.7))
;     stopeclipse = min(where(timearr gt 96.3))
;     xarr = [xarr[0:starteclipse],xarr[stopeclipse:n_elements(timearr)-1]]
;     yarr = [yarr[0:starteclipse],yarr[stopeclipse:n_elements(timearr)-1]]
;     fluxerr = [fluxerr[0:starteclipse],fluxerr[stopeclipse:n_elements(timearr)-1]]
;     flux = [flux[0:starteclipse],flux[stopeclipse:n_elements(timearr)-1]]
;     corrflux = [corrflux[0:starteclipse],corrflux[stopeclipse:n_elements(timearr)-1]]
;     corrfluxerr = [corrfluxerr[0:starteclipse],corrfluxerr[stopeclipse:n_elements(timearr)-1]]
;     bmjd = [bmjd[0:starteclipse],bmjd[stopeclipse:n_elements(timearr)-1]]
;     utcs = [utcs[0:starteclipse],utcs[stopeclipse:n_elements(timearr)-1]]
;     timearr = [timearr[0:starteclipse],timearr[stopeclipse:n_elements(timearr)-1]]
;
;     starteclipse= max(where(timearr lt 8.2))
;     stopeclipse = min(where(timearr gt 11.7))
;     xarr = [xarr[0:starteclipse],xarr[stopeclipse:n_elements(timearr)-1]]
;     yarr = [yarr[0:starteclipse],yarr[stopeclipse:n_elements(timearr)-1]]
;     fluxerr = [fluxerr[0:starteclipse],fluxerr[stopeclipse:n_elements(timearr)-1]]
;     flux = [flux[0:starteclipse],flux[stopeclipse:n_elements(timearr)-1]]
;     corrflux = [corrflux[0:starteclipse],corrflux[stopeclipse:n_elements(timearr)-1]]
;     corrfluxerr = [corrfluxerr[0:starteclipse],corrfluxerr[stopeclipse:n_elements(timearr)-1]]
;     bmjd = [bmjd[0:starteclipse],bmjd[stopeclipse:n_elements(timearr)-1]]
;     utcs = [utcs[0:starteclipse],utcs[stopeclipse:n_elements(timearr)-1]]
;     timearr = [timearr[0:starteclipse],timearr[stopeclipse:n_elements(timearr)-1]]
; endif

print, 'ntime, flux', n_elements(timearr), n_elements(flux), n_elements(xarr)
;  jun = plot(timearr, flux)

;can I find the breaks in time?
  count = 0
  print, 'starting breaks', timearr(10) - timearr(9)
  breakpt = fltarr(9)
  for ti = 0, n_elements(timearr) - 5 do begin
     if timearr(ti+1) - timearr(ti) gt 1. then begin
        print, 'found a break', timearr(ti)
        breakpt[count] = ti
        count = count + 1
     endif
  endfor
  print, count, 'count'
  breakpt = breakpt[0:count-1]
  print, 'breakpt', breakpt

  ;make the last element equal to the end of the array
  if count gt 0 then begin
     breakpt = [breakpt,n_elements(timearr)-1] 
  endif else begin
     breakpt = n_elements(timearr) - 1
  endelse

  nsegments = count +1
  print, 'nsegments', nsegments
  print, 'breakpt', breakpt



;test this
;  print,'testing timearr',  timearr[breakpt[0]], timearr[breakpt[0]+1], timearr[breakpt[1]], timearr[breakpt[1] + 1], timearr[breakpt[2]], timearr[breakpt[2]+1]


;make individual segments
;at the same time cut the first 'timestart' time
for nseg = 1, nsegments do begin

 if nseg eq 1 then begin
     x = timearr[0+timestart :breakpt[0]]
     y =  flux[0+timestart :breakpt[0]]
     yerr = fluxerr[0+timestart :breakpt[0]]
     xcenarr = xarr[0+timestart :breakpt[0]]
     ycenarr = yarr[0+timestart :breakpt[0]]
     amx = median(xarr[0+timestart :breakpt[0]])         
     amy = median(yarr[0+timestart :breakpt[0]])        
     adx = xarr[0+timestart :breakpt[0]] - amx           
     ady = yarr[0+timestart :breakpt[0]] - amy    
     ;nfits = n_elements(x) / nframes
     bmjdarr =  bmjd[0+timestart :breakpt[0]]
     utcsarr= utcs[0+timestart :breakpt[0]]
  endif
 if nseg eq 2 then begin
     x = timearr[breakpt[0] +timestart :breakpt[1]]
     y = flux[breakpt[0] +timestart :breakpt[1]]
     yerr = fluxerr[breakpt[0] +timestart :breakpt[1]]
     xcenarr = xarr[breakpt[0] +timestart :breakpt[1]]
     ycenarr = yarr[breakpt[0] +timestart :breakpt[1]]
     amx = median(xcenarr)         
     amy = median(ycenarr)        
     adx = xcenarr - amx           
     ady = ycenarr - amy    
     ;nfits = n_elements(x) / nframes
     bmjdarr = bmjd[breakpt[0] +timestart :breakpt[1]]
     utcsarr=utcs[breakpt[0] +timestart :breakpt[1]]
  endif
  if nseg eq 3 then begin
     x = timearr[breakpt[1] +timestart :breakpt[2]]
     y = flux[breakpt[1] +timestart :breakpt[2]]
     yerr = fluxerr[breakpt[1] +timestart :breakpt[2]]
     xcenarr = xarr[breakpt[1] +timestart :breakpt[2]]
     ycenarr = yarr[breakpt[1] +timestart :breakpt[2]]
     amx = median(xcenarr)         
     amy = median(ycenarr)        
     adx = xcenarr - amx           
     ady = ycenarr - amy    
     ;nfits = n_elements(x) / nframes
     bmjdarr = bmjd[breakpt[1] +timestart :breakpt[2]]
     utcsarr=utcs[breakpt[1] +timestart :breakpt[2]]
  endif
 if nseg eq 4 then begin
    x = timearr[breakpt[2] +timestart :breakpt[3]]
     y = flux[breakpt[2] +timestart :breakpt[3]]
     yerr = fluxerr[breakpt[2] +timestart :breakpt[3]]
     xcenarr = xarr[breakpt[2] +timestart :breakpt[3]]
     ycenarr = yarr[breakpt[2] +timestart :breakpt[3]]
     amx = median(xcenarr)         
     amy = median(ycenarr)        
     adx = xcenarr - amx           
     ady = ycenarr - amy    
     ;nfits = n_elements(x) / nframes
     bmjdarr = bmjd[breakpt[2] +timestart :breakpt[3]]
     utcsarr=utcs[breakpt[2] +timestart :breakpt[3]]
  endif
 if nseg eq 5 then begin
     x = timearr[breakpt[3] +timestart :breakpt[4]]
     y = flux[breakpt[3] +timestart :breakpt[4]]
     yerr = fluxerr[breakpt[3] +timestart :breakpt[4]]
     xcenarr = xarr[breakpt[3] +timestart :breakpt[4]]
     ycenarr = yarr[breakpt[3] +timestart :breakpt[4]]
     amx = median(xcenarr)         
     amy = median(ycenarr)        
     adx = xcenarr - amx           
     ady = ycenarr - amy    
     ;nfits = n_elements(x) / nframes
     bmjdarr= bmjd[breakpt[3] +timestart :breakpt[4]]
     utcsarr=utcs[breakpt[3] +timestart :breakpt[4]]
  endif
 if nseg eq 6 then begin
     x = timearr[breakpt[4] +timestart :breakpt[5]]
     y = flux[breakpt[4] +timestart :breakpt[5]]
     yerr = fluxerr[breakpt[4] +timestart :breakpt[5]]
     xcenarr = xarr[breakpt[4] +timestart :breakpt[5]]
     ycenarr = yarr[breakpt[4] +timestart :breakpt[5]]
     amx = median(xcenarr)         
     amy = median(ycenarr)        
     adx = xcenarr - amx           
     ady = ycenarr - amy    
     ;nfits = n_elements(x) / nframes
     bmjdarr= bmjd[breakpt[4] +timestart :breakpt[5]]
     utcsarr=utcs[breakpt[4] +timestart :breakpt[5]]
  endif
 if nseg eq 7 then begin
     x = timearr[breakpt[5] +timestart :breakpt[6]]
     y = flux[breakpt[5] +timestart :breakpt[6]]
     yerr = fluxerr[breakpt[5] +timestart :breakpt[6]]
     xcenarr = xarr[breakpt[5] +timestart :breakpt[6]]
     ycenarr = yarr[breakpt[5] +timestart :breakpt[6]]
     amx = median(xcenarr)         
     amy = median(ycenarr)        
     adx = xcenarr - amx           
     ady = ycenarr - amy    
     ;nfits = n_elements(x) / nframes
     bmjdarr= bmjd[breakpt[5] +timestart :breakpt[6]]
     utcsarr=utcs[breakpt[5] +timestart :breakpt[6]]
  endif
 if nseg eq 8 then begin
     x = timearr[breakpt[6] +timestart :*]
     y = flux[breakpt[6] +timestart :*]
     yerr = fluxerr[breakpt[6] +timestart :*]
     xcenarr = xarr[breakpt[6] +timestart :*]
     ycenarr = yarr[breakpt[6] +timestart :*]
     amx = median(xcenarr)         
     amy = median(ycenarr)        
     adx = xcenarr - amx           
     ady = ycenarr - amy    
     ;nfits = n_elements(x) / nframes
     bmjdarr= bmjd[breakpt[6] +timestart :*]
     utcsarr=utcs[breakpt[6] +timestart :*]
  endif


  print, 'n x, y, yerr', n_elements(x), n_elements(y), n_elements(yerr), n_elements(flux), nfits
;the plots!
 ; p1 = errorplot(x, y, yerr, '1-s',  xtitle = 'Time(hours)', ytitle = ' Flux', title = 'Data Fit', errorbar_capsize = .05 , sym_size = 0.4, sym_filled = 1.,yrange = [1.4, 1.5]) ;, xrange = [2,4]
  ;p1 = plot(x, y,  '1s',  xtitle = 'Time(hours)', ytitle = ' Flux',  sym_size = 0.4, sym_filled = 1.,yrange = fyrange, title = titlename);, xrange = [2,4]

  
;;--------------------------------------------------------------------------
;;--------------------------------------------------------------------------
;now work on the data fitting with a polynomial for comparison:
;can I get rid of pixel phase effect?

; Initial guesses	
  if keyword_set(sine) then begin  ;use a fitting function with a sine curve or not
     known_period = 84.593976; hours  304538.3 ;seconds
     pa0 = [median(y), 1.,1.,1.,1.,1.,0.001,known_period*(2*!Pi),0.9]
     func = 'fpa1_xfunc3'
  endif else begin
     pa0 = [median(y), 1.,1.,1.,1.,1.]
     func = 'fpa1_xfunc2'

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

  ;is there some way to overplot the y vs. time plot with the model?
  ;XXX not working
  model_y = pa[0] * (1+  pa[2] * ady + pa[3] * adx * ady  + pa[5] * ady * ady)

  
  ;plot residuals
  sub = y - model_fit
 ; f1 = plot(x,   model  , '6b1.',/overplot)
  ;f2 = plot(x, sub+add, '1sr',/overplot , sym_size = 0.4, sym_filled = 1.)

  print, 'nsub', n_elements(sub)

;then put it all back together again
  if nseg eq 1 then begin
     full_xarr = xcenarr
     full_yarr = ycenarr
     full_flux = y
     full_fluxerr = yerr
     full_bmjd = bmjdarr
     full_utcs = utcsarr
     full_model = model_fit
     full_model_y = model_y
     full_sub = sub
     full_time = x
  endif else begin
     full_xarr = [full_xarr, xcenarr]
     full_yarr = [full_yarr, ycenarr]
     full_flux = [full_flux,y]
     full_fluxerr = [full_fluxerr, yerr]
     full_bmjd = [full_bmjd,bmjdarr]
     full_utcs = [full_utcs, utcsarr]
     full_model = [full_model,model_fit]
     full_model_y = [full_model_y, model_y]
     full_sub = [full_sub,sub]
     full_time = [full_time, x]

  endelse
print, 'nfull', n_elements(full_xarr)

endfor                          ;for each of the nseg segments
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
 ;binning after fitting
  if keyword_set(binning) then begin
     
     nfits = n_elements(full_xarr) / bin_level

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
     bin_model_y = dblarr(nfits)

     print, 'nfits', nfits
     help, nfits
     for si = 0L, long(nfits) - 1L do begin
      ;  print, 'working on si', si, n_elements(y), si*nframes, si*nframes + (nframes - 1)
        idata = full_flux[si*nframes:si*nframes + (nframes - 1)]
        idataerr = full_fluxerr[si*nframes:si*nframes + (nframes - 1)]
        bin_flux[si] = mean(idata,/nan)
        bin_fluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        ;bin_corrflux[si] = mean(corrflux[si*nframes:si*nframes + (nframes - 1)],/nan)
        ;bin_corrfluxerr[si] = sqrt(total((corrfluxerr[si*nframes:si*nframes + (nframes - 1)])^2))/ (n_elements(corrfluxerr[si*nframes:si*nframes + (nframes - 1)]))
        bin_time[si] = mean(full_time[si*nframes:si*nframes + (nframes - 1)])
        bin_xcen[si] = mean(full_xarr[si*nframes:si*nframes + (nframes - 1)])
        bin_ycen[si] = mean(full_yarr[si*nframes:si*nframes + (nframes - 1)])
        bin_bmjd[si] = mean(full_bmjd[si*nframes:si*nframes + (nframes - 1)])
        bin_utcs[si] = mean(full_utcs[si*nframes:si*nframes + (nframes - 1)])
        bin_sub[si] = mean(full_sub[si*nframes:si*nframes + (nframes - 1)])
        bin_model[si] =  mean(full_model[si*nframes:si*nframes + (nframes - 1)])
        bin_model_y[si] =  mean(full_model_y[si*nframes:si*nframes + (nframes - 1)])
     endfor                     ;for each fits image
     
  print, 'nbin', n_elements(bin_time)
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
     st2 = plot(bin_time, bin_model_y,'1sb',/overplot, sym_size = 0.2, sym_filled = 1.)
     print, 'bin_model_y', mean(bin_model_y)
  ;-----------------
  ;x vs time
     st1 = plot(bin_time, bin_xcen,'1o', xtitle = 'Time from transit (hours)',ytitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4, yrange =xyrange, title = titlename) ;
  ;-----------------
  ;x vs y
    ; st2 = plot(bin_xcen, bin_ycen,'1o', ytitle = 'Y pix',xtitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4, title = titlename, yrange = yrange, xrange = xrange) ;
   ;-----------------

     x =  phase_time             ;/ max(timearr)        ; fix to 1 phase
     y = bin_flux               ;/ fluxarr[1]          ;normalize
     yerr = bin_fluxerr         ;/ fluxarr[1]    ;normalize
     print, 'n x, y, yerr', n_elements(x), n_elements(y), n_elements(yerr)
     p1 = plot(x, y,  '1s',  xtitle = 'Orbital Phase', ytitle = ' Flux',  sym_size = 0.2, sym_filled = 1.,yrange = fyrange, title = titlename);, xrange = [-0.7,0.7]) ;, xrange = [2,4]
     f1 = plot(x,bin_model,  '1sb',/overplot , sym_size = 0.2, sym_filled = 1.) ;-0.0027
     f2 = plot(x, bin_sub+add, '1sr',/overplot , sym_size = 0.2, sym_filled = 1.)

     ;p1.save, "/Users/jkrick/irac_warm/hd209458/fluxvsphase_"+titlename +".png"
     ;st.save, "/Users/jkrick/irac_warm/hd209458/yvstime_"+titlename +".png"
     ;st1.save,  "/Users/jkrick/irac_warm/hd209458/xvstime_"+titlename +".png"
     ;st2.save,  "/Users/jkrick/irac_warm/hd209458/xvsy_"+titlename +".png"
     
  endif                         ;binning

;save, /all, filename = '/Users/jkrick/irac_warm/hd209458/hd209458_plot_'+titlename+'_bad.sav'

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
