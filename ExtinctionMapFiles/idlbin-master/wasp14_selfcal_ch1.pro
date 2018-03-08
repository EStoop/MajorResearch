pro wasp14_selfcal_ch1, binning=binning, sine=sine

  aorname = ['r45426688',  'r45428224',  'r45428480',  'r45428736',  'r45428992'] ;ch1

;  AORwasp14 = replicate({wasp14ob, ra:0D,dec:0D,xcen:dblarr(nfits),ycen:dblarr(nfits),flux:dblarr(nfits),fluxerr:dblarr(nfits), corrflux:dblarr(nfits), corrfluxerr:dblarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', bmjdarr:dblarr(nfits), utcsarr:dblarr(nfits), bkgd:dblarr(nfits), bkgderr:dblarr(nfits)},n_elements(aorname))

  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp14/wasp14_phot_ch1.sav'
  period = 2.519961 ;days
  timestart = 0

;set phase = 0.0 in the middle of the transit
 ; for a = 0, n_elements(aorname) - 1 do  begin
 ;     AORwasp14[a].bmjdarr = ((AORwasp14[a].bmjdarr - utmjd_center) / period) 
 ; endfor
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
  xarr = AORwasp14.xcen
  yarr = AORwasp14.ycen
  timearr = AORwasp14.timearr
  fluxerr = AORwasp14.fluxerr
  flux = AORwasp14.flux
  corrflux = AORwasp14.corrflux
  corrfluxerr = AORwasp14.corrfluxerr
  bmjd = AORwasp14.bmjdarr
  utcs = AORwasp14.utcsarr

;get rid of outliers
 print, ' fluxarr', n_elements(flux)
 
 ; the position outliers
;  good = where(AORwasp14[a].xcen lt mean(AORwasp14[a].xcen)
;  +3.0*stddev(AORwasp14[a].xcen) and AORwasp14[a].xcen gt
;  mean(AORwasp14[a].xcen) -3.0*stddev(AORwasp14[a].xcen) and
;  AORwasp14[a].ycen lt mean(AORwasp14[a].ycen)
;  +3.0*stddev(AORwasp14[a].ycen) and AORwasp14[a].ycen gt
;  mean(AORwasp14[a].ycen) - 3.0*stddev(AORwasp14[a].ycen),ngood_pmap,
;                            complement=bad) 

;q = where(xarr gt 15.4 and xarr lt 16.2, nq)
;print,' off x', xarr(q)
;print, 'y', yarr(q)
;print, 'time' ,  (timearr(q) - timearr(0))/60./60.
;print, 'flux', flux(q)

 print, 'positions', mean(xarr) + 1.0*stddev(xarr),  mean(xarr) -1.0*stddev(xarr),  mean(yarr) +1.0*stddev(yarr), mean(yarr) - 1.0*stddev(yarr)

  good = where(xarr lt mean(xarr) + 5.*stddev(xarr)  and xarr gt mean(xarr) -5.0*stddev(xarr) and yarr lt mean(yarr) +5.0*stddev(yarr) and yarr gt mean(yarr) - 5.0*stddev(yarr),ngood_pmap, complement=bad) 
  print, 'bad position',n_elements(bad), n_elements(good)

;  xarr = AORwasp14[a].xcen[good]
;  yarr = AORwasp14[a].ycen[good]
;  timearr = AORwasp14[a].timearr[good]
;  fluxerr = AORwasp14[a].fluxerr[good]
;  flux = AORwasp14[a].flux[good]
;  corrflux = AORwasp14[a].corrflux[good]
;  corrfluxerr = AORwasp14[a].corrfluxerr[good]
;  bmjd = AORwasp14[a].bmjdarr[good]
;  utcs = AORwasp14[a].utcsarr[good]
 
  xarr = xarr[good]
  yarr = yarr[good]
  timearr = timearr[good]
  fluxerr = fluxerr[good]
  flux = flux[good]
  corrflux =corrflux[good]
  corrfluxerr = corrfluxerr[good]
  bmjd = bmjd[good]
  utcs = utcs[good]
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
  print, 'timearr', min(timearr), max(timearr)

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------
     ;add the transit as a break in the light curve
  starttransit = max(where(timearr lt 3.0)) ;3.0
  stoptransit = min(where(timearr gt 7.0)) ;6.5
  xarr = [xarr[0:starttransit],xarr[stoptransit:n_elements(timearr)-1]]
  yarr = [yarr[0:starttransit],yarr[stoptransit:n_elements(timearr)-1]]
  fluxerr = [fluxerr[0:starttransit],fluxerr[stoptransit:n_elements(timearr)-1]]
  flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
  corrflux = [corrflux[0:starttransit],corrflux[stoptransit:n_elements(timearr)-1]]
  corrfluxerr = [corrfluxerr[0:starttransit],corrfluxerr[stoptransit:n_elements(timearr)-1]]
  bmjd = [bmjd[0:starttransit],bmjd[stoptransit:n_elements(timearr)-1]]
  utcs = [utcs[0:starttransit],utcs[stoptransit:n_elements(timearr)-1]]
  timearr = [timearr[0:starttransit],timearr[stoptransit:n_elements(timearr)-1]]
  
  starttransit = max(where(timearr lt 31.5)) ;3.0
  stoptransit = min(where(timearr gt 35.)) ;6.5
  xarr = [xarr[0:starttransit],xarr[stoptransit:n_elements(timearr)-1]]
  yarr = [yarr[0:starttransit],yarr[stoptransit:n_elements(timearr)-1]]
  fluxerr = [fluxerr[0:starttransit],fluxerr[stoptransit:n_elements(timearr)-1]]
  flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
  corrflux = [corrflux[0:starttransit],corrflux[stoptransit:n_elements(timearr)-1]]
  corrfluxerr = [corrfluxerr[0:starttransit],corrfluxerr[stoptransit:n_elements(timearr)-1]]
  bmjd = [bmjd[0:starttransit],bmjd[stoptransit:n_elements(timearr)-1]]
  utcs = [utcs[0:starttransit],utcs[stoptransit:n_elements(timearr)-1]]
  timearr = [timearr[0:starttransit],timearr[stoptransit:n_elements(timearr)-1]]

 starttransit = max(where(timearr lt 12.5)) ;3.0
  stoptransit = min(where(timearr gt 13.2)) ;6.5
  xarr = [xarr[0:starttransit],xarr[stoptransit:n_elements(timearr)-1]]
  yarr = [yarr[0:starttransit],yarr[stoptransit:n_elements(timearr)-1]]
  fluxerr = [fluxerr[0:starttransit],fluxerr[stoptransit:n_elements(timearr)-1]]
  flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
  corrflux = [corrflux[0:starttransit],corrflux[stoptransit:n_elements(timearr)-1]]
  corrfluxerr = [corrfluxerr[0:starttransit],corrfluxerr[stoptransit:n_elements(timearr)-1]]
  bmjd = [bmjd[0:starttransit],bmjd[stoptransit:n_elements(timearr)-1]]
  utcs = [utcs[0:starttransit],utcs[stoptransit:n_elements(timearr)-1]]
  timearr = [timearr[0:starttransit],timearr[stoptransit:n_elements(timearr)-1]]

;really the jumps in pcrspeakup observing to move back to sweet spot
 starttransit = max(where(timearr lt 39.0)) ;3.0
  stoptransit = min(where(timearr gt 39.2)) ;6.5
  xarr = [xarr[0:starttransit],xarr[stoptransit:n_elements(timearr)-1]]
  yarr = [yarr[0:starttransit],yarr[stoptransit:n_elements(timearr)-1]]
  fluxerr = [fluxerr[0:starttransit],fluxerr[stoptransit:n_elements(timearr)-1]]
  flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
  corrflux = [corrflux[0:starttransit],corrflux[stoptransit:n_elements(timearr)-1]]
  corrfluxerr = [corrfluxerr[0:starttransit],corrfluxerr[stoptransit:n_elements(timearr)-1]]
  bmjd = [bmjd[0:starttransit],bmjd[stoptransit:n_elements(timearr)-1]]
  utcs = [utcs[0:starttransit],utcs[stoptransit:n_elements(timearr)-1]]
  timearr = [timearr[0:starttransit],timearr[stoptransit:n_elements(timearr)-1]]

 starttransit = max(where(timearr lt 52.0)) ;3.0
  stoptransit = min(where(timearr gt 53.)) ;6.5
  xarr = [xarr[0:starttransit],xarr[stoptransit:n_elements(timearr)-1]]
  yarr = [yarr[0:starttransit],yarr[stoptransit:n_elements(timearr)-1]]
  fluxerr = [fluxerr[0:starttransit],fluxerr[stoptransit:n_elements(timearr)-1]]
  flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
  corrflux = [corrflux[0:starttransit],corrflux[stoptransit:n_elements(timearr)-1]]
  corrfluxerr = [corrfluxerr[0:starttransit],corrfluxerr[stoptransit:n_elements(timearr)-1]]
  bmjd = [bmjd[0:starttransit],bmjd[stoptransit:n_elements(timearr)-1]]
  utcs = [utcs[0:starttransit],utcs[stoptransit:n_elements(timearr)-1]]
  timearr = [timearr[0:starttransit],timearr[stoptransit:n_elements(timearr)-1]]

 starttransit = max(where(timearr lt 57.0)) ;3.0
  stoptransit = min(where(timearr gt 62.)) ;6.5
  xarr = [xarr[0:starttransit],xarr[stoptransit:n_elements(timearr)-1]]
  yarr = [yarr[0:starttransit],yarr[stoptransit:n_elements(timearr)-1]]
  fluxerr = [fluxerr[0:starttransit],fluxerr[stoptransit:n_elements(timearr)-1]]
  flux = [flux[0:starttransit],flux[stoptransit:n_elements(timearr)-1]]
  corrflux = [corrflux[0:starttransit],corrflux[stoptransit:n_elements(timearr)-1]]
  corrfluxerr = [corrfluxerr[0:starttransit],corrfluxerr[stoptransit:n_elements(timearr)-1]]
  bmjd = [bmjd[0:starttransit],bmjd[stoptransit:n_elements(timearr)-1]]
  utcs = [utcs[0:starttransit],utcs[stoptransit:n_elements(timearr)-1]]
  timearr = [timearr[0:starttransit],timearr[stoptransit:n_elements(timearr)-1]]

;can I find the breaks in time?
  count = 0
  print, 'starting breaks', timearr(10) - timearr(9)
  breakpt = fltarr(9)
  for ti = 0, n_elements(timearr) - 5 do begin
     if timearr(ti+1) - timearr(ti) gt 0.1 then begin
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

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

;now work on the data fitting with a polynomial for comparison:
;can I get rid of pixel phase effect?

;     x = timearr
;     y =flux
;     yerr = fluxerr
;     xcenarr = xarr
;     ycenarr = yarr
;     amx = median(xcenarr)         
;     amy = median(ycenarr)        
;     adx = xcenarr - amx           
;     ady = ycenarr - amy    
;     ;nfits = n_elements(x) / nframes
;     bmjdarr= bmjd
;     utcsarr=utcs
 ;    nseg = 1

; Initial guesses	
  if keyword_set(sine) then begin  ;use a fitting function with a sine curve or not
     known_period = period*24.; hours 
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
endfor  ;for nseg

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

  bin_level = 63L
  nframes = bin_level
  nfits = long(n_elements(flux)) / nframes
  print, 'nframes, nfits', nframes, nfits
  help, nfits

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
  tt = 0  ;XXX
  bin_time = bin_time -tt

;set time in orbital phase from 0 at transit to 1 at next transit
  phase_time = bin_time / (period*24.)

  ;y vs. time
 ;    st = plot(bin_time, bin_ycen,'1o',xtitle = 'Time from transit (hours)',ytitle = 'Y pix',   color = 'black', sym_filled = 1, sym_size = 0.4, yrange = yyrange, title = titlename) ;, axis_style = 1) 
;tjunk = plot(bin_time,bin_ycen,/current, xrange = st.xrange, yrange = yyrange, nodata = 1)
;topx = axis('X', location = [0, 14.8],  target = tjunk, title = 'Time from transit (hours)')
;     st2 = plot(bin_time, bin_model_y+14.5,'1sb',/overplot, sym_size = 0.2, sym_filled = 1.)
     print, 'bin_model_y', mean(bin_model_y)
  ;-----------------
  ;x vs time
 ;    st1 = plot(bin_time, bin_xcen,'1o', xtitle = 'Time from transit (hours)',ytitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4, yrange =xyrange, title = titlename) ;
  ;-----------------
  ;x vs y
    ; st2 = plot(bin_xcen, bin_ycen,'1o', ytitle = 'Y pix',xtitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4, title = titlename, yrange = yrange, xrange = xrange) ;
   ;-----------------

     x =  bin_time             ;/ max(timearr)        ; fix to 1 phase
     y = bin_flux               ;/ fluxarr[1]          ;normalize
     yerr = bin_fluxerr         ;/ fluxarr[1]    ;normalize
     print, 'n x, y, yerr', n_elements(x), n_elements(y), n_elements(yerr)
     p1 = plot(x, y,  '1s',  xtitle = 'Time(hrs)', ytitle = ' Flux',  sym_size = 0.2, sym_filled = 1., title = 'Wasp 14 ch1 selfcal', xrange = [0,65], yrange = [0.093, 0.1]) ;, xrange = [2,4]
     f1 = plot(x,bin_model,  '1sb',/overplot , sym_size = 0.2, sym_filled = 1.) ;-0.0027
     f2 = plot(x, bin_sub+0.094, '1sr',/overplot , sym_size = 0.2, sym_filled = 1.)

     
     
  endif                         ;binning
  selfcal_timearr = full_time
  selfcal_flux = full_flux
  save, selfcal_timearr, filename='/Users/jkrick/irac_warm/pcrs_planets/wasp14/selfcal_timearr_ch1.sav'
  save, selfcal_flux, filename='/Users/jkrick/irac_warm/pcrs_planets/wasp14/selfcal_flux_ch1.sav'
  save, bin_time, filename='/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_time_ch1.sav'
  save, bin_sub, filename='/Users/jkrick/irac_warm/pcrs_planets/wasp14/bin_sub_ch1.sav'

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

