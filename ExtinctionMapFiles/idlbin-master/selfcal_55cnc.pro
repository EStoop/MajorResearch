pro selfcal_55cnc, binning=binning, sine=sine
 
;  AOR55cnc = replicate({ob55cnc, ra:0D,dec:0D,xcen:dblarr(nfits),ycen:dblarr(nfits),flux:dblarr(nfits),fluxerr:dblarr(nfits), corrflux:dblarr(nfits), corrfluxerr:dblarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', bmjdarr:dblarr(nfits), utcsarr:dblarr(nfits), bkgd:dblarr(nfits), bkgderr:dblarr(nfits)},n_elements(aorname)) ;

 ; restore, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/55cnc_phot_ch2_39524608.sav'
 ; aorname_55cnc  = ['r39524608'] ;ch2  39524608 transit
 
 restore, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/55cnc_phot_ch2_3_3_7.sav'
  aorname_55cnc = ['r43981312','r43981568','r43981824','r43981056'] ;ch2 eclipses

  utmjd_center = [double(55944.25889), double(55947.20505),double(55949.41467), double(55957.51661)]  

  duration = 105.7              ;in minutes
  duration = duration /60./24.  ; in days
  period = 0.73654 ;days
  

;set time = 0.5 in the middle of the eclipse
;     bin_time = bin_time -tt
  for a = 0, n_elements(aorname_55cnc) - 1 do  begin
     test =  ((AOR55cnc[a].bmjdarr - utmjd_center[a]) / period) + 0.5
     
     ;if a eq 1 then begin
        ;for t = 0, 100*63,63  do begin
        ;   print, t/63, ' ', AOR55cnc[1].bmjdarr[t], ' ',test[t], ' ', utmjd_center[a], format = '(F0,A, F0, A, F0, A, F0)'
       ; endfor
     ;endif

     ;AOR55cnc[a].bmjdarr = ((AOR55cnc[a].bmjdarr - utmjd_center[a]) / period) + 0.5
     bmjd =  ((AOR55cnc[a].bmjdarr - utmjd_center[a]) / period) + 0.5
  endfor

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

for a = 0,  n_elements(aorname_55cnc) - 1 do begin

   print, 'time test', AOR55cnc[a].timearr[0:100] - 1.0113223E9

;   startnum = [57600, 0, 56000, 56000]
 startnum = [0, 0, 0, 0]

; for phase set time = 0.5 in the middle of the eclipse
   bmjd =  ((AOR55cnc[a].bmjdarr[startnum(a):*] - utmjd_center[a]) / period) + 0.5

   print, ' fluxarr first', n_elements(AOR55cnc[a].flux)
  xarr = AOR55cnc[a].xcen[startnum(a):*]
  yarr = AOR55cnc[a].ycen[startnum(a):*]
  timearr = AOR55cnc[a].timearr[startnum(a):*]
  fluxerr = AOR55cnc[a].fluxerr[startnum(a):*]
  flux = AOR55cnc[a].flux[startnum(a):*]
  corrflux = AOR55cnc[a].corrflux[startnum(a):*]
  corrfluxerr = AOR55cnc[a].corrfluxerr[startnum(a):*]
  ;bmjd = AOR55cnc[a].bmjdarr[startnum(a):*]
  utcs = AOR55cnc[a].utcsarr[startnum(a):*]

;get rid of outliers
 print, ' fluxarr', n_elements(flux)
 print, 'positions', mean(xarr) + 3.0*stddev(xarr),  mean(xarr) -3.0*stddev(xarr),  mean(yarr) +3.0*stddev(yarr), mean(yarr) - 3.0*stddev(yarr)

  good = where(xarr lt mean(xarr) + 3.0*stddev(xarr) and xarr gt mean(xarr) -3.0*stddev(xarr) and yarr lt mean(yarr) +3.0*stddev(yarr) and yarr gt mean(yarr) - 3.0*stddev(yarr),ngood_pmap, complement=bad) 
  print, 'bad position',n_elements(bad), n_elements(good)
 
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
;now work on the data fitting with a polynomial for comparison:
;can I get rid of pixel phase effect?

    x = timearr
    y =flux
    yerr = fluxerr
    xcenarr = xarr
    ycenarr = yarr
    amx = median(xcenarr)         
    amy = median(ycenarr)        
    adx = xcenarr - amx           
    ady = ycenarr - amy    
    ;nfits = n_elements(x) / nframes
    bmjdarr= bmjd
    utcsarr=utcs
    nseg = 1

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
  sub = y - model_fit + pa[0]  ; add back in the overall level in Jy

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
;endfor  ;for nseg

  f = strcompress('/Users/jkrick/irac_warm/pcrs_planets/55cnc/selfcal_bmjd_' + string(a) + '.sav',/remove_all)
  save, full_bmjd, filename = f
  f = strcompress('/Users/jkrick/irac_warm/pcrs_planets/55cnc/selfcal_sub_' + string(a) + '.sav',/remove_all)
  save, full_sub, filename = f
  f = strcompress('/Users/jkrick/irac_warm/pcrs_planets/55cnc/selfcal_fluxerr_' + string(a) + '.sav',/remove_all)
  save, full_fluxerr, filename = f

;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------

  bin_level = 100.*63L
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
;     st = plot(bin_time, bin_ycen,'1o',xtitle = 'Time from transit (hours)',ytitle = 'Y pix',   color = 'black', sym_filled = 1, sym_size = 0.4, yrange = yyrange, title = titlename) ;, axis_style = 1) 
;tjunk = plot(bin_time,bin_ycen,/current, xrange = st.xrange, yrange = yyrange, nodata = 1)
;topx = axis('X', location = [0, 14.8],  target = tjunk, title = 'Time from transit (hours)')
;     st2 = plot(bin_time, bin_model_y+14.5,'1sb',/overplot, sym_size = 0.2, sym_filled = 1.)
     print, 'bin_model_y', mean(bin_model_y)
  ;-----------------
  ;x vs time
  ;   st1 = plot(bin_time, bin_xcen,'1o', xtitle = 'Time from transit (hours)',ytitle = 'X pix',  color = 'black', sym_filled = 1, sym_size = 0.4) ;
  ;-----------------
  ;x vs y
 ;    st2 = plot(bin_time, bin_ycen,'1o', ytitle = 'Y pix',xtitle = 'Time',  color = 'black', sym_filled = 1, sym_size = 0.4) ;
   ;-----------------

     x =  bin_bmjd             ;/ max(timearr)        ; fix to 1 phase
     y = bin_flux             ;/ fluxarr[1]          ;normalize
     print, 'n x, y, yerr', n_elements(x), n_elements(y), n_elements(yerr)
;     p1 = plot(x, y,  '1o',  xtitle = 'Orbital Phase', ytitle = ' Flux',  sym_size = 0.2, sym_filled = 1.,title = titlename, xrange = [0.32,0.62])
;     f1 = plot(x,bin_model,  '1ob',/overplot , sym_size = 0.2, sym_filled = 1.) ;-0.0027
;     f2 = plot(x, bin_sub, '1or',/overplot , sym_size = 0.2, sym_filled = 1.)

     y = (bin_sub )/mean(bin_sub)
     yerr = bin_fluxerr/ mean(bin_sub)
     o = errorplot(x, y,yerr,'1o',sym_size = 0.7, sym_filled = 1, xtitle = 'Orbital Phase', ytitle = 'Normalized corrected binned flux', xrange = [0.30, 0.70])
     
  endif                         ;binning
  selfcal_timearr = full_time
  selfcal_flux = full_flux
  
  save, /all, filename='/Users/jkrick/irac_warm/pcrs_planets/55cnc/selfcal.sav'
;  save, selfcal_timearr, filename='/Users/jkrick/irac_warm/pcrs_planets/55cnc/selfcal_timearr.sav'
;  save, selfcal_flux, filename='/Users/jkrick/irac_warm/pcrs_planets/55cnc/selfcal_flux.sav'
;  save, bin_time, filename='/Users/jkrick/irac_warm/pcrs_planets/55cnc/bin_time.sav'
;  save, bin_sub, filename='/Users/jkrick/irac_warm/pcrs_planets/55cnc/bin_sub.sav'
endfor

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
