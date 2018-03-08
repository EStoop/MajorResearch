pro hd7924_datafit, binning = binning
  
  restore,  '/Users/jkrick/irac_warm/snapshots/hd7924/hd7924_corr.sav'


  if keyword_set(binning) then begin
     sa = 0 
;try binning the fluxes by subarray image
     nfits =  17; 1700.;3400                
     nframes =  12600.;126.;63.
     bin_snaps = replicate({hd7924ob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' '},1)
     
      for si = 0, nfits - 1 do begin
         ;print, 'si', si, n_elements(snapshots[0].flux), si*nframes, si*nframes + (nframes - 1)
         idata = snapshots[sa].flux[si*nframes:si*nframes + (nframes - 1)]
         idataerr = snapshots[sa].fluxerr[si*nframes:si*nframes + (nframes - 1)]
         binned_flux = mean(idata,/nan)
         binned_fluxerr =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
         binned_corrflux = mean(snapshots[sa].corrflux[si*nframes:si*nframes + (nframes - 1)],/nan)
         binned_corrfluxerr = sqrt(total((snapshots[sa].corrfluxerr[si*nframes:si*nframes + (nframes - 1)])^2))/ (n_elements(snapshots[sa].corrfluxerr[si*nframes:si*nframes + (nframes - 1)]))
         binned_timearr = mean(snapshots[sa].timearr[si*nframes:si*nframes + (nframes - 1)])
         binned_xcen = mean(snapshots[sa].xcen[si*nframes:si*nframes + (nframes - 1)])
         binned_ycen = mean(snapshots[sa].ycen[si*nframes:si*nframes + (nframes - 1)])
         
         bin_snaps[sa].flux[si] = binned_flux
         bin_snaps[sa].fluxerr[si] = binned_fluxerr
         bin_snaps[sa].corrflux[si] = binned_corrflux
         bin_snaps[sa].corrfluxerr[si] = binned_corrfluxerr
         bin_snaps[sa].timearr[si] = binned_timearr
         bin_snaps[sa].xcen[si] = binned_xcen
         bin_snaps[sa].ycen[si] = binned_ycen
         
      endfor
      bin_snaps[sa].sclktime_0 = snapshots[sa].sclktime_0
   
endif   ;binning

if not keyword_set(binning) then begin
   bin_snaps = snapshots
endif



  xarr = bin_snaps[0].xcen
  yarr = bin_snaps[0].ycen
  timearr = bin_snaps[0].timearr
  fluxarr = bin_snaps[0].flux
  fluxerrarr = bin_snaps[0].fluxerr
  corrflux = bin_snaps[0].corrflux
  corrfluxerr = bin_snaps[0].corrfluxerr

;get rid of outliers in the pointing
;  good = where(xarr gt 15.0 and xarr lt 15.15 and yarr gt 14.85 and yarr lt 15.1 );
;  xarr = xarr[good]
;  yarr = yarr[good]
;  timearr = timearr[good]
;  fluxerrarr = fluxerrarr[good]
;  fluxarr = fluxarr[good]
;  corrflux = corrflux[good]
;  corrfluxerr = corrfluxerr[good]


   good_pmap = where(finite(bin_snaps[0].corrflux) gt 0 and finite(bin_snaps[0].corrfluxerr) gt 0,ngood_pmap, complement=bad_pmap)
   ;print, 'bad pmap', bin_snaps[a].corrfluxerr(bad_pmap)
  xarr = xarr[good_pmap]
  yarr = yarr[good_pmap]
  timearr = timearr[good_pmap]
  fluxerrarr = fluxerrarr[good_pmap]
  fluxarr = fluxarr[good_pmap]
  corrflux = corrflux[good_pmap]
  corrfluxerr = corrfluxerr[good_pmap]


  x = timearr ;/ max(timearr)        ; fix to 1 phase
  y = fluxarr ;/ fluxarr[1]          ;normalize
  yerr = fluxerrarr ;/ fluxarr[1]    ;normalize
  print, 'n x, y, yerr', n_elements(x), n_elements(y), n_elements(yerr)
;the plots!
  p1 = errorplot(x, y/median(y), yerr/median(y), '1-s',  xtitle = 'Time(hours)', ytitle = 'Normalized Flux', title = 'Data Fit', yrange = [0.998, 1.002] ,errorbar_capsize = .05 , sym_size = 0.4, sym_filled = 1.) ;, xrange = [2,4]
     
 ;--------------------------------------------------------------------------------------------------
;can I get rid of pixel phase effect?
; Initial guess	
  pa0 = [median(y), 1.,1.,1.,1.,1.,0.006,20.618,0.5]
;  pa0 = [median(y), 1.,1.,1.,1.,1.]
  amx = median(xarr)            ;xa[agptr])
  amy = median(yarr)            ; ya[agptr])	
  adx = xarr - amx              ; xa[agptr] - amx
  ady = yarr - amy              ; ya[agptr] - amy 
  
  parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.D,0.D]}, n_elements(pa0))
                                ;limit the range of the sin cuve phase
  ;parinfo[8].fixed = 1 
  parinfo[7].fixed = 1          
  ;parinfo[6].limited[0] = 1
  ;parinfo[6].limits[0] = 0.0
  

  afargs = {FLUX:y, DX:adx, DY:ady, T:x, ERR:yerr}
  pa = mpfit('fpa1_xfunc3', pa0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg, parinfo = parinfo) ;,/quiet)
  print, 'status', status
  print, errmsg
  achi = achi / adof
  print, 'reduced chi squared', achi
  
     ;t = text(0.2, 0.965, 'chi sq' + string(achi), /data, /current)
     
     
  ;look at correlations amongs parameters
  ;PCOR = COV * 0
  ;n = n_elements(spa)
  ;FOR i = 0, n-1 DO FOR j = 0, n-1 DO PCOR[i,j] = COV[i,j]/sqrt(COV[i,i]*COV[j,j])
  ;print, 'pcor', pcor
  
  ;;f1 = plot(x,   pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + $
  ;;                        pa[5] * ady * ady )  , '6b1.',/overplot)
  
  ;plot residuals
  ;sub = y - (pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + $
 ;                     pa[5] * ady * ady ) )
  sub = y - (pa[0] * (1+ pa[1] * adx + pa[2] * ady + pa[3] * adx * ady + pa[4] * adx * adx + $
                             pa[5] * ady * ady ) + pa[6]*sin(x/pa[7] + pa[8]))
  f2 = errorplot(x, sub+1.447, yerr, '1-sr',/overplot ,errorbar_capsize = .05 , sym_size = 0.4, sym_filled = 1.)
  
 ;compare with the light curve from using Jim's pixel phase map.
    
  delta_sub = bin_snaps[0].corrflux - sub 
  print,'delta_sub',  mean(delta_sub)
  fdelta = plot(timearr, delta_sub, '1-*', xtitle = 'Time(Hours)', ytitle = 'Delta Flux')

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


;fit residuals to know what the errors look like.
;     plothist, sub, xhist, yhist,bin=0.0005, /noprint,/noplot
;     r0 = barplot(xhist, yhist, fill_color = 'green', xrange = [-0.01, 0.01], xtitle = 'Light Curve', ytitle = 'Number')
    
;     noise = fltarr(n_elements(yhist)) + 1.0 ;  
 ;;    by = where(yhist lt 100.)
;     noise(by) = 100.
;     ;print, 'noise', noise
;     start = [0.0001, 0.01, 10000.]
;     result= MPFITFUN('mygauss',xhist,yhist, noise, start) 
;     r1 = plot( xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.),  '3',/overplot)
    
;     t_sig = text(0.02, 200, string(result(1)),'3',/data)
