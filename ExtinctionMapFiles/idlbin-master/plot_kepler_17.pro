pro plot_kepler_17, planetname, bin_level, apradius, chname, selfcal=selfcal, errorbars = errorbars, phaseplot = phaseplot, fit_eclipse = fit_eclipse

;COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycenp, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr

;get all the necessary saved info/photometry
  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  intended_phase = planetinfo[planetname, 'intended_phase']
  stareaor = planetinfo[planetname, 'stareaor']
  delta_red = -0.02
  delta_grey =  0.01
  delta_blue = 0.02  
  delta_cyan = 0.03
  delta_green =- 0.01
  planetname_final = planetname
  dirname = strcompress(basedir + planetname +'/')
  print, 'dirname', dirname
  startaor = 0
  stopaor = 1


  ;making things look better on the plot
  ;leave space for the selfcal plot, or not.
 ; if keyword_set(selfcal) then delta_green = delta_green else delta_red = delta_green

  rawnorm = 0.017717

  ;need to look for all endings of the pixphasecorr file saved
  filename =strcompress(dirname +'pixphasecorr_ch'+chname+'_'+string(apradius)+'_*.sav',/remove_all)
  savfiles = FILE_SEARCH(filename,count=nfiles)
  print, 'restoring',nfiles , 'files'
;  print, 'flux_np', flux_np
  for f = 0, nfiles -1 do begin
     restore, savfiles[f]

     ;which input version am I working with?  xy, np, fw
     ending = strmid(savfiles[f], 5, 2, /reverse_offset)
     if ending eq 'fw' then plotcolor='cyan'
     if ending eq 'np' then plotcolor = 'blue'
     if ending eq 'xy' then plotcolor = 'red'

     pmapcorr = 1
     numberarr = findgen(n_elements(flux_marr))
     h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
     print, 'omin', om, 'nh', n_elements(h)
        
        
;mean together the flux values in each phase bin
     bin_bmjd = dblarr(n_elements(h))
     bin_corrflux = dblarr(n_elements(h))
     bin_corrfluxerr = dblarr(n_elements(h))
     bin_flux_m = dblarr(n_elements(h))
     bin_fluxerr_m = dblarr(n_elements(h))
     bin_flux_np = dblarr(n_elements(h))
     bin_fluxerr_np = dblarr(n_elements(h))
     bin_flux = dblarr(n_elements(h))
     bin_fluxerr = dblarr(n_elements(h))
     bin_time = dblarr(n_elements(h))
     bin_time_0 = dblarr(n_elements(h))
     bin_phase = dblarr(n_elements(h))
     bin_flux_fwhm = dblarr(n_elements(h))
     bin_fluxerr_fwhm= dblarr(n_elements(h))
     c = 0
     for j = 0L, n_elements(h) - 1 do begin
           
;;get rid of the bins with no values and low numbers, meaning low overlap
        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
 
;commenting this out for the case that I have no corrfluxes          
;           if finite(corrfluxarr[ri[ri[j]]]) gt 0 and finite(corrfluxarr[ri[ri[j+1]-1]]) gt 0 then begin
;              meanclip, corrfluxarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;              bin_corrflux[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
;           endif else begin
;              bin_corrflux[c] = alog10(-1)
;           endelse
           
;           icorrdataerr = corrfluxerrarr[ri[ri[j]:ri[j+1]-1]]
;           bin_corrfluxerr[c] =   sqrt(total(icorrdataerr^2))/ (n_elements(icorrdataerr))
           
           meanclip, flux_marr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_flux_m[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           ifluxmerr = fluxerr_marr[ri[ri[j]:ri[j+1]-1]]
           bin_fluxerr_m[c] =   sqrt(total(ifluxmerr^2))/ (n_elements(ifluxmerr))
           
;              meanclip, flux_np[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;              bin_flux_np[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
;              ifluxnperr = fluxerr_np[ri[ri[j]:ri[j+1]-1]]
;              bin_fluxerr_np[c] =   sqrt(total(ifluxnperr^2))/ (n_elements(ifluxnperr))
;           print, flux[ri[ri[j]:ri[j+1]-1]]
           if total(finite(flux[ri[ri[j]:ri[j+1]-1]] )) lt 0.9 then begin ; entire bin is NAN's
              bin_flux[c] = !VALUES.F_NAN 
           endif else begin
              meanclip, flux[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
              bin_flux[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           endelse 
           
           idataerr = fluxerr[ri[ri[j]:ri[j+1]-1]]
           bin_fluxerr[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
           
;              meanclip, flux_fwhm[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;              bin_flux_fwhm[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
;              idataerr = fluxerr_fwhm[ri[ri[j]:ri[j+1]-1]]
;              bin_fluxerr_fwhm[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
           
           meanclip, timearr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_time[c] = meanx  ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip, time_0[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_time_0[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip, phasearr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_phase[c] = meanx ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           meanclip, bmjdarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_bmjd[c] = meanx  ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
           
           c = c + 1
        endif
     endfor
     
;     bin_corrflux =bin_corrflux[0:c-1]
;     bin_corrfluxerr = bin_corrfluxerr[0:c-1]
     bin_flux_m = bin_flux_m[0:c-1]
     bin_flux_np = bin_flux_np[0:c-1]
     bin_flux = bin_flux[0:c-1]
     bin_time = bin_time[0:c-1]
     bin_time_0 =  bin_time_0[0:c-1]
     bin_phase = bin_phase[0:c-1]

     offp = where(bin_phase lt 0)
     bin_phase[offp] = bin_phase[offp] + 1.0

     bin_bmjd = bin_bmjd[0:c-1]
     bin_fluxerr_m = bin_fluxerr_m[0:c-1]
     bin_fluxerr_np = bin_fluxerr_np[0:c-1]
     bin_fluxerr = bin_fluxerr[0:c-1]
     bin_flux_fwhm = bin_flux_fwhm[0:c-1]
     bin_fluxerr_fwhm = bin_fluxerr_fwhm[0:c-1]
     
;if working on a secondary, then want to plot the phase around 0.5,
;and not split it around 0.
;     if intended_phase eq 0.5 then begin
;        pa = where(bin_phase lt 0.0,pacount)
;        if pacount gt 0 then bin_phase[pa] = bin_phase[pa] + 1.0
;     endif
     
     
;test the levels
     print, 'mean raw black flux', mean(bin_flux_m / median(bin_flux_m),/nan)
;     print, 'mean corr gray flux', mean((bin_corrflux /median( bin_corrflux)),/nan)
     print, 'mean pixphasecorr flux', ending, mean(bin_flux,/nan)
;     print, 'mean np blue flux', mean(bin_flux_np,/nan)
;     print, 'mean fwhm cyan flux', mean(bin_flux_fwhm,/nan)
;     print, 'phase', bin_phase
     
;print out levels so that I can use them for TAP (or exofast I suppose)
                                ;    for exofast = 0, n_elements(bin_time) - 1 do begin
;        if finite(bin_flux(exofast)) gt 0 then begin
;           print, bin_bmjd(exofast) + 2400000.5D, ' ',bin_corrflux(exofast), ' ',bin_corrfluxerr(exofast), format = '(F0, A,F0, A, F0)'
;           print, bin_bmjd(exofast) - 56081.26 , ' ',(bin_flux_np(exofast) / median(bin_flux_np)) - 0.0005, format = '(F0, A,F0)'
                                ;a guess at mid-transit
;        endif
                                ;endfor
     outfilename =  strcompress(dirname +'phot_ch'+chname+'_TAP_nn.ascii',/remove_all)
     openw, outlun, outfilename,/GET_LUN
     endval = fix(0.3*n_elements(bin_flux))
     normfactor = median(bin_flux[0:endval])
     print, 'n_elements(bin_phase)',n_elements(bin_phase)
     for te = 0, n_elements(bin_phase) -1 do begin
        printf, outlun, bin_bmjd(te) , ' ', bin_flux(te)/normfactor,  format = '(F0, A,F0)'
     endfor
     close, outlun
;-----------------

;-----------------
;plot the results
     ;middle one is out of wack.
     bad = where(bin_phase lt 1E-3)
     bin_phase(bad) = 0.5
     print, 'bin_phase', bin_phase


     p1 = errorplot(bin_phase, bin_flux,bin_fluxerr,$ ;median(bin_flux_m)
                    '1s', sym_size = 0.5, sym_filled = 1,xrange = [0.40, 0.60],$
                    color = 'black', xtitle = 'Phase', ytitle = 'Normalized Corrected Flux',$ ; title = planetname, $
                    name = 'raw flux',$;  yrange =[0.97, 1.03], $
                    errorbar_capsize = 0.1) ; xstyle = 1,axis_style = 1,

  ;overplot a simple trapezoid based on the TAP output
     ; this just lets me control the plot more than TAP

     xp = findgen(100) / 100.
     yp = fltarr(100) + 1.0
     ax = where(xp gt 0.475 and xp lt 0.53)
     yp(ax) = 0.996
;     p2 = plot(xp, yp, color = 'blue', thick = 3, overplot = p1)
     


     ;read in the TAP model
     readcol, '/Users/jkrick/irac_warm/pcrs_planets/kepler-17/TAPmcmc_20141020_1309/ascii_phased_model.ascii', MJD, hours, model
     p3 = plot(bin_phase, model, color = 'cyan', thick = 2, overplot = p1)

  endfor
end

function fit_eclipse, xphase, ynorm, ynormerr, t1, dt, t3, d, delta, plotcolor
     print, 'inside fit_eclipse, working on ', plotcolor

     ;fit the eclipse with a trapezoid
     ;start with the red curve in phase plot only
     ;starting value educated guesses
;     t1 = 0.486; 1.                    ; hrs start of ingress
;     dt = 0.005;.5                    ; hrs duration of ingress (or egress)
;     t3 =0.515                 ; hrs start of egress
;     d = 0.001

  MESSAGE, /RESET
  ;remove indices with nans
  badnan = where(finite(ynormerr) lt 1 or finite(ynorm) lt 1 , badcount)
  if badcount gt 0 and badcount ne n_elements(ynorm) then begin
     remove, badnan, ynorm
     remove, badnan, ynormerr
     remove, badnan, xphase
  endif

     pa0 = [mean(ynorm,/nan), d,t1,t3,dt] ;f0, d, t1, t3, dt
     afargs = {FLUXfit:ynorm, Xfit:xphase, YERRfit:ynormerr}
     pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},n_elements(pa0))
     pi(4).limited(0) = 1
     pi(4).limits(0) = 1E-4
     pi(3).limited(0) = 1
     pi(3).limits(0) = 0.505
     pi(2).limited(0) = 1
     pi(2).limits(0) = xphase(2)
     pi(1).limited(0) = 1
     pi(2).limits(0) = 0.
     pa = mpfit('eclipse_trapezoid', pa0, FUNCTARGS=afargs, PERROR=spa, BESTNORM=achi, DOF=adof, COVAR = COV, status = status, errmsg = errmsg, parinfo = pi) ;, parinfo = parinfo) ;,/quiet)
     print, 'status', status
     print, 'errmsg', errmsg
     achi = achi / adof
     print, 'reduced chi squared', achi
     print, 'depth of eclipse', pa(1)
          
   ;plotting overlay
     flat = where(xphase le pa(2), flatcount)
     ;catch a bug: if there is no initial flat part
     if flatcount lt 1 then begin
        print, 'inside if statement'
        flatcount = 1
        flat = 1
     endif

;     print, 'xphase', xphase
;     print, 'flatcount', flatcount , pa(2)
     x1 = xphase[flat]
     y1 = fltarr(flatcount) + pa(0) + delta ;  pa(0)
    ; print, 'x,y', x1[15], y1[15]
 ;;    f1 = plot(x1,   y1  , '2',color = plotcolor, /overplot)
     
     ingress = where(xphase ge pa(2) and xphase le (pa(2) + pa(4)))
     x2 = fltarr(2)
     x2(0) = pa(2)
     x2(1) = pa(2) + pa(4)
     y2 = fltarr(2)
     y2(0) = pa(0) + delta
     y2(1) = pa(0) + delta - pa[1]
 ;;    f2 = plot(x2, y2, '2',color = plotcolor,/overplot)
     
; In eclipse
     eclipse = where(xphase ge (pa[2] + pa[4]) and xphase le pa[3], count)
     y3 = fltarr(count) + pa(0) + delta - pa[1] ; pa[0] - pa[1]
     x3 = xphase[eclipse]
;;     f3 = plot(x3, y3, '2',color = plotcolor,/overplot)
     
; Leaving eclipse
     egress = where(xphase ge pa[3] and xphase le (pa[3] + pa[4]), count)
;   y4 = y[egress] + pa[1] * ((x[egress] - pa[3]) / pa[4] - 1.0)
;   x4 = x[egress]
     x4 = fltarr(2)
     x4(0) = pa(3)
     x4(1) = pa(3) + pa(4)
     y4 = fltarr(2)
     y4(0) = pa(0) + delta - pa[1]
     y4(1) = pa(0) + delta
     
;;     f4 = plot(x4, y4, '2',color = plotcolor,/overplot)
     
     flat = where(xphase ge (pa[3] + pa[4]), flatcount)
     x5 = xphase[flat]
     y5 = fltarr(flatcount) + pa(0) + delta ; pa(0)
;;     f5 = plot(x5,   y5  , '2',color = plotcolor,/overplot)

return, 0

end



function eclipse_trapezoid, p, FLUXfit=fluxfit, Xfit=xfit, yERRfit=yerrfit
; p is parameters, f0, d, t1, t3, dt
; t is time
; Start of ingress is t1, end of ingress is t2, start of egress is t3, end of egress is t4
; assume ingress and egress have same durations, dt 
   n = n_elements(xfit)
   model = fltarr(n) + p[0]

   t2 = p[2] + p[4]
   t4 = p[3] + p[4]

; Out of eclipse
    ptr = where(xfit le p[2], count)
    if (count gt 0) then model[ptr] = p[0]

    ptr = where(xfit gt t4, count)
    if (count gt 0) then model[ptr] = p[0]

; Beginning eclipse
   ptr = where(xfit gt p[2] and xfit le t2, count)
   if (count gt 0) then $
       model[ptr] = model[ptr] - p[1] * (1.0 - (t2 - xfit[ptr]) / p[4])

; In eclipse
   ptr = where(xfit gt t2 and xfit le p[3], count)
   if (count gt 0) then model[ptr] = model[ptr] - p[1]

; Leaving eclipse
   ptr = where(xfit gt p[3] and xfit le t4, count)
   if (count gt 0) then $
       model[ptr] = model[ptr] + p[1] * ((xfit[ptr] - p[3]) / p[4] - 1.0)

 model = (fluxfit - model) / yerrfit

return, model
end
