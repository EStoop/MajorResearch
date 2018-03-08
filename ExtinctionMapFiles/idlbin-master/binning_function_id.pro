function binning_function_id, a,bin_level, pmapcorr, ch, set_nbins = set_nbins, n_nbins = n_nbins

  common bin_block
  ;-------------------------------------------------------------------------
  ;;setup

;;  print, 'starting binning', ch
;;keys =['ra', 'dec', 'pid', 'campaign', 'ch','xcen', 'ycen', 'timearr', 'aor', 'bmjdarr','exptime']

  timearr = [planethash[aorname(a),'timearr']]
;;  fluxarr = [planethash[aorname(a),'flux']]
;;  fluxerrarr = [planethash[aorname(a),'fluxerr']]
;;  corrfluxarr=[planethash[aorname(a),'corrflux']]
;;  corrflux_darr = [planethash[aorname(a),'corrflux_d']]
;;  corrfluxerrarr = [planethash[aorname(a),'corrfluxerr']]
  xarr = [planethash[aorname(a),'xcen']]
  yarr = [planethash[aorname(a),'ycen']]
;;  bkgd = [planethash[aorname(a),'bkgd']]
  bmjd = [planethash[aorname(a),'bmjdarr']]
;;  np = [planethash[aorname(a),'np']]
;;  npcentroids = planethash[aorname(a),'npcentroids']
;;  phase = [planethash[aorname(a),'phase']]
  ;;fix a small bug
;;  bp = where(phase gt 1.0, bpcount)
;;  if bpcount gt 0 then phase(bp) = phase(bp) - 1.
;;  xfwhmarr = [planethash[aorname(a),'xfwhm']]
;;  yfwhmarr = [planethash[aorname(a),'yfwhm']]

  ;-------------------------------------------------------------------------
  ;;remove outliers 
  ;;are there non-nan corrfluxes aka sweet spot
  goodpmap = where(xarr lt mean(xarr,/nan) + 2.5*stddev(xarr,/nan) and xarr gt mean(xarr,/nan) -2.5*stddev(xarr,/nan) and xarr lt mean(xarr,/nan) +3.0*stddev(yarr,/nan) and yarr gt mean(yarr,/nan) - 3.0*stddev(yarr,/nan)   ,ngood_pmap, complement=badpmap) ;and finite(corrfluxarr) gt 0
  ;;just care if it made the pixel
  good = where(xarr lt mean(xarr,/nan) + 2.5*stddev(xarr,/nan) and xarr gt mean(xarr,/nan) -2.5*stddev(xarr,/nan) and xarr lt mean(xarr,/nan) +3.0*stddev(yarr,/nan) and yarr gt mean(yarr,/nan) - 3.0*stddev(yarr,/nan) ,ngood, complement=bad)
     

  print, 'bad, good, percentage ',n_elements(bad), n_elements(good), n_elements(bad)/(float(n_elements(bad)) + n_elements(good))
  print, 'badp ',n_elements(badpmap), n_elements(goodpmap)
  xarr = xarr[good]
  yarr = yarr[good]
  timearr = timearr[good]
;;  flux = fluxarr[good]
;;  fluxerr = fluxerrarr[good]
;;  corrflux = corrfluxarr[good]
;;  corrflux_d = corrflux_darr[good]
;;  corrfluxerr = corrfluxerrarr[good]
  bmjdarr = bmjd[good]
;;  bkgdarr = bkgd[good]
;;  phasearr = phase[good]
;;  nparr = np[good]
;;  npcentarr = npcentroids[good]
 ;; xfwhm = xfwhmarr[good]
 ;; yfwhm = yfwhmarr[good]
  
  ;;and a second set for those that are in the sweet spot
  xarrp = xarr[goodpmap]
  yarrp = yarr[goodpmap]
  timearrp = timearr[goodpmap]
;;  fluxp = fluxarr[goodpmap]
;;  fluxerrp= fluxerrarr[goodpmap]
;;  corrfluxp = corrfluxarr[goodpmap]
;;  corrflux_dp = corrflux_darr[goodpmap]
;;  corrfluxerrp = corrfluxerrarr[goodpmap]
  bmjdarrp = bmjd[goodpmap]
;;  bkgdarrp = bkgd[goodpmap]
;;  phasearrp = phase[goodpmap]
;;  nparrp = np[goodpmap]
;;  npcentarrp = npcentroids[goodpmap]
;;  xfwhmp = xfwhmarr[goodpmap]
;;  yfwhmp = yfwhmarr[goodpmap]
  
  
  ;-------------------------------------------------------------------------
  ;; binning
  numberarr = findgen(n_elements(xarr))
  if keyword_set(set_nbins) then begin
     h = histogram(numberarr, OMIN=om, nbins = n_nbins, reverse_indices = ri)
  endif else begin
     h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
  endelse
;  print, 'omin', om, 'nh', n_elements(h)
  
  numberarrp = findgen(n_elements(xarrp))
  if keyword_set(set_nbins) then begin
     hp = histogram(numberarrp, OMIN=omp,nbins =n_nbins, reverse_indices = rip)
  endif else begin
     hp = histogram(numberarrp, OMIN=omp, binsize = bin_level, reverse_indices = rip)
  endelse
;  print, 'ominp', omp, 'nhp', n_elements(hp), n_elements(xarrp)
  ;;--------
  ;;setup arrays
  ;;mean together the flux values in each phase bin
  bin_flux = dblarr(n_elements(h))
  bin_fluxerr = bin_flux
  bin_corrflux= bin_flux
  bin_corrflux_d = bin_flux
  bin_corrfluxerr= bin_flux
  bin_ncorr = bin_flux
  bin_timearr = bin_flux
  bin_bmjdarr = bin_flux
  bin_bkgd = bin_flux
  bin_bkgderr = bin_flux
  bin_xcen = bin_flux
  bin_ycen = bin_flux
  bin_phase = bin_flux
  bin_centerpix = bin_flux
  bin_np = bin_flux
  bin_npcent = bin_flux
  bin_xfwhm = bin_flux
  bin_yfwhm = bin_flux
  
  bin_fluxp = dblarr(n_elements(hp))
  bin_fluxerrp = bin_fluxp
  bin_corrfluxp= bin_fluxp
  bin_corrflux_dp = bin_fluxp
  bin_corrfluxerrp= bin_fluxp
  bin_ncorrp = bin_fluxp
  bin_timearrp = bin_fluxp
  bin_bmjdarrp = bin_fluxp
  bin_bkgdp = bin_fluxp
  bin_bkgderrp = bin_fluxp
  bin_xcenp = bin_fluxp
  bin_ycenp= bin_fluxp
  bin_phasep = bin_fluxp
  bin_centerpixp = bin_fluxp
  bin_nparrp = bin_fluxp
  bin_npcentarrp = bin_fluxp
  bin_xfwhmp = bin_fluxp
  bin_yfwhmp = bin_fluxp
  
  ;; now for the actual binning 
  c = 0
  for j = 0L, n_elements(h) - 1 do begin
     
     ;;get rid of the bins with no values and low numbers, meaning
     ;;low overlap
     if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
        
        meanclip, xarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
        bin_xcen[c] = meanx     ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
        meanclip, yarr[ri[ri[j]:ri[j+1]-1]], meany, sigmay
        bin_ycen[c] = meany     ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
;;        meanclip, xfwhm[ri[ri[j]:ri[j+1]-1]], meanxfwhm, sigmaxfwhm
;;        bin_xfwhm[c] = meanxfwhm ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
;;        meanclip, yfwhm[ri[ri[j]:ri[j+1]-1]], meanyfwhm, sigmayfwhm
;;        bin_yfwhm[c] = meanyfwhm ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
;;        meanclip, bkgdarr[ri[ri[j]:ri[j+1]-1]], meansky, sigmasky
;;        bin_bkgd[c] = meansky   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
;;        meanclip, flux[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
;;        bin_flux[c] = meanflux  ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
;;        idataerr = fluxerr[ri[ri[j]:ri[j+1]-1]]
;;        bin_fluxerr[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        ;;meanclip, fluxerr[ri[ri[j]:ri[j+1]-1]], meanfluxerr, sigmafluxerr
        ;;bin_fluxerr[c] = sigmafluxerr
        
;;        meanclip, nparr[ri[ri[j]:ri[j+1]-1]], meannp, sigmanp
;;        bin_np[c] = meannp      ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
;;        meanclip, npcentarr[ri[ri[j]:ri[j+1]-1]], meannpcent, sigmanpcent
;;        bin_npcent[c] =  meannpcent ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        
;;        junk = where(finite(corrflux[ri[ri[j]:ri[j+1]-1]]) gt 0,ngood)
;;        bin_ncorr[c] = ngood
        
        ;;meanclip, timearr[ri[ri[j]:ri[j+1]-1]], meantimearr, sigmatimearr
        meantimearr = median(timearr[ri[ri[j]:ri[j+1]-1]])
        bin_timearr[c]=meantimearr
        
        ;; meanclip, phasearr[ri[ri[j]:ri[j+1]-1]], meanphasearr, sigmaphasearr
;;        meanphasearr = mean( phasearr[ri[ri[j]:ri[j+1]-1]],/nan)
        ;;having trouble with the boundary between 0.5 and -0.5
        ;;if a eq 0 then print, 'diff', meanphasearr - bin_phase[c-1]
        ;;if meanphasearr - bin_phase[c-1] lt -0.1 and meanphasearr - bin_phase[c-1] gt -0.5 then meanphasearr = mean(abs(phasearr[ri[ri[j]:ri[j+1]-1]]),/nan)
;;        bin_phase[c]= meanphasearr
        
;           if a eq 0 then begin
;              print, 'phase', phasearr[ri[ri[j]:ri[j+1]-1]]
;              print, 'meanphase', meanphasearr, mean(abs(phasearr[ri[ri[j]:ri[j+1]-1]]),/nan), abs(meanphasearr - bin_phase[c-1])
;           endif
        
        meanbmjdarr = mean( bmjdarr[ri[ri[j]:ri[j+1]-1]],/nan)
        bin_bmjdarr[c]= meanbmjdarr
        
        
        c = c + 1
     endif
  endfor
  
  bin_xcen = bin_xcen[0:c-1]
  bin_ycen = bin_ycen[0:c-1]
;;  bin_bkgd = bin_bkgd[0:c-1]
;;  bin_flux = bin_flux[0:c-1]
;;  bin_fluxerr = bin_fluxerr[0:c-1]
  ;;bin_corrflux = bin_corrflux[0:c-1]
  bin_timearr = bin_timearr[0:c-1]
  bin_bmjdarr = bin_bmjdarr[0:c-1]
  ;;bin_corrfluxerr = bin_corrfluxerr[0:c-1]
;;  bin_phase = bin_phase[0:c-1]
;;  bin_ncorr = bin_ncorr[0:c-1]
;;  bin_np = bin_np[0:c-1]
;;  bin_npcent = bin_npcent[0:c-1]
;;  bin_xfwhm = bin_xfwhm[0:c-1]
;;  bin_yfwhm = bin_yfwhm[0:c-1]
  
  ;;bin_centerpix = bin_centerpix[0:c-1]
  ;;print, 'bin_xcen', bin_xcen
  ;;bin_bkgderr = bin_bkgderr[0:c-1]
     
  ;;---------------------------------------------
  ;;do it again for the sweet spot points.
  ;;xxx should clean this up and make it a function
  cp = 0
  for j = 0L, n_elements(hp) - 1 do begin
     
;get rid of the bins with no values and low numbers, meaning low overlap
     if (rip[j+1] gt rip[j] + 2)  then begin ;require 3 elements in the bin
;           print, 'binning together', n_elements(numberarr[rip[rip[j]:rip[j+1]-1]])
                                ;print, 'binning', numberarr[rip[rip[j]:rip[j+1]-1]]
        
        meanclip, xarrp[rip[rip[j]:rip[j+1]-1]], meanx, sigmax
        bin_xcenp[cp] = meanx   ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])
        
        meanclip, yarrp[rip[rip[j]:rip[j+1]-1]], meany, sigmay
        bin_ycenp[cp] = meany   ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])
        
;;        meanclip, xfwhmp[rip[rip[j]:rip[j+1]-1]], meanxfwhm, sigmaxfwhm
;;        bin_xfwhmp[cp] = meanxfwhm ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])
        
;;        meanclip, yfwhmp[rip[rip[j]:rip[j+1]-1]], meanyfwhm, sigmayfwhm
;;        bin_yfwhmp[cp] = meanyfwhm ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])
        
        ;;meanclip, centerpixarrp[ri[ri[j]:ri[j+1]-1]], meancenterpix, sigmacenterpix
        ;;bin_centerpixp[cp]= meancenterpix
        
;;        meanclip, bkgdarrp[rip[rip[j]:rip[j+1]-1]], meansky, sigmasky
;;        bin_bkgdp[cp] = meansky ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])
        
;;        meanclip, fluxp[rip[rip[j]:rip[j+1]-1]], meanflux, sigmaflux1
;;        bin_fluxp[cp] = meanflux ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])
        
;;        meanclip, nparrp[rip[rip[j]:rip[j+1]-1]], meannp, sigmanp
;;        bin_nparrp[cp] = meannp ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])
        
;;        meanclip, npcentarrp[rip[rip[j]:rip[j+1]-1]], meannpcent, sigmanpcent
;;        bin_npcentarrp[cp] = meannpcent ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])
        
;;        junk = where(finite(corrfluxp[rip[rip[j]:rip[j+1]-1]]) gt 0,ngood)
 ;;       bin_ncorrp[cp] = ngood
        
        meanclip, timearrp[rip[rip[j]:rip[j+1]-1]], meantimearr, sigmatimearr
        bin_timearrp[cp]=meantimearr
        
;           meanclip, phasearrp[rip[rip[j]:rip[j+1]-1]], meanphasearr, sigmabmjdarr
;;        meanphasearr = mean(phasearrp[rip[rip[j]:rip[j+1]-1]],/nan)
;           if a eq 24 then print, meanphasearr, bin_phasep[cp-1],  meanphasearr - bin_phasep[cp-1], phasearrp[rip[rip[j]:rip[j+1]-1]]
;;        if bin_phasep[cp-1] gt 0.48 and (meanphasearr - bin_phasep[cp-1]) lt -0.1 $
;;           and (meanphasearr - bin_phasep[cp-1]) gt -0.9 then meanphasearr = mean(abs(phasearrp[rip[rip[j]:rip[j+1]-1]]),/nan)
;           if a eq 24 then print, 'final meanphase', meanphasearr
;;        bin_phasep[cp]= meanphasearr
        
        meanbmjdarr = mean( bmjdarrp[rip[rip[j]:rip[j+1]-1]],/nan)
        bin_bmjdarrp[cp]= meanbmjdarr
        
                                ;xxxx this could change
                                ;ripght now it is just the scatter in the bins
;;        icorrdataerr = corrfluxerrp[rip[rip[j]:rip[j+1]-1]]
;;        icorrdata = corrfluxp[rip[rip[j]:rip[j+1]-1]]
;;        bin_corrfluxerrp[cp] =  sqrt(total(icorrdataerr^2))/ (n_elements(icorrdataerr))
        
        
                                ;can only compute means if there are values in there
;;        if pmapcorr eq 1 then begin
;;           meanclip, corrfluxp[rip[rip[j]:rip[j+1]-1]], meancorrflux, sigmacorrflux
;              meancorrflux = mean(corrflux[rip[rip[j]:rip[j+1]-1]],/nan)
;;           bin_corrfluxp[cp] = meancorrflux
                                ;bin_corrfluxerrp[cp] = sigmacorrflux
;;           meanclip, corrflux_dp[rip[rip[j]:rip[j+1]-1]], meancorrflux_d, sigmacorrflux_d
;;           bin_corrflux_dp[cp] = meancorrflux_d
;;        endif
        
;;        idataerr = fluxerrp[rip[rip[j]:rip[j+1]-1]]
;;        bin_fluxerrp[cp] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        
;           meanclip, corrfluxerrp[rip[rip[j]:rip[j+1]-1]], meancorrfluxerr, sigmacorrfluxerr
;           bin_corrfluxerrp[cp] = sigmacorrfluxerr
;           meanclip, fluxerrp[rip[rip[j]:rip[j+1]-1]], meanfluxerr, sigmafluxerr
;           bin_fluxerrp[cp] = sigmafluxerr
        
;           meanclip, bkgderrarr[rip[rip[j]:rip[j+1]-1]], meanbkgderrarr, sigmabkgderrarr
                                ;          bin_bkgderr[cp] = meanbkgderrarr
        
        cp = cp + 1
                                ;print, 'testing', j, phasearr[ri[ri[j]:ri[j+1]-1]]
     endif
  endfor
  
  bin_xcenp = bin_xcenp[0:cp-1]
  bin_ycenp = bin_ycenp[0:cp-1]
;;  bin_bkgdp = bin_bkgdp[0:cp-1]
;;  bin_fluxp = bin_fluxp[0:cp-1]
;;  bin_fluxerrp = bin_fluxerrp[0:cp-1]
;;  bin_corrfluxp = bin_corrfluxp[0:cp-1]
;;  bin_corrflux_dp = bin_corrflux_dp[0:cp-1]
  bin_timearrp = bin_timearrp[0:cp-1]
  bin_bmjdarrp = bin_bmjdarrp[0:cp-1]
;;  bin_corrfluxerrp = bin_corrfluxerrp[0:cp-1]
;;  bin_phasep = bin_phasep[0:cp-1]
;;  bin_ncorrp = bin_ncorrp[0:cp-1]
;;  bin_nparrp = bin_nparrp[0:cp-1]
;;  bin_npcentarrp = bin_npcentarrp[0:cp-1]
;;  bin_xfwhmp = bin_xfwhmp[0:cp-1]
;;  bin_yfwhmp = bin_yfwhmp[0:cp-1]
;
;     bin_centerpixp = bin_centerpixp[0:cp-1]
;  bin_bkgderrp = bin_bkgderrp[0:cp-1]
  
  
  
  
  return, 0
  
end


                               ;can only compute means if there are values in there
;           if pmapcorr eq 1 then begin
;              meanclip, corrflux[ri[ri[j]:ri[j+1]-1]], meancorrflux, sigmacorrflux
;              meancorrflux = mean(corrflux[ri[ri[j]:ri[j+1]-1]],/nan)
;              bin_corrflux[c] = meancorrflux
;           print, 'finished corrfluxx'
;           endif
        
                                ;xxxx this could change
                                ;right now it is just the scatter in the bins
;           meanclip, corrfluxerr[ri[ri[j]:ri[j+1]-1]], meancorrfluxerr, sigmacorrfluxerr
;           bin_corrfluxerr[c] = sigmacorrfluxerr
        
        
;           meanclip, centerpixarr[ri[ri[j]:ri[j+1]-1]], meancenterpix, sigmacenterpix
;           bin_centerpix[c]= meancenterpix
        
;            meanclip, phasearr[ri[ri[j]:ri[j+1]-1]], meanphase, sigmaphasearr
;            bin_phase[c]= meanphase
        
        
        
;           meanclip, bkgderrarr[ri[ri[j]:ri[j+1]-1]], meanbkgderrarr, sigmabkgderrarr
                                ;          bin_bkgderr[c] = meanbkgderrarr
        
 
