pro plot_hatp17_ch2

  restore, '/Users/jkrick/irac_warm/pcrs_planets/hatp17/hatp17_phot_ch2.sav'

;AORs = replicate({planetob, ra:0D,dec:0D,xcen:dblarr(nfits),ycen:dblarr(nfits),flux:dblarr(nfits),fluxerr:dblarr(nfits), corrflux:dblarr(nfits), corrfluxerr:dblarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', bmjdarr:dblarr(nfits), utcsarr:dblarr(nfits)},n_elements(aorname))


  a = plot(AORs.xcen, AORs.ycen, '1.', xrange = [14.5, 15.5], yrange = [14.5, 15.5], xtitle = 'x pix', ytitle = 'y pix')
  b = plot((AORs.timearr- AORs.sclktime_0)/60./60., AORs.xcen, '1.', yrange =  [14.8, 15.2], xtitle = 'time(hrs)', ytitle = 'x pix')
  c = plot((AORs.timearr- AORs.sclktime_0)/60./60., AORs.ycen, '1.', yrange =  [14.8, 15.2], xtitle = 'time(hrs)', ytitle = 'y pix');
;  d = plot((AORs.timearr- AORs.sclktime_0)/60./60., AORs.corrflux, '1.', yrange = [0.06, 0.064], xtitle = 'time(hrs)', ytitle = 'pmap corrected flux')
  
;practice binning
  bin_level = 63
  nframes = bin_level
  nfits = long(n_elements(AORs.corrflux)) / nframes
  print, 'nframes, nfits', nframes, nfits
  
  bin_corrflux= dblarr(nfits)
  bin_corrfluxerr= dblarr(nfits)
  bin_timearr = dblarr(nfits)
  bin_bmjdarr = dblarr(nfits)

  for si = 0L, long(nfits) - 1L do begin
                                ;print, 'working on si', si, n_elements(corrflux), si*nframes, si*nframes + (nframes - 1)
     idata = AORs.corrflux[si*nframes:si*nframes + (nframes - 1)]
;       idata2 = corrflux2[si*nframes:si*nframes + (nframes - 1)]
     id = AORs.flux[si*nframes:si*nframes + (nframes - 1)]
     idataerr = AORs.corrfluxerr[si*nframes:si*nframes + (nframes - 1)]
 
     bin_corrflux[si] = mean(idata,/nan)
     bin_corrfluxerr[si] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
     bin_timearr[si]= mean(AORs.timearr[si*nframes:si*nframes + (nframes - 1)])
     bin_bmjdarr[si]= mean(AORs.bmjdarr[si*nframes:si*nframes + (nframes - 1)])
  endfor                        ;for each fits image

;   d = plot((bin_timearr- bin_timearr[0])/60./60., bin_corrflux, '1.', yrange = [0.0622, 0.063], xtitle = 'time(hrs)', ytitle = 'pmap corrected flux')
   d = plot(bin_bmjdarr, bin_corrflux, '1.', yrange = [0.0622, 0.063], xtitle = 'Time (bmjd)', ytitle = 'pmap corrected flux')

  
end
