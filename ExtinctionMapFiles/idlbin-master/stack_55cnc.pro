pro stack_55cnc
  aorname_55cnc = ['r43981312','r43981568','r43981824','r43981056'] ;ch2
  utmjd_center = [double(55944.25889), double(55947.20505),double(55949.41467), double(55957.51661)]  
  duration = 105.7              ;in minutes
  duration = duration /60./24.  ; in days
  period = 0.73654 ;days

;chain together all 4 AORs

  for a = 0, n_elements(aorname_55cnc) - 1 do  begin

     fb = strcompress('/Users/jkrick/irac_warm/pcrs_planets/55cnc/selfcal_bmjd_' + string(a) + '.sav',/remove_all)
     fs = strcompress('/Users/jkrick/irac_warm/pcrs_planets/55cnc/selfcal_sub_' + string(a) + '.sav',/remove_all)
     fe = strcompress('/Users/jkrick/irac_warm/pcrs_planets/55cnc/selfcal_fluxerr_' + string(a) + '.sav',/remove_all)

     restore,fb
     restore, fs
     restore, fe

      if a eq 0 then begin
         phasearr = full_bmjd
         print, 'mean nele', a, n_elements(full_bmjd), n_elements(full_sub), n_elements(full_fluxerr)
         fluxarr = full_sub / mean(full_sub)
         fluxerrarr = full_fluxerr / mean(full_sub)
      endif else begin
         phasearr = [phasearr, full_bmjd]
         print, 'mean nele', a, n_elements(full_bmjd), n_elements(full_sub), n_elements(full_fluxerr)
         add = full_sub/mean(full_sub)
         adderr = full_fluxerr / mean(full_sub)
         fluxarr = [fluxarr, add]
         fluxerrarr = [fluxerrarr, adderr]
      endelse


;      phasearr = [AOR55cnc[0].bmjdarr, AOR55cnc[1].bmjdarr, AOR55cnc[2].bmjdarr, AOR55cnc[3].bmjdarr]
;      fluxarr =  [AOR55cnc[0].flux, AOR55cnc[1].flux, AOR55cnc[2].flux, AOR55cnc[3].flux]
;      fluxerrarr =  [AOR55cnc[0].fluxerr, AOR55cnc[1].fluxerr, AOR55cnc[2].fluxerr, AOR55cnc[3].fluxerr]
;      corrfluxarr =  [AOR55cnc[0].corrflux, AOR55cnc[1].corrflux, AOR55cnc[2].corrflux, AOR55cnc[3].corrflux]
;      corrfluxerrarr =  [AOR55cnc[0].corrfluxerr, AOR55cnc[1].corrfluxerr, AOR55cnc[2].corrfluxerr, AOR55cnc[3].corrfluxerr]
     print, 'n phases', n_elements(phasearr), n_elements(fluxarr), n_elements(fluxerrarr)
      
  endfor

  ;quick check on the fluxes and their errors
  plothist, fluxarr, xhist, yhist, /noplot, bin = 0.001
  pt = plot(xhist, yhist)
  plothist, fluxerrarr, xhist, yhist, /noplot, bin = 0.001
  pt2 = plot(xhist, yhist)

  print, 'n flux', n_elements(fluxarr), n_elements(fluxerrarr)
;histogram phase by 0.01 as a technique to stack

  h = histogram(phasearr, OMIN=om, binsize = 0.008, reverse_indices = ri)
  print, 'omin', om, 'nh', n_elements(h)

;mean together the flux values in each phase bin
  binflux = dblarr(n_elements(h))
  binphase = dblarr(n_elements(h))
  binfluxerr= dblarr(n_elements(h))
  c = 0
  for j = 0L, n_elements(h) - 1 do begin

;get rid of the bins with no values and low numbers, meaning low overlap
     if (ri[j+1] gt ri[j]) and ( n_elements(fluxarr[ri[ri[j]:ri[j+1]-1]]) gt 37000) then begin
        print, 'binning together', n_elements(fluxarr[ri[ri[j]:ri[j+1]-1]])
        meanclip, fluxarr[ri[ri[j]:ri[j+1]-1]], mean, sigma
        binflux[c] = mean; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])
        binphase[c] = phasearr[ri[ri[j]]]
        n = n_elements(fluxarr[ri[ri[j]:ri[j+1]-1]])
        binfluxerr[c] = sigma/sqrt(n); (mean(fluxerrarr[ri[ri[j]:ri[j+1]-1]])) / sqrt(n)
        c = c + 1
        ;print, 'testing', j, phasearr[ri[ri[j]:ri[j+1]-1]]
     endif
  endfor

  binflux = binflux[0:c-1]
  binphase = binphase[0:c-1]
  binfluxerr = binfluxerr[0:c-1]

;plot the output
;  p = plot(binphase, binflux,  '1o',sym_filled  = 1, xtitle = 'Phase', ytitle = 'Normalized Corrected Flux')
  p1 = errorplot(binphase, binflux, binfluxerr, '1o',sym_filled  = 1, sym_size = 0.8, xtitle = 'Phase', ytitle = 'Normalized Corrected Flux', title = '4 stacked light curves', xrange = [0.32, 0.62], yrange = [0.9997, 1.0002], name = 'simple reduction')

;save, binflux, filename = '/Users/jkrick/idlbin/binflux.sav'
;save, binphase, filename = '/Users/jkrick/idlbin/binphase.sav'


;ok, what does Demory's global fit look like overplotted here?
x = [0.3, 0.455, 0.46, 0.54, 0.545, 0.7]
y = [1.0, 1.0, 0.99988, 0.99988, 1.0, 1.0]
p2= plot(x, y, color = 'red',/overplot, name = 'Demory et al. model')

;l = legend(target = [p1,p2], position = [0.55,0.9998],/data)

;p2.save, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/stacked_selfcal_3_11_15.png'
end
