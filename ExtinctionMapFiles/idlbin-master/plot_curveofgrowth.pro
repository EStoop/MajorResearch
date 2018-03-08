pro plot_curveofgrowth, sigma = sigma, selfcal = selfcal, nn = nn
  planetname = ['wasp14', 'wasp62', 'hatp26']
  colorarr = ['red', 'black', 'blue']
  for pn = 0, n_elements(planetname) - 1 do begin
     fn = '/Users/jkrick/irac_warm/pcrs_planets/'+planetname(pn)+'/curveofgrowth.sav' 
     print, fn
     restore, filename = fn
     normfactor = meanarr(4)/stdarr(4)
     if keyword_set(sigma) then begin
        normfactor = stdarr(4)
        p = plot(aps1, (stdarr)/normfactor, '1s-', xtitle = 'aperture size', ytitle = 'Sigma', $
                 name = planetname(pn), sym_filled = 1, thick = 2, yrange = [0.9, 1.1], xrange = [1.5,3.5], $
                 color = colorarr(pn),/overplot)
        l = legend(TARGET=planetname, POSITION=[3.0,0.98], /DATA, /AUTO_TEXT_COLOR)
  
     endif else begin
        normfactor = meanarr(4)/stdarr(4)
        p = plot(aps1, (meanarr/stdarr)/normfactor, '1s-', xtitle = 'aperture size', ytitle = 'SNR', $
                 name = planetname(pn), sym_filled = 1, thick = 2, yrange = [0.9, 1.05], xrange = [1.5,3.5], $
                 color = colorarr(pn), /overplot)
        l = legend(TARGET=planetname, POSITION=[2.8,0.94], /DATA, /AUTO_TEXT_COLOR)
  
     endelse
  endfor  ;n_elements(planetname)

 
;-------------
;now to check the selfcal corrected fluxes
  if keyword_set(selfcal) then begin
     p.close
     aps1 = [ 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25]
     meanarr = fltarr(n_elements(aps1))
     stdarr = fltarr(n_elements(aps1))
     for ap = 0, n_elements(aps1) -1 do begin
        filename=strcompress('/Users/jkrick/irac_warm/pcrs_planets/wasp14/selfcal'+string(aps1(ap))+'.sav',/remove_all)
        print, filename
        restore, filename
        print, 'y217', y(217)
        meanclip, y, meanflux, sigmacorr, clipsig = 3
        
        meanarr(ap) = meanflux
        stdarr(ap) = sigmacorr
     endfor
     
     print, 'sigma', stdarr

     if keyword_set(sigma) then begin
        normfactor = stdarr(4)
        s = plot(aps1, (stdarr)/normfactor, '1s-', xtitle = 'aperture size', ytitle = 'Sigma', $
                 sym_filled = 1, thick = 2, yrange = [0.9, 1.1], xrange = [1.5,3.5], $
                 color = 'black',/overplot)
        
     endif else begin
        normfactor = meanarr(4)/stdarr(4)
        s = plot(aps1, (meanarr/stdarr)/normfactor, '1s-', xtitle = 'aperture size', ytitle = 'SNR', $
                 sym_filled = 1, thick = 2, yrange = [0.9, 1.05], xrange = [1.5,3.5], $
                 color = 'black', /overplot)
        
     endelse
  endif  ;selfcal



;-------------
;now to check the selfcal corrected fluxes
  if keyword_set(nn) then begin
     p.close
     aps1 = [ 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25]
     meanarr = fltarr(n_elements(aps1))
     stdarr = fltarr(n_elements(aps1))
     for ap = 0, n_elements(aps1) -1 do begin
        filename =strcompress('/Users/jkrick/irac_warm/pcrs_planets/wasp14/pixphasecorr_ch2_r45428992'+string(aps1(ap))+'.sav',/remove_all)
        print, filename
        restore, filename
        print, 'flux_np', flux_np(217)
        meanclip, flux_np, meanflux, sigmacorr, clipsig = 3
        
        meanarr(ap) = meanflux
        stdarr(ap) = sigmacorr
     endfor
     
     print, 'sigma', stdarr

     if keyword_set(sigma) then begin
        normfactor = stdarr(4)
        s = plot(aps1, (stdarr)/normfactor, '1s-', xtitle = 'aperture size', ytitle = 'Sigma', $
                 sym_filled = 1, thick = 2, yrange = [0.9, 1.1], xrange = [1.5,3.5], $
                 color = 'black',/overplot)
        
     endif else begin
        normfactor = meanarr(4)/stdarr(4)
        s = plot(aps1, (meanarr/stdarr)/normfactor, '1s-', xtitle = 'aperture size', ytitle = 'SNR', $
                 sym_filled = 1, thick = 2, yrange = [0.9, 1.05], xrange = [1.5,3.5], $
                 color = 'black', /overplot)
        
     endelse
  endif  ;selfcal


end

