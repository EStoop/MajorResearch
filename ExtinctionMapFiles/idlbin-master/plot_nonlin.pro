pro plot_nonlin, planetname, chname, apradius, breatheap = breatheap
;example call:
;plot_nonlin, ['HD209458', 'pmap_star_ch2','WASP-14b', 'WASP-62b', 'HD158460','HAT-P-8', 'pmap_star_ch2'],'2', 2.25
;plot_nonlin, ['HD209458', 'pmap_star_ch2','HAT-P-8', 'HD158460', 'pmap_star_ch2', 'WASP-14b'],'2', 2.25
;plot_nonlin, ['HD209458','HAT-P-8', 'HD158460', 'WASP-62b', 'WASP-14b'],'2', 2.25
;'pmap_star_ch2', 'WASP-14b','WASP-62b', 'HD209458','HAT-P-8', 'HD158460'], '2'
 
 colorarr = ['magenta',  'cyan','black', 'red','gray',  'Maroon', 'blue']
  for pl = 0,  n_elements(planetname) - 1 do begin
     undefine, xcen
     print, 'working on planet', planetname(pl)
  ;read in the relevant photometry
     planetinfo = create_planetinfo()
     if chname eq '2'  then aorname= planetinfo[planetname(pl), 'aorname_ch2'] else aorname = planetinfo[planetname(pl), 'aorname_ch1'] 
     basedir = planetinfo[planetname(pl), 'basedir']
     dirname = strcompress(basedir + planetname(pl) +'/')
    if keyword_set(breatheap) then begin
       savename = strcompress(dirname + planetname(pl) +'_phot_ch'+chname+'_varap.sav')
    endif else begin
       savename = strcompress(dirname + planetname(pl) +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
    endelse
    print, 'savename', savename

     restore, savename
     
                                ;put all the relevant information from the aors together and sort one for exptime
     if planetname(pl) eq 'pmap_star_ch2' then begin
                                ;track only the 0.1s frametimes
        if pl eq 0 then begin
           for a = 0, n_elements(aorname) -1 do begin
              if planethash[aorname(a),'exptime'] lt 0.2 then begin
                 if n_elements(xcen) eq 0 then begin ; hasn't been defined yet, aka, the first aor with this exptime
                    xcen= planethash[aorname(a),'xcen']
                    ycen= planethash[aorname(a),'ycen']
                    flux = planethash[aorname(a),'flux']
                    fluxerr = planethash[aorname(a),'fluxerr']
                    corrflux = planethash[aorname(a),'corrflux']
                    peakpixDN = planethash[aorname(a),'peakpixDN']
                    np = planethash[aorname(a),'npcentroids']
;                    framedly =  planethash[aorname(a),'framedly']
                 endif else begin
                    xcen= [xcen, planethash[aorname(a),'xcen']]
                    ycen= [ycen, planethash[aorname(a),'ycen']]
                    flux= [flux, planethash[aorname(a),'flux']]
                    fluxerr= [fluxerr, planethash[aorname(a),'fluxerr']]
                    corrflux= [corrflux, planethash[aorname(a),'corrflux']]
                    peakpixDN = [peakpixDN, planethash[aorname(a),'peakpixDN']]
                    np = [np, planethash[aorname(a),'npcentroids']]
;                    framedly =  [framedly, planethash[aorname(a),'framedly']]

                endelse
              endif
           endfor
        endif                              
;track only the 0.4s frametimes
        if pl eq 3 then begin
           for a = 0, n_elements(aorname) -1 do begin
              if planethash[aorname(a),'exptime'] gt 0.2 then begin
                 if n_elements(xcen) eq 0 then begin ; hasn't been defined yet, aka, the first aor with this exptime
                    xcen= planethash[aorname(a),'xcen']
                    ycen= planethash[aorname(a),'ycen']
                    flux = planethash[aorname(a),'flux']
                    fluxerr = planethash[aorname(a),'fluxerr']
                    corrflux = planethash[aorname(a),'corrflux']
                    peakpixDN = planethash[aorname(a),'peakpixDN']
                    np = planethash[aorname(a),'npcentroids']
;                    framedly =  planethash[aorname(a),'framedly']

                 endif else begin
                    xcen= [xcen, planethash[aorname(a),'xcen']]
                    ycen= [ycen, planethash[aorname(a),'ycen']]
                    flux= [flux, planethash[aorname(a),'flux']]
                    fluxerr= [fluxerr, planethash[aorname(a),'fluxerr']]
                    corrflux= [corrflux, planethash[aorname(a),'corrflux']]
                    peakpixDN = [peakpixDN, planethash[aorname(a),'peakpixDN']]
                    np = [np, planethash[aorname(a),'npcentroids']]
;                    framedly =  [framedly, planethash[aorname(a),'framedly']]

                 endelse
              endif
           endfor
        endif                                


     endif else begin  ; not working with the pmap_star_ch2
        
        for a = 0, n_elements(aorname) -1 do begin
           if a eq 0 then begin
              xcen= planethash[aorname(a),'xcen']
              ycen= planethash[aorname(a),'ycen']
              flux = planethash[aorname(a),'flux']
              fluxerr = planethash[aorname(a),'fluxerr']
              corrflux = planethash[aorname(a),'corrflux']
              peakpixDN = planethash[aorname(a),'peakpixDN']
              np = planethash[aorname(a),'npcentroids']
;              framedly =  planethash[aorname(a),'framedly']
           endif else begin
              xcen= [xcen, planethash[aorname(a),'xcen']]
              ycen= [ycen, planethash[aorname(a),'ycen']]
              flux= [flux, planethash[aorname(a),'flux']]
              fluxerr= [fluxerr, planethash[aorname(a),'fluxerr']]
             corrflux= [corrflux, planethash[aorname(a),'corrflux']]
              peakpixDN = [peakpixDN, planethash[aorname(a),'peakpixDN']]
              np = [np, planethash[aorname(a),'npcentroids']]
;              framedly =  [framedly,planethash[aorname(a),'framedly']]

           endelse
           
        endfor
     endelse
     
;fix frame delay to have a value for each subframe 
; starts with one value per fits image, since that is all we have from the headers
;set the rest to zero?
;     fdarr = fltarr(n_elements(framedly) * 63)
;     for nfd = 0, n_elements(framedly) - 1 do begin
;        fdarr[nfd*63] = framedly[nfd]
;     endfor

                                ;what is the peakpix in DN
 ;    plothist, peakpixDN, xhist, yhist, /noplot,/nan
 ;    ppp = plot(xhist, yhist, overplot = ppp, color = colorarr[pl], xrange = [0, 2E4], yrange = [0, 400])
 ;    bz = where(peakpixDN eq 0, zerocount)
 ;    print, 'n elements', n_elements(xcen), zerocount
                                ;calculate distance from sweet spot
     xsweet = 15.12
     ysweet = 15.1  ; 15.0
     xdist = (xcen - xsweet)
     ydist = (ycen- ysweet)
     dsweet = sqrt(xdist^2 + ydist^2)
     gain = flux / corrflux  ; back out what the gain is at the position of the observation
     ;get rid of some nan's
     badgain = where(corrflux eq 0, badcount)
     print, 'number of bad corrfluxes', badcount, n_elements(gain)

                                ;bin as a function of dsweet (maybe in bins of 0.015 pixels)
     ;need to sort on gain first  XXX
     ds = sort(dsweet)
     dsweet = dsweet(ds)
     flux = flux(ds)
     np = np(ds)
     gain = gain(ds)
     fluxerr = fluxerr(ds)
     corrflux = corrflux(ds)
     peakpixDN = peakpixDN(ds)
     print, 'gain', gain[10:100]
;     framedly = framedly(ds)
  ;plot one of the flux/dsweets just to see what that looks like
  ; does it have the shape I expect it to?
 ;   if pl eq 0 then  p2 = plot(dsweet, flux, '1s', sym_size = 0.1, color = colorarr(pl), xrange = [0, 0.5], $
 ;              xtitle = 'Distance from Sweet Spot (Pix)', ytitle = 'Flux', yrange = [0.05, 0.07])


 ;    numberarr = findgen(n_elements(dsweet))
;     bin_number = 30
     bin_level = 0.015          ; 0.005
 ;    h = histogram(numberarr, OMIN=om, nbins  = bin_number, reverse_indices = ri)
;     h = histogram(gain, OMIN=om, binsize  = bin_level, reverse_indices = ri);, max = 0.5)
     h = histogram(dsweet, OMIN=om, binsize  = bin_level, reverse_indices = ri, max = 0.5)
     print, 'n_elements(h)', n_elements(h)
     bin_flux = dblarr(n_elements(h))
     bin_dsweet = bin_flux
     bin_peakpix = bin_flux
     bin_np = bin_flux
     bin_fd = bin_flux
     bin_gain = bin_flux
     bin_fluxerr = bin_flux
     c = 0
     for j = 0L, n_elements(h) - 1 do begin
;        print, 'testing ri', ri[j+1], ri[j] + 2
                                ;get rid of the bins with no values and low numbers, meaning low overlap
        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
;           print, 'dsweet', dsweet[ri[ri[j]:ri[j+1]-1]]
           meanclip_jk, dsweet[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_dsweet[c] = meanx    
           
           meanclip_jk, flux[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_flux[c] = meanx    
           idataerr = fluxerr[ri[ri[j]:ri[j+1]-1]]
           bin_fluxerr[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))

           meanclip_jk, peakpixDN[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_peakpix[c] = meanx    

           meanclip_jk, np[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_np[c] = meanx    

           meanclip_jk, gain[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_gain[c] = meanx    
;           meanclip_jk, framedly[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
;           bin_fd[c] = meanx 
   
           c = c +1
        endif
     endfor
     print, 'c', c
     bin_flux = bin_flux[0:c-1]
     bin_dsweet = bin_dsweet[0:c-1]
     bin_peakpix = bin_peakpix[0:c-1]
     bin_np = bin_np[0:c-1]
     bin_gain = bin_gain[0:c-1]
     bin_fluxerr = bin_fluxerr[0:c-1]
;     bin_fd = bin_fd[0:c-1]
;------------------------------------
;plot
                               
     normdist = 14
     normflux = bin_flux(normdist)
     print, 'normflux', normflux
     print, 'bin_dsweet', bin_dsweet
     print, 'bin_np', bin_np
     print, 'bin_flux', bin_flux
     p = errorplot(bin_dsweet, bin_flux/normflux, bin_fluxerr, '1s', sym_size = 1,   $
                   sym_filled = 1,xtitle = 'Distance from Sweet Spot (Pix)', ytitle = 'F/ F(0.3)',$
                   xrange = [0, 0.5], yrange = [0.99, 1.01], color = colorarr(pl), overplot = p,$
                  errorbar_color = colorarr(pl))
;     pt = text(0.01, 0.991 + pl*0.001, planetname(pl) + string(fix(bin_peakpix(0))), /current, /data, color = colorarr(pl))
;     a = where(finite(flux) lt 1, nancount)
;     print, 'nancount', nancount, n_elements(flux)

     p2 = plot(dsweet, peakpixDN / mean(peakpixDN[0:100]), '1s', sym_size = 0.1,   sym_filled = 1, xtitle = 'Dsweet', ytitle = 'PeakPix (DN)',$
              xrange = [0, 0.3], color = colorarr(pl), overplot = p2, yrange = [0.6, 1.1])
     pt = text(0.01, 0.62 + pl*0.02, planetname(pl) + string(fix(bin_peakpix(0))), /current, /data, color = colorarr(pl))

;try a 3D plot with distance from sweet spot, normalized flux, and np
     ;what have I done by binning on one variable?

;     frametimearr = ['aqua', 'aqua ', 'black', ' orange', 'black', ' orange', 'black']
     frametimearr = ['aqua', 'aqua ', 'black', 'black', ' orange', 'black', ' orange']
;     frametimearr = ['aqua', 'black', ' orange', 'black',  'black']
 ;    p3d = plot3d(bin_dsweet,bin_np, bin_flux/normflux, 'o', /sym_filled, xrange = [0, 0.3], yrange = [4, 6], $ ;
 ;                 zrange = [0.985, 1.005], color = frametimearr(pl), overplot = p3d, axis_style = 2,$; margin =[0.2, 0.3, 0.3, 0], $
 ;                 XMINOR=0, YMINOR=0, ZMINOR=0, $
 ;                 DEPTH_CUE=[0, 2], /PERSPECTIVE, $
 ;                 SHADOW_COLOR=frametimearr(pl), $
 ;                 XY_SHADOW=1, YZ_SHADOW=1, XZ_SHADOW=1, $
 ;                 XTITLE='Dsweet', YTITLE='NP',ztitle = 'F/F(0)',/current)
 ;    ax = p3d.AXES
 ;    ax[2].HIDE = 1
 ;    ax[6].HIDE = 1
 ;    ax[7].HIDE = 1
 ;    if pl eq n_elements(planetname) -1 then zaxis = axis('Z', location = [0.0,6], title = 'F/F(0)')

;;     pt = text(0.1, 4.2, 0.992 + pl*0.002, ftarr(pl), /current, /data, color = colorarr(pl))
;;    xplot3d, bin_dsweet,bin_np, bin_flux/normflux,/overplot, color = colorarr(pl), XTITLE='Dsweet', YTITLE='NP', Ztitle = 'F/F(d = 0.3)',$
;;              xrange = [0, 0.3], yrange = [4,6], zrange = [0.98, 1.01]
     
;     pe = plot(bin_dsweet, bin_np, '1s', overplot = pe)



  endfor                        ; for each planet pl
;     pt = text(0.2, 4.1, .992, '0.1s', /current, /data, color = 'orange')
;     pt = text(0.2, 4.25, .992, '0.4s', /current, /data, color = 'aqua')
;     pt = text(0.2, 4.4, .992, '2.0s', /current, /data, color = 'black')


end
