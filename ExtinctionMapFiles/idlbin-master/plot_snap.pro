pro plot_snap, planetname, bin_level

;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
  chname = planetinfo[planetname, 'chname']
  ra_ref = planetinfo[planetname, 'ra']
  dec_ref = planetinfo[planetname, 'dec']
  aorname = planetinfo[planetname, 'aorname']
  basedir = planetinfo[planetname, 'basedir']
  utmjd_center =  planetinfo[planetname, 'utmjd_center']
  transit_duration =  planetinfo[planetname, 'transit_duration']
  period =  planetinfo[planetname, 'period']
  intended_phase = planetinfo[planetname, 'intended_phase']
  stareaor = planetinfo[planetname, 'stareaor']
  plot_norm= planetinfo[planetname, 'plot_norm']
  plot_corrnorm = planetinfo[planetname, 'plot_corrnorm']
  
  dirname = strcompress(basedir + planetname +'/')
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'.sav')
  restore, savefilename

 colorarr = ['blue', 'red', 'deep_pink','fuchsia', 'magenta', 'medium_purple','medium_orchid', 'orchid', 'violet', 'plum', 'thistle', 'pink', 'orange_red', 'light_pink', 'rosy_brown','pale_violet_red',  'chocolate', 'saddle_brown', 'maroon', 'hot_pink', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red', 'deep_pink','fuchsia', 'magenta', 'medium_purple','medium_orchid', 'orchid', 'violet', 'plum', 'thistle', 'pink', 'orange_red', 'light_pink', 'rosy_brown','pale_violet_red',  'chocolate', 'saddle_brown', 'maroon', 'hot_pink', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki','black', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua' ]


;-------------------------------------------------------------------------------------
;try binning the fluxes by subarray image
 for a = 0,  n_elements(aorname) - 1 do begin
     print, 'working on AOR', a, '   ', aorname(a)

     timearr = [planethash[aorname(a),'timearr']]
     fluxarr = [ planethash[aorname(a),'flux']]
     fluxerrarr = [ planethash[aorname(a),'fluxerr']]
     corrfluxarr = [ planethash[aorname(a),'corrflux']]
     corrfluxerrarr = [planethash[aorname(a),'corrfluxerr']]
     xarr = [ planethash[aorname(a),'xcen']]
     yarr = [planethash[aorname(a),'ycen']]
     bkgd = [ planethash[aorname(a),'bkgd']]
     bmjd = [ planethash[aorname(a),'bmjdarr']]
     np = [ planethash[aorname(a),'np']]
     phase = [ planethash[aorname(a),'phase']]

     ;check if I should be using pmap corr or not
     ncorr = where(finite(corrfluxarr) gt 0)
     ;if 20% of the values are correctable than go with the pmap corr 
     print, 'nflux, ncorr, ', n_elements(fluxarr), n_elements(ncorr)
     if n_elements(ncorr) gt 0.2*n_elements(fluxarr) then pmapcorr = 1 else pmapcorr = 0

     ;remove outliers, 
;     if pmapcorr eq 1 then begin 
     goodpmap = where(xarr lt mean(xarr) + 2.5*stddev(xarr) and xarr gt mean(xarr) -2.5*stddev(xarr) and xarr lt mean(xarr) +3.0*stddev(yarr) and yarr gt mean(yarr) - 3.0*stddev(yarr) and finite(corrfluxarr) gt 0  ,ngood_pmap, complement=badpmap) ;
 ;    endif else begin
     good = where(xarr lt mean(xarr) + 2.5*stddev(xarr) and xarr gt mean(xarr) -2.5*stddev(xarr) and xarr lt mean(xarr) +3.0*stddev(yarr) and yarr gt mean(yarr) - 3.0*stddev(yarr) ,ngood, complement=bad)
 ;    endelse
     

     print, 'bad ',n_elements(bad), n_elements(good)
     print, 'badp ',n_elements(badpmap), n_elements(goodpmap)

     xarr = xarr[good]
     yarr = yarr[good]
     timearr = timearr[good]
     flux = fluxarr[good]
     fluxerr = fluxerrarr[good]
     corrflux = corrfluxarr[good]
     corrfluxerr = corrfluxerrarr[good]
     bmjdarr = bmjd[good]
     bkgdarr = bkgd[good]
     phasearr = phase[good]
     nparr = np[good]
     ;print, 'test phase early ', phasearr[0:10]

;and a second set for those that are in the sweet spot
     xarrp = xarr[goodpmap]
     yarrp = yarr[goodpmap]
     timearrp = timearr[goodpmap]
     fluxp = fluxarr[goodpmap]
     fluxerrp= fluxerrarr[goodpmap]
     corrfluxp = corrfluxarr[goodpmap]
     corrfluxerrp = corrfluxerrarr[goodpmap]
     bmjdarrp = bmjd[goodpmap]
     bkgdarrp = bkgd[goodpmap]
     phasearrp = phase[goodpmap]
     nparrp = np[goodpmap]
     phasearrp = phase[goodpmap]

; binning
     numberarr = findgen(n_elements(xarr))
     h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
     ;print, 'omin', om, 'nh', n_elements(h)

     numberarrp = findgen(n_elements(xarrp))
     hp = histogram(numberarrp, OMIN=omp, binsize = bin_level, reverse_indices = rip)
     ;print, 'ominp', omp, 'nhp', n_elements(hp), n_elements(xarrp)


;mean together the flux values in each phase bin
     bin_flux = dblarr(n_elements(h))
     bin_fluxerr = dblarr(n_elements(h))
     bin_corrflux= dblarr(n_elements(h))
     bin_corrfluxerr= dblarr(n_elements(h))
     bin_ncorr = dblarr(n_elements(h))
     bin_timearr = dblarr(n_elements(h))
     bin_bmjdarr = dblarr(n_elements(h))
     bin_bkgd = dblarr(n_elements(h))
     bin_bkgderr = dblarr(n_elements(h))
     bin_xcen = dblarr(n_elements(h))
     bin_ycen = dblarr(n_elements(h))
     bin_phase = dblarr(n_elements(h))
     bin_centerpix = dblarr(n_elements(h))
     bin_np = dblarr(n_elements(h))

     bin_fluxp = dblarr(n_elements(hp))
     bin_fluxerrp = dblarr(n_elements(hp))
     bin_corrfluxp= dblarr(n_elements(hp))
     bin_corrfluxerrp= dblarr(n_elements(hp))
     bin_ncorrp = dblarr(n_elements(hp))
     bin_timearrp = dblarr(n_elements(hp))
     bin_bmjdarrp = dblarr(n_elements(hp))
     bin_bkgdp = dblarr(n_elements(hp))
     bin_bkgderrp = dblarr(n_elements(hp))
     bin_xcenp = dblarr(n_elements(hp))
     bin_ycenp= dblarr(n_elements(hp))
     bin_phasep = dblarr(n_elements(hp))
     bin_centerpixp = dblarr(n_elements(hp))
     bin_nparrp = dblarr(n_elements(hp))

     c = 0
     for j = 0L, n_elements(h) - 1 do begin

;get rid of the bins with no values and low numbers, meaning low overlap
        if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
        ;print, 'debug inside if loop'
           meanclip, xarr[ri[ri[j]:ri[j+1]-1]], meanx, sigmax
           bin_xcen[c] = meanx   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

           meanclip, yarr[ri[ri[j]:ri[j+1]-1]], meany, sigmay
           bin_ycen[c] = meany   ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

           meanclip, bkgdarr[ri[ri[j]:ri[j+1]-1]], meansky, sigmasky
           bin_bkgd[c] = meansky ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

           meanclip, flux[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
           bin_flux[c] = meanflux ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

           idataerr = fluxerr[ri[ri[j]:ri[j+1]-1]]
           bin_fluxerr[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
;           meanclip, fluxerr[ri[ri[j]:ri[j+1]-1]], meanfluxerr, sigmafluxerr
;           bin_fluxerr[c] = sigmafluxerr

           meanclip, nparr[ri[ri[j]:ri[j+1]-1]], meannp, sigmanp
           bin_np[c] = meannp ; mean(fluxarr[ri[ri[j]:ri[j+1]-1]])

           junk = where(finite(corrflux[ri[ri[j]:ri[j+1]-1]]) gt 0,ngood)
           bin_ncorr[c] = ngood

           meanclip, timearr[ri[ri[j]:ri[j+1]-1]], meantimearr, sigmatimearr
           bin_timearr[c]=meantimearr
           
           meanclip, phasearr[ri[ri[j]:ri[j+1]-1]], meanphasearr, sigmaphasearr
           ;print, 'meanphasearr', meanphasearr
           bin_phase[c]= meanphasearr

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


           c = c + 1
        endif
     endfor
     
     bin_xcen = bin_xcen[0:c-1]
     bin_ycen = bin_ycen[0:c-1]
     bin_bkgd = bin_bkgd[0:c-1]
     bin_flux = bin_flux[0:c-1]
     bin_fluxerr = bin_fluxerr[0:c-1]
;     bin_corrflux = bin_corrflux[0:c-1]
     bin_timearr = bin_timearr[0:c-1]
;     bin_bmjdarr = bin_bmjdarr[0:c-1]
;     bin_corrfluxerr = bin_corrfluxerr[0:c-1]
     bin_phase = bin_phase[0:c-1]
     bin_ncorr = bin_ncorr[0:c-1]
     bin_np = bin_np[0:c-1]
;     bin_centerpix = bin_centerpix[0:c-1]
;     print, 'bin_phase', bin_phase
;  bin_bkgderr = bin_bkgderr[0:c-1]
     
     ;---------------------------------------------
     ;do it again for the sweet spot points.
     ;xxx should clean this up and make it a function
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

;           meanclip, centerpixarrp[ri[ri[j]:ri[j+1]-1]], meancenterpix, sigmacenterpix
;           bin_centerpixp[cp]= meancenterpix

           meanclip, bkgdarrp[rip[rip[j]:rip[j+1]-1]], meansky, sigmasky
           bin_bkgdp[cp] = meansky ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])

           meanclip, fluxp[rip[rip[j]:rip[j+1]-1]], meanflux, sigmaflux1
           bin_fluxp[cp] = meanflux ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])

           meanclip, nparrp[rip[rip[j]:rip[j+1]-1]], meannp, sigmanp
           bin_nparrp[cp] = meannp ; mean(fluxarr[rip[rip[j]:rip[j+1]-1]])

           junk = where(finite(corrfluxp[rip[rip[j]:rip[j+1]-1]]) gt 0,ngood)
           bin_ncorrp[cp] = ngood
           ;can only compute means if there are values in there
           if pmapcorr eq 1 then begin
              meanclip, corrfluxp[rip[rip[j]:rip[j+1]-1]], meancorrflux, sigmacorrflux
;              meancorrflux = mean(corrflux[rip[rip[j]:rip[j+1]-1]],/nan)
              bin_corrfluxp[cp] = meancorrflux
           endif

           meanclip, timearrp[rip[rip[j]:rip[j+1]-1]], meantimearr, sigmatimearr
           bin_timearrp[cp]=meantimearr
           
           meanclip, phasearrp[rip[rip[j]:rip[j+1]-1]], meanphasearr, sigmabmjdarr
           bin_phasep[cp]= meanphasearr

           ;xxxx this could change
           ;ripght now it is just the scatter in the bins
           icorrdataerr = corrfluxerrp[rip[rip[j]:rip[j+1]-1]]
           bin_corrfluxerrp[cp] =   sqrt(total(icorrdataerr^2))/ (n_elements(icorrdataerr))
           idataerr = fluxerrp[rip[rip[j]:rip[j+1]-1]]
           bin_fluxerrp[cp] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))

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
     bin_bkgdp = bin_bkgdp[0:cp-1]
     bin_fluxp = bin_fluxp[0:cp-1]
     bin_fluxerrp = bin_fluxerrp[0:cp-1]
     bin_corrfluxp = bin_corrfluxp[0:cp-1]
     bin_timearrp = bin_timearrp[0:cp-1]
;     bin_bmjdarrp = bin_bmjdarr[0:cp-1]
     bin_corrfluxerrp = bin_corrfluxerrp[0:cp-1]
     bin_phasep = bin_phasep[0:cp-1]
     bin_ncorrp = bin_ncorrp[0:cp-1]
     bin_nparrp = bin_nparrp[0:cp-1]
;     bin_centerpixp = bin_centerpixp[0:cp-1]
;  bin_bkgderrp = bin_bkgderrp[0:cp-1]
 ;-------------------------------------


;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------



;
;the plotting
;-------------------------------------------------------------------------------------
;-------------------------------------------------------------------------------------


     print, 'test phase', bin_phase(0)
     print, 'test normflux', bin_flux(0)/plot_norm
     
     setxrange = [0.,1.0]       ;[0.45, 0.55]
     
     if a eq 0 then begin       ; for the first AOR
        print, 'debug first aor plotting'

        pp = plot(bin_phase, bin_xcen, '1s', sym_size = 0.3,   sym_filled = 1, title = planetname, $
                  color = colorarr[a], xtitle = 'Orbital Phase', ytitle = 'X position', $
                  xrange = setxrange, yrange = [14.5, 15.5])
        
        pq = plot(bin_phase, bin_ycen, '1s', sym_size = 0.3,   sym_filled = 1, color = colorarr[a], $
                  xtitle = 'Orbital Phase', ytitle = 'Y position', title = planetname, $
                  xrange = setxrange, yrange = [14.5, 15.5])

        if keyword_set(errorbars) then begin
           pr = errorplot(bin_phase, bin_flux/plot_norm,bin_fluxerr/plot_norm,  '1s', sym_size = 0.3,   $
                          sym_filled = 1,  color = 'black',  xtitle = 'Orbital Phase', $
                          ytitle = 'Normalized Flux', title = planetname, yrange = [0.985,1.005], $
                          xrange = setxrange) 
        endif else begin
           pr = plot(bin_phase, bin_flux/plot_norm, '1s', sym_size = 0.3,   sym_filled = 1,  $
                     color ='black',  xtitle = 'Orbital Phase', ytitle = 'Normalized Flux', $
                     title = planetname, yrange = [0.985,1.005], xrange = setxrange) 
        endelse                 ;/errorbars
        
        if pmapcorr eq 1 then begin
           print, 'inside pmapcorr eq 1', median( (bin_corrfluxp/plot_corrnorm) - 0.01)
           if keyword_set(errorbars) then begin
              pr = errorplot(bin_phasep, (bin_corrfluxp/plot_corrnorm) - 0.01,  $
                             bin_corrfluxerrp/plot_corrnorm ,/overplot, '1s', sym_size = 0.3, $
                             sym_filled = 1, color = colorarr[a])
           endif else begin
              pr = plot(bin_phasep, (bin_corrfluxp/plot_corrnorm) - 0.01 ,/overplot, '1s', $
                        sym_size = 0.3,   sym_filled = 1, color =  colorarr[a])
           endelse              ;/errorbars
           
        endif                   ;enough pmap corrections
        
        ps= plot(bin_phase, bin_np, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], $
                 xtitle = 'Orbital Phase', ytitle = 'Noise Pixel', title = planetname,$
                 xrange = setxrange, yrange = [2., 10.])
        
        pt = plot(bin_phase, bin_bkgd, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], $
                  xtitle = 'Orbital Phase', ytitle = 'Background', title = planetname, $
                  xrange =setxrange)

     endif                      ; if a = 0

     if (a gt 0) and (a le stareaor) then begin

           pp.window.SetCurrent
           pp = plot(bin_phase, bin_xcen, '1s', sym_size = 0.3,   sym_filled = 1,color = colorarr[a],  /overplot,/current)
           pq.window.SetCurrent
           pq = plot(bin_phase, bin_ycen, '1s', sym_size = 0.3,   sym_filled = 1, color = colorarr[a], /overplot,/current)
           pr.window.SetCurrent
           pr = plot(bin_phase, bin_flux/plot_norm , '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], /overplot,/current) 
           if pmapcorr eq 1 then begin
              print, ' a gt 0, a lt stareaor and pmapcorr eq 1'
              pr = plot(bin_phasep, (bin_corrfluxp/plot_corrnorm)- 0.01,/overplot, 's1', sym_size = 0.3,   sym_filled = 1, color = 'black',/current)
           endif

           ps.window.SetCurrent
           ps = plot(bin_phase, bin_np, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], /overplot,/current) 
           pt.window.SetCurrent
           pt = plot(bin_phase, bin_bkgd, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], /overplot,/current)

        endif



        if a gt stareaor then begin
           print, 'inside a gt stareaor'
           pp.window.SetCurrent
           pp = plot(bin_phase, bin_xcen, '1s', sym_size = 0.3,   sym_filled = 1,color = colorarr[a],$
                     /overplot,/current)

           pq.window.SetCurrent
           pq = plot(bin_phase, bin_ycen, '1s', sym_size = 0.3,   sym_filled = 1, color = colorarr[a],$
                     /overplot,/current)

           pr.window.SetCurrent
           if keyword_set(errorbars) then begin
              pr = errorplot(bin_phase, bin_flux/(plot_norm) , bin_fluxerr/plot_norm, '1s', $
                        sym_size = 0.3,   sym_filled = 1,  color = 'black', /overplot,/current) 
           endif else begin
              pr = plot(bin_phase, bin_flux/(plot_norm) , '1s', sym_size = 0.3,   sym_filled = 1,  $
                        color = 'black', /overplot,/current) 
           endelse


           if pmapcorr eq 1 then begin
              print, 'inside pmapcorr eq 1', median( (bin_corrfluxp/plot_corrnorm) - 0.01)

               if keyword_set(errorbars) then begin
                  pr = errorplot(bin_phasep, (bin_corrfluxp/plot_corrnorm)- 0.01,$
                            bin_corrfluxerrp/ plot_corrnorm, '1s', sym_size = 0.3,   $
                            sym_filled = 1, color = colorarr[a],/overplot,/current)
               endif else begin
                  pr = plot(bin_phasep, (bin_corrfluxp/plot_corrnorm) -0.01, '1s', $
                            sym_size = 0.3,   sym_filled = 1, color = colorarr[a],/overplot,/current)
                endelse
            endif

           ps.window.SetCurrent
           ps = plot(bin_phase, bin_np, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a],$
                     /overplot,/current) 

           pt.window.SetCurrent
           pt = plot(bin_phase, bin_bkgd, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], $
                     /overplot,/current) 


        endif ; a gt stareaor = is a staring AOR

        ps.save, dirname +'binnp_phase.png'
        pt.save, dirname +'binbkg_phase.png'

        if keyword_set(selfcal) then begin
         ;  restore, strcompress(dirname + 'selfcal.sav')          
          ; print, 'test selfcal', y[0:10], 'x', x[0:10]
          ; pr.window.SetCurrent
           ;pr = errorplot(x, y-0.01,yerr,'1bo',sym_size = 0.2, sym_filled = 1, errorbar_color = 'black',color = 'black',/overplot,/current)
           ;pr = plot(x, y -0.01, '1s', sym_size = 0.3,   sym_filled = 1, color = 'black',/overplot,/current)
        endif



  endfor  ; end for all aors

 
;

  ;save, /all, filename='/Users/jkrick/irac_warm/snapshots/snaps_fit.sav'

end


function bin_data

;but how do I return all those variables cleanly?
;need to make some sort of structure.

return, 0
end

