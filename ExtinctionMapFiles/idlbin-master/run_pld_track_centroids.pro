pro run_pld_track_centroids

  

  ;;eventually want this to be a function called by track_centroids
  
  savenames = [ '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_01.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_02.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_03.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_04.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_05.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_06.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_07.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_08.sav']

  newsavenames = [ '/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_01_pld.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_02_pld.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_03_pld.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_04_pld.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_05_pld.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_06_pld.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_07_pld.sav','/Users/jkrick/Library/Mobile Documents/com~apple~CloudDocs/centroids_save/track_centroids_pixval_08_pld.sav']

  extra={xthick: 2, ythick:2, margin: 0.2, sym_size: 0.2,   sym_filled: 1, xstyle:1}
  
  starts = 0
  stops =  0 ; n_elements(savenames) - 1
  totalaorcount = 0
  for s = starts, stops do begin
     TIC
     print, 'restoring ', newsavenames(s)
     restore, newsavenames(s)
     aorname = planethash.keys()
     
     
     ;;planetinfo = create_planetinfo()
     ;;aorname= planetinfo[planetname, 'aorname_ch2']
     ;; utmjd_center = planetinfo[planetname, 'utmjd_center']
     ;;period = planetinfo[planetname, 'period']
     ;;basedir = planetinfo[planetname, 'basedir']    
     
     
     for a = 0, n_elements(aorname) - 1 do begin
        pixgrid = planethash[aorname(a),'piarr']
        flux = planethash[aorname(a),'flux']
        sigma_flux = planethash[aorname(a),'fluxerr']
        bmjd = planethash[aorname(a),'bmjdarr']
        xcen = planethash[aorname(a), 'xcen']
        ycen = planethash[aorname(a), 'ycen']
        min_dur = planethash[aorname(a), 'min_dur']
        starname = planethash[aorname(a), 'starname']
        nparr = planethash[aorname(a), 'npcentroids']
        ;print, 'min_dur', min_dur
        if min_dur[0] gt 120. then begin  ;; don't want to work on the short ones
           print, 'correcting AOR ',a, ' ', aorname(a)
           meanpix = mean(pixgrid,/nan) 
           bad = where(finite(pixgrid) lt 1,nbad)
           
           if (meanpix ne 0.) and (nbad lt 1) then begin ; if pixgrid is zeros, PLD will crash
              help, pixgrid
              help, flux
              help, sigma_flux
              corrected_flux = PLD_CORRECT_WITHFIT(pixgrid,flux,sigma_flux,CORR_UNC=corr_unc)
            
              ;;quick plot of the results
              meanclip, corrected_flux, meancorr, meansigmacorr
              ;;p2 = errorplot(bmjd, corrected_flux/meancorr, corr_unc/meancorr, '1s', sym_size = 0.4, sym_filled = 1, $
              ;;               ytitle = "PLD corrected flux", xtitle = 'BMJD', yrange =[0.9, 1.1], title = aorname(a))
              
              ;; binning
              numberarr = findgen(n_elements(bmjd))
              h = histogram(numberarr, OMIN=om, binsize = 63, reverse_indices = ri)
              
              bin_pldflux=dblarr(n_elements(h))
              bin_bmjd = bin_pldflux
              bin_corrunc = bin_pldflux
              bin_xcen = bin_pldflux
              bin_ycen = bin_pldflux
              bin_flux = bin_pldflux
              bin_np = bin_pldflux
             
              c = 0
              for j = 0L, n_elements(h) - 1 do begin
                 
                 ;;get rid of the bins with no values and low numbers, meaning
                 ;;low overlap
                 if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
                    
                    meanclip, corrected_flux[ri[ri[j]:ri[j+1]-1]], meancorrflux, sigmacorrflux
                    bin_pldflux[c] = meancorrflux
                    
                    meanclip, flux[ri[ri[j]:ri[j+1]-1]], meanflux, sigmaflux
                    bin_flux[c] = meanflux

                    idataerr = corr_unc[ri[ri[j]:ri[j+1]-1]]
                    bin_corrunc[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
                    
                    meanbmjdarr = mean( bmjd[ri[ri[j]:ri[j+1]-1]],/nan)
                    bin_bmjd[c]= meanbmjdarr
                    
                    meanxcen = mean( xcen[ri[ri[j]:ri[j+1]-1]],/nan)
                    bin_xcen[c]= meanxcen
                    
                    meanycen = mean( ycen[ri[ri[j]:ri[j+1]-1]],/nan)
                    bin_ycen[c]= meanycen
                    
                    meannp = mean( nparr[ri[ri[j]:ri[j+1]-1]],/nan)
                    bin_np[c]= meannp

                    c = c + 1
                 endif
              endfor
              bin_pldflux = bin_pldflux[0:c-1]
              bin_corrunc = bin_corrunc[0:c-1]
              bin_bmjd = bin_bmjd[0:c-1]
              bin_xcen = bin_xcen[0:c-1]
              bin_ycen = bin_ycen[0:c-1]
              bin_flux = bin_flux[0:c-1]
              bin_np = bin_np[0:c-1]
              
              normfactor = median(bin_pldflux)
              plot_norm = median(bin_flux)
              

              ;;set up some plots to be saved and not displayed
              plotx = (bin_bmjd - bin_bmjd[0]) * 24.
              pp = plot(plotx, bin_xcen, '1s',  title = aorname(a) + '  ' + starname, $
                         ytitle = 'X position', position =[0.2, 0.75, 0.9, 0.91], ytickinterval = 0.1, $
                        xshowtext = 0, ytickformat = '(F10.2)', dimensions = [600, 900], _extra = extra, yminor = 0,$
                        buffer = 1)
                        
              
              pp= plot(plotx, bin_ycen, '1s',  $
                       ytitle = 'Y position',  position = [0.2, 0.59, 0.9, 0.75],/current,  ytickinterval = 0.1,$
                       xshowtext = 0,ytickformat = '(F10.2)', _extra = extra, yminor = 0, buffer = 1)
              
              pp= plot(plotx, bin_np, '1s', $
                       ytitle = 'Noise Pixel',  position = [0.2, 0.42, 0.9, 0.58] , /current, $
                       xshowtext = 0,ytickformat = '(F10.1)', _extra = extra, ytickinterval = 1.0, yminor = 0, yrange = [4.0, 8.0],$
                       buffer = 1) ;,$ title = planetname,, ymajor = 4;xrange = setxrange) position =  [0.2, 0.50, 0.9, 0.63]
              
              pp = plot(plotx, bin_flux/plot_norm, '1s',  ytitle = 'Raw Flux',xshowtext = 0,$
                        position =[0.2, 0.25, 0.9, 0.41], /current, _extra = extra, $;ytickinterval = 0.002, $
                        yminor = 0, buffer = 1)


              pp = plot(plotx, bin_pldflux/normfactor, '1s',$
                             ytitle = "PLD Flux", xtitle = 'Time in hours',/current,_extra = extra, $
                             position = [0.2, 0.08, 0.9, 0.24],buffer = 1 );, ytickinterval = 0.002 )
                             
              ;;pp = errorplot(plotx, bin_pldflux/normfactor, bin_corrunc/normfactor,'1s',$
              ;;               ytitle = "PLD Flux", xtitle = 'Time in hours',/current,_extra = extra, $
              ;;               position = [0.2, 0.08, 0.9, 0.24],buffer = 1 ) ;, ytickinterval = 0.002 )
              
              plotname = strcompress('/Users/jkrick/external/irac_warm/trending/plots/multi_' + string(aorname(a)) + '.png',/remove_all)
              print, plotname
              pp.save,plotname

              
           endif else begin
              ;;PLD won't work on these because the SVDC will
              ;;barf at non-finite values so set the corrected values
              ;;to be non-finite
              corrected_flux = flux / alog10(-1)
           endelse

           ;;for the pre-AORs set their corrected fluxes to zeros
           corrected_flux = flux - flux  
        endif
        ;;now I need to get these pld corrected fluxes into the database
        ;;alter the planethash
        planethash[aorname(a),'pld_flux'] = corrected_flux

     endfor
     
     save, planethash, filename=newsavenames(s)

     
  endfor
  


  TOC


  
end
