pro run_pld_correct_withfit, planetname, apradius, chname

  
  ;specific to HD75289 for now
  ;restore, '/Users/jkrick/external/irac_warm/HD75289/HD75289_phot_ch2_2.25000_160126_pi.sav'
  
  ;aorname = ['r52584192', 'r52616448', 'r52617728', 'r52619008', 'r52620288', 'r52621568',  'r52608256', 'r52616704', 'r52617984', 'r52619264', 'r52620544', 'r52621824', 'r52615680', 'r52616960', 'r52618240', 'r52619520', 'r52620800', 'r52622080', 'r52615936', 'r52617216', 'r52618496', 'r52619776', 'r52621056', 'r52622336', 'r52616192', 'r52617472', 'r52618752', 'r52620032', 'r52621312', 'r52622592']
 
 ; utmjd_center = 57398.541D
 ; period = 3.5091651
 ;restore, '/Users/jkrick/external/irac_warm/stauffer/EPIC204117263/EPIC204117263_phot_ch2_2.25000_160126.sav'
 ; aorname = [ 'r63203840'] 
 ; utmjd_center = double(57915.3073228)
 ; period =54.

  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2']
  utmjd_center = planetinfo[planetname, 'utmjd_center']
  period = planetinfo[planetname, 'period']
  basedir = planetinfo[planetname, 'basedir']

  dirname = strcompress(basedir + planetname +'/')                ;;+'/hybrid_pmap_nn/')
  if chname eq '2' then  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_160126.sav',/remove_all)
  if chname eq '1' then  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_150722.sav',/remove_all) 
  restore, savefilename

  
  phasearr = fltarr(n_elements(aorname))
  fluxerrarr=phasearr
  fluxarr=phasearr
  corrfluxerrarr=phasearr
  corrfluxarr=phasearr


  
  
  for a = 0,  n_elements(aorname) - 1 do begin
     pixgrid = planethash[aorname(a),'pixvals']
     flux = planethash[aorname(a),'flux']
     sigma_flux = planethash[aorname(a),'fluxerr']
     bmjd = planethash[aorname(a),'bmjdarr']
     xcen = planethash[aorname(a), 'xcen']
     ycen = planethash[aorname(a), 'ycen']
     
     corrected_flux = PLD_CORRECT_WITHFIT(pixgrid,flux,sigma_flux,CORR_UNC=corr_unc)

     ;;--------------------
     ;;set up for some testing plots
     bmjd_dist =bmjd - utmjd_center ; how many UTC away from the transit center
     
     phase =( bmjd_dist / period )- fix(bmjd_dist/period)
     phase = phase + (phase lt -0.5 and phase ge -1.0)
     phase = phase- (phase gt 0.5 and phase le 1.0)
     
     ;;p1 = plot(phase, flux, '1s', xtitle = 'phase', ytitle = 'Flux',sym_size= 0.2,   sym_filled= 1, overplot = p1, $
     ;;          xrange = [-0.5, 0.5], yrange = [1.50, 1.58])
     ;;p1 = plot(phase, corrected_flux, '1s', sym_size= 0.2,   sym_filled= 1, color = 'red', overplot = p1)


     ;;make some meanclip arrays
     meanclip, phase, meanphase, sigmaphase
     meanclip, flux, meanflux, sigmaflux
     meanclip, sigma_flux, meansigma_flux, sigmasigma
     meanclip, corrected_flux, meancorr, meansigmacorr
     meanclip, corr_unc, meancorr_unc, sigmacorr_unc
     phasearr[a] = meanphase
     fluxerrarr[a] = meansigma_flux
     fluxarr[a] = meanflux
     corrfluxerrarr[a] = meancorr_unc
     corrfluxarr[a] = meancorr
  endfor

  ;;now plot
  ;;p2 = errorplot(phasearr, fluxarr, fluxerrarr, '1s',   sym_filled= 1, overplot = p1, $
  ;;               xrange = [-0.5, 0.5], yrange = [1.52, 1.56], title = 'HD75289 PLD', xtitle = 'Phase', ytitle = 'Flux')
  ;;p2 = errorplot(phasearr, corrfluxarr, corrfluxerrarr, '1s', color = 'red',errorbar_color = 'red',sym_filled= 1, overplot = p2)

  ;;help, corrfluxarr
  ;;print, corrfluxarr[0:23]

  
  p2 = errorplot(bmjd, corrected_flux/meancorr, corr_unc/meancorr, '1s', sym_size = 0.4, sym_filled = 1, $
                 ytitle = "PLD corrected flux", xtitle = 'BMJD', yrange =[0.9, 1.1], title = planetname)

  ;;ok try binning
   ;; binning
  numberarr = findgen(n_elements(bmjd))
  ;;if keyword_set(set_nbins) then begin
  ;;   h = histogram(numberarr, OMIN=om, nbins = n_nbins, reverse_indices = ri)
  ;;endif else begin
  h = histogram(numberarr, OMIN=om, binsize = 63, reverse_indices = ri)
 ;; endelse
  bin_pldflux=dblarr(n_elements(h))
  bin_bmjd = bin_pldflux
  bin_corrunc = bin_pldflux
  bin_xcen = bin_pldflux
  bin_ycen = bin_pldflux
  
  c = 0
  for j = 0L, n_elements(h) - 1 do begin
     
     ;;get rid of the bins with no values and low numbers, meaning
     ;;low overlap
     if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin

        meanclip, corrected_flux[ri[ri[j]:ri[j+1]-1]], meancorrflux, sigmacorrflux
        bin_pldflux[c] = meancorrflux
         
        idataerr = corr_unc[ri[ri[j]:ri[j+1]-1]]
        bin_corrunc[c] =   sqrt(total(idataerr^2))/ (n_elements(idataerr))
        
        meanbmjdarr = mean( bmjd[ri[ri[j]:ri[j+1]-1]],/nan)
        bin_bmjd[c]= meanbmjdarr
        
        meanxcen = mean( xcen[ri[ri[j]:ri[j+1]-1]],/nan)
        bin_xcen[c]= meanxcen
        
        meanycen = mean( ycen[ri[ri[j]:ri[j+1]-1]],/nan)
        bin_ycen[c]= meanycen

        c = c + 1
     endif
  endfor
   bin_pldflux = bin_pldflux[0:c-1]
   bin_corrunc = bin_corrunc[0:c-1]
   bin_bmjd = bin_bmjd[0:c-1]
   bin_xcen = bin_xcen[0:c-1]
   bin_ycen = bin_ycen[0:c-1]
   
   normfactor = median(bin_pldflux)

   p3 = errorplot(bin_bmjd, bin_pldflux/normfactor, bin_corrunc/normfactor,'1s', sym_size = 0.4, sym_filled = 1, $
                 ytitle = "PLD corrected binned flux", xtitle = 'BMJD', yrange =[0.97, 1.05], title = planetname)


   ;;and finally write out a text file
   textname = strcompress(dirname + planetname +'_data.txt',/remove_all)
   openw, unit, textname,/get_lun

   for countprint = 0, n_elements(bin_bmjd) - 1 do begin
      printf, unit,bin_bmjd(countprint),  bin_xcen(countprint), bin_ycen(countprint), bin_pldflux(countprint), bin_corrunc(countprint)
   endfor
   free_lun, unit
   
end
