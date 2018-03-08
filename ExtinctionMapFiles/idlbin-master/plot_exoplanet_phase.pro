pro plot_exoplanet_phase, planetname, bin_level, apradius, chname, function_fit = function_fit, set_nbins = set_nbins, single_bin = single_bin, all_plots = all_plots

  COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycenp, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr, bin_xfwhm, bin_yfwhm,  bin_corrflux_dp

  if planetname eq 'WASP-14b' or planetname eq 'HD209458' then begin
     colorarr = ['gray', 'gray','gray','gray','gray','burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
  endif else begin
 colorarr = ['burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']  
;     colorarr = ['blue', 'red','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple',  'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua' ,'blue', 'red','black','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo']
;
  endelse


;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  stareaor = planetinfo[planetname, 'stareaor']
  print, 'stareaor', stareaor
  exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
  exosystem = strmid(planetname, 0, 8 )+ ' b' ;'HD 209458 b' ;
  
  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  stareaor = planetinfo[planetname, 'stareaor']
  exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
  exosystem = strmid(planetname, 0, 8 )+ ' b' ;'HD 209458 b' ;
;exosystem = planetname
  if planetname eq 'WASP-13b' then exosystem = 'WASP-13 b'
  if planetname eq 'WASP-14b' then exosystem = 'WASP-14 b'
  if planetname eq 'WASP-15b' then exosystem = 'WASP-15 b'
  if planetname eq 'WASP-16b' then exosystem = 'WASP-16 b'
  if planetname eq 'WASP-38b' then exosystem = 'WASP-38 b'
  if planetname eq 'WASP-62b' then exosystem = 'WASP-62 b'
  if planetname eq 'WASP-52b' then exosystem = 'WASP-52 b'
  if planetname eq 'HAT-P-22' then exosystem = 'HAT-P-22 b'
  if planetname eq 'GJ1214' then exosystem = 'GJ 1214 b'
  if planetname eq '55Cnc' then exosystem = '55 Cnc e'
  if planetname eq 'HD209458' then exosystem = 'HD 209458 b'
  if planetname eq 'XO3' then exosystem = 'XO-3 b' 

  print, exosystem, 'exosystem'
  if planetname eq 'WASP-52b' then teq_p = 1315
  if planetname eq 'HD 7924 b' then begin
     inclination = 85.
     rp_rstar = 0.001
  endif
  
  if chname eq '2' then lambdaname  = '4.5'
  if chname eq '1' then lambdaname  = '3.6'
  get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                     TEQ_P=1315,TEFF_STAR=teff_star,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA=lambdaname,$
                     INCLINATION=inclination,MJD_TRANSIT=mjd_transit,P_ORBIT=p_orbit,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                     DISTANCE=distance,ECC=ecc,T14=t14,F36=f36,F45=f45,FP_FSTAR0=fp_fstar0,/verbose
  if ra lt 400 then begin       ; that means get_exoplanet_data actually found the target
;   ra_ref = ra*15.       ; comes in hours!;
;   dec_ref = dec
     utmjd_center = mjd_transit
     period = p_orbit
  endif else begin
                                ;warning these could be junk as well
;   ra_ref = planetinfo[planetname, 'ra']
;   dec_ref = planetinfo[planetname, 'dec']
     utmjd_center = planetinfo[planetname, 'utmjd_center']
     period = planetinfo[planetname, 'period']
  endelse

  if planetname eq 'WASP-14b' then begin
     ;;ok, but actually use the Wong et al. updated parameters
     p_orbit = 2.24376543       ;days
     mjd_transit = 2456034.21290 - 2400000.5
     ecc = 0.0828
     omega = 252.11             ;degrees
  endif

;---------------
  ;;read in the photometry save file
  dirname = strcompress(basedir + planetname +'/')                                                 
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_150226.sav',/remove_all) ;
  print, 'restoring ', savefilename
  restore, savefilename
  print, 'aorname', aorname(0)
  time_0 = (planethash[aorname(0),'timearr']) 
  time_0 = time_0(0)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;setup some arrays to hold binned values
  bin_corrfluxfinal = fltarr(n_elements(aorname))
  bin_phasefinal = fltarr(n_elements(aorname))
  bin_corrfluxerrfinal = fltarr(n_elements(aorname))

 ;for debugging: skip some AORs
  startaor =   0                 ;  n_elements(aorname) -29
  stopaor =   n_elements(aorname) - 1
  mean_corrfluxarr = fltarr(stopaor - startaor + 1,/nozero)
  mean_dsweetarr = mean_corrfluxarr
  stddev_corrfluxarr = mean_corrfluxarr
  stddev_dsweetarr = mean_corrfluxarr
 
;--------------------------------------------
;pick the normalization
;--------------------------------------------
  phase0 = planethash[aorname(0),'phase']
  corrflux0 = planethash[aorname(0),'corrflux_d']
  flux0 = planethash[aorname(0),'flux']
;  se = where(phase0 lt -.012 );and phase0 lt 0.51)
  se = where(phase0 gt 0.47 and phase0 lt 0.6)
  plot_corrnorm =  mean(corrflux0(se),/nan)
  plot_norm = mean(flux0(se), /nan) 
  addoffset = 0.0
  snap_addoffset = 0.0
  if n_elements(aorname) gt stareaor + 4 then begin
     ;pick a different normalization for the snapshots
     count = 0
     for a = startaor, stopaor do begin
        phase = planethash[aorname(a),'phase']
        if mean(phase,/nan) gt 0.45 and mean(phase,/nan) lt 0.51 and count eq 0 then sec_corrflux = [planethash[aorname(a),'corrflux_d']]
        if mean(phase,/nan) gt 0.45 and mean(phase,/nan) lt 0.51 and count gt 0 then sec_corrflux = [sec_corrflux, planethash[aorname(a),'corrflux_d']]
     endfor
     snap_corrnorm = mean(sec_corrflux, /nan)


     ;;except I think it is additive since I think it is a latent
     ;;issue, and not multiplicative
     ;;if planetname eq 'WASP-14b' then 
     snap_addoffset = (plot_corrnorm - snap_corrnorm) / plot_corrnorm
  endif
print, 'snap_addoffset', snap_addoffset

;--------------------------------------------


  for a = startaor,stopaor do begin
     print, '------------------------------------------------------'
     print, 'working on AOR', a, '   ', aorname(a), startaor

     ;check if I should be using pmap corr or not
     ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
     ;if 20% of the values are correctable than go with the pmap corr 
     print, '0.2nflux, ncorr, ', 0.2*n_elements([planethash[aorname(a),'flux']]), corrcount
     if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
     print, 'pmapcorr', pmapcorr

     if keyword_set(single_bin) then begin     
        ;;make one data point per AOR
        bin_level = n_elements([planethash[aorname(a),'flux']]) 
     endif

     if keyword_set(set_nbins) then begin  ; fix the number of bins instead of the number in each bin
        if a ge stareaor then begin
           ;;variable binning
           junkpar = binning_function(a, bin_level, pmapcorr,chname, /set_nbins, n_nbins = 4)
        endif else begin
           
           junkpar = binning_function(a, bin_level, pmapcorr,chname, /set_nbins, n_nbins = 95) ; wasp14 frametime specific
        endelse
     endif else begin 
        junkpar = binning_function(a, bin_level, pmapcorr,chname)
     endelse

     ;use time degraded corrfluxes
     bin_corrfluxp = bin_corrflux_dp

     ;fake the eccentricities effect on phase
     if planetname eq 'XO3' then bin_phasep = bin_phasep - 0.07
     

;------------------------------------------------------
     ;;normalize for each AOR
     ;;remove this for snapshots
     plot_corrnorm =  mean(bin_corrfluxp,/nan)
     plot_norm = mean(bin_flux, /nan) 


;------------------------------------------------------
     ;;put all the variables together and save them for later
         if a eq startaor then begin
        bin_all_phasep = bin_phasep
        bin_all_corrfluxp = bin_corrfluxp
        bin_all_fluxerrp = bin_fluxerrp
        bin_all_corrfluxerrp = bin_corrfluxerrp
        bin_all_time = bin_timearr
        bin_all_xcen = bin_xcen
        bin_all_ycen = bin_ycen
     endif else begin
        bin_all_phasep = [bin_all_phasep, bin_phasep]
        bin_all_corrfluxp = [bin_all_corrfluxp, bin_corrfluxp]
        bin_all_corrfluxerrp = [bin_all_corrfluxerrp, bin_corrfluxerrp]
        bin_all_fluxerrp = [bin_all_fluxerrp, bin_fluxerrp]
        bin_all_time = [bin_all_time,bin_timearr]
        bin_all_xcen = [bin_all_xcen, bin_xcen]
        bin_all_ycen = [bin_all_ycen, bin_ycen]
     endelse

;------------------------------------------------------
;now try plotting
;------------------------------------------------------
     setyrange =  [0.996, 1.005] ;[0.996, 1.004];
     if planetname eq 'HD158460' then setyrange = [0.997, 1.003];[0.999, 1.001]
     setynormfluxrange = [0.985, 1.005] ;[0.97, 1.005]
     setxrange = [0.40, 0.66];
     extra={ sym_size: 0.2, sym_filled: 1, xrange: setxrange, color: colorarr[a], title:planetname}

     if keyword_set(all_plots) then begin
        print, 'min, max phase', min(bin_phase), max(bin_phase)
        pp = plot(bin_phase, bin_xcen, '1s', xtitle = 'Orbital Phase', ytitle = 'X position', _extra = extra, overplot = pp)          
        pq = plot(bin_phase, bin_ycen, '1s', xtitle = 'Orbital Phase', ytitle = 'Y position', _extra = extra, overplot = pq)
        ps= plot(bin_phase, bin_npcent, '1s',xtitle = 'Orbital Phase', ytitle = 'Noise Pixel', _extra = extra, overplot = ps)
        pt = plot(bin_phase, bin_bkgd, '1s', xtitle = 'Orbital Phase', ytitle = 'Background', _extra = extra, overplot = pt)
        pxf = plot(bin_phase, bin_xfwhm, '1s', xtitle = 'Orbital Phase', ytitle = 'XFWHM', _extra = extra, $
                   yrange = [1.9, 2.3], overplot = pxf)
        pyf = plot(bin_phase, bin_yfwhm, '1s',xtitle = 'Orbital Phase', ytitle = 'YFWHM', _extra = extra, $
                   yrange = [1.9, 2.3], overplot = pyf)
        pr = errorplot(bin_phase, bin_flux/plot_norm, bin_fluxerr/plot_norm,  '1s', xtitle = 'Orbital Phase', $
                       ytitle = 'Normalized Flux', yrange = setynormfluxrange, sym_size = 0.7, sym_filled = 1, $
                       xrange = setxrange, color = colorarr[a], errorbar_color = colorarr[a], overplot = pr)
        
        
        
     endif   ;;keyword all_plots
     ;;work on the corrflux snap normalization being different than
     ;;the stare normalization

     if a gt stareaor then begin
        ;plot_corrnorm =  snap_corrnorm
        addoffset = snap_addoffset
     endif

     if a eq 3 then begin
        corrflux3 = planethash[aorname(a),'corrflux_d']
        phase3 = planethash[aorname(a),'phase']
        goodt = where(phase3 lt -.012)
        plot_corrnorm = mean(corrflux3(goodt),/nan)
     endif

     if pmapcorr eq 1 then begin
        print, 'plot_corrnorm', plot_corrnorm
;        print, 'corrflux', bin_corrfluxp/plot_corrnorm
;        print, 'bin_phasep', bin_phasep
        pu = errorplot(bin_phasep, bin_corrfluxp/plot_corrnorm + addoffset, $
                       bin_corrfluxerrp/plot_corrnorm,  '1s', sym_size = 0.7,  $
                       symbol = plotsym, sym_filled = 1,color =colorarr[a] ,xtitle = 'Orbital Phase',$
                       errorbar_color =  colorarr[a], title = planetname, ytitle = 'Pmap Corrected Flux', $
                       yrange = setyrange, xrange = setxrange, overplot =pu) 
        if keyword_set(function_fit) then begin
           pu.xshowtext = 0
           pu.position =  [0.2,0.35,0.9,0.9]
        endif
        
     endif
     
     
;-------------------------------------
     
     if a eq 0 and planetname eq 'HD209458' then begin
        print, 'plotting Lewis curve'
                                ;plot curves from Lewis
        restore, '/Users/jkrick/irac_warm/hd209458/hd209_ch2_for_jessica.sav'
                                ;need to get time in phase
        time = time / period    ; now in phase
        
                                ;need to bin
        numberarr = findgen(n_elements(time))
        h = histogram(numberarr, OMIN=om, binsize = bin_level, reverse_indices = ri)
        print, 'omin', om, 'nh', n_elements(h)
        bin_corr_flux = dblarr(n_elements(h))
        bin_time = dblarr(n_elements(h))
        c = 0
        for j = 0L, n_elements(h) - 1 do begin
           
           if (ri[j+1] gt ri[j] + 2)  then begin ;require 3 elements in the bin
              
              meanclip_jk, time[ri[ri[j]:ri[j+1]-1]], meant, sigmat
              bin_time[c] = meant
              
              meanclip_jk, corr_flux[ri[ri[j]:ri[j+1]-1]], meanc, sigmac
              bin_corr_flux[c] = meanc
              
              c = c + 1
              
           endif
        endfor
        
        bin_time = bin_time[0:c-1]
        bin_corr_flux = bin_corr_flux[0:c-1]
        pu = plot(bin_time, bin_corr_flux , '1s', sym_size =  0.3, sym_filled = 1, color ='black', overplot = pu)
     endif
      
  endfor                        ;binning for each AOR
  
  ;;print out the median x and y centers
  print, 'bin_all xcen', median(bin_all_xcen), mean(bin_all_xcen, /nan)
  print, 'bin_all ycen', median(bin_all_ycen), mean(bin_all_ycen, /nan)

  
;------------------------------------
;fit a mandel and agol function and overplot it
     if keyword_set(function_fit) then begin
        modelfilename = strcompress(dirname + planetname +'_model_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
        ;;z = function_fit_lightcurve(planetname, bin_all_phasep, bin_all_corrfluxp, bin_all_corrfluxerrp,modelfilename)
        z = function_fit_lightcurve(planetname,bin_phasep, bin_corrfluxp, bin_corrfluxerrp,modelfilename)
        print, 'bin_phasep', bin_phasep
        restore, modelfilename
        pu = plot(phase, rel_flux, color = 'Sky Blue', thick = 4, overplot = pu)
        print, 'rel_flux', rel_flux
        print, 'phase', phase


        if planetname eq 'WASP-14b' then begin
           ;;overplot the Wong et al. fit to the staring mode data
           readcol, '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/wong/ch2.dat', w_phase, w_model
           w_model = w_model + 0.0016
           pu = plot(w_phase, w_model , thick = 4, overplot = pu, color = 'green') ;+ 0.0016
           
           ;;make some residual plots
           ;;wong model
           wong_close = fltarr(n_elements(bin_all_phasep))
           for nbp = 0, n_elements(bin_all_phasep) -1 do begin
              wong_close[nbp] = closest(w_phase, bin_all_phasep(nbp))
           endfor
           wmodel_snaps = w_model(wong_close)
           resid_wong = ((bin_all_corrfluxp/plot_corrnorm)) - wmodel_snaps + 0.0018 ; .0025
           pun = plot(bin_all_phasep, resid_wong, '1s', sym_size = 0.5,   sym_filled = 1, $
                      color = 'green', yrange = [-0.004, 0.004],xtitle = 'Orbital Phase',$
                      position = [0.2, 0.15, 0.9, 0.35],/current, ytitle = 'Residuals', $
                      xrange = pu.xrange,  ytickinterval = 0.003)
        endif

           ;;ca08 model
        ca08_close = fltarr(n_elements(bin_all_phasep))
        for nbp = 0, n_elements(bin_all_phasep) -1 do begin
           ca08_close[nbp] = closest(phase, bin_all_phasep(nbp))
        endfor
        model_snaps = rel_flux(ca08_close)
        resid_ca08 = ((bin_all_corrfluxp/plot_corrnorm)) - model_snaps + 0.0018; 0.004
        ;;pun = plot(bin_all_phasep, resid_ca08, '1s', sym_size = 0.5,   sym_filled = 1, $
        ;;           color = 'sky blue',  overplot = pun)

     endif

     ;;print out some simple facts for others to play with.
     outfilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_150226.txt',/remove_all)
     openw, outlun, outfilename, /GET_LUN
     for a = 0, n_elements(aorname) - 1 do begin
        aphase = planethash[aorname(a),'phase']
        cde = (planethash[aorname(a),'corrfluxerr'])
        
        for i = 0, n_elements(aphase) - 1 do begin
           printf, outlun, (planethash[aorname(a),'bmjdarr'])(i), (planethash[aorname(a),'phase'])(i), $
                   (planethash[aorname(a),'corrflux_d'])(i), cde(i),$
                   format = '(D, F10.6, F10.6,F10.6)'
        endfor
        
        
     endfor
     free_lun, outlun
     
end


