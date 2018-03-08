pro plot_exoplanet_time, planetname, bin_level, apradius, chname, line_fit = line_fit, montecarlo = montecarlo, set_nbins = set_nbins, single_bin = single_bin, all_plots = all_plots

  COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycenp, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr, bin_xfwhm, bin_yfwhm,  bin_corrflux_dp

  if planetname eq 'WASP-14b' or planetname eq 'HD209458' then begin
     colorarr = ['gray', 'gray','gray','gray','gray','burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
  endif else begin
   
     colorarr = ['blue', 'red','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple',  'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua' ,'blue', 'red','black','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo']
;
  endelse

ra = 500
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
  ;;ra = planetinfo[planetname, 'ra']
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
   ra_ref = planetinfo[planetname, 'ra']
   dec_ref = planetinfo[planetname, 'dec']
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
  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_150226_bcdsdcorr.sav',/remove_all) ;
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
  startaor =  0                 ;  n_elements(aorname) -29
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
  se = where(phase0 gt 0.47 and phase0 lt 0.51, secount)
  if secount gt 0 then begin
     print, 'using phase to normalize'
     plot_corrnorm =  mean(corrflux0(se),/nan)
     plot_norm = mean(flux0(se), /nan) 
  endif else begin
     corrflux2 = [corrflux0,planethash[aorname(1),'corrflux_d'], planethash[aorname(stopaor),'corrflux_d'] ]
     flux2 = [flux0, planethash[aorname(1),'flux'],planethash[aorname(stopaor),'flux'] ]
     plot_corrnorm = mean(corrflux2, /nan)
     plot_norm = mean(flux2, /nan)
  endelse

  if n_elements(aorname) gt stareaor + 1 then begin
     ;pick a different normalization for the snapshots
     count = 0
     for a = startaor, stopaor do begin
        phase = planethash[aorname(a),'phase']
        if mean(phase,/nan) gt 0.45 and mean(phase,/nan) lt 0.51 and count eq 0 then sec_corrflux = [planethash[aorname(a),'corrflux_d']]
        if mean(phase,/nan) gt 0.45 and mean(phase,/nan) lt 0.51 and count gt 0 then sec_corrflux = [sec_corrflux, planethash[aorname(a),'corrflux_d']]
     endfor
     snap_corrnorm = mean(sec_corrflux, /nan)
     ;if planetname eq 'WASP-14b' then 
     snap_addoffset = (plot_corrnorm - snap_corrnorm) / plot_corrnorm

  endif


;--------------------------------------------
;setup plot for calstar
  for a = startaor,stopaor do begin
     ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
     if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
     if pmapcorr gt 0 then begin
        timea = planethash[aorname(a),'timearr']
        if a eq startaor then timeazero = timea(0)
        timea = (timea - timeazero)/60./60.
        corrfluxa = planethash[aorname(a),'corrflux_d']
        
        pu = plot(timea, corrfluxa + 0.011,  '1s', sym_size = 0.2,  $
                  symbol = plotsym, sym_filled = 1,color =colorarr[a] ,xtitle = 'Time(Hours)',$
                  title = planetname, ytitle = 'Pmap Corrected Flux', $
                  yrange = [0.98, 1.02], overplot =pu) 
     endif
  endfor


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

;------------------------------------------------------
     ;;put all the variables together and save them for later
     if a eq startaor then begin
        bin_all_phasep = bin_phasep
        bin_all_corrfluxp = bin_corrfluxp
        bin_all_fluxerrp = bin_fluxerrp
        bin_all_corrfluxerrp = bin_corrfluxerrp
        bin_all_time = bin_timearr
     endif else begin
        bin_all_phasep = [bin_all_phasep, bin_phasep]
        bin_all_corrfluxp = [bin_all_corrfluxp, bin_corrfluxp]
        bin_all_corrfluxerrp = [bin_all_corrfluxerrp, bin_corrfluxerrp]
        bin_all_fluxerrp = [bin_all_fluxerrp, bin_fluxerrp]
        bin_all_time = [bin_all_time,bin_timearr]
     endelse

;------------------------------------------------------
;now try plotting
;------------------------------------------------------
     setyrange = [0.990, 1.005]
     if planetname eq 'HD158460' then setyrange = [0.998, 1.002];[0.999, 1.001]
     setynormfluxrange = [0.98, 1.05] ;[0.97, 1.005]
     extra={ sym_size: 0.2, sym_filled: 1,color: colorarr[a], title:planetname}
     xarr = (bin_timearr - time_0)/60./60.

     if keyword_set(all_plots) then begin
        pp = plot(xarr, bin_xcen, '1s', xtitle = 'Time(Hours)', ytitle = 'X position', _extra = extra, overplot = pp)          
        pq = plot(xarr, bin_ycen, '1s', xtitle = 'Time(Hours)', ytitle = 'Y position', _extra = extra, overplot = pq)
        ps= plot(xarr, bin_npcent, '1s',xtitle = 'Time(Hours)', ytitle = 'Noise Pixel', _extra = extra, overplot = ps)
        pt = plot(xarr, bin_bkgd, '1s', xtitle = 'Time(Hours)', ytitle = 'Background', _extra = extra, overplot = pt)
        pxf = plot(xarr, bin_xfwhm, '1s', xtitle = 'Time(Hours)', ytitle = 'XFWHM', _extra = extra, $
                   yrange = [1.9, 2.3], overplot = pxf)
        pyf = plot(xarr, bin_yfwhm, '1s',xtitle = 'Time(Hours)', ytitle = 'YFWHM', _extra = extra, $
                   yrange = [1.9, 2.3], overplot = pyf)
        pr = errorplot(xarr, bin_flux/plot_norm, bin_fluxerr/plot_norm,  '1s', xtitle = 'Time(Hours)', $
                       ytitle = 'Normalized Flux', yrange = setynormfluxrange, sym_size = 0.7, sym_filled = 1, $
                       color = colorarr[a], errorbar_color = colorarr[a], overplot = pr)

        
     endif   ;;keyword all_plots

     ;;work on the corrflux snap normalization being different than
     ;;the stare normalization

     if a gt stareaor then begin
        ;plot_corrnorm = snap_corrnorm
        addoffset = snap_addoffset
     endif

     if pmapcorr eq 1 then begin
        timea = planethash[aorname(a),'timearr']
        if a eq startaor then timeazero = timea(0)
        timea = (timea - timeazero)/60./60.
        corrfluxa = planethash[aorname(a),'corrflux_d']

;        pu = plot(timea, corrfluxa + 0.01,  '1s', sym_size = 0.2,  $
;                       symbol = plotsym, sym_filled = 1,color =colorarr[a] ,xtitle = 'Time(Hours)',$
;                        title = planetname, ytitle = 'Pmap Corrected Flux', $
;                       yrange = [0.98, 1.02], overplot =pu) 
        pu = plot(xarr, bin_corrfluxp/plot_corrnorm + snap_addoffset, $
                         '1s', sym_size = 1.0,  $
                       symbol = plotsym, sym_filled = 1,color = 'black',$
                       errorbar_color =  'black', overplot =pu) 
;        pu = errorplot(xarr, bin_corrfluxp/plot_corrnorm + snap_addoffset, $
;                       bin_corrfluxerrp/plot_corrnorm,  '1s', sym_size = 0.7,  $
;                       symbol = plotsym, sym_filled = 1,color =colorarr[a] ,$
;                       errorbar_color =  colorarr[a], title = planetname,  $
;                       yrange = setyrange, overplot =pu) 
        
     endif
          
  
      
  endfor                        ;binning for each AOR
  
  ;;------------------------------------------------------

  if keyword_set(line_fit) then begin
     
     ;;reset time to hours starting at the first observation
     bin_all_time = (bin_all_time - time_0)/60./60.

     ;;need to do some sigma clipping rejection
     meanclip, bin_all_corrfluxp, meanval, sigmaval, subs = goodnum
     print, 'n', n_elements(bin_all_corrfluxp), n_elements(goodnum)
     bin_all_time = bin_all_time(goodnum)
     bin_all_corrfluxp = bin_all_corrfluxp(goodnum)
     bin_all_corrfluxerrp = bin_all_corrfluxerrp(goodnum)
     bin_all_fluxerrp = bin_all_fluxerrp(goodnum)
     
     start = [-1E-3,meanval] ;median(bin_all_corrfluxp)
     startflat = [1.0]
     result= MPFITFUN('linear',bin_all_time, bin_all_corrfluxp/plot_corrnorm, $  ;median(bin_all_corrfluxp)
                      bin_all_corrfluxerrp/plot_corrnorm, $ ;median(bin_all_corrfluxp)
                      start, perror = perror, bestnorm = bestnorm)    
        
     pu = plot(bin_all_time, result(0)*(bin_all_time) + result(1),  color = 'black', overplot = pu, linestyle = 2, thick = 3)
        
     DOF     = N_ELEMENTS(bin_all_time) - N_ELEMENTS(result) ; deg of freedom
     PCERROR = PERROR * SQRT(BESTNORM / DOF)                 ; scaled uncertainties
     print, 'linear result', result, 'pcerror', pcerror, ' perror', perror, 'bestnorm', bestnorm, 'dof', DOF, sqrt(bestnorm/dof)

     ;;t1 = text(0.2, 1.001, sigfig(result(0), 3, /scientific) + '    ' + sigfig(bestnorm/DOF, 5), /data, color = 'black', font_style = 'bold', target = pu) ;pcerror(0)
     
     ;;result2 = MPFITFUN('flatline',bin_all_time,bin_all_corrfluxp/median(bin_all_corrfluxp), bin_all_corrfluxerrp/median(bin_all_corrfluxp), $
;                        startflat, perror = perror, bestnorm = bestnorm)    
     ;;pu = plot(bin_all_time, 0*bin_all_time + result2(0),  color = 'blue', overplot = pu, thick = 3, linestyle = 2)
     ;;tt1 = text(0.2, 1.0015, '0    ' + sigfig(bestnorm/DOF, 5), /data, color = 'blue', font_style = 'bold', target = pu)
     
     
  endif

;----------------------------------------------------
;Monte Carlo
;randomize the time variables, re-measure the slopes.  
; What is the distribution of slopes measured from S such randomizations?
  
  if keyword_set(montecarlo) then begin
     MC = 1000                  ; how many implementations should I use?
     slopearr = dindgen(MC)
     for m = 0, MC - 1 do begin
;use Fisher-Yates shuffle to randomize the time array
        sbin_all_time = fisher_yates_shuffle(bin_all_time)
           
;re-fit
        l_result= MPFITFUN('linear',sbin_all_time, bin_all_corrfluxp/median(bin_all_corrfluxp), $
                           bin_all_corrfluxerrp/median(bin_all_corrfluxp), $
                           start, perror = perror, bestnorm = bestnorm,/quiet)    
        f_result= MPFITFUN('flatline',sbin_all_time,bin_all_corrfluxp/median(bin_all_corrfluxp), $
                           bin_all_corrfluxerrp/median(bin_all_corrfluxp), $
                           startflat, perror = perror, bestnorm = bestnorm,/quiet) 
           
;keep an array of slopes
        slopearr[m] = l_result(0)
     endfor
        
;now plot a distribution of the slopes
     plothist, slopearr, xhist, yhist,  bin=1E-7, /noprint,/noplot
     titlename = string(MC) + ' Monte Carlo realizations'

         
     ph = barplot(xhist, yhist,  xtitle = 'Slope', ytitle = 'Number', fill_color = 'sky_blue' ,$ ; title = titlename, $
                  xrange = [-2.5E-6, 2.5E-6], yrange = [0, 60])
     ph = plot(intarr(200) + result(0), indgen(200), linestyle = 2, thick = 2, overplot = ph)
        
        
;XXXtake this one step further, what fraction of the slopes are closer to the mean?
;fit a gaussian
     start = [0.000, 1E-6, MC]
     junkerr = fltarr(n_elements(xhist)) + 1.0
     g_result= MPFITFUN('mygauss',xhist, yhist,junkerr, start, perror = perror, bestnorm = bestnorm)    
        ;ph = plot(xhist, g_result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - g_result(0))/g_result(1))^2.), color = 'black', overplot = ph)
;        ph = plot(intarr(200) - g_result(1) , indgen(200), thick = 2, overplot = ph)

        ;; this works because the mean is zero slope
     g = where(slopearr lt -1*(result(0)) and slopearr gt result(0), gcount) 
     gm = where(xhist gt -1*(g_result(1)) and xhist lt g_result(1),  complement = outside) 
     ph = barplot(xhist(outside), yhist(outside), fill_color = 'gray', overplot = ph)

     print, 'fraction with slopes closer to zero',100*(gcount / float(MC)), gcount
  endif
  


end


