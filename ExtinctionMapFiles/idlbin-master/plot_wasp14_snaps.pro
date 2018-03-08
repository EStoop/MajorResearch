pro plot_wasp14_snaps, bin_level, time_plot = time_plot, position_plot=position_plot, unbinned_bkgplot=unbinned_bkgplot, normalize = normalize, wong_model = wong_model, secondary = secondary, transit = transit, printout = printout

  COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycenp, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr, bin_xfwhm, bin_yfwhm, bin_corrflux_dp

  colorarr = ['gray', 'gray','gray','gray','gray','burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']



  exosystem = 'WASP-14 b'
  planetname = 'WASP-14b'
;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  stareaor = planetinfo[planetname, 'stareaor']
  exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
  
  lambdaname  = '4.5'
  get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                     TEQ_P=1315,TEFF_STAR=teff_star,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA=lambdaname,$
                     INCLINATION=inclination,MJD_TRANSIT=mjd_transit,P_ORBIT=p_orbit,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                     DISTANCE=distance,ECC=ecc,T14=t14,F36=f36,F45=f45,FP_FSTAR0=fp_fstar0,/verbose



  dirname = strcompress(basedir + planetname +'/')                                                            ;+'/hybrid_pmap_nn/')
;  savefilename = strcompress(dirname + planetname +'_phot_ch2_2.25000_150226_bcdsdcorr.sav',/remove_all)
  savefilename = strcompress(dirname + planetname +'_phot_ch2_2.25000_160126.sav',/remove_all) 
  restore, savefilename

 
  ;;ok, but actually use the Wong et al. updated parameters
  p_orbit = 2.24376543          ;days
  mjd_transit = 2456034.21290 - 2400000.5
  ecc = 0.0828
  omega = 252.11                ;degrees
  


  startaor =  0;0                 ;  n_elements(aorname) -29
  stopaor =   n_elements(aorname) - 1
  time_0 = (planethash[aorname(startaor),'timearr']) 
  time_0 = time_0(0)
;  stareaor = 5

;--------------------------------------------
;pick the normalization
;--------------------------------------------
  phase0 = planethash[aorname(0),'phase']
  corrflux0 = planethash[aorname(0),'corrflux_d']
  flux0 = planethash[aorname(0),'flux']
  se = where(phase0 gt 0.47 and phase0 lt 0.51)
  plot_corrnorm =  mean(corrflux0(se),/nan)
  plot_norm = mean(flux0(se), /nan) 
  addoffset = 0.0

  if n_elements(aorname) gt stareaor + 1 then begin
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
     snap_addoffset = (plot_corrnorm - snap_corrnorm) / plot_corrnorm
  endif

;;------------------------------------------------------------------------

  if keyword_set(secondary) then begin
     ;;make a plot of a zoom in on the secondary eclipse region with
     ;;stares, snaps, and a model from blecic et al to compare
     ;;measured depths.
     
     for a = 0, stopaor do begin
        print, '------------------------------------------------------'
        print, 'working on AOR', a, '   ', aorname(a), startaor, colorarr[a]
        ;;check if I should be using pmap corr or not
        ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
        ;;if 20% of the values are correctable than go with the pmap corr 
        if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
        print, '0.2nflux, ncorr, ', 0.2*n_elements([planethash[aorname(a),'flux']]), corrcount, pmapcorr

        if a ge stareaor then begin
           ;;variable binning
           junkpar = binning_function(a, bin_level, pmapcorr,/set_nbins, n_nbins = 4*3)
           addoffset = snap_addoffset
           ss = 0.7
        endif else begin
           junkpar = binning_function(a, bin_level, pmapcorr,/set_nbins, n_nbins = 95*3)
           ss = 0.5
        endelse
     

        print, 'mean bin_corrflux_dp', mean(bin_corrflux_dp)
        psec = errorplot(bin_phasep, bin_corrflux_dp/plot_corrnorm + addoffset, bin_corrfluxerrp/plot_corrnorm, '1s', $
                    sym_size = ss, sym_filled = 1, errorbar_color =colorarr[a],$
                    color = colorarr[a], yrange = [0.990, 1.008], xrange = [0.4, 0.6], overplot = psec, $
                        xtitle = 'Orbital Phase', ytitle = 'Normalized Corrected Flux')
        ;;unwrap the negative ones
        psec = errorplot(bin_phasep + 1.0, bin_corrflux_dp/plot_corrnorm + addoffset, bin_corrfluxerrp/plot_corrnorm,$
                    '1s', sym_size = ss, sym_filled = 1, errorbar_color =colorarr[a],$
                    color = colorarr[a], overplot = psec)


        if a eq startaor then begin
        ;;line for the second secondary depth as measured in the literature

        xp = [0.47, 0.51]
        yw1 = 1. - ([0.22, 0.22] / 100.)
        yw2 = 1. - ([0.24, 0.24] / 100.)
        yb1 = 1. - ([0.224, 0.224] / 100.)
;        psec = plot(xp, yw1, color = 'gray', thick = 3, overplot = psec)
;        psec = plot(xp, yw2, color = 'gray', thick = 3, overplot = psec)
;        psec = errorplot([0.46], [0.99776],[0.00018], '1D', color = 'black', overplot = psec, $
;                         sym_filled = 1, sym_size = 1.5)
;        stext = text(.41, .99776, 'Blecic et al. 13',/data, vertical_alignment = 0.5, font_style = 'Bold')
;        psec = plot([0.4, 0.6], [1.0, 1.0], color = 'gray', thick = 2, overplot = psec)
        
        ;;try arrows instead of lines
;        a1 = arrow([0.50, 0.5], [1.0,1. - (0.1997 / 100.)] , color = 'black', /data, thick = 2, /current)
;        a1 = arrow([0.47, 0.47], [1.0,1. - (0.2249 / 100.)] , color = 'black', /data, thick = 2, /current)

        ;;an unbinnned version
        phase = planethash[aorname(a),'phase']
        corrflux = planethash[aorname(a),'corrflux_d']
        corrfluxerr = planethash[aorname(a),'corrfluxerr']
        good = where(finite(corrflux) gt 0)
        phase = phase(good)
        corrflux = corrflux(good)
        corrfluxerr = corrfluxerr(good)

        neg = where(phase lt 0)
        phase(neg) = phase(neg) + 1.0

;        psec = plot(phase, corrflux/normcorr + corroffset, '1s', sym_size = 0.5, sym_filled = 1, color = colorarr[a], $
;                 yrange = [0.980, 1.02], xrange = [0.4, 0.6], overplot = psec)
     endif

     endfor  ; for each AOR


  endif


;;------------------------------------------------------------------------

  if keyword_set(transit) then begin
     ;;make a plot of a zoom in on the transit  region with
     ;;stares, snaps, and a model from blecic et al to compare
     ;;measured depths.
     
     for a = 0, stopaor do begin
        print, '------------------------------------------------------'
        print, 'working on AOR', a, '   ', aorname(a), startaor, colorarr[a]
        ;;check if I should be using pmap corr or not
        ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
        ;;if 20% of the values are correctable than go with the pmap corr 
        if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
        print, '0.2nflux, ncorr, ', 0.2*n_elements([planethash[aorname(a),'flux']]), corrcount, pmapcorr

        if a ge stareaor then begin
           ;;variable binning
           junkpar = binning_function(a, bin_level, pmapcorr,/set_nbins, n_nbins = 4*3)
           addoffset = snap_addoffset
           ss = 0.7
        endif else begin
           ss = 0.5
           junkpar = binning_function(a, bin_level, pmapcorr,/set_nbins, n_nbins = 95*3)
        endelse

        ptran = errorplot(bin_phasep, bin_corrflux_dp/plot_corrnorm + addoffset, bin_corrfluxerrp/plot_corrnorm, '1s', $
                    sym_size = ss, sym_filled = 1, errorbar_color =colorarr[a],$
                    color = colorarr[a], yrange = [0.990, 1.008], xrange = [-0.1, 0.1], overplot = ptran, $
                        xtitle = 'Orbital Phase', ytitle = 'Normalized Corrected Flux')
        ;;unwrap the negative ones
        ptran = errorplot(bin_phasep + 1.0, bin_corrflux_dp/plot_corrnorm + addoffset, bin_corrfluxerrp/plot_corrnorm,$
                    '1s', sym_size = ss, sym_filled = 1, errorbar_color =colorarr[a],$
                    color = colorarr[a], overplot = ptran)


        if a eq stopaor then begin
        ;;line for the second secondary depth as measured in the literature

           xp = [-0.015, 0.02]
           yw1 = 1. - [0.00887, 0.00887] 
           yj1 = 1. - [0.0102, 0.0102] 
;           ptran = plot(xp, yw1, color = 'gray', thick = 3, overplot = ptran)
;           ;ptran = plot(xp, yj1, color = 'black', thick = 3, overplot = ptran, linestyle = 2)
;           ptran = errorplot([-0.02], [0.9898], [.00119], '1D', color = 'black',  overplot = ptran, $
;                            sym_filled = 1.0, sym_size = 1.5)
;           ttext = text(-.07, .9898, 'Joshi et al. 09',/data, vertical_alignment = 0.5, font_style = 'Bold')
        endif

     endfor                     ; for each AOR


  endif

;;------------------------------------------------------------------------

  if keyword_set(time_plot) then begin
     nphasesub = 63. ; 763.0 number of phases between MJD and the snaps
     for a = 0, stopaor do begin
        print, '------------------------------------------------------'
        print, 'working on AOR', a, '   ', aorname(a), startaor, colorarr[a]
                                ;check if I should be using pmap corr or not
        ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
                                ;if 20% of the values are correctable than go with the pmap corr 
        print, '0.2nflux, ncorr, ', 0.2*n_elements([planethash[aorname(a),'flux']]), corrcount
        if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
        print, 'pmapcorr', pmapcorr
        
        if a gt stareaor then begin
           ;;variable binning
           junkpar = binning_function(a, bin_level, pmapcorr,/set_nbins, n_nbins = 4)
           addoffset = snap_addoffset
           ss = 0.7
        endif else begin
           ss = 0.5
           junkpar = binning_function(a, bin_level, pmapcorr,/set_nbins, n_nbins = 95)
           addoffset = 0.0
        endelse
        ;;use the time degraded corrfluxes
        bin_corrfluxp = bin_corrflux_dp
        
        ;;save some arrays for later
        CASE 1 of 
           ;;take the  stare AORs and make one array of phase and corrflux
           ;;and corrfluxerr
           
           (a eq 0):BEGIN
              stare_phasearr = bin_phasep
              stare_corrfluxarr = bin_corrfluxp
              stare_corrfluxerrarr = bin_corrfluxerrp
              stare_colorarr = 'gray'
           END
           
           (a gt 0) and (a lt stareaor): BEGIN
              stare_phasearr = [stare_phasearr,bin_phasep]
              stare_corrfluxarr = [stare_corrfluxarr,bin_corrfluxp]
              stare_corrfluxerrarr = [stare_corrfluxerrarr,bin_corrfluxerrp]
              stare_colorarr = [stare_colorarr, 'gray']
           end
           
           (a eq stareaor):BEGIN
              
              ;;make stare_phasearr all positive
              stare_phasearr = stare_phasearr + 1.0
              
              ;;then duplicate that array with phasexN          
              ;;there are 170 periods in 1 years
              n_phases = 170
              ext_phasearr = fltarr(n_elements(stare_phasearr) *n_phases)
              ext_corrfluxarr = ext_phasearr
              ext_corrfluxerrarr = ext_phasearr
              for i = 1, n_phases  do begin
                 ext_phasearr[(i - 1)*n_elements(stare_phasearr):(i*n_elements(stare_phasearr)) - 1] = stare_phasearr + (i-1)
                 ext_corrfluxarr[(i - 1)*n_elements(stare_phasearr):(i*n_elements(stare_phasearr)) - 1] = stare_corrfluxarr
                 ext_corrfluxerrarr[(i - 1)*n_elements(stare_phasearr):(i*n_elements(stare_phasearr)) - 1] = stare_corrfluxerrarr
              endfor
              ;stare_norm = mean(planethash[aorname(0),'corrflux'],/nan)
              pt = plot(ext_phasearr, ext_corrfluxarr/plot_corrnorm + addoffset, '1s', sym_size = 0.2,   $
                        sym_filled = 1, color = 'grey', yrange = [0.989, 1.005], xrange = [0,10]);[157,169]
              
              
              ;;need a new phasing for the snapshots
              ;;phase =( bmjd_dist / period )
              snap_phasearr = ( (bin_bmjdarr-mjd_transit) / p_orbit )
              snap_corrfluxarr = bin_corrfluxp
              snap_corrfluxerrarr = bin_corrfluxerrp
              snap_colorarr = colorarr(a)
              
              ;;overplot the snaps
;              snap_norm =  stare_norm
              pt = plot(snap_phasearr - nphasesub, snap_corrfluxarr/plot_corrnorm + addoffset, '1s', $
                        sym_size = 0.5,   sym_filled = 1, color = colorarr[a], overplot = pt,$
                        xtitle = 'Number of Periods', ytitle = 'Normalized Corrected Flux')
              
              
           end
           
           (a gt stareaor):BEGIN
              
              snap_phasearr = [snap_phasearr,( (bin_bmjdarr-mjd_transit) / p_orbit )]
                            
              pt = plot(( (bin_bmjdarr-mjd_transit) / p_orbit )- nphasesub, bin_corrfluxp/plot_corrnorm + addoffset,$
                        '1s', sym_size = 0.5,   sym_filled = 1, color = colorarr[a], overplot = pt)
              
              snap_corrfluxarr = [snap_corrfluxarr,bin_corrfluxp]
              snap_corrfluxerrarr = [snap_corrfluxerrarr,bin_corrfluxerrp]
              snap_colorarr = [snap_colorarr, colorarr(a)]
           end
        endcase
     endfor
     pt.save,'/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/snaps_paper/allfluxes_period_ch2_zoom1_3_7.eps'
  endif                         ; time_plot
  
;;-------------------------------------------------
  
  if keyword_set(printout) then begin
  
     ;;print out some simple facts for others to play with.
     openw, outlun, '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/wasp14_snap_latenterr_150226.txt', /GET_LUN
     for a = 5, n_elements(aorname) - 1 do begin
        aphase = planethash[aorname(a),'phase']
        cde = (planethash[aorname(a),'corrfluxerr'])

        if (aphase(0) lt 0.03) and (aphase(0) gt -0.03) then begin
           print, 'adding error', a, ' ',  aphase(0), colorarr(a)
           cde = cde + 0.001  ; add error due to likely latent where we know the latent as an effect
        endif

        if (aphase(0) lt -.455) then begin; and ((planethash[aorname(a),'corrflux_d'])(0) gt .0030) then begin
           print,' adding error', a, ' ', aphase(0), colorarr(a), (planethash[aorname(a),'corrflux_d'])(0) 
           cde = cde + 0.001  ; add error due to likely latent where we know the latent as an effect
        endif

           for i = 0, n_elements(aphase) - 1 do begin
              printf, outlun, (planethash[aorname(a),'bmjdarr'])(i), (planethash[aorname(a),'phase'])(i), $
                      (planethash[aorname(a),'corrflux_d'])(i), cde(i),$
                      format = '(D, F10.6, F10.6,F10.6)'
;              printf, outlun, (planethash[aorname(a),'bmjdarr'])(i), (planethash[aorname(a),'phase'])(i), $
;                      (planethash[aorname(a),'corrflux_d'])(i), (planethash[aorname(a),'corrfluxerr'])(i),$
;                      format = '(D, F10.6, F10.6,F10.6)'
           endfor

  
     endfor
     free_lun, outlun
;----
     openw, outlun, '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/wasp14_stare_150226.txt', /GET_LUN
     for a = 0, 4 do begin
        for i = 0, n_elements(planethash[aorname(a),'phase']) - 1 do begin
           printf, outlun, (planethash[aorname(a),'bmjdarr'])(i), (planethash[aorname(a),'phase'])(i), $
                   (planethash[aorname(a),'corrflux_d'])(i), (planethash[aorname(a),'corrfluxerr'])(i),$
                   format = '(D, F10.5, F10.6,F10.6)'
        endfor
  
     endfor
     free_lun, outlun
  
endif

;;-------------------------------------------------
  
  
  if keyword_set(position_plot) then begin
     fits_read, '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_121120.fits', pmapdata, pmapheader
     c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, axis_style = 0, aspect_ratio = 1, xrange = [0,500], yrange = [0,500])
     
     for a = stareaor ,  stopaor  do begin
        xcen500 = 500.* ((planethash[aorname(a),'xcen']) - 14.5)
        ycen500 = 500.* ((planethash[aorname(a),'ycen']) - 14.5)
        an = plot(xcen500, ycen500, '1s', sym_size = 0.05,   sym_filled = 1, color = colorarr[a],/overplot)
     endfor
;     xs = 500.*(xsweet - 14.5)
;     ys = 500.*(ysweet - 14.5)
;   an = text(xs, ys, 'X', /data, alignment = 0.5, vertical_alignment = 0.5)
     an.save, dirname+'position_ch2.png'
  endif 

;;-------------------------------------------------


  if keyword_set(unbinned_bkgplot) then begin

     for a = stareaor,n_elements(aorname) -1 do begin ; 0, n_elements(aorname) - 1 do begin
        if a eq stareaor then begin
           ay = plot( planethash[aorname(a),'phase'] , planethash[aorname(a),'bkgd'],'1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], ytitle = 'Bkgd', xtitle = 'Phase') ;,title = planetname
        endif else begin
           ay = plot( planethash[aorname(a),'phase'] , planethash[aorname(a),'bkgd'],'1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/overplot)
        endelse
        
     endfor
     ay.save, dirname +'bkgd_time_ch2.png'
  endif
;;-------------------------------------------------

  if keyword_set(normalize) then begin
     normval = fltarr(n_elements(aorname))
     time_0 = (planethash[aorname(0),'timearr']) 
     time_0 = time_0(0)

     for a = 0, n_elements(aorname) - 1 do begin


        ;;----------------------
        ;;make a correction for flux degradation over time
        
        ;;0.05% per year or 0.0005*(time in
        ;;years since start of observation)
        deltatime = planethash[aorname(a),'timearr'] - time_0 ; in seconds
        deltatime = deltatime / 60./60./24./365.
        degrade = .0005*deltatime
        planethash[aorname(a),'corrflux'] =  planethash[aorname(a),'corrflux'] + (planethash[aorname(a),'corrflux'] * degrade)
        
        ;;----------------------


        ;;try with binning
        ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
        if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
        if a ge stareaor then begin
           ;;variable binning
           junkpar = binning_function(a, bin_level, pmapcorr, '2', /set_nbins, n_nbins = 4)
        endif else begin
           junkpar = binning_function(a, bin_level, pmapcorr,'2', /set_nbins, n_nbins = 95)
        endelse

        normval[a] =mean(planethash[aorname(a),'corrflux'],/nan) ; mean(bin_corrfluxp);
;        x = planethash[aorname(a),'phase']
;        y =  planethash[aorname(a),'corrflux']
;        g = where(finite(y) gt 0 )
;        x = x(g)
;        y =  y(g)
;        pnn = plot(x, y, '1s', sym_size=0.2, sym_filled =1, overplot = pnn

        if a lt stareaor then c = 'grey' else c = 'blue'
       
        pnn = plot(bin_phasep, bin_corrfluxp, '1s', sym_size=0.2, sym_filled =1, overplot = pnn, yrange = [0.056, 0.06])
     endfor
     print, 'normval', normval
     plothist, normval, xhist, yhist, bin = 0.0001, /noplot
;     print, 'xhist', xhist
     pn = barplot(xhist, yhist, fill_color = 'green', xrange = [0.056, 0.06]) 

endif

if keyword_set(wong_model) then begin
   readcol, '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/wong/ch2.dat', phase, model
   pw = plot(phase, model, xrange = [-0.7, 0.7])
   pw = plot(findgen(10)/ 5 - 1.0, fltarr(10) + 1.00, overplot = pw)
   pw = plot(findgen(10)/ 5 - 1.0, fltarr(10) + .9985, overplot = pw)
   pw = plot(findgen(10)/ 5 - 1.0, fltarr(10) + 1.0005, overplot = pw)

endif


end



;     if a eq startaor then begin ; for the first AOR
                                ;set the normalization values from the
                                ;medians of the first AOR...is close enough
;        plot_corrnorm = median(bin_corrfluxp)
;        plot_norm = median(bin_flux)
        
;        if pmapcorr eq 1 then begin
;           print, 'inside pmapcorr eq 1', median(bin_corrfluxp), median(bin_flux)
           
;           pu = errorplot(((bin_timearr - time_0)/60./60.) , (bin_corrfluxp/plot_corrnorm) + corroffset,$ ;- 5.E3)/24.
;                          bin_corrfluxerrp/plot_corrnorm, '1s' , sym_size = 0.2,   sym_filled = 1 ,$
;                          color = colorarr[a] ,xtitle = 'Time (Hours)', errorbar_color =  colorarr[a], $
;                          title = planetname, ytitle = 'Pmap Corrected Flux') ;
           
;           endif  ; if pmapcorr eq 1
;     endif        ; if a = startaor
     
;     if (a gt 0) and (a le stareaor) then begin
;        print, 'inside a gt 0 a le stareaor', a
        
;        if pmapcorr eq 1 then begin
;           pu =  errorplot(((bin_timearr - time_0)/60./60.) , (bin_corrfluxp/plot_corrnorm) + corroffset, $ ;- 5.5E3)/24.
;                           bin_corrfluxerrp/plot_corrnorm, '1s' , sym_size = 0.2,   sym_filled = 1, $
;                          color = colorarr[a], errorbar_color =  colorarr[a], overplot = pu, /current)
;        endif
        
;     endif 
     

;     if a gt stareaor then begin
;        print, 'inside a gt stareaor', a
        
;        if pmapcorr eq 1 then begin
;           pu =  errorplot(((bin_timearr - time_0)/60./60.) , (bin_corrfluxp/plot_corrnorm) + corroffset, $ ;- 5.5E3)/24.
;                           bin_corrfluxerrp/plot_corrnorm, '1s' , sym_size = 0.2,   sym_filled = 1, $
;                           color = colorarr[a], errorbar_color =  colorarr[a], overplot = pu)
           
;        endif
        
;     endif
