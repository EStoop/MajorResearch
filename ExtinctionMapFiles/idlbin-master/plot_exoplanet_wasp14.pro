pro plot_exoplanet, planetname, bin_level, apradius, chname, phaseplot = phaseplot, sclkplot = sclkplot, timeplot = timeplot, Dsweet = dsweet, selfcal=selfcal, centerpixplot = centerpixplot, errorbars = errorbars, unbinned_xplot = unbinned_xplot, unbinned_yplot = unbinned_yplot, unbinned_fluxplot=unbinned_fluxplot, unbinned_npplot = unbinned_npplot, unbinned_bkgplot = unbinned_bkgplot, Position_plot = position_plot, DN_plot = DN_plot, unbinned_xyfwhm = unbinned_xyfwhm, line_fit = line_fit, montecarlo = montecarlo, function_fit = function_fit, set_nbins = set_nbins, all_plots = all_plots
;example call plot_exoplanet, 'wasp15', 2*63L

  COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycenp, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr, bin_xfwhm, bin_yfwhm,  bin_corrflux_dp

  if planetname eq 'WASP-14b' then begin
     colorarr = ['gray', 'gray','gray','gray','gray','burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
  endif else begin
   

 colorarr = ['black','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia'] ;'burlywood'
;     colorarr = ['black', 'red','blue','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple',  'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua' ,'blue', 'red','black','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua' ,'blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua']
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
  if planetname eq 'WASP-80b' then exosystem = 'WASP-80 b'
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
  if planetname eq 'simul_XO3' then exosystem = 'XO-3 b' 

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
;  if ra lt 400 then begin       ; that means get_exoplanet_data actually found the target
;   ra_ref = ra*15.       ; comes in hours!;
;   dec_ref = dec
;     utmjd_center = mjd_transit
;     period = p_orbit
;  endif else begin
     
                                ;warning these could be junk as well
;   ra_ref = planetinfo[planetname, 'ra']
;   dec_ref = planetinfo[planetname, 'dec']
     utmjd_center = planetinfo[planetname, 'utmjd_center']
     period = planetinfo[planetname, 'period']
;  endelse

  if planetname eq 'WASP-14b' then begin
     ;;ok, but actually use the Wong et al. updated parameters
     p_orbit = 2.24376543       ;days
     mjd_transit = 2456034.21290 - 2400000.5
     ecc = 0.0828
     omega = 252.11             ;degrees
  endif


print, 'using UTMJD', utmjd_center
;---------------
  
  dirname = strcompress(basedir + planetname +'/')                                                            ;+'/hybrid_pmap_nn/')
;  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'.sav',/remove_all) ;
  if chname eq '2' then  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_150723_newtime.sav',/remove_all) ;'_150226_bcdnosdcorr.sav'
  if chname eq '1' then  savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_150722_newtime.sav',/remove_all) ;140716
  print, 'restoring ', savefilename
  restore, savefilename
  print, 'aorname', aorname(0)
  time_0 = (planethash[aorname(0),'timearr']) 
  time_0 = time_0(0)
;  plothist, planethash[aorname(0),'npcentroids'], xhist, yhist, /noplot, bin = 0.01
;  testplot = plot(xhist, yhist,  xtitle = 'NP', ytitle = 'Number', thick =3, color = 'blue')

  ;print, 'testing phase', (planethash[aorname(0),'phase'] )[0:10], (planethash[aorname(0),'phase'] )[600:610]

 

;for debugging: skip some AORs
  startaor =  0                 ;  n_elements(aorname) -29
  stopaor =    n_elements(aorname) - 1


  xsweet = 15.12
  ysweet = 15.00
 
;-------------------------------------------------------
;make a bunch of plots before binning
;-------------------------------------------------------

;this one overlays the positions on the pmap contours
  if keyword_set(position_plot) then begin
    ;; if chname eq '2' then fits_read, '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_121120.fits', pmapdata, pmapheader
    ;; if chname eq '1' then fits_read, '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch1_500x500_0043_120409.fits', pmapdata, pmapheader
     if chname eq '2' then fits_read, '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_100x100_ch2_r2p25_s3_7_0p1s_x4_150723.fits', pmapdata, pmapheader
     if chname eq '1' then fits_read, '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_500x500_ch1_r2p25_s3_7_sd_0p4s_sdark_150722.fits', pmapdata, pmapheader
     c = contour(pmapdata, /fill, n_levels = 15, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', title = planetname, aspect_ratio = 1, xrange = [0,100], yrange = [0,150], axis_style = 0) ;21
      if planethash[aorname(0),'exptime'] lt 10. then begin
        xss = 14.5
        yss = 14.5
     endif else begin
        xss = 23.5
        yss = 232.5
     endelse
;  print, 'naor', n_elements(aorname)
     for a = startaor ,  stopaor  do begin
;     print, 'testing aors', a, colorarr[a]
;      if planethash[aorname(a),'exptime'] gt 0.2 then begin
         ;scolor = 'blue'
;         scolor = 'red'

;         xcen500 = 500.* ((planethash[aorname(a),'xcen']) - 14.5)
;         ycen500 = 500.* ((planethash[aorname(a),'ycen']) - 14.5)
;         an = plot(xcen500, ycen500, '1s', sym_size = 0.1,   sym_filled = 1, color = colorarr(a),/overplot)
;      endif
;      if planethash[aorname(a),'exptime'] gt 0.2 then begin
         ;scolor = 'blue'
;        scolor = 'red'
        print, 'mean position', mean(planethash[aorname(a),'xcen']), mean(planethash[aorname(a),'ycen'])
        xcen500 = 100.* ((planethash[aorname(a),'xcen']) - xss)
        ycen500 = 100.* ((planethash[aorname(a),'ycen']) - yss)
        ;limit = where(indgen(n_elements(planethash[aorname(a),'xcen'])) mod 100 eq 0) 
        an = plot(xcen500, ycen500, '1s', sym_size = 0.3,   sym_filled = 1, color = colorarr[a],/overplot)
;;        an = plot(xcen500(limit), ycen500(limit), '1s', sym_size = 0.3,   sym_filled = 1, color = colorarr[a],/overplot)
;      endif

;     npmax = [4.8, 37., 4.5]
;     testnp = where((planethash[aorname(a),'np']) gt npmax(a) , complement = nospike)
;     xcen500 = ((planethash[aorname(a),'xcen'])[nospike] )
;     ycen500 = ((planethash[aorname(a),'ycen'])[nospike])
;     print, 'mean np', mean((planethash[aorname(a),'np'])[nospike])
;     an = plot(xcen500, ycen500, '1s', sym_size = 0.2,   sym_filled = 1, color = 'blue', xtitle = 'xcen', ytitle = 'ycen', xrange = [14.9, 15.2], yrange = [14.75, 15.25]) 
                                ;    xcen500 = ((planethash[aorname(a),'xcen'])[testnp] )
                                ;    ycen500 = ((planethash[aorname(a),'ycen'])[testnp])
                                ;    print, 'mean np', mean((planethash[aorname(a),'np'])[testnp])
                                ;    an = plot(xcen500, ycen500, '1s', sym_size = 0.4,   sym_filled = 1, color = 'red',/overplot)
      
; ;    print, 'xcen500', xcen500
; ;    print, 'ycen500', ycen500
     endfor
     xs = 500.*(xsweet - xss)
     ys = 500.*(ysweet - yss)
;   an = text(xs, ys, 'X', /data, alignment = 0.5, vertical_alignment = 0.5)
;   an.save, dirname+'position_ch'+chname+'.png'
  endif 
;GOTO, Jumpend

;-----

  ;print, 'min', min((planetob[a].timearr - planetob[0].timearr(0))/60./60.)
 ; z = pp_multiplot(multi_layout=[1,3], global_xtitle='Time (hrs) ')
  if keyword_set(unbinned_xplot) then begin
     for a =  startaor,stopaor do begin
        if a eq startaor then begin
           xmin = 14.0          ; mean(planethash[aorname(a),'xcen'])-0.25
           xmax = 15.5          ; mean(planethash[aorname(a),'xcen'])+0.25
                                ;print, 'xmin.xmax', xmin, xmax

           am = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(a),'timearr'])(0))/60./60., planethash[aorname(a),'xcen'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], ytitle = 'X pix', title = planetname, xtitle = 'Time(hrs)', xrange = [0.,8], yrange = [14.9, 15.2])
        endif else begin
           am = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'xcen'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;           am = plot( (planethash[aorname(a),'bmjdarr'] - (planethash[aorname(a),'bmjdarr'])(0)), planethash[aorname(a),'xcen'],'6r1s-', sym_size = 0.3,   sym_filled = 1, color = colorarr[a],yrange = [14.47, 14.62], xrange = [0.1501,0.15018], ytitle = 'xcen');,/overplot)
        endelse
        
     endfor
     am.save, dirname + 'x_time_ch'+chname+'.png'
  endif

;------
  if keyword_set(unbinned_yplot) then begin
     for a = startaor,stopaor do begin ; 0, n_elements(aorname) - 1 do begin
        if a eq startaor then begin
           ymin = 14.9          ; mean(planethash[aorname(a),'ycen'])-0.25
           ymax = 16.0          ;mean(planethash[aorname(a),'ycen'])+0.25
                                ;print, 'ymin.ymax', ymin, ymax
           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(a),'timearr'])(0))/60./60., planethash[aorname(a),'ycen'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], ytitle = 'Y pix', title = planetname, xtitle = 'Time(hrs)', xrange = [0.,8], yrange = [14.9, 15.2])
        endif else begin
           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'ycen'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;           ay = plot( (planethash[aorname(a),'bmjdarr'] - (planethash[aorname(a),'bmjdarr'])(0)), planethash[aorname(a),'ycen'],'6r1s-', sym_size = 0.3,   sym_filled = 1, color = colorarr[a] , yrange = [14.3, 14.45], xrange = [0.1501,0.15018], ytitle = 'ycen');,/overplot)
        endelse
        
     endfor
     ay.save, dirname +'y_time_ch'+chname+'.png'
  endif
 
;------
  if keyword_set(unbinned_fluxplot) then begin
     unormfactor = median(planethash[aorname(0),'corrflux'])

     for a =0, n_elements(aorname) - 1 do begin
        ;;check if I should be using pmap corr or not
        ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
        ;;if 20% of the values are correctable than go with the pmap corr 
        print, '0.2nflux, ncorr, ', 0.2*n_elements([planethash[aorname(a),'flux']]), corrcount
        if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
        if a eq 0 then begin
                                ;       ymin = 14.9             ; mean(planethash[aorname(a),'ycen'])-0.25
                                ;       ymax = 16.0             ;mean(planethash[aorname(a),'ycen'])+0.25
                                ;print, 'ymin.ymax', ymin, ymax
           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., $
                      (planethash[aorname(a),'corrflux'])/unormfactor,'1s', sym_size = 0.1,   sym_filled = 1, $
                      color = colorarr[a], ytitle = 'Flux', xtitle = 'Time(hrs)',title = planetname,$
                      yrange = [0.99, 1.005], xrange = [0,190])
           ;ay = plot( (planethash[aorname(a),'phase']), (planethash[aorname(a),'flux'])/plot_corrnorm - 0.001,'1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], ytitle = 'Raw Flux', xtitle = 'Phase', yrange = [0.985, 1.015], xrange = [-0.5, 0.5]) ;,title = planetname, xtitle = 'Time(hrs)'
        endif else begin
           if pmapcorr gt 0 then begin
              ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., $
                         (planethash[aorname(a),'corrflux'])/unormfactor,'1s', sym_size = 0.1,   sym_filled = 1,$
                         color = colorarr[a], overplot = ay)
                                ; ay = plot(
                                ; (planethash[aorname(a),'phase']),
                                ; (planethash[aorname(a),'flux'])/plot_corrnorm - 0.001,'1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/overplot)
           endif
        endelse
        
     endfor
;     ay.xrange =[0,54]
;     ay.yrange = [0.0168, 0.0183]
     ay.save, dirname +'raw_flux_time_ch'+chname+'.png'
  endif

;------
  if keyword_set(unbinned_bkgplot) then begin

     for a = 0,n_elements(aorname) -1 do begin ; 0, n_elements(aorname) - 1 do begin
        if a eq 0 then begin
           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'bkgd'],'1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], ytitle = 'Bkgd', xtitle = 'Time(hrs)') ;,title = planetname
        endif else begin
           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'bkgd'],'1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/overplot)
        endelse
        
     endfor
     ay.save, dirname +'bkgd_time_ch'+chname+'.png'
  endif

;------
  if keyword_set(unbinned_xyfwhm) then begin

     for a = startaor,stopaor do begin ; 0, n_elements(aorname) - 1 do begin
        if a eq startaor then begin
           ax = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(startaor),'timearr'])(0))/60./60., planethash[aorname(a),'xfwhm'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], ytitle = 'XFWHM', xtitle = 'Time(hrs)', yrange = [1.9, 2.2]) ;,title = planetname
           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(startaor),'timearr'])(0))/60./60., planethash[aorname(a),'yfwhm'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], ytitle = 'YFWHM', xtitle = 'Time(hrs)', yrange = [1.9, 2.2]) ;,title = planetname
        endif else begin
           ax = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(startaor),'timearr'])(0))/60./60., planethash[aorname(a),'xfwhm'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], overplot = ax)
           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(startaor),'timearr'])(0))/60./60., planethash[aorname(a),'yfwhm'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], overplot = ay)
        endelse
        
     endfor
     ay.save, dirname +'xyfwhm_time_ch'+chname+'.png'
  endif

;------
  if keyword_set(unbinned_npplot) then begin

;     for a = 0 ,n_elements(aorname) -1 do begin ; 0, n_elements(aorname) - 1 do begin
;        if a eq 0 then begin
;           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(a),'timearr'])(0))/60./60., planethash[aorname(a),'np'],'1s', sym_size = 0.1,   sym_fil;led = 1, color = colorarr[a], ytitle = 'NP',title = planetname, xtitle = 'Time(hrs)')
;        endif else begin
;           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'np'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], /overplot)
;           ay = plot( (planethash[aorname(a),'bmjdarr'] - (planethash[aorname(a),'bmjdarr'])(0)), planethash[aorname(a),'np'],'6r1s-', sym_size = 0.3,   sym_filled = 1, color = colorarr[a], xrange = [0.1501,0.15018], yrange = [7.3,8.5], ytitle = 'NP');,/overplot)
;        endelse
        
;     endfor
;     ay.save, dirname +'np_time_ch'+chname+'.png'

;compare to np coming from box_centroider
     for a = 0 , n_elements(aorname) - 1 do begin
        print, 'should be plotting npcent'
        if a eq 0 then begin
           ayc = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(a),'timearr'])(0))/60./60., planethash[aorname(a),'npcentroids'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], ytitle = 'NP',title = planetname, xrange = [0,8],yrange = [4,8], xtitle = 'Time(hrs)')
        endif else begin
           ayc = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., planethash[aorname(a),'npcentroids'],'1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a], /overplot)
;           ay = plot( (planethash[aorname(a),'bmjdarr'] - (planethash[aorname(a),'bmjdarr'])(0)), planethash[aorname(a),'np'],'6r1s-', sym_size = 0.3,   sym_filled = 1, color = colorarr[a], xrange = [0.1501,0.15018], yrange = [7.3,8.5], ytitle = 'NP');,/overplot)
        endelse
        
     endfor
  endif
;------
  if keyword_set(DN_plot) then begin
     fluxconv = 0.1469 ;MJy/sr / DN/s;
;     exptime = planethash[aorname(a),'exptime']
     exptime = 1.92 ;XXX
     gain = 3.71
     nch = 2 ; XXX
     pxscal1 = [-1.22334117768332D, -1.21641835430637D, -1.22673962032422D, -1.2244968325831D]
     pxscal1 = abs(pxscal1)
     pxscal2 = [1.22328355209902D, 1.21585676679388D, 1.22298117494211D, 1.21904995758086D]
     pscale2 = pxscal1[nch-1] * pxscal2[nch-1]
     scale = pscale2 * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06

   for a = 0,n_elements(aorname) -1 do begin ; 0, n_elements(aorname) - 1 do begin

      mjy = planethash[aorname(a),'flux']
      mjypersr = mjy/ scale
      DNpers = mjypersr / fluxconv
      DN = DNpers *exptime
      e = DN * gain
        if a eq 0 then begin
           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., e,'1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], ytitle = 'Electrons', xtitle = 'Time(hrs)', yrange = [7000, 8E4]) ;,title = planetname
        endif else begin
           ay = plot( (planethash[aorname(a),'timearr'] - (planethash[aorname(0),'timearr'])(0))/60./60., e,'1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/overplot)
        endelse
        
     endfor
    ay.save, dirname +'e_time_ch'+chname+'.png'

  endif

;------
  
;------
     
     
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;goto, jumpend

;setup some arrays to hold binned values
  bin_corrfluxfinal = fltarr(n_elements(aorname))
  bin_phasefinal = fltarr(n_elements(aorname))
  bin_corrfluxerrfinal = fltarr(n_elements(aorname))

  mean_corrfluxarr = fltarr(stopaor - startaor + 1,/nozero)
  mean_dsweetarr = mean_corrfluxarr
  stddev_corrfluxarr = mean_corrfluxarr
  stddev_dsweetarr = mean_corrfluxarr
 

;pick the normalization
  phase0 = planethash[aorname(0),'phase']
  corrflux0 = planethash[aorname(0),'corrflux_d']
  se = where(phase0 gt 0.47 and phase0 lt 0.51)
  
  plot_corrnorm =  mean(corrflux0(se),/nan)
  
  for a = startaor,stopaor do begin
     print, '------------------------------------------------------'
     print, 'working on AOR', a, '   ', aorname(a), startaor
     ;divide plots by exposure time
     ;if planethash[aorname(a),'exptime'] lt 0.2 then scolor = 'blue' else scolor = 'red'

 
     ;check if I should be using pmap corr or not
     ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
     ;if 20% of the values are correctable than go with the pmap corr 
     print, '0.2nflux, ncorr, ', 0.2*n_elements([planethash[aorname(a),'flux']]), corrcount
     if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
     print, 'pmapcorr', pmapcorr
; bin within each AOR only, don't want to bin across boundaries.
;     junk =planethash[aorname(a),'corrflux'] 
;     print, 'corrflux before binning', junk[0:100]
     
     ;make one data point per AOR
     ;;bin_level = n_elements([planethash[aorname(a),'flux']]) 

     if keyword_set(set_nbins) then begin
        if a ge stareaor then begin
           ;;variable binning
           junkpar = binning_function(a, bin_level, pmapcorr,chname, /set_nbins, n_nbins = 4)
        endif else begin
           
           junkpar = binning_function(a, bin_level, pmapcorr,chname, /set_nbins, n_nbins = 95)
        endelse
     endif else begin 
        junkpar = binning_function(a, bin_level, pmapcorr,chname)
     endelse

     ;use time degraded corrfluxes
     bin_corrfluxp = bin_corrflux_dp


;     print, 'testing phase', (planethash[aorname(a),'phase'] )[0:10], (planethash[aorname(a),'phase'] )[600:610]
;     print, 'testing x', bin_xcen[0:10]
;     print, 'testing y', bin_ycen[0:10]
;     print, 'testing npcent', bin_npcent[0:10]
;     print, 'testing error bars', colorarr[a], bin_corrfluxerrp
;------------------------------------------------------
     ;possible comparison statistic
     ;what is the distribution of standard deviations among the corrfluxes?

     meanclip_jk, planethash[aorname(a),'corrflux'] , mean_corrflux, stddev_corrflux
     stddev_corrfluxarr[a - startaor] = stddev_corrflux
;------------------------------------------------------
     ;save some variables for later
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


;now try plotting

     if keyword_set(dsweet) then begin
        ;measure mean distance from the sweet spot
        xdist = (planethash[aorname(a),'xcen'] - xsweet)
        ydist = (planethash[aorname(a),'ycen'] - ysweet)
        meanclip_jk, ydist, meany, stddevy
        meanclip_jk, xdist, meanx, stddevx

        dsweet = sqrt(xdist^2 + ydist^2)
        meanclip_jk, dsweet, mean_dsweet, stddev_dsweet
        
        meanclip_jk, planethash[aorname(a),'corrflux'] , mean_corrflux, stddev_corrflux

        theta = 180./!Pi*ATAN(ydist, xdist) ; in degrees
        
        ;ahhh quadrants
;        if xdist gt 0 and ydist gt 0 then theta = theta
;        if mean(xdist) ge 0 and mean(ydist) lt 0 then theta = theta + 360
;        if mean(xdist) lt 0 then theta = theta +180

       meanclip_jk, theta, mean_theta, stddev_theta

        if a eq startaor then begin 
           normfactor = mean_corrflux
;           sweetplot = errorplot([mean_theta],[mean_corrflux/normfactor],  [stddev_theta], [stddev_corrflux/normfactor],$
;                                 '1s', sym_size = 0.4,   sym_filled = 1,  xtitle = 'Theta', ytitle = 'Corrflux', title = planetname,$
;                                color = colorarr[a],  yrange = [0.99, 1.01]) 
;           stddevplot = plot([mean_theta],[stddev_corrflux/normfactor], '1s', sym_size = 2.0, color = colorarr[a],  $
;                             sym_filled = 1,  xtitle = 'Theta', ytitle = 'Stddev Corrflux', title = planetname) 
           dtheta = plot([mean_theta], [stddevx],$
                                 '1s', sym_size = 1.0,   sym_filled = 1,  xtitle = 'Theta', ytitle = 'Stddev X centroids', title = planetname,$
                                color = colorarr[a]) 
        endif else begin
;           sweetplot.window.SetCurrent
;           sweetplot = errorplot([mean_theta],[mean_corrflux/normfactor],  [stddev_theta], [stddev_corrflux/normfactor],$
;                                 '1s', sym_size = 0.4,  sym_filled = 1, color = colorarr[a],/overplot) 
;           stddevplot.window.SetCurrent
;           stddevplot = plot([mean_theta],[stddev_corrflux/normfactor], '1s', sym_size = 2.0, color = colorarr[a],  $
;                             sym_filled = 1,/overplot) 
;           dtheta.window.SetCurrent
           dtheta = plot([mean_theta], [stddevx],$
                                 '1s', sym_size = 1.0,   sym_filled = 1, color = colorarr[a],/overplot) 
        endelse


     endif

;------------------------------------------------------


     if keyword_set(sclkplot) then begin

        sclk_0 = planethash[aorname(a),'sclktime_0']  ; in seconds
        
        ;turn bmjd into sclk
        elapsed_days= bin_bmjdarr(n_elements(bin_bmjdarr)-1) - bin_bmjdarr(0)
        elapsed_secs = elapsed_days*24.*60.*60.
        bin_sclk = dindgen(n_elements(bin_bmjdarr))*elapsed_secs / n_elements(bin_bmjdarr)
        bin_sclk = bin_sclk + sclk_0

        print, 'testing sclk', bin_sclk[0:10], bin_sclk(n_elements(bin_sclk)-1)
        print, 'showing bmjd', bin_bmjdarr[0:10], bin_bmjdarr(n_elements(bin_bmjdarr)-1)

        ;psn= plot(bin_sclk, bin_np, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], $
        ;            xtitle = 'Sclk_time', ytitle = 'Noise Pixel', title = planetname)
        psx = plot(bin_sclk, bin_xcen, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], $
                    xtitle = 'Sclk_time', ytitle = 'Xcen', title = planetname)
        psy= plot(bin_sclk, bin_ycen, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], $
                    xtitle = 'Sclk_time', ytitle = 'Ycenl', title = planetname)
     endif

;------------------------------------------------------
     ;make plot symbold based on the size of the AOR
;     if n_elements([planethash[aorname(a),'flux']])  lt 1000 then plotsym = 's' else plotsym = '*'
     plotsym = 's'
     setyrange = [0.99, 1.005]
     if planetname eq 'WASP-80b' then setyrange = [0.96, 1.02];[0.999, 1.001]
     if keyword_set(phaseplot) then begin ;make the plot as a function of phase
        ;print, ' phase', bin_phase
        print, 'corrflux',  (bin_corrfluxp/plot_corrnorm)
        setxrange =  [-0.5, 0.5] 
        corrnormoffset =  0; 0.02
        corroffset = 1.0
        if planetname eq 'WASP14-b' then corroffset =  1.000;0.0015
        if planetname eq 'HD209458' then corroffset = 0.001
       
        setynormfluxrange = [0.99, 1.005];[0.97, 1.005]

 ;       cphase = where(bin_phase lt 0)
 ;       bin_phase(cphase) = bin_phase(cphase) + 1.0


        if a eq startaor then begin    ; for the first AOR
                                ;set the normalization values from the
                                ;medians of the first AOR...is close enough
;           print, 'beginning AOR bin phase', bin_phasep
;           print, 'beginning AOR bin xcen', bin_xcen
;           print, 'beginning AOR bin ycen', bin_ycen
           plot_norm = median(bin_flux)  ; only happens for the first AOR
           print, 'plot_corrnorm', plot_corrnorm, bin_corrfluxp[0]
           if keyword_set(all_plots) then begin
              pp = plot(bin_phase, bin_xcen, '1s', sym_size = 0.2,   sym_filled = 1, title = planetname, $
                        color = colorarr[a], xtitle = 'Orbital Phase', ytitle = 'X position', $
                        xrange = setxrange)

              pq = plot(bin_phase, bin_ycen, '1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], $
                        xtitle = 'Orbital Phase', ytitle = 'Y position', title = planetname, $
                        xrange = setxrange)

              if keyword_set(errorbars) then begin
                 pr = errorplot(bin_phase, bin_flux/plot_norm,bin_fluxerr/plot_norm,  '1s', sym_size = 0.2,   $
                                sym_filled = 1,  color =colorarr[a],  xtitle = 'Orbital Phase', $
                                ytitle = 'Normalized Flux', title = planetname, yrange = setynormfluxrange, $
                                xrange = setxrange )
              endif else begin
                                ;             print, 'testing bin_flux', bin_flux[1:10]/plot_norm
                 pr = plot(bin_phase, bin_flux/plot_norm, '1s', sym_size = 0.2,   sym_filled = 1,  $
                        color = colorarr[a],  xtitle = 'Orbital Phase', ytitle = 'Normalized Flux', $
                        title = planetname, yrange = setynormfluxrange, xrange = setxrange) 
              endelse           ;/errorbars


              if pmapcorr eq 1 then begin
                 print, 'inside pmapcorr eq 1', median( (bin_corrfluxp/plot_corrnorm) ), median(bin_flux)
              if keyword_set(errorbars) then begin
                 pr = errorplot(bin_phasep, (bin_corrfluxp/plot_corrnorm) -corrnormoffset,  $
                        bin_corrfluxerrp/plot_corrnorm ,/overplot, '1s', sym_size = 0.2, $
                        sym_filled = 1, color = colorarr[a])
              endif else begin
                 print, 'plotting corrfluxp',  (bin_corrfluxp(1)/plot_corrnorm) -corrnormoffset
                 pr = plot(bin_phasep, (bin_corrfluxp/plot_corrnorm) -corrnormoffset ,overplot = pr, '1s', $
                           sym_size = 0.2,   sym_filled = 1, color = colorarr[a]);
              endelse           ;/errorbars
                 
              endif             ;enough pmap corrections
              
              ps= plot(bin_phase, bin_npcent, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], $
                       xtitle = 'Orbital Phase', ytitle = 'Noise Pixel', title = planetname,$
                       xrange = setxrange)
              
              pt = plot(bin_phase, bin_bkgd, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], $
                        xtitle = 'Orbital Phase', ytitle = 'Background', title = planetname, $
                        xrange =setxrange)
              
              pxf = plot(bin_phase, bin_xfwhm, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], $
                         xtitle = 'Orbital Phase', ytitle = 'XFWHM', title = planetname, $
                         xrange =setxrange, yrange = [1.9, 2.3])
              
              pyf = plot(bin_phase, bin_yfwhm, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], $
                         xtitle = 'Orbital Phase', ytitle = 'YFWHM', title = planetname, $
                         xrange =setxrange, yrange = [1.9, 2.3])
           endif

                                ;setup a plot for just the snaps
;           if n_elements(aorname) gt 10 and
           if pmapcorr eq 1 then begin
;              pu = errorplot(bin_phase, (bin_flux/plot_norm) + corroffset, bin_fluxerr/plot_norm,  sym_size = 0.7,  $
;                        symbol = plotsym, sym_filled = 1,color =scolor ,xtitle = 'Orbital Phase', errorbar_color =  scolor, $
;                        title = planetname, ytitle = 'Flux',  yrange
;                        = [0.995, 1.02], xrange = setxrange) ;
              if a le stareaor then begin
                 scolor = colorarr[a]; 'gray' 
                 if planetname eq 'WASP-14b' then corroffset =  1.000; 0.0015    
                 if planetname eq 'HD158460' then corroffset = 0.0
              endif else begin
                 scolor = colorarr[a]
                 if planetname eq 'WASP-14b' then corroffset =  1.000; 0.0015
                 if planetname eq 'HD209458' then corroffset = 0.001
              endelse
              print, 'corroffset', corroffset
              pu = errorplot(bin_phasep, (bin_corrfluxp/plot_corrnorm) * corroffset, $
                             bin_corrfluxerrp/plot_corrnorm,  '1s', sym_size = 0.7,  $
                             symbol = plotsym, sym_filled = 1,color =scolor ,xtitle = 'Orbital Phase',$
                             errorbar_color =  scolor, title = planetname, ytitle = 'Pmap Corrected Flux', $
                             yrange = setyrange, xrange = setxrange);, position =[0.2,0.35,0.9,0.9]) ;
              if keyword_set(function_fit) then begin
                 pu.xshowtext = 0
                 pu.position =  [0.2,0.35,0.9,0.9]
              endif

           endif

        endif                   ; if a = 0

      ;  test_med(a) = median(bin_corrfluxp)
        print,'median bin_corrfluxp ',  median(bin_corrfluxp), median(bin_flux)
;        print, 'a', a
;-------------------------------------
        if (a gt 0) and (a lt stareaor) then begin
           if planetname eq 'WASP-14b' then corroffset =  1.000; 0.0015 ;
           print, 'inside a gt 0 a le stareaor', a, corroffset
           print, 'just before pmapcorr',mean(bin_flux/plot_norm)

           if pmapcorr eq 1 then begin
              print, 'inside pmapcorr agt 0 a le stareaor'
              print, 'mean levle', plot_corrnorm, mean((bin_corrfluxp/plot_corrnorm)-corrnormoffset)
;              pr = plot(bin_phasep, (bin_corrfluxp/plot_corrnorm)-corrnormoffset, /overplot, 's1', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/current)
           endif

           if keyword_set(all_plots) then begin
              pp.window.SetCurrent
              pp = plot(bin_phase, bin_xcen, '1s', sym_size = 0.2,   sym_filled = 1,color = colorarr[a],  /overplot,/current)
              pq.window.SetCurrent
              pq = plot(bin_phase, bin_ycen, '1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], /overplot,/current)
;           pr.window.SetCurrent
;           pr = plot(bin_phase, bin_flux/plot_norm , '1s', sym_size = 0.2,   sym_filled = 1,  color = 'green', /overplot,/current,/nodata)
              ps.window.SetCurrent
              ps = plot(bin_phase, bin_npcent, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], /overplot,/current) ;, xtitle = 'Orbital Phase', ytitle = 'Normalized Flux',) 
              
              pt.window.SetCurrent
              pt = plot(bin_phase, bin_bkgd, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a],/overplot,/current) 
              

              pxf = plot(bin_phase, bin_xfwhm, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], overplot = pxf)
              
              pyf = plot(bin_phase, bin_yfwhm, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], overplot = pyf)
              
           if planetname eq 'WASP-15b' then begin
              ;get a quick calculation of depth of eclipse
              ec = where(bin_phase ge 0.49 and bin_phase le 0.51)
              ef = where(bin_phase ge 0.525 or bin_phase le 0.47)
              delta = 1 - (mean(bin_corrfluxp(ec))) / (mean(bin_corrfluxp(ef)))
              print, 'estimated depth', delta, mean(bin_corrfluxp(ef)), mean(bin_corrfluxp(ec))
           endif
        endif
           
           if pmapcorr eq 1 then begin
              pu.window.SetCurrent
;              pu =  errorplot(bin_phase, (bin_flux/plot_norm) + corroffset, bin_fluxerr/plot_norm, sym_size = 0.7,  $
;                              symbol = plotsym, sym_filled = 1,color
;                              = scolor, errorbar_color =  scolor,
;                              /overplot, /current)
              
              pu = errorplot(bin_phasep, (bin_corrfluxp/plot_corrnorm) * corroffset , bin_corrfluxerrp/plot_corrnorm,$
                             '1s', sym_size = 0.7, symbol = plotsym, sym_filled = 1,color =colorarr[a], $
                             errorbar_color = colorarr[a], overplot = pu) ;
           endif
           
        endif  ; a gt 0 a lt stareor
        


        if a gt stareaor then begin
           if planetname eq 'WASP-14b' then corroffset =  1.0015; 0.0015
           print, 'inside a ge stareaor', a, corroffset
           if keyword_set(all_plots) then begin
              pp.window.SetCurrent
              pp = plot(bin_phase, bin_xcen, '1s', sym_size = 0.2,   sym_filled = 1,color = colorarr[a],$
                        /overplot,/current)
              
              pq.window.SetCurrent
              pq = plot(bin_phase, bin_ycen, '1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],$
                        /overplot,/current)
              
              pr.window.SetCurrent
              if keyword_set(errorbars) then begin
                 pr = errorplot(bin_phase, bin_flux/(plot_norm) , bin_fluxerr/plot_norm, '1s', $
                                sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], /overplot,/current) 
              endif else begin
                 pr = plot(bin_phase, bin_flux/(plot_norm) , '1s', sym_size = 0.2,   sym_filled = 1,  $
                           color = colorarr[a], /overplot,/current)
              endelse
              
              ps.window.SetCurrent
              ps = plot(bin_phase, bin_npcent, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a],$
                        /overplot,/current) 

              pt.window.SetCurrent
              pt = plot(bin_phase, bin_bkgd, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], $
                        /overplot,/current) 
              
              pxf = plot(bin_phase, bin_xfwhm, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], overplot = pxf)
              
              pyf = plot(bin_phase, bin_yfwhm, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], overplot = pyf)
           endif

           if pmapcorr eq 1 then begin
              print, 'inside pmapcorr eq 1', median( (bin_corrfluxp/plot_corrnorm) +corroffset)
              
;               if keyword_set(errorbars) then begin
;                  pr = errorplot(bin_phasep, (bin_corrfluxp/plot_corrnorm)- corrnormoffset,$
;                            bin_corrfluxerrp/ plot_corrnorm, '1s', sym_size = 0.2,   $
;                            sym_filled = 1, color = colorarr[a],/overplot,/current)
;               endif else begin
;                  pr = plot(bin_phasep, (bin_corrfluxp/plot_corrnorm) - corrnormoffset, '1s', $
;                            sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/overplot,/current)
;               endelse
              
              
                                ;setup a plot for just the snaps
              if n_elements(aorname) gt 10 then begin
                 
                 pu = errorplot(bin_phasep, (bin_corrfluxp/plot_corrnorm) * corroffset , bin_corrfluxerrp/plot_corrnorm, '1s', $
                                sym_size = 0.7, symbol = plotsym, sym_filled = 1,color =colorarr[a] , $
                                errorbar_color = colorarr[a], overplot = pu ) ;
              endif
           endif
           
           
        endif
        
        
;        ps.save, dirname +'binnp_phase_ch'+chname+'.png'
;        pt.save, dirname +'binbkg_phase_ch'+chname+'.png'
 ;       pu.save, dirname + 'binsnaps_phase_ch'+chname+'.png'
;        pr.save, dirname + 'binflux_phase_ch'+chname+'.png'

        if keyword_set(selfcal) then begin
         ;  restore, strcompress(dirname + 'selfcal.sav')          
          ; print, 'test selfcal', y[0:10], 'x', x[0:10]
          ; pr.window.SetCurrent
           ;pr = errorplot(x, y-0.01,yerr,'1bo',sym_size = 0.2, sym_filled = 1, errorbar_color = 'black',color = 'black',/overplot,/current)
           ;pr = plot(x, y -0.01, '1s', sym_size = 0.3,   sym_filled = 1, color = 'black',/overplot,/current)
        endif




        if a eq 0 and planetname eq 'HD209458' then begin
           print, 'plotting Lewis curve'
           ;plot curves from Lewis
           restore, '/Users/jkrick/irac_warm/hd209458/hd209_ch2_for_jessica.sav'
           ;need to get time in phase
           time = time / period ; now in phase

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

           ;pr.window.SetCurrent
           pu = plot(bin_time, bin_corr_flux , '1s', sym_size =  0.3, sym_filled = 1, color ='grey', overplot = pu)
        endif

        

     endif 

;------------------------------------------------------
;     print, 'n_elements time', n_elements(bin_timearr), n_elements(bin_xcen), n_elements(bin_phase)
     if keyword_set(timeplot) then begin ;make the plot as a function of time

        corroffset = 0.0000
        print, 'making plot as a function of time'
        corrnormoffset = 0;0.02
        setynormfluxrange = [0.98, 1.01];[0.97, 1.005]
        if a eq startaor then begin ; for the first AOR
                                ;set the normalization values from the
                                ;medians of the first AOR...is close enough
           plot_corrnorm = median(bin_corrfluxp)
           print, 'plot_corrnorm', plot_corrnorm, mean(bin_corrfluxp,/nan), corroffset
           plot_norm = median(bin_flux)
           pp = plot((bin_timearr - time_0)/60./60., bin_xcen, '1s', sym_size = 0.3,   sym_filled = 1, $;title = planetname, $
                     color = colorarr[a], xtitle = 'Time(hrs)', ytitle = 'X position', aspect_ratio = 0.0, margin = 0.2);, $
                     ;xrange = setxrange)
                     

           pqy = plot((bin_timearr - time_0)/60./60., bin_ycen, '1s', sym_size = 0.3, sym_filled = 1, color = colorarr[a], $
                     xtitle = 'Time(hrs)', ytitle = 'Y position', aspect_ratio = 0.0, margin = 0.2);, $, title = planetname
                     ;xrange = setxrange)

 ;          if keyword_set(errorbars) then begin
 ;             pr = errorplot((bin_timearr - time_0)/60./60., bin_flux/plot_norm,bin_fluxerr/plot_norm,  '1s', sym_size = 0.3,  $
 ;                            sym_filled = 1,  color =colorarr[a],  xtitle = 'Time(hrs)', $
 ;                            ytitle = 'Normalized Flux' ,yrange = setynormfluxrange, $ ;, title = planetname,
 ;                            xrange = setxrange , aspect_ratio = 0.0, margin = 0.2)
 ;          endif else begin
 ;             pr = plot((bin_timearr - time_0)/60./60., bin_flux/plot_norm, '1s', sym_size = 0.3,   sym_filled = 1,  $
 ;                       color = colorarr[a],  xtitle = 'Time(hrs)', ytitle = 'Normalized Flux', $ ;title = planetname, 
 ;                       yrange = setynormfluxrange, aspect_ratio = 0.0, margin = 0.2);, xrange = setxrange);,/nodata) 
 ;          endelse  ;/errorbars

           if pmapcorr eq 1 then begin
              print, 'inside pmapcorr eq 1', median(bin_corrfluxp), median(bin_flux)
              if keyword_set(errorbars) then begin
                 pr = errorplot((bin_timearrp - time_0)/60./60., (bin_corrfluxp/plot_corrnorm) -corrnormoffset,  $
                        bin_corrfluxerrp/plot_corrnorm , '1s', sym_size = 0.3, $
                        sym_filled = 1, color = colorarr[a], xtitle = 'Time(hrs)', $
                                ytitle = 'Normalized Corrected Flux', title = planetname,aspect_ratio = 0.0, margin = 0.2)
              endif else begin
                 pr = plot((bin_timearr - time_0)/60./60., (bin_flux/plot_corrnorm) -corrnormoffset,  '1s', $
                           sym_size = 0.3,   sym_filled = 1, color = colorarr[a],  xtitle = 'Time(hrs)', $
                           ytitle = 'Normalized Flux', title = planetname, $
                           yrange = setynormfluxrange, aspect_ratio = 0.0, margin = 0.2)
              endelse           ;/errorbars

           endif  ;enough pmap corrections

           ps= plot((bin_timearr - time_0)/60./60., bin_npcent, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], $
                    xtitle = 'Time(hrs)', ytitle = 'Noise Pixel',  aspect_ratio = 0.0, margin = 0.2);,$ title = planetname,
                    ;xrange = setxrange)

           pt = plot((bin_timearr - time_0)/60./60., bin_bkgd, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], $
                     xtitle = 'Time(hrs)', ytitle = 'Background',  aspect_ratio = 0.0, margin = 0.2);, $ title = planetname,
                    ; xrange =setxrange)

           pxf = plot((bin_timearr - time_0)/60./60., bin_xfwhm, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], $
                     xtitle = 'Time(hrs)', ytitle = 'XFWHM',  aspect_ratio = 0.0, margin = 0.2);, $ title = planetname,

           pyf = plot((bin_timearr - time_0)/60./60., bin_yfwhm, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], $
                     xtitle = 'Time(hrs)', ytitle = 'YFWHM',  aspect_ratio = 0.0, margin = 0.2);, $ title = planetname,

           ;setup a plot for just the snaps
           if n_elements(aorname) gt 6 and pmapcorr eq 1 then begin
              print, '((bin_timearr - time_0)/60./60.) - 2.E3', ((bin_timearr - time_0)/60./60.), ((bin_timearr - time_0)/60./60.) - 2.E3
              help, bin_timearr
              help, time_0
;              pu = plot(((bin_timearr - time_0)/60./60.) , (bin_corrfluxp/median(bin_all_corrfluxp)) + corroffset,$ ;- 5.E3)/24.
;                              sym_size = 1.0,  $
;                             symbol = plotsym, sym_filled = 1,color ='black' ,$
;                              overplot = ay) ;
              pu = errorplot(((bin_timearr - time_0)/60./60.) , (bin_corrfluxp/plot_corrnorm) + corroffset,$ ;- 5.E3)/24.
                             bin_corrfluxerrp/plot_corrnorm,  sym_size = 0.7,  $
                             symbol = plotsym, sym_filled = 1,color =colorarr[a] ,xtitle = 'Time (Hours)', $
                             errorbar_color =  colorarr[a], title = planetname, ytitle = 'Pmap Corrected Flux',  $
                             yrange = setyrange, xrange = [-5, 190]) ;

           endif



        endif                   ; if a = 0

      ;  test_med(a) = median(bin_corrfluxp)
        print,'median bin_corrfluxp ',  median(bin_corrfluxp), median(bin_flux)
        if (a gt 0) and (a le stareaor) then begin
           print, 'inside a gt 0 a le stareaor', a
           pp = plot((bin_timearr - time_0)/60./60., bin_xcen, '1s', sym_size = 0.3,   sym_filled = 1,color = colorarr[a],  overplot=pp)
           pqy = plot((bin_timearr - time_0)/60./60., bin_ycen, '1s', sym_size = 0.3,   sym_filled = 1, color = colorarr[a], overplot=pq)
;           pr.window.SetCurrent
;           pr = plot((bin_timearr - time_0)/60./60., bin_flux/plot_norm , '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], /overplot,/current,/nodata)

           ps = plot((bin_timearr - time_0)/60./60., bin_npcent, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], overplot=ps) ;, xtitle = 'Time(hrs)', ytitle = 'Normalized Flux',) 

           pt = plot((bin_timearr - time_0)/60./60., bin_bkgd, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a],overplot=pt) 

           if pmapcorr eq 1 then begin
;                  pu.window.SetCurrent
;                  pu =  errorplot(((bin_timearr - time_0)/60./60.) - 5.E3, (bin_flux/plot_norm) + corroffset, bin_fluxerr/plot_norm, $
;                                  sym_size = 0.7,symbol = plotsym, sym_filled = 1,color = scolor, errorbar_color =  scolor, $
;                                  overplot = pu, /current)
              ;;print, 'time, corrflux', (((bin_timearr - time_0)/60./60.) - 5.5E3)/24., (bin_corrfluxp/plot_corrnorm) + corroffset
;              pu =  plot(((bin_timearr - time_0)/60./60.) , (bin_corrfluxp/median(bin_all_corrfluxp)) + corroffset, $ ;- 5.5E3)/24.
;                               sym_size = 1.0, symbol = plotsym, sym_filled = 1,$
;                              color = 'black', overplot = ay,
;                              /current)  ;pu

              ;;pr = plot((bin_timearr - time_0)/60./60., (bin_corrfluxp/plot_corrnorm)-corrnormoffset, overplot=pr, 's1', sym_size = 0.2,   sym_filled = 1, color = colorarr[a])
              pr =  errorplot(((bin_timearrp - time_0)/60./60.) , (bin_corrfluxp/plot_corrnorm) + corroffset, $ ;- 5.5E3)/24.
                              bin_corrfluxerrp/plot_corrnorm, sym_size = 0.7, symbol = plotsym, sym_filled = 1,$
                              color = colorarr[a], errorbar_color =  colorarr[a], overplot = pr)  ;pu

           endif

        endif



        if a gt stareaor then begin
           print, 'inside a gt stareaor', a
           pp = plot((bin_timearr - time_0)/60./60., bin_xcen, '1s', sym_size = 0.3,   sym_filled = 1,color = colorarr[a],$
                     overplot=pp)

           pq.window.SetCurrent
           pq = plot((bin_timearr - time_0)/60./60., bin_ycen, '1s', sym_size = 0.3,   sym_filled = 1, color = colorarr[a],$
                     overplot=pq)

           if keyword_set(errorbars) then begin
;           pr.window.SetCurrent
;              pr = errorplot((bin_timearr - time_0)/60./60., bin_flux/(plot_norm) , bin_fluxerr/plot_norm, '1s', $
;                        sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], /overplot,/current) 
           endif else begin
;              pr = plot((bin_timearr - time_0)/60./60., bin_flux/(plot_norm) , '1s', sym_size = 0.3,   sym_filled = 1,  $
;                        color = colorarr[a], /overplot,/current,/nodata)
           endelse


           if pmapcorr eq 1 then begin
              print, 'inside pmapcorr eq 1', median( (bin_corrfluxp/plot_corrnorm) +corroffset)

               if keyword_set(errorbars) then begin
                  pr = errorplot((bin_timearr - time_0)/60./60., (bin_flux/plot_corrnorm)- corrnormoffset,$
                            bin_corrfluxerrp/ plot_corrnorm, '1s', sym_size = 0.3,   $
                            sym_filled = 1, color = colorarr[a],overplot=pr)
               endif else begin
                  pr = plot((bin_timearr - time_0)/60./60., (bin_flux/plot_corrnorm) - corrnormoffset, '1s', $
                            sym_size = 0.3,   sym_filled = 1, color = colorarr[a],overplot=pr)
               endelse
               
               if pmapcorr eq 1 then begin
;                  pu =  plot(((bin_timearr - time_0)/60./60.) , (bin_corrfluxp/median(bin_all_corrfluxp)) + corroffset, $ ;- 5.5E3)/24.
;                                  sym_size = 1.0, symbol = plotsym, $
;                                  sym_filled = 1,color = 'black',  overplot = ay) ;pu
                 pu =  errorplot(((bin_timearr - time_0)/60./60.) , (bin_corrfluxp/plot_corrnorm) + corroffset, $ ;- 5.5E3)/24.
                                  bin_corrfluxerrp/plot_corrnorm, sym_size = 0.7, symbol = plotsym, $
                                  sym_filled = 1,color = colorarr[a], errorbar_color =  colorarr[a], overplot = pu) ;pu

               endif
               
            endif

           ps = plot((bin_timearr - time_0)/60./60., bin_npcent, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a],$
                     overplot=ps) 

           pt = plot((bin_timearr - time_0)/60./60., bin_bkgd, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], $
                    overplot=pt) 
;           pu =  plot((bin_timearr - time_0)/60./60., (bin_corrfluxp/plot_corrnorm) + corroffset, '1s', sym_size = 0.2,  $
;                      sym_filled = 1,color = colorarr[a], overplot = pu, /current)
           pxf = plot((bin_timearr - time_0)/60./60., bin_xfwhm, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], $
                     overplot = pxf, /current)

           pyf = plot((bin_timearr - time_0)/60./60., bin_yfwhm, '1s', sym_size = 0.3,   sym_filled = 1,  color = colorarr[a], $
                     overplot = pyf, /current)

        endif


     endif

     if keyword_set(centerpixplot) then begin
        if a eq 0 then begin
           cplot = plot((bin_timearr - time_0)/60./60., bin_centerpix, '6rs1', sym_size = 0.2, sym_filled = 1, xtitle =  'Time (hrs)', ytitle = 'mean Pixel value', title = planetname) 
        endif else begin
           cplot = plot((bin_timearr - time_0)/60./60., bin_centerpix, '6rs1', sym_size = 0.2, sym_filled = 1,/overplot) 
        endelse
     endif

;keep these for fitting later
;     (bin_timearr - bin_timearr(0))/60./60.final(a) = bin_phasep
;     bin_corrfluxfinal(a) = (bin_corrfluxp/plot_corrnorm) 
;    bin_corrfluxerrfinal(a) = bin_corrfluxerrp / plot_corrnorm

  endfor                        ;binning for each AOR


  if keyword_set(line_fit) then begin
     
     ;reset time to hours starting at the first observation
     bin_all_time = (bin_all_time - time_0)/60./60.

     ;need to do some sigma clipping rejection
     meanclip, bin_all_corrfluxp, meanval, sigmaval, subs = goodnum
     print, 'n', n_elements(bin_all_corrfluxp), n_elements(goodnum)
     bin_all_time = bin_all_time(goodnum)
     bin_all_corrfluxp = bin_all_corrfluxp(goodnum)
     bin_all_corrfluxerrp = bin_all_corrfluxerrp(goodnum)
     bin_all_fluxerrp = bin_all_fluxerrp(goodnum)

     start = [-1E-3,median(bin_all_corrfluxp)]
     startflat = [1.0]
     result= MPFITFUN('linear',bin_all_time, bin_all_corrfluxp/median(bin_all_corrfluxp), $
                      bin_all_corrfluxerrp/median(bin_all_corrfluxp), $
                      start, perror = perror, bestnorm = bestnorm)    
        
     pu = plot(bin_all_time, result(0)*(bin_all_time) + result(1),  color = 'black', overplot = pu, linestyle = 2, thick = 3)
        
     DOF     = N_ELEMENTS(bin_all_time) - N_ELEMENTS(result) ; deg of freedom
     PCERROR = PERROR * SQRT(BESTNORM / DOF)                 ; scaled uncertainties
     print, 'linear result', result, 'pcerror', pcerror, ' perror', perror, 'bestnorm', bestnorm, 'dof', DOF, sqrt(bestnorm/dof)
;     t1 = text(0.2, 1.001, sigfig(result(0), 3, /scientific) + '    ' + sigfig(bestnorm/DOF, 5), /data, color = 'black', font_style = 'bold', target = pu) ;pcerror(0)
     
;     result2 = MPFITFUN('flatline',bin_all_time,bin_all_corrfluxp/median(bin_all_corrfluxp), bin_all_corrfluxerrp/median(bin_all_corrfluxp), $
;                        startflat, perror = perror, bestnorm = bestnorm)    
;     pu = plot(bin_all_time, 0*bin_all_time + result2(0),  color = 'blue', overplot = pu, thick = 3, linestyle = 2)
;     tt1 = text(0.2, 1.0015, '0    ' + sigfig(bestnorm/DOF, 5), /data, color = 'blue', font_style = 'bold', target = pu)
     
     
  endif

;----------------------------------------------------
;Monte Carlo
;randomize the time variables, re-measure the slopes.  
; What is the distribution of slopes measured from S such randomizations?
  
     if keyword_set(montecarlo) then begin
        MC = 1000               ; how many implementations should I use?
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

;------------------------------------
;fit a mandel and agol function and overplot it
     if keyword_set(function_fit) then begin
        modelfilename = strcompress(dirname + planetname +'_model_ch'+chname+'_'+string(apradius)+'.sav',/remove_all)
;        z = function_fit_lightcurve('WASP-14b', bin_all_phasep, bin_all_corrfluxp, bin_all_corrfluxerrp,modelfilename)
        
        restore, modelfilename
        pu = plot(phase, rel_flux, color = 'Sky Blue', thick = 4, overplot = pu)
        print, 'rel_flux', rel_flux


        ;overplot the Wong et al. fit to the staring mode data
        readcol, '/Users/jkrick/irac_warm/pcrs_planets/WASP-14b/wong/ch2.dat', w_phase, w_model
        pu = plot(w_phase, w_model + 0.0016, thick = 4, overplot = pu, color = 'green')

        ;;make some residual plots
        ;;wong model
        wong_close = fltarr(n_elements(bin_all_phasep))
        for nbp = 0, n_elements(bin_all_phasep) -1 do begin
           wong_close[nbp] = closest(w_phase, bin_all_phasep(nbp))
        endfor
        wmodel_snaps = w_model(wong_close)
        resid_wong = ((bin_all_corrfluxp/plot_corrnorm)) - wmodel_snaps + 0.0025
        pun = plot(bin_all_phasep, resid_wong, '1s', sym_size = 0.5,   sym_filled = 1, $
                   color = 'green', yrange = [-0.004, 0.004],xtitle = 'Orbital Phase',$
                  position = [0.2, 0.15, 0.9, 0.35],/current, ytitle = 'Residuals', $
                   xrange = pu.xrange,  ytickinterval = 0.003)

        ;;ca08 model
        ca08_close = fltarr(n_elements(bin_all_phasep))
        for nbp = 0, n_elements(bin_all_phasep) -1 do begin
           ca08_close[nbp] = closest(phase, bin_all_phasep(nbp))
        endfor
        model_snaps = rel_flux(ca08_close)
        resid_ca08 = ((bin_all_corrfluxp/plot_corrnorm)) - model_snaps +0.004
        pun = plot(bin_all_phasep, resid_ca08, '1s', sym_size = 0.5,   sym_filled = 1, $
                   color = 'sky blue',  overplot = pun)


     endif


  savename = strcompress(dirname + planetname +'_plot_ch'+chname+'.sav')
  save, /all, filename=savename
  print, 'saving planethash', savename


jumpend: print, 'ending '
end


;what is the distrubition of standard deviations in the corrfluxes?
;plothist, stddev_corrfluxarr, xhist, yhist, /autobin, /noplot,/nan
;di = plot(xhist, yhist, thick = 2, xtitle = 'standard deviation of corrflux', ytitle = 'Number', title = planetname, yrange = [0,20], xrange = [2.1E-4, 2.6E-4])

 ;if keyword_set(phaseplot) then  pr.save, dirname +'binflux_phase.png' else pl.save, dirname +'binflux_time.png'

 ;test

;;;;;;;-----------------------------------------------------------------
;overplot predicted light curve
;XXX change this to fitting next
;  exosystemname = "WASP-14 b"
;;  exoplanet_light_curve, t, rel_flux, exosystem = exosystemname,/verbose, /tphase
;  fit_exoplanet_light_curve, exosystemname, apradius
;  strput, savefilename, 'model',69
;  print, 'savefilename', savefilename
;  restore, savefilename
;  pu.window.SetCurrent
;  pu = plot(phase, rel_flux,'4', color = 'Sky Blue', /overplot,/current)
;  pu.save,dirname + 'snap_bin_corrflux_fit.png'
;;;;;;;-----------------------------------------------------------------




;  pp.save, dirname +'binxcen_time_ch'+chname+'.png'
;  pq.save, dirname +'binycen_time_ch'+chname+'.png'
;  pt.save, dirname +'binbkg_time_ch'+chname+'.png'
;  ps.save, dirname + 'binnp_time_ch'+chname+'.png'
;  pr.save, dirname + 'binflux_time_ch'+chname+'.png'
;  pxf.save, dirname + 'binxfwhm_time_ch'+chname+'.png'
;  pyf.save, dirname + 'binyfwhm_time_ch'+chname+'.png'
  
