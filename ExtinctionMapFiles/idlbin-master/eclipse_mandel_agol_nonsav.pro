pro eclipse_mandel_agol_nonsav, planetname
;planetname = 'XO3'


 colorarr = ['burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']

;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  
  exoplanet_data_file = '/Users/jkrick/idlbin/exoplanets.csv'
  exosystem = strmid(planetname, 0, 8 )+ ' b' ;'HD 209458 b' ;
  if planetname eq 'XO3' then exosystem = 'XO-3 b'
  lambdaname  = '4.5'           ; '3.6'
  get_exoplanet_data,EXOSYSTEM=exosystem,MSINI=msini,MSTAR=mstar,TRANSIT_DEPTH=transit_depth,RP_RSTAR=rp_rstar,AR_SEMIMAJ=ar_semimaj,$
                     TEQ_P=1315,TEFF_STAR=teff_star,SECONDARY_DEPTH=secondary_depth,SECONDARY_LAMBDA=lambdaname,$
                     INCLINATION=inclination,MJD_TRANSIT=mjd_transit,P_ORBIT=p_orbit,EXODATA=exodata,RA=ra,DEC=dec,VMAG=vmag,$
                     DISTANCE=distance,ECC=ecc,T14=t14,F36=f36,F45=f45,FP_FSTAR0=fp_fstar0;,/verbose
;;  print, 'p_orbit, mjd_transit', p_orbit, mjd_transit
  dirname = strcompress(basedir + planetname +'/')       
  s_file = ['Stevenson-XO3b2-Visit01-CorrectedLC.txt','Stevenson-XO3b2-Visit02-CorrectedLC.txt','Stevenson-XO3b2-Visit09-CorrectedLC.txt','Stevenson-XO3b2-Visit10-CorrectedLC.txt']

;    time_0 = (planethash[aorname(0),'timearr']) 
;  time_0 = time_0(0)

  modelfilename = strcompress(dirname + planetname +'_model_bliss.sav',/remove_all)

  startaor = 0
  stopaor = n_elements(s_file) - 1 
  for a = startaor,stopaor, 1 do begin
     print, '------------------------------------------------------'
     print, 'working on AOR', a, '   ', s_file(a)
     ;;unbinned stuff

     ;;read in the stevenson files
     readcol, strcompress(dirname + s_file(a),/remove_all), bmjd, corrflux, corrfluxerr, format = '(D0, D0, D0)'
 
     ;;get phasing out of the way here
     bmjd_dist = bmjd - mjd_transit ; how many UTC away from the transit center
     phase =( bmjd_dist / p_orbit )- fix(bmjd_dist/p_orbit)
;;     phase = phase + (phase lt -0.5 and phase ge -1.0)
;;     phase = phase - (phase gt 0.5 and phase le 1.0)

    ;;not actually binning
     bin_corrfluxp = corrflux
     bin_phasep = phase
     bin_corrfluxerrp = corrfluxerr
     print, n_elements(bin_corrfluxp), n_elements(bin_phasep), n_elements(bin_corrfluxerrp)
     

     ;fake that this is a transit/ not eclipse
     bin_phasep = bin_phasep -0.665; 0.595; 0.585

     ;pick the normalization
     se = where(bin_phasep gt -0.015 and bin_phasep lt 0.015)
     plot_corrnorm =  mean(bin_corrfluxp(se),/nan)
    ;; print, 'plot_corrnorm', plot_corrnorm, mean(bin_corrfluxp)

     pu = errorplot(bin_phasep, (bin_corrfluxp/plot_corrnorm), $
                             bin_corrfluxerrp/plot_corrnorm,  '1s', sym_size = 0.3,  $
                              sym_filled = 1,color =colorarr[a] ,xtitle = 'Orbital Phase',$
                             errorbar_color =  colorarr[a], title = s_file(a), ytitle = 'Corrected Flux')

    
     test = function_fit_lightcurve(planetname, bin_phasep, bin_corrfluxp/plot_corrnorm,  bin_corrfluxerrp/plot_corrnorm, modelfilename)
     p2 = plot(findgen(20) - 5, fltarr(20) + 1.0 + 0.001875, overplot = pu , xrange = [-0.08, 0.08]) ; 0.00158
     plotname = strcompress(dirname + planetname +'_stevenson_'+strmid(s_file(a),16, 7) +'.png',/remove_all) 
     p2.save, plotname

     ;write output to be used by TAP
;;     outfilename =  strcompress(dirname +'phot_ch'+chname+'_TAP_'+aorname(a)+'_nobin_bliss.ascii',/remove_all)
;;     openw, outlun, outfilename,/GET_LUN
;;     for te = 0, n_elements(bin_bmjdarrp) -1 do begin ;n_elements(bin_phasep) -1
;;        printf, outlun,  bin_bmjdarrp(te) , ' ', bin_corrfluxp(te)/plot_corrnorm,  bin_corrfluxerrp(te)/plot_corrnorm, format = '(F0, A,F0, A, F0)'
 ;;       if finite(corrflux(te)) then printf, outlun, bmjd(te), ' ', corrflux(te)/plot_corrnorm, ' ', corrfluxerr(te)/plot_corrnorm, format = '(F0, A,F0, A, F0)'
;;     endfor
;;     free_lun, outlun
     
  endfor



end
