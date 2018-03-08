pro delta_superdark, planetname
  
  COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycenp, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr, bin_xfwhm, bin_yfwhm,  bin_corrflux_dp

  if planetname eq 'WASP-14b' or planetname eq 'HD209458' then begin
     colorarr = ['gray', 'gray','gray','gray','gray','burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
  endif else begin
   
     colorarr = ['blue', 'red','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple',  'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua' ,'blue', 'red','black','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo']
;
  endelse

  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  stareaor = planetinfo[planetname, 'stareaor']
  dirname = strcompress(basedir + planetname +'/')     

  ;;-------------------------------------
  ;;read in superdark photometry
  ;;-------------------------------------
                                      
  savefilename = strcompress(dirname + planetname +'_phot_ch2_'+string(2.25)+'_150226_bcdsdcorr.sav',/remove_all) ;
  restore, savefilename
  startaor = stareaor
  for a = startaor, n_elements(aorname) - 1 do begin
          
     ;;check if I should be using pmap corr or not
     ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
     ;;if 20% of the values are correctable than go with the pmap corr 
     if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
     print, 'a, pmapcorr', a, ' ', pmapcorr

     ;;make one data point per AOR, do binning
     bin_level = n_elements([planethash[aorname(a),'flux']]) 
     junkpar = binning_function(a, bin_level, pmapcorr,'2')

     ;use time degraded corrfluxes
     bin_corrfluxp = bin_corrflux_dp

     ;;put all the variables together and save them for later
     if a eq startaor then begin
        bin_sdcorr_phasep = bin_phasep
        bin_sdcorr_corrfluxp = bin_corrfluxp
        bin_sdcorr_corrfluxerrp = bin_corrfluxerrp
        bin_sdcorr_time = bin_timearr
     endif else begin
        bin_sdcorr_phasep = [bin_sdcorr_phasep, bin_phasep]
        bin_sdcorr_corrfluxp = [bin_sdcorr_corrfluxp, bin_corrfluxp]
        bin_sdcorr_corrfluxerrp = [bin_sdcorr_corrfluxerrp, bin_corrfluxerrp]
        bin_sdcorr_time = [bin_sdcorr_time,bin_timearr]
     endelse

  endfor ; all aors

  ;;-------------------------------------
  ;;read in *no*superdark photometry
  ;;-------------------------------------
                                      
  savefilename = strcompress(dirname + planetname +'_phot_ch2_'+string(2.25)+'_150226_bcdnosdcorr.sav',/remove_all) ;
  restore, savefilename

  for a = startaor, n_elements(aorname) - 1 do begin
          
     ;;make one data point per AOR, do binning
     bin_level = n_elements([planethash[aorname(a),'flux']]) 
     junkpar = binning_function(a, bin_level, pmapcorr,'2')

     ;use time degraded corrfluxes
     bin_corrfluxp = bin_corrflux_dp

     ;;put all the variables together and save them for later
     if a eq startaor then begin
        bin_nosdcorr_phasep = bin_phasep
        bin_nosdcorr_corrfluxp = bin_corrfluxp
        bin_nosdcorr_corrfluxerrp = bin_corrfluxerrp
        bin_nosdcorr_time = bin_timearr
     endif else begin
        bin_nosdcorr_phasep = [bin_nosdcorr_phasep, bin_phasep]
        bin_nosdcorr_corrfluxp = [bin_nosdcorr_corrfluxp, bin_corrfluxp]
        bin_nosdcorr_corrfluxerrp = [bin_nosdcorr_corrfluxerrp, bin_corrfluxerrp]
        bin_nosdcorr_time = [bin_nosdcorr_time,bin_timearr]
     endelse

  endfor ; all aors
  normfactor = median(bin_nosdcorr_corrfluxp)
  p = plot(bin_nosdcorr_phasep, bin_nosdcorr_corrfluxp / normfactor, '1s', color = 'green',$
           sym_filled =1, yrange = [0.990, 1.005], xrange = [-0.5, 0.5])
  p2 =  plot(bin_sdcorr_phasep, bin_sdcorr_corrfluxp / normfactor, '1s', color = 'blue',$
           sym_filled =1, overplot = p)

  ;;-------------------------------------
  ;;look at deltacorrflux
  ;;-------------------------------------
  delta = (bin_nosdcorr_corrfluxp/normfactor) - (bin_sdcorr_corrfluxp / normfactor)
  plothist, delta, xhist, yhist, bin = 5E-5, /noplot,/noprint
  ph = barplot(xhist, yhist,  xtitle = 'Delta Corrflux', ytitle = 'Number', fill_color = 'sky_blue', $
               title = planetname, xrange = [-1E-3, 1E-3] ) 
;  print, 'n delta', n_elements(delta), n_elements(aorname)
;  for a = startaor, n_elements(aorname) - 1 do print, a, ' ',  aorname(a), ' ', colorarr(a), ' ', delta(a), format = '(I10, A,A, A, A, A, G)'
  

end
