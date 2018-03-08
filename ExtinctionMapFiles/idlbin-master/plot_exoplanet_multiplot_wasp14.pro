pro plot_exoplanet_multiplot, planetname, bin_level, apradius, chname,  timeplot = timeplot
;example call plot_exoplanet, 'wasp15', 2*63L

  COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycenp, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr, bin_xfwhm, bin_yfwhm

  if planetname eq 'WASP-14b' then begin
     colorarr = ['gray', 'gray','gray','gray','gray','burlywood','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']
   
  endif else begin
     colorarr = ['black','sandy_brown', 'rosy_brown','saddle_brown', 'brown', 'maroon', 'firebrick', 'crimson', 'salmon', 'orange_red', 'dark_orange', 'orange', 'goldenrod', 'gold', 'yellow','khaki', 'green_yellow', 'lime', 'lime_green', 'green', 'dark_green', 'olive', 'olive_drab', 'sea_green', 'light_green', 'medium_spring_green', 'medium_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'cyan', 'light_sky_blue', 'dodger_blue', 'steel_blue', 'blue', 'dark_blue', 'indigo', 'medium_slate_blue', 'purple', 'blue_violet', 'dark_orchid', 'orchid', 'pink', 'pale_violet_red', 'deep_pink', 'fuchsia']

;colorarr = ['red','blue','green','grey','purple', 'deep_pink','fuchsia', 'magenta', 'medium_purple','medium_orchid', 'orchid', 'violet', 'plum', 'thistle', 'pink', 'orange_red', 'light_pink', 'rosy_brown','pale_violet_red',  'chocolate', 'saddle_brown', 'maroon', 'hot_pink', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red', 'deep_pink','fuchsia', 'magenta', 'medium_purple','medium_orchid', 'orchid', 'violet', 'plum', 'thistle', 'pink', 'orange_red', 'light_pink', 'rosy_brown','pale_violet_red',  'chocolate', 'saddle_brown', 'maroon', 'hot_pink', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua' ]
;
  endelse


;run code to read in all the input planet parameters
  planetinfo = create_planetinfo()
  if chname eq '2' then aorname= planetinfo[planetname, 'aorname_ch2'] else aorname = planetinfo[planetname, 'aorname_ch1'] 
  basedir = planetinfo[planetname, 'basedir']
  period =  planetinfo[planetname, 'period']
  stareaor = planetinfo[planetname, 'stareaor']
;for debugging: skip some AORs
  startaor = 5; 0                  ;5
  stopaor =   n_elements(aorname) - 1
  print, 'stopaor', stopaor
;  plot_norm= planetinfo[planetname, 'plot_norm']
;  plot_corrnorm = planetinfo[planetname, 'plot_corrnorm']
  
  dirname = strcompress(basedir + planetname +'/');+'/hybrid_pmap_nn/')
  if chname eq '2' then savefilename = strcompress(dirname + planetname +'_phot_ch2_2.25000_160126.sav',/remove_all) 
  ;;savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_150723.sav',/remove_all)
  if chname eq '1' then savefilename = strcompress(dirname + planetname +'_phot_ch'+chname+'_'+string(apradius)+'_140716.sav',/remove_all)
  print, 'restoring ', savefilename
  restore, savefilename
  print, 'aorname', aorname(0)
  time_0 = (planethash[aorname(startaor),'timearr']) 
  time_0 = time_0(0)
  print, 'time_0', time_0
;  plothist, planethash[aorname(0),'npcentroids'], xhist, yhist, /noplot, bin = 0.01
;  testplot = plot(xhist, yhist,  xtitle = 'NP', ytitle = 'Number', thick =3, color = 'blue')

  ;print, 'testing phase', (planethash[aorname(0),'phase'] )[0:10], (planethash[aorname(0),'phase'] )[600:610]

 xsweet = 15.12
 ysweet = 15.00
 

;setup some arrays to hold binned values
  bin_corrfluxfinal = fltarr(n_elements(aorname))
  bin_phasefinal = fltarr(n_elements(aorname))
  bin_corrfluxerrfinal = fltarr(n_elements(aorname))

  mean_corrfluxarr = fltarr(stopaor - startaor + 1,/nozero)
  mean_dsweetarr = mean_corrfluxarr
  stddev_corrfluxarr = mean_corrfluxarr
  stddev_dsweetarr = mean_corrfluxarr
 
  for a = startaor,stopaor do begin
     print, '------------------------------------------------------'
     print, 'working on AOR', a, '   ', aorname(a)

     ;check if I should be using pmap corr or not
     ncorr = where(finite([ planethash[aorname(a),'corrflux']]) gt 0, corrcount,/L64)
     ;if 20% of the values are correctable than go with the pmap corr 
     print, '0.2nflux, ncorr, ', 0.2*n_elements([planethash[aorname(a),'flux']]), corrcount
     if corrcount gt 0.2*n_elements([planethash[aorname(a),'flux']]) then pmapcorr = 1 else pmapcorr = 0
     print, 'pmapcorr', pmapcorr
;;     peakpixDN = planethash[aorname(a), 'peakpixDN']
;;     print, 'median peakpixDN', median(peakpixDN)

; bin within each AOR only, don't want to bin across boundaries.
     junkpar = binning_function(a, bin_level, pmapcorr)
;     print, 'testing phase', (planethash[aorname(a),'phase'] )[0:10], (planethash[aorname(a),'phase'] )[600:610]
;     print, 'testing x', bin_xcen[0:10]
;     print, 'testing y', bin_ycen[0:10]
;     print, 'testing npcent', bin_npcent[0:10]
;------------------------------------------------------
     ;possible comparison statistic
     ;what is the distribution of standard deviations among the corrfluxes?

     meanclip_jk, planethash[aorname(a),'corrflux'] , mean_corrflux, stddev_corrflux
     stddev_corrfluxarr[a- startaor] = stddev_corrflux
;------------------------------------------------------
;now try plotting
     xcen = planethash[aorname(a),'xcen']
     ycen = planethash[aorname(a),'ycen']
     flux = planethash[aorname(a),'flux']
;     print, 'first', xcen[350:360], ycen[350:360], flux[350:360]
;     print, 'second', xcen[6395:6405], ycen[6395:6405], flux[6395:6405]
;------------------------------------------------------
     print, 'n_elements time', n_elements(bin_timearr), n_elements(bin_xcen), n_elements(bin_phase)
     if keyword_set(timeplot) then begin ;make the plot as a function of time
        print, 'making plot as a function of time'
        plotx = (bin_timearr-  time_0)/60./60./24.
      
        plotx2 =(bin_timearr - bin_timearr(0))/60./60.
;;        help, bin_timearr
;;        print,bin_timearr
        if planetname eq 'WASP-14b' then begin
;;           plotx = plotx2
;;           setxrange = [0,359]
        endif
        
        if plotx(0) gt 300. then begin
           print, 'large plotx', plotx(0)
           plotx = plotx - (7800./24.)
           setxrange = [0, 900/24.]
        endif
        ending = 'time'
        
     endif else begin
        plotx =  bin_phase ;;-1.09
        setxrange =    [-0.5,0.5];[0.5,0.65]; [0.45,0.67];
        ending = 'phase'
     endelse
     
     
     extra={xthick: 2, ythick:2, margin: 0.2, sym_size: 0.2,   sym_filled: 1, xstyle:1}
     corrnormoffset = 0         ;0.02
     corroffset = 0.001
     setynormfluxrange = [0.98, 1.01] ;[0.97, 1.005]
;     setxrange = [-1, 18]

     ;;normalize each AOR seperately
     ;plot_norm =  median(bin_flux)
     ;plot_corrnorm = mean(bin_corrfluxp,/nan)
     ;bkgd_norm =  mean(bin_bkgd,/nan)

     if a eq startaor then begin ; for the first AOR
                                ;set the normalization values from the
                                ;medians of the first AOR...is close enough
        plot_norm =  median(bin_flux)
        plot_corrnorm = mean(bin_corrfluxp,/nan)
        bkgd_norm =  mean(bin_bkgd,/nan)
;           print, 'bin_xcen', bin_xcen
        print, 'plot_corrnorm', plot_corrnorm, mean(bin_corrfluxp)
        pp = plot(plotx, bin_xcen, '1s',  $ ;title = planetname, $
                  color = colorarr[a], ytitle = 'X position', position = [0.2,0.78,0.9,0.91], ytickinterval = 0.1, $
                  xshowtext = 0, ytickformat = '(F10.2)', dimensions = [600, 900], _extra = extra, yminor = 0,$
                  xrange = setxrange, yrange = [15.0,15.4]);, title = planetname ) ;, $, ymajor = 4 [22.85, 23.2]
        
                                ;turn off refreshing to make this quicker hopefully
;           pp.Refresh, /Disable
        
        
        pq= plot(plotx, bin_ycen, '1s',  color = colorarr[a], $
                 ytitle = 'Y position',  position = [0.2, 0.64, 0.9, 0.77],/current,  ytickinterval = 0.1,$
                 xshowtext = 0,ytickformat = '(F10.2)', _extra = extra, yminor = 0,$
                 xrange = setxrange, yrange = [14.9,15.5]) ;, $, title = planetname , ymajor = 4
                                ;xrange = setxrange); [230.85,231.15]
        
         pxy = plot(plotx, bin_xfwhm, '1s', color = 'blue', $
                   ytitle = 'X & Y FWHM',  position = [0.2, 0.50, 0.9, 0.63], /current, $
                   xshowtext = 0,ytickformat = '(F10.2)', _extra = extra, ytickinterval = 0.1, yminor = 0,$
                   xrange = setxrange, yrange = [1.9, 2.2]) ;
        pxy = plot(plotx, bin_yfwhm, '1s', color = 'black', $
                   overplot = pxy, _extra = extra,$
                   xrange = setxrange) ;
       
         
        ps= plot(plotx, bin_npcent, '1s', color = colorarr[a], $
                 ytitle = 'Noise Pixel',  position = [0.2, 0.36, 0.9, 0.49] , /current, $
                 xshowtext = 0,ytickformat = '(F10.1)', _extra = extra, ytickinterval = 1.0, yminor = 0, yrange = [4.0, 7.0],$
                 xrange = setxrange) ;,$ title = planetname,, ymajor = 4;xrange = setxrange) position =  [0.2, 0.50, 0.9, 0.63]
        
        pt = plot(plotx, bin_bkgd/ bkgd_norm, '1s' , color = colorarr[a], $
                  ytitle = 'Norm. Bkgd',  margin = 0.2,position = [0.2, 0.22, 0.9, 0.35] , /current, xshowtext = 0,$
                  ytickformat = '(F10.2)', _extra = extra, ytickinterval = .2, yminor = 0,$ 
                  xrange = setxrange, yrange = [0.70, 1.4])  ;, $ title = planetname,ymajor = 4, position = [0.2, 0.36, 0.9, 0.49]
        
        
        pr = plot(plotx, bin_flux/plot_norm, '1s',  $
                  color = colorarr[a],   ytitle = 'Norm. Flux', xtitle = 'Phase',$ 
                   position =[0.2, 0.08, 0.9, 0.21], /current, _extra = extra, ytickinterval = 0.01, yminor = 0,$
                  xrange = setxrange, yrange = [0.98, 1.01]) ;position =   [0.2, 0.22, 0.9, 0.35]
        
 
     endif                      ; if a = startaor
     
                                ;  test_med(a) = median(bin_corrfluxp)
     print,'median bin_corrfluxp ',  median(bin_corrfluxp), median(bin_flux)
     if (a gt 0) and (a le stareaor) then begin
        print, 'inside a gt 0 a le stareaor', a
;           pp.window.SetCurrent
        pp = plot(plotx, bin_xcen, '1s', sym_size = 0.2,   sym_filled = 1,color = colorarr[a],  overplot = pp,/current)
;           pq.window.SetCurrent
        pq = plot(plotx, bin_ycen, '1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a], overplot = pq,/current)
;           pr.window.SetCurrent
        pr = plot(plotx, bin_flux/plot_norm , '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], overplot = pr,/current)
                                ;          if pmapcorr eq 1 then begin
;              pr = plot((bin_timearr - time_0)/60./60., (bin_corrfluxp/plot_corrnorm)-corrnormoffset, /overplot, 's1', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/current,overplot = pr)
;           endif
        
;           ps.window.SetCurrent
        ps = plot(plotx, bin_npcent, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], overplot = ps,/current) ;, xtitle = 'Time(hrs)', ytitle = 'Normalized Flux',) 
        
;           pt.window.SetCurrent
        pt = plot(plotx, bin_bkgd/ bkgd_norm, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], overplot = pt,/current) 
        
;;        pxy = plot(plotx, bin_xfwhm, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], overplot = pxy,/current)
;;        pxy = plot(plotx, bin_yfwhm, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], overplot = pxy,/current)
        
        
     endif
     
     
     
     if a gt stareaor then begin
        print, 'inside a gt stareaor', a
        pp = plot(plotx, bin_xcen, '1s', sym_size = 0.2,   sym_filled = 1,color = colorarr[a],$
                  overplot = pp, layout = [1,5,1])
        
        
        pq = plot(plotx, bin_ycen, '1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],$
                  overplot = pq, layout = [1,5,2])
        
        pr = plot(plotx, bin_flux/(plot_norm) , '1s', sym_size = 0.2,   sym_filled = 1,  $
                  color = colorarr[a], overplot = pr)
        
        
        ps = plot(plotx, bin_npcent, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a],$
                  overplot = ps) 
        
        pt = plot(plotx, bin_bkgd/bkgd_norm, '1s', sym_size = 0.2,   sym_filled = 1,  color = colorarr[a], $
                  overplot = pt) 
        pxy = plot(plotx, bin_xfwhm, '1s',  sym_size = 0.2,   sym_filled = 1, color ='blue', overplot = pxy)    ;
        pxy = plot(plotx, bin_yfwhm, '1s', sym_size = 0.2,   sym_filled = 1, color = 'black', overplot = pxy)    ;
        
        
     endif
     
;        pp.Refresh
;        pp.save, dirname +'binxcen_time_ch'+chname+'.eps'
;        pq.save, dirname +'binycen_time_ch'+chname+'.eps'
;        pt.save, dirname +'binbkg_time_ch'+chname+'.eps'
;        ps.save, dirname + 'binnp_time_ch'+chname+'.eps'
;        pr.save, dirname + 'binflux_time_ch'+chname+'.eps'

;       

 
  endfor                        ;binning for each AOR

  if keyword_set(timeplot) then begin
     if planetname eq 'WASP-14b' then begin
        pp.xrange = [0, 800/24.]
        pq.xrange = [0, 800/24.]
         pt.xrange = [0, 800/24.]
         ps.xrange = [0, 800/24.]
         pr.xrange = [0, 800/24.]
         pxy.xrange = [0, 800/24.]
         pp = plot([550/24., 520/24.], [min(pp.yrange), max(pp.yrange)], linestyle = 4, overplot = pp)
         pq = plot([550/24., 520/24.], [min(pq.yrange), max(pq.yrange)], linestyle = 4, overplot = pq)
         pt = plot([550/24., 520/24.], [min(pt.yrange), max(pt.yrange)], linestyle = 4, overplot = pt)
         ps = plot([550/24., 520/24.], [min(ps.yrange), max(ps.yrange)], linestyle = 4, overplot = ps)
         pr = plot([550/24., 520/24.], [min(pr.yrange), max(pr.yrange)], linestyle = 4, overplot = pr)
         pxy = plot([550/24., 520/24.], [min(pxy.yrange), max(pxy.yrange)], linestyle = 4, overplot = pxy)
         pr.xtickname = ['0', '5', '10', '15','20', '350', '355'] ; account for the time shifting
     endif
     
;;      pr.xtitle = 'Time(Days)'
pr.xtitle = 'Time(hours)'
  endif

 ;; pp.save, dirname +'multiplot_ch'+chname+'_'+ending+'.eps'
end

