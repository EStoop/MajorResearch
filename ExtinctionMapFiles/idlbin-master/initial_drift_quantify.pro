pro initial_drift_quantify, bin_level, make_plot = make_plot

  COMMON bin_block, aorname, planethash, bin_xcen, bin_ycen, bin_bkgd, bin_flux, bin_fluxerr,  bin_timearr, bin_phase, bin_ncorr,bin_np, bin_npcent, bin_xcenp, bin_ycep, bin_bkgdp, bin_fluxp, bin_fluxerrp,  bin_corrfluxp,  bin_timearrp, bin_corrfluxerrp,  bin_phasep,  bin_ncorrp, bin_nparrp, bin_npcentarrp, bin_bmjdarr, bin_xfwhm, bin_yfwhm,  bin_corrflux_dp


  aorname = [ 'r50678528',  'r50665984',  'r51824640',  'r51816192',  'r51836416',  'r51833088', 'r51843584',  'r51819776',   'r51832320',  'r51820544',  'r51844096',  'r51832576',  'r51837184',  'r51828480',  'r51840000',  'r51823360', 'r51840256',  'r51841536',  'r51842816',  'r51823872',  'r51831040',  'r51840512',  'r51834880',  'r51835136',  'r51820800',  'r51835904',  'r51827968',  'r51816704',  'r51834368',  'r51829504', 'r51816448',  'r51826176',  'r51822848',  'r51832064',  'r51838976',  'r51815936',   'r51835648',  'r51828224',  'r51831552',  'r51829248',  'r52364288', 'r52364544',  'r52364800',  'r52365312',  'r50651392',  'r50653440',   'r50678016',  'r50676480',   'r50653184',  'r50652928',   'r51831296',  'r51818496', 'r50654976',  'r50652672', 'r51827200',  'r51842304',  'r50655232',  'r50652416', 'r51818752',  'r51830016',  'r51817728',  'r51816960',  'r47037184',  'r47029248',  'r51817216',  'r51821568', 'r51837440',  'r51830528',  'r47030784',  'r47053824',  'r51837952',  'r51826688',    'r51836672',  'r51833600',  'r51833344',  'r51838464',  'r51825152',  'r51825920',  'r51832832',  'r51838720', 'r51837696',  'r51823104',  'r51834624',  'r51826432',   'r51819520',  'r51824384', 'r53521152',  'r53521664',  'r51843072', 'r51821056', 'r53522176',  'r53521920',  'r51836160',  'r51827456',    'r53522688',  'r53522432',    'r54315776',  'r54315520','r50651904',  'r50650112',  'r53523200',  'r53522944',   'r50653696',  'r50655488'];, 'r51843328','r51839488'] ;'r51822080',  'r51834112','r49711360',  'r50653184',  'r50652928',  'r50654976',  'r50652672',  'r50655232',  'r50652416',  'r51883264','r51884032', 'r51841024', 'r51842048', 'r50651136',  'r50656000','r50651648', 'r50655744',    'r50649856',  'r50656768',  'r50650624',  'r50656256', 
;wasp-64 both AOR  'r51830784', 'r51842560',

colorarr = ['blue', 'red','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple',  'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua' ,'blue', 'red','black','green','black','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo']
;
  ;;read in the saved centroids from all AORs
  restore, '/Users/jkrick/external/initial_drift/centroids.sav'
  extra={ sym_size: 0.1, sym_filled: 1, xrange: [0,6], yrange: [0.975, 1.01]};[14.5, 20.0]};[20.0, 26.5]}

  ;;read in the stored pitch angle table from Carl
  restore, '/Users/jkrick/external/initial_drift/ctab_data.6378.ecsv.sav'

  spreadfactor = 0.1
  startaor = 1; 54 ;67 ;44
  stopaor =  n_elements(aorname) -1 ; 53; 104; 113; 89; 66
  delta_x = fltarr(n_elements(aorname))
  delta_y = delta_x
  deltapitcharr = fltarr(n_elements(aorname)/2.) 
  c = 0
  ai = 0 ; only the inital AORs
  for a = startaor, stopaor, 2 do begin; only the long ones
     print, 'working on ', aorname(a), ' ', a
     xarr = planethash[aorname(a), 'xcen']
     yarr = planethash[aorname(a), 'ycen']
     timearr = planethash[aorname(a), 'timearr'] 
     bmjdarr = planethash[aorname(a), 'bmjdarr'] 

     ;;binning
     junkpar = binning_function_id(a, bin_level, 1,'2')

     starttime = where((bin_timearr - bin_timearr(0)) / 60./60. ge 0);  .6
     print, 'n_ele', n_elements(bin_xcen), n_elements(starttime)
     xarr_n = bin_xcen(starttime)
     yarr_n = bin_ycen(starttime)
     timearr_n = bin_timearr(starttime)

     xarr_n = xarr_n / xarr_n(0)
     yarr_n = yarr_n / yarr_n(0)
     timearr_n = (timearr_n - timearr_n(0))/60./60.

     ;;try to figure out what the delta pitch angle is for this
     ;;observation
     b = a - 1                  ; want to look at the previous AORs start SCLK time.
     pretimearr = planethash[aorname(b), 'timearr'] 
     ;;find that same time in the pitch angle table.
     match = min(where(time gt pretimearr(0) - 10 and time lt pretimearr(0) + 10))
     pitchtime = time(match) ; should be the start time
     deltapitch = pitchb(match) - pitchb(match - 360) 
     print, 'deltapitch', deltapitch
     if deltapitch lt 12 then plotcolor = 'blue'
     if deltapitch lt 24 and deltapitch ge 12 then plotcolor = 'cyan'
     if deltapitch ge 24 then plotcolor = 'red'
     deltapitcharr[ai] = deltapitch
     ai++

     pnx = plot(timearr_n, xarr_n , '1', xtitle = 'Time (Hrs)', title = 'Current Obs =30min pre AOR',$
               ytitle = 'Normalized X position', _extra = extra, color = plotcolor, overplot = pnx)  ; colorarr[a]
     pny = plot(timearr_n, yarr_n, '1', xtitle = 'Time (Hrs)', title = 'Current Obs =30min pre AOR',$
               ytitle = 'Normalized Y position', _extra = extra, color = plotcolor, overplot = pny) ; colorarr[a]


     ;(bin_timearr - bin_timearr(0)) / 60./60.

     ;;normalize to the first value and then plot
     ;;xarr_n = xarr / xarr(0)
     ;;yarr_n = yarr/ yarr(0)
     ;;print, 'normalized', xarr_n[0:10], xarr_n[0:10]
     ;;try gaussian smoothing
     ;;xarr_gauss = GAUSS_SMOOTH(xarr_n,5,/nan) 
     ;;yarr_gauss = GAUSS_SMOOTH(yarr_n,5,/nan)
     ;;print, 'gauss', xarr_gauss[10:20], xarr_gauss[10:20]

     ;;pnx = plot((timearr - timearr(0)) / 60./60., xarr_gauss , '1s', xtitle = 'Time (Hrs)', $
     ;;          ytitle = 'Normalized X position', _extra = extra, color = plotcolor, overplot = pnx)
     ;;pny = plot((timearr - timearr(0)) / 60./60., yarr_gauss, '1s', xtitle = 'Time (Hrs)', $
     ;;          ytitle = 'Normalized Y position', _extra = extra, color = plotcolor, overplot = pny)


     ;;track the first 2 minutes and then 2 minutes at 4 hours
     howlong = 30; 240  ;after 4 hours
     firsttime = where((timearr - timearr(0)) / 60. lt 5 and (timearr - timearr(0)) / 60. ge 0,nfirst)
     lasttime = where((timearr - timearr(0)) / 60. lt howlong and (timearr - timearr(0)) / 60. ge (howlong - 5), nlast)
     if nlast gt 0 then begin  ;;otherwise less than 4 hours in duration
        x_mean_first = mean(xarr(firsttime))
        x_mean_last = mean(xarr(lasttime))
        y_mean_first = mean(yarr(firsttime))
        y_mean_last = mean(yarr(lasttime))

        delta_x[c] = x_mean_last - x_mean_first
        delta_y[c] = y_mean_last - y_mean_first
        c++
     endif

     ;;plot them all together at first.
     if keyword_set(make_plot) then begin

        
        if a eq startaor then begin ; plot the first one regardless
           p_timearr = (timearr - timearr(0)) / 60./60.
           
           ;;try gaussian smoothing
           xarr_gauss = GAUSS_SMOOTH(xarr,5,/nan) + a*spreadfactor
           yarr_gauss = GAUSS_SMOOTH(yarr,5,/nan)+ a*spreadfactor
        
           px = plot(p_timearr, xarr_gauss, '1s', xtitle = 'Time (Hrs)', ytitle = 'X position', _extra = extra, color = plotcolor)
           py = plot(p_timearr, yarr_gauss, '1s', xtitle = 'Time (Hrs)', ytitle = 'Y position', _extra = extra, color = plotcolor)
           
        endif else begin
           print, colorarr(a), ' ',colorarr(a-1), planethash[aorname(a), 'ra'], $
                  planethash[aorname(a-1), 'ra'], format = '( A, A, A, F, F)'
           if abs(planethash[aorname(a), 'ra'] - planethash[aorname(a-1), 'ra']) lt 0.01 then begin
              print, 'same ra as previous AOR', colorarr(a - 1)
              ;;plot them with the same color
              t0 = planethash[aorname(a-1), 'timearr'] 
              x0 = planethash[aorname(a-1), 'xcen'] 
              y0 = planethash[aorname(a-1), 'ycen'] 
              p_timearr = [p_timearr, (timearr - t0(0)) / 60./60.]
              ;;try gaussian smoothing
              xarr_gauss = [xarr_gauss, GAUSS_SMOOTH(xarr,5,/nan) + (a - 1)*spreadfactor]
              yarr_gauss = [yarr_gauss, GAUSS_SMOOTH(yarr,5,/nan) + (a - 1)*spreadfactor]
              
              px = plot(p_timearr, xarr_gauss , '1s',  _extra = extra, color = plotcolor,overplot = px) ; colorarr[a-1]
              py = plot(p_timearr, yarr_gauss , '1s',  _extra = extra, color = plotcolor,overplot = py) ; colorarr[a-1]
           endif else begin
              ;;plot them with a different color
              p_timearr = (timearr - timearr(0)) / 60./60.
              
              ;;try gaussian smoothing
              xarr_gauss = GAUSS_SMOOTH(xarr,5,/nan) + a*spreadfactor
              yarr_gauss = GAUSS_SMOOTH(yarr,5,/nan) + a*spreadfactor
;           print, 'xarr_gauss new AOR', xarr_gauss[100:110]
;
              px = plot(p_timearr, xarr_gauss, '1s', _extra = extra, color = plotcolor,overplot = px)
              py = plot(p_timearr, yarr_gauss, '1s', _extra = extra, color = plotcolor,overplot = py)
           endelse
        endelse
     endif                      ; keyword set make_plotxs
     
  endfor                        ; for each AOR

  ;;annotate the normalized plots
  ;;pnx = plot(findgen(7), fltarr(7) + 1.00 - 0.06, overplot = pnx)
  pny = plot(findgen(3) + 4.0, fltarr(3) + 1.00 - 0.006, overplot = pny, color = 'black')
  pny = plot(findgen(3) + 4.0, fltarr(3) + 1.00 + 0.006, overplot = pny, color = 'black')
  pny = plot([4.0, 4.0], [0.9, 0.994], overplot = pny, color = 'black')
  pny = plot([4.0, 4.0], [1.006, 2.0], overplot = pny, color = 'black')

  pnx = plot(findgen(3) + 4.0, fltarr(3) + 1.00 - 0.006, overplot = pnx, color = 'black')
  pnx = plot(findgen(3) + 4.0, fltarr(3) + 1.00 + 0.006, overplot = pnx, color = 'black')
  pnx = plot([4.0, 4.0], [0.9, 0.994], overplot = pnx, color = 'black')
  pnx = plot([4.0, 4.0], [1.006, 2.0], overplot = pnx, color = 'black')

  print, 'there are ', c, ' AORs'
  delta_x = delta_x[0:c-1]
  delta_y = delta_y[0:c-1]
  plothist, delta_x, xxhist, xyhist, bin = 0.06, /noprint, /noplot
  pdx = barplot(xxhist, xyhist,  fill_color = 'blue', xtitle = 'Delta X position after 4 hours',ytitle = 'Number', xrange = [-0.35, 0.15])

  plothist, delta_y, yxhist, yyhist, bin = 0.06, /noprint, /noplot
  pdx = barplot(yxhist, yyhist,  fill_color = 'blue', xtitle = 'Delta Y position after 4 hours',ytitle = 'Number', xrange = [-0.35, 0.15])

  ;;plot delta pitch vs. delta position
  dpx = plot( deltapitcharr,delta_x, '1o',ytitle = 'Delta X position after 0.5 hours', xtitle = 'Pitch angle change  30m',$
            sym_filled = 1)
  dpy = plot(deltapitcharr, delta_y, '1o', ytitle = 'Delta Y position after 0.5 hours', xtitle = 'Pitch angle change 30m',$
            sym_filled = 1)

  print, 'deltapitcharr', deltapitcharr

; an interesting set to plot on its own:
;  aorname = [ 'r50651392',  'r50653440',  'r50653184',  'r50652928',  'r50654976',  'r50652672',  'r50655232',  'r50652416']
;  for a = 0, n_elements(aorname) - 1 do begin
;     xarr = planethash[aorname(a), 'xcen']
;     yarr = planethash[aorname(a), 'ycen']
;     timearr = planethash[aorname(a), 'timearr'] 
;     help, timearr
;     help, xarr
;     help, yarr
;     p1 = plot(timearr, xarr, '1s',sym_size = 0.1,sym_filled = 1, overplot = p1)
;     p2 = plot(timearr, yarr, '1s', _extra, overplot = p2)
;  endfor


end
