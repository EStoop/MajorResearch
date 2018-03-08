pro plot_peakpix

colorarr = ['blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'tomato', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple','light_sea_green', 'teal', 'cadet_blue', 'aquamarine', 'dark_turquoise', 'aqua','blue', 'red','black','green','grey','purple', 'deep_pink', 'magenta', 'medium_purple', 'orchid', 'thistle', 'pink', 'orange_red', 'rosy_brown',  'chocolate', 'saddle_brown', 'maroon', 'dark_orange', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'khaki', 'tomato','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green','blue', 'red','deep_pink', 'magenta', 'medium_purple']

  planetname = 'upsandb'
  chname = '2'
  apradius = 2.25
  planetinfo = create_planetinfo()
  aorname= planetinfo[planetname, 'aorname_ch2'] 
  basedir = planetinfo[planetname, 'basedir']
  dirname = strcompress(basedir + planetname +'/')
  savename = strcompress(dirname + planetname+'_phot_ch'+chname+'_'+string(apradius)+'_150226.sav',/remove_all)
  
  restore, savename
  xsweet = 15.12
  ysweet = 15.1                 ; 15.0

  for a = 0, n_elements(aorname) -1 do begin
     xcen= planethash[aorname(a),'xcen']
     ycen= planethash[aorname(a),'ycen']
     flux = planethash[aorname(a),'flux']
     fluxerr = planethash[aorname(a),'fluxerr']
     corrflux = planethash[aorname(a),'corrflux']
     peakpixDN = planethash[aorname(a),'peakpixDN']
     np = planethash[aorname(a),'npcentroids']

     xdist = (xcen - xsweet)
     ydist = (ycen- ysweet)
     dsweet = sqrt(xdist^2 + ydist^2)
     gain = flux / corrflux  ; back out what the gain is at the position of the observation

     p = plot(dsweet, peakpixDN, '1s', sym_size = 0.1,   sym_filled = 1, xtitle = 'Dsweet', ytitle = 'PeakPix (DN)',$
              xrange = [0, 0.4], color = colorarr(a), overplot = p, yrange = [4000, 12000])

 ;    p2 = plot(gain, flux / mean(flux,/nan), '1s', sym_size = 0.1,   sym_filled = 1, xtitle = 'Gain', ytitle = 'Normalized Flux ',$
 ;              color = colorarr(a), overplot = p2, xrange = [0.995, 1.02], yrange = [0.985, 1.015]) ;yrange = [8000, 16000], 
  endfor  ; for each aor

  
end
