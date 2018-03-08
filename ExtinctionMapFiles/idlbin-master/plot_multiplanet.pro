pro plot_multiplanet

  chname = '2'
  planetname = ['WASP-13b', 'WASP-38b', 'HAT-P-22', 'WASP-15b', 'WASP-16b', 'WASP-62b']
  apradius = [2.25, 2.25, 2.25, 1.75, 2.25, 2.5]
  colorarr = ['blue', 'red','black','lime green','aqua','purple' ]
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_121120.fits', pmapdata, pmapheader
;
;  chname = '1'
;  planetname = [ 'HAT-P-22', 'WASP-16b', 'WASP-62b']
;  apradius = [1.75, 2.0, 2.25]
;  colorarr = ['dark orange', 'pale green','magenta' ]
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch1_500x500_0043_120409.fits', pmapdata, pmapheader

;for debugging: skip some AORs
  startaor = 1
  stopaor =  1

  xsweet = 15.12
  ysweet = 15.00
  planetinfo = create_planetinfo()

  c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', aspect_ratio = 1, xrange = [0,500], yrange = [0,500],xshowtext = 0, yshowtext = 0)
  undefine, pmapdata

  for p = 0, n_elements(planetname) - 1 do begin
     print, 'working on planet', planetname(p)
;run code to read in all the input planet parameters
     aorname= planetinfo[planetname(p), 'aorname_ch2'] 
;     aorname = planetinfo[planetname(p), 'aorname_ch1'] 
     print, 'aorname', aorname(0)
     basedir = planetinfo[planetname(p), 'basedir']
     dirname = strcompress(basedir + planetname(p) +'/') 
     savefilename = strcompress(dirname + planetname(p) +'_phot_ch'+chname+'_'+string(apradius(p))+'.sav',/remove_all)
     print, 'restoring ', savefilename
     restore, savefilename
     

;  print, 'naor', n_elements(aorname)
     for a = startaor ,  stopaor  do begin
;     print, 'testing aors', a, colorarr[a]
        print, 'aorname ', aorname(a)
        xcen500 = 500.* ((planethash[aorname(a),'xcen']) - 14.5)
        ycen500 = 500.* ((planethash[aorname(a),'ycen']) - 14.5)
        an = plot(xcen500, ycen500, '1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[p],/overplot)
        
       endfor
                                ;need some way of cleaning up to let this keep going through multiple planets
;     print, 'memory', memory()
;     undefine, planethash, xcen500, ycen500
;     print, 'memory', memory()
;     help
  endfor                        ; each planetname

;legend
  an = text(20, 450, 'WASP-13b', /data, color = colorarr[0], /overplot, font_style = 1)
  an = text(20, 420, 'WASP-15b', /data, color = colorarr[3], /overplot, font_style = 1)
  an = text(20, 390, 'WASP-16b', /data, color = colorarr[4], /overplot, font_style = 1)
  an = text(20, 360, 'WASP-38b', /data, color = colorarr[1], /overplot, font_style = 1)
  an = text(20, 330, 'WASP-62b', /data, color = colorarr[5], /overplot, font_style = 1)
  an = text(20, 300, 'HAT-P-22', /data, color = colorarr[2], /overplot, font_style = 1)
  

;  an = text(20, 450, 'WASP-16b', /data, color = colorarr[1], /overplot, font_style = 1)
;  an = text(20, 420, 'WASP-62b', /data, color = colorarr[2], /overplot, font_style = 1)
;  an = text(20, 390, 'HAT-P-22', /data, color = colorarr[0], /overplot, font_style = 1)


  an.save, '/Users/jkrick/irac_warm/pcrs_planets/secondaries_paper/position_ch'+chname+'.png'

end
