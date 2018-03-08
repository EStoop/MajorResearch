pro plot_pmap_cals
   colorarr = ['deep_pink', 'magenta', 'medium_purple', 'hot_pink', 'light_pink', 'rosy_brown', 'chocolate', 'saddle_brown', 'maroon', 'orange_red', 'dark_orange', 'peach_puff', 'pale_goldenrod','red',  'aquamarine', 'teal', 'steel_blue', 'dodger_blue', 'dark_blue', 'indigo','dark_slate_blue', 'blue_violet', 'purple','dim_grey', 'slate_grey', 'dark_slate_grey', 'khaki', 'light_cyan', 'lavender','gold', 'green_yellow', 'lime', 'green', 'olive_drab', 'pale_green', 'spring_green']

   aorname_pmap = ['IRAC029900/bcd/0045217536', 'IRAC029900/bcd/0045217792', 'IRAC030000/bcd/0045255424', 'IRAC030000/bcd/0045255680','IRAC030100/bcd/0045286656', 'IRAC030100/bcd/0045286912','IRAC030200/bcd/0045381376', 'IRAC030200/bcd/0045381632','IRAC030300/bcd/0045423872', 'IRAC030300/bcd/0045424128','IRAC030400/bcd/0045528576', 'IRAC030400/bcd/0045528832','IRAC030500/bcd/0045567232', 'IRAC030500/bcd/0045567488']
   restore, '/Users/jkrick/irac_warm/pmap/pmap_cals.sav'

;print, 'test', pmap[1].xcen(34), pmap[2].ycen(200)
;now plot them all on the image 
   fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p1s_x4_500x500_0043_120124/pmap_ch2_0p1s_x4_500x500_0043_120124.fits', pmapdata, pmapheader

 c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', aspect_ratio = 1, xrange = [0,500], yrange = [0,500], axis_style = 0)


   for a = 0, n_elements(aorname_pmap) - 1 do begin
      xcen500 = 500.* (pmap[a].xcen - 14.5)
      ycen500 = 500.* (pmap[a].ycen - 14.5)
;      an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
      an = text(xcen500[1], ycen500[1], '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5, color = colorarr[a] , /data)
   endfor

;-------------------------------------------------------------------------
;add the other pcrs peakup exoplanet observations to date:

  restore, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/55cnc_phot_ch2.sav'
  aorname_55cnc = ['r43981312','r43981568','r43981824','r43981056'] ;ch2

  for a = 0, n_elements(aorname_55cnc) - 1 do begin
     xcen500 = 500.* (AOR55cnc[a].xcen - 14.5)
     ycen500 = 500.* (AOR55cnc[a].ycen - 14.5)
;     an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14],/overplot)
     an = text(xcen500[1], ycen500[1],  '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5, color = colorarr[a+14],/data)
  endfor
;-----
  restore, '/Users/jkrick/irac_warm/pcrs_planets/hatp17/hatp17_phot_ch2.sav'
  xcen500 = 500.* (AORhatp17[0].xcen - 14.5)
  ycen500 = 500.* (AORhatp17[0].ycen - 14.5)
;  an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+17],/overplot)
  an = text(xcen500[1], ycen500[1],  '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5, color = colorarr[a+17],/data)
;-----
  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp33/wasp33_phot_ch2.sav'
  aorname_wasp33 = ['r45383424', 'r45384448', 'r45384704'] ;ch2

  for a = 0, n_elements(aorname_wasp33) - 1 do begin
     xcen500 = 500.* (AORwasp33[a].xcen - 14.5)
     ycen500 = 500.* (AORwasp33[a].ycen - 14.5)
;    an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+20],/overplot)
    an = text(xcen500[1], ycen500[1],  '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5, color = colorarr[a+20],/data)
  endfor
;-----
  restore, '/Users/jkrick/irac_warm/pcrs_planets/hatp26/hatp26_phot_ch2.sav'
  xcen500 = 500.* (AORhatp26[0].xcen - 14.5)
  ycen500 = 500.* (AORhatp26[0].ycen - 14.5)
;  an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+23],/overplot)
  an = text(xcen500[1], ycen500[1], '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5, color = colorarr[a+23],/data)
 ;-----
  restore, '/Users/jkrick/irac_warm/pcrs_planets/hd97658/hd97658_phot_ch2.sav'
  xcen500 = 500.* (AORhd97658[0].xcen - 14.5)
  ycen500 = 500.* (AORhd97658[0].ycen - 14.5)
;  an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+23],/overplot)
  an = text(xcen500[1], ycen500[1], '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5, color = colorarr[a+23],/data)
;-----
  restore, '/Users/jkrick/irac_warm/pcrs_planets/hip57274/hip57274_phot_ch2.sav'
  xcen500 = 500.* (AORhip57274[0].xcen - 14.5)
  ycen500 = 500.* (AORhip57274[0].ycen - 14.5)
;  an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+25],/overplot)
  an = text(xcen500[1], ycen500[1], '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5, color = colorarr[a+25],/data)
;-----
  restore, '/Users/jkrick/irac_warm/pcrs_planets/wasp39/wasp39_phot_ch2.sav'
  xcen500 = 500.* (AORwasp39[0].xcen - 14.5)
  ycen500 = 500.* (AORwasp39[0].ycen - 14.5)
;  an = plot(xcen500, ycen500, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+26],/overplot)
  an = text(xcen500[1], ycen500[1],  '$\diamondsuit$', alignment = 0.5, vertical_alignment = 0.5, color = colorarr[a+26],/data)
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;x vs time and y vs time

;  for a = 0, n_elements(aorname_pmap) - 1 do begin
;    if a eq 0 then begin
;       am = plot( (pmap[a].timearr - pmap[a].timearr(0))/60./60., pmap[a].xcen, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],xtitle = 'Time(hrs)', ytitle = 'X pix', yrange = [14.5, 15.5], xrange = [0,12])
;    endif else begin
;       am = plot((pmap[a].timearr - pmap[a].timearr(0))/60./60.,pmap[a].xcen,  '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;    endelse
;
;  endfor
;  
;  for a = 0, n_elements(aorname_55cnc) - 1 do begin
;     am = plot( (AOR55cnc[a].timearr - AOR55cnc[a].timearr(0))/60./60.,AOR55cnc[a].xcen, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14],/overplot)
;  endfor
;
;  am= plot((AORhatp17[0].timearr - AORhatp17[0].timearr(0))/60./60.,AORhatp17[0].xcen,  '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+17],/overplot)
;
;  for a = 0, n_elements(aorname_wasp33) - 1 do begin
;     am = plot( (AORwasp33[a].timearr - AORwasp33[a].timearr(0))/60./60., AORwasp33[a].xcen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+20],/overplot)
;  endfor
;;--------
;
;  for a = 0, n_elements(aorname_pmap) - 1 do begin
;    if a eq 0 then begin
;       am = plot( (pmap[a].timearr - pmap[a].timearr(0))/60./60., pmap[a].ycen, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],xtitle = 'Time(hrs)', ytitle = 'Y pix', yrange = [14.5, 15.5], xrange = [0,12])
;    endif else begin
;       am = plot((pmap[a].timearr - pmap[a].timearr(0))/60./60.,pmap[a].ycen,  '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a],/overplot)
;    endelse
;
;  endfor
;  
;  for a = 0, n_elements(aorname_55cnc) - 1 do begin
;     am = plot( (AOR55cnc[a].timearr - AOR55cnc[a].timearr(0))/60./60., AOR55cnc[a].ycen,'6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+14],/overplot)
;  endfor
;
;  am= plot( (AORhatp17[0].timearr - AORhatp17[0].timearr(0))/60./60.,AORhatp17[0].ycen, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+17],/overplot)
;
;  for a = 0, n_elements(aorname_wasp33) - 1 do begin
;     am = plot((AORwasp33[a].timearr - AORwasp33[a].timearr(0))/60./60., AORwasp33[a].ycen, '6r1s', sym_size = 0.1,   sym_filled = 1, color = colorarr[a+20],/overplot)
;  endfor


end
