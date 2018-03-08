pro hd158460_snap
  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' PINK',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'TAN', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'WHEAT'    , 'WHITE_SMOKE'   , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]

  aorname = ['0044497408','0044497152','0044497664','0044497920','0044498176','0044498432','0044498688','0044498944','0044499200','0044499456']

m = pp_multiplot(multi_layout=[n_elements(aorname),2],global_xtitle='Test X axis title',global_ytitle='Test Y axis title')


  for a = 0, n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC029100/bcd/' + string(aorname(a)) 
     CD, dir                    ; change directories to the correct AOR directory
     command  =  " ls *bcd_fp.fits > /Users/jkrick/irac_warm/snapshots/cbcdlist.txt"
     spawn, command
     
     readcol,'/Users/jkrick/irac_warm/snapshots/cbcdlist.txt',fitsname, format = 'A', /silent
     yarr = fltarr(n_elements(fitsname)*63)
     xarr = fltarr(n_elements(fitsname)*63)
     timearr = fltarr(n_elements(fitsname)*63)
     fluxarr = fltarr(n_elements(fitsname)*63)
     fluxerrarr = fltarr(n_elements(fitsname)*63)
     corrfluxarr = fltarr(n_elements(fitsname)*63)
     for i =0.D, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
                                ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

        for nt = 0, 63 - 1 do begin
           timearr[i*63+nt] = ((sclk_obs - sclk_0) + 0.4*nt)/60./60.
           timearr[i*63+nt] = ((sclk_obs - sclk_0) + 0.4*nt)/60./60.
        endfor

        if i eq 0 then begin
           ra_ref = sxpar(header, 'RA_REF')
           dec_ref = sxpar(header, 'DEC_REF')
           print, 'ra, dec', ra_ref, dec_ref
        endif
        
        get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [3],ra = ra_ref, dec = dec_ref, /silent    

        ;correct for pixel phase effect
        ;corrflux = pixel_phase_correct_gauss(abcdflux,x_center,y_center,2,  '3_12_20')

        if i eq 0 then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
  ;         corrfluxarr = corrflux[1:*]
        endif else begin
           xarr = [xarr, x_center[1:*]]
           yarr = [yarr, y_center[1:*]]
           fluxarr = [fluxarr, abcdflux[1:*]]
           fluxerrarr = [fluxerrarr, fs[1:*]]
 ;          corrfluxarr = [corrfluxarr, corrflux[1:*]]
        endelse
        
     endfor                     ; for each fits file in the AOR
     
     if a eq 0 then begin
        xy = plot(xarr, yarr, '6r1.', xtitle = 'X (pixel)',ytitle = 'Y (pixel)',xrange = [14.5,15.5], yrange = [14.5,15.5],  title = 'HD158460_snap',  color = colorarr[a], aspect_ratio = 1, xthick = 1.5, ythick =1.5) ;
                               
     endif
     
     if a gt 0  then begin
        xy.Select
        xy = plot( xarr, yarr, '6b1.',/overplot, /current,  color = colorarr[a])
     endif
     
;make the timearr 
;     timearr = (findgen(n_elements(xarr)) * 0.4) / 60./60.  ;cheating way
     
           if a eq 0 then begin
        st = plot(timearr, yarr,'6r1.', yrange = [14.8, 15.1],xtitle = 'Time(hours)',ytitle = 'Y (pixel)',  title = 'HD158460_stare',  color = colorarr[a], xthick = 1.5, ythick =1.5) ;
        st2 = plot(timearr, xarr,'6r1.',yrange = [15.05, 15.35], xtitle = 'Time(hours)',ytitle = 'X (pixel)',  title = 'HD158460_stare',  color = colorarr[a], ystyle = 1, xthick = 1.5, ythick =1.5) ;
        ;st3 = errorplot(timearr, fluxarr,fluxerrarr,'6r1.', xtitle = 'Time(hours)',ytitle = 'Flux',  title = 'HD158460_stare') ;
        st3 = plot(timearr, fluxarr,'6r1.', yrange = [1.05, 1.07],xtitle = 'Time(hours)',ytitle = 'Flux',  title = 'HD158460_stare',  color = colorarr[a], xthick = 1.5, ythick =1.5) ;
        ;st4 = plot(timearr, corrfluxarr,'6r1.', yrange = [1.01, 1.08],xtitle = 'Time(hours)',ytitle = 'Corrected Flux',  title = 'HD158460_stare') ;
     endif
     
     if a gt 0 then begin
        st.Select
        st = plot( timearr, yarr, '6b1.',/overplot,/current,  color = colorarr[a])
        st2.Select
        st2 = plot( timearr, xarr, '6b1.',/overplot, /current,  color = colorarr[a])
        st3.Select
        ;st3 = errorplot(timearr, fluxarr-0.2,fluxerrarr,'6r1.',/overplot, /current)
        st3 = plot(timearr, fluxarr,'6b1.',/overplot, /current,  color = colorarr[a])
  ;     st4.Select
        ;st4 = errorplot(timearr, fluxarr-0.2,fluxerrarr,'6r1.',/overplot, /current)
  ;      st4 = plot(timearr, corrfluxarr-0.02,'6b1.',/overplot, /current)

     endif
        
  endfor                        ;for each AOR
  
  ;mode_x =  14.7456 
  ;mode_y =  15.0683
  xsweet = 15.120
  ysweet = 15.085
  xy.Select
;  box_x = [mode_x-0.1, mode_x-0.1, mode_x + 0.1, mode_x + 0.1, mode_x -0.1]
;  box_y = [mode_y-0.1, mode_y +0.1, mode_y +0.1, mode_y - 0.1,mode_y -0.1]
;  line3 = polyline(box_x, box_y, thick = 3, color = !color.black,/data)
  
  box_x = [xsweet-0.1, xsweet-0.1, xsweet + 0.1, xsweet + 0.1, xsweet -0.1]
  box_y = [ysweet-0.1, ysweet +0.1, ysweet +0.1, ysweet - 0.1,ysweet -0.1]
  line4 = polyline(box_x, box_y, thick = 2, color = !color.black,/data)
  
; save these files
  st.Save , '/Users/jkrick/irac_warm/snapshots/yvstime_snap.png'
  st2.Save ,  '/Users/jkrick/irac_warm/snapshots/xvstime_snap.png'
  st3.Save, '/Users/jkrick/irac_warm/snapshots/fluxvstime_snap.png'
;  st4.Save, '/Users/jkrick/irac_warm/snapshots/corrfluxvstime.png'
  xy.Save ,  '/Users/jkrick/irac_warm/snapshots/xy_snap.png'


;save these variables
save, /ALL, filename = '/Users/jkrick/irac_warm/snapshots/HD158460_snap.sav'
end
