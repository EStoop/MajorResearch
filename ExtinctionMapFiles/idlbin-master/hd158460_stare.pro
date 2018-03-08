pro hd158460_stare
  
  aorname = ['r42051584','r42506496']
  
  for a = 0,1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/snapshots/' + string(aorname(a)) +'/ch2/bcd'
     CD, dir                    ; change directories to the correct AOR directory
                                ;first = 5
     command  =  " ls *_bcd.fits > /Users/jkrick/irac_warm/snapshots/cbcdlist.txt"
     spawn, command
     
     readcol,'/Users/jkrick/irac_warm/snapshots/cbcdlist.txt',fitsname, format = 'A', /silent
     yarr = fltarr(n_elements(fitsname)*64)
     xarr = fltarr(n_elements(fitsname)*64)
     timearr = fltarr(n_elements(fitsname)*64)
     fluxarr = fltarr(n_elements(fitsname)*64)
     fluxerrarr = fltarr(n_elements(fitsname)*64)
     corrfluxarr = fltarr(n_elements(fitsname)*64)
     for i =0.D, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
                                ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

        for nt = 0, 64 - 1 do begin
           timearr[i*64+nt] = ((sclk_obs - sclk_0) + 0.4*nt)/60./60.
           timearr[i*64+nt] = ((sclk_obs - sclk_0) + 0.4*nt)/60./60.

        endfor

        if i eq 0 then begin
           ra_ref = sxpar(header, 'RA_REF')
           dec_ref = sxpar(header, 'DEC_REF')
           print, 'ra, dec', ra_ref, dec_ref
        endif
        
        get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [3],ra = ra_ref, dec = dec_ref, /silent    

        ;correct for pixel phase effect

        ;corrflux = correct_pixel_phase(2,x_center,y_center,abcdflux,/SUB)
        corrflux = pixel_phase_correct_gauss(abcdflux,x_center,y_center,2,  '3_12_20')

        if i eq 0 then begin
           xarr = x_center
           yarr = y_center
           fluxarr = abcdflux
           fluxerrarr = fs
           corrfluxarr = corrflux
        endif else begin
           xarr = [xarr, x_center]
           yarr = [yarr, y_center]
           fluxarr = [fluxarr, abcdflux]
           fluxerrarr = [fluxerrarr, fs]
           corrfluxarr = [corrfluxarr, corrflux]
        endelse
        
     endfor                     ; for each fits file in the AOR
     
     
     if a eq 0 then begin
        xy = plot(xarr, yarr, '6r1.', xtitle = 'X pix',ytitle = 'Y pix',xrange = [14.0,16.5], yrange = [14.0,16.5],  title = 'HD158460_stare') ;
                               
     endif
     
     if a gt 0  then begin
        xy.Select
        xy = plot( xarr, yarr, '6b1.',/overplot, /current)
     endif
     
;make the timearr 
;     timearr = (findgen(n_elements(xarr)) * 0.4) / 60./60.  ;cheating way
     
     
     if a eq 0 then begin
        st = plot(timearr, yarr,'6r1.', yrange = [14.5, 15.2],xtitle = 'Time(hours)',ytitle = 'Y pix',  title = 'HD158460_stare') ;
        st2 = plot(timearr, xarr,'6r1.',yrange = [14.5, 15.2], xtitle = 'Time(hours)',ytitle = 'X pix',  title = 'HD158460_stare') ;
        ;st3 = errorplot(timearr, fluxarr,fluxerrarr,'6r1.', xtitle = 'Time(hours)',ytitle = 'Flux',  title = 'HD158460_stare') ;
        st3 = plot(timearr, fluxarr,'6r1.', yrange = [1.01, 1.08],xtitle = 'Time(hours)',ytitle = 'Flux',  title = 'HD158460_stare') ;
        st4 = plot(timearr, corrfluxarr,'6r1.', yrange = [1.01, 1.08],xtitle = 'Time(hours)',ytitle = 'Corrected Flux',  title = 'HD158460_stare') ;
     endif
     
     if a gt 0 then begin
        st.Select
        st = plot( timearr, yarr, '6b1.',/overplot,/current)
        st2.Select
        st2 = plot( timearr, xarr- 0.2, '6b1.',/overplot, /current)
        st3.Select
        ;st3 = errorplot(timearr, fluxarr-0.2,fluxerrarr,'6r1.',/overplot, /current)
        st3 = plot(timearr, fluxarr-0.02,'6b1.',/overplot, /current)
       st4.Select
        ;st4 = errorplot(timearr, fluxarr-0.2,fluxerrarr,'6r1.',/overplot, /current)
        st4 = plot(timearr, corrfluxarr-0.02,'6b1.',/overplot, /current)

     endif
     
     
     
  endfor                        ;for each AOR
  
  ;mode_x =  14.7456 
  ;mode_y =  15.0683
  ;xsweet = 15.120
  ;ysweet = 15.085
  ;xy.Select
  ;box_x = [mode_x-0.1, mode_x-0.1, mode_x + 0.1, mode_x + 0.1, mode_x -0.1]
  ;box_y = [mode_y-0.1, mode_y +0.1, mode_y +0.1, mode_y - 0.1,mode_y -0.1]
  ;line3 = polyline(box_x, box_y, thick = 3, color = !color.black,/data)
  
  ;box_x = [xsweet-0.1, xsweet-0.1, xsweet + 0.1, xsweet + 0.1, xsweet -0.1]
  ;box_y = [ysweet-0.1, ysweet +0.1, ysweet +0.1, ysweet - 0.1,ysweet -0.1]
;line4 = polyline(box_x, box_y, thick = 3, color = !color.black,/data)
  
; save these files
  st.Save , '/Users/jkrick/irac_warm/snapshots/yvstime.png'
  st2.Save ,  '/Users/jkrick/irac_warm/snapshots/xvstime.png'
  st3.Save, '/Users/jkrick/irac_warm/snapshots/fluxvstime.png'
  st4.Save, '/Users/jkrick/irac_warm/snapshots/corrfluxvstime.png'
  xy.Save ,  '/Users/jkrick/irac_warm/snapshots/xy.png'


;save these variables
save, /ALL, filename = '/Users/jkrick/irac_warm/snapshots/HD158460_stare.sav'
end
