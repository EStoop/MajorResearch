pro hd158460_stare
  
  aorname = ['r42506496','r42051584']
  
  for a = 0,1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/snapshots/' + string(aorname(a)) +'/ch2/bcd'
     CD, dir                    ; change directories to the correct AOR directory
                                ;first = 5
     command  =  " ls *_bcd.fits > /Users/jkrick/irac_warm/snapshots/cbcdlist.txt"
     spawn, command
     command2 =  "ls *bunc.fits > /Users/jkrick/irac_warm/snapshots/bunclist.txt"
     spawn, command2

     readcol,'/Users/jkrick/irac_warm/snapshots/cbcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/snapshots/bunclist.txt',buncname, format = 'A', /silent
     
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
        
        ;get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [3],ra = ra_ref, dec = dec_ref, /silent    

      ;print, 'working on ', fitsname
      fits_read, fitsname(i), im, h
      fits_read, buncname(i), unc, hunc
     ; fits_read, covname, covdata, covheader
      get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                   x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                   x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                   xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                   xfwhm, yfwhm, /WARM
      ;I think flux comes out of this in electrons

        ;correct for pixel phase effect
      ;flux input needs to be in Jy
      ;choose 3 pixel aperture, 3-7 pixel background
      mjyflux = electron_to_mjy( f[*,9],sxpar(h,'PXSCAL2'),  sxpar(h, 'GAIN') , sxpar(h, 'EXPTIME') , sxpar(h,'FLUXCONV'))
      jyflux = mjyflux / 1E3
      mjyfluxerr = electron_to_mjy( fs[*,9],sxpar(h,'PXSCAL2'),  sxpar(h, 'GAIN') , sxpar(h, 'EXPTIME') , sxpar(h,'FLUXCONV')); 
      corrflux = pixel_phase_correct_gauss(jyflux,x3,y3,2,  '3_12_20',/subarray)
     corrflux = corrflux * 1E6
       ; print, 'ncorrflux', n_elements(corrflux)
        if i eq 0 then begin
           xarr = x3[1:*]
           yarr = y3[1:*]
           fluxarr = mjyflux[1:*]*1E3; f[1:*,1]
           fluxerrarr = mjyfluxerr[1:*]*1E3
           corrfluxarr = corrflux[1:*]
        endif else begin
           xarr = [xarr, x3[1:*]]
           yarr = [yarr, y3[1:*]]
           fluxarr = [fluxarr, mjyflux[1:*]*1E3]
           fluxerrarr = [fluxerrarr, mjyfluxerr[1:*]*1E3]
           corrfluxarr = [corrfluxarr, corrflux[1:*]]
        endelse
        
     endfor                     ; for each fits file in the AOR
     
     
     if a eq 0 then begin
        xy = plot(xarr, yarr, '6r1.', xtitle = 'X pix',ytitle = 'Y pix',xrange = [22.5,23.5], yrange = [230.5,231.5],  title = 'HD158460_stare') ;
                               
     endif
     
     if a gt 0  then begin
        xy.Select
        xy = plot( xarr, yarr, '6b1.',/overplot, /current)
     endif
     
;make the timearr 
;     timearr = (findgen(n_elements(xarr)) * 0.4) / 60./60.  ;cheating way
     
     
     if a eq 0 then begin
        st = plot(timearr, yarr,'6r1.', yrange = [230.5,231.5],xtitle = 'Time(hours)',ytitle = 'Y pix',  title = 'HD158460_stare') ;
        st2 = plot(timearr, xarr,'6r1.',yrange = [22.5,23.5], xtitle = 'Time(hours)',ytitle = 'X pix',  title = 'HD158460_stare') ;
        ;st3 = errorplot(timearr, fluxarr,fluxerrarr,'6r1.', xtitle = 'Time(hours)',ytitle = 'Flux',  title = 'HD158460_stare') ;
        st3 = plot(timearr, fluxarr,'6r1.', yrange = [4.0, 4.3],xtitle = 'Time(hours)',ytitle = 'Flux',  title = 'HD158460_stare') ;
        st4 = plot(timearr, corrfluxarr,'6r1.', yrange = [4.0, 4.15],xtitle = 'Time(hours)',ytitle = 'Corrected Flux',  title = 'HD158460_stare') ;
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
     
     if a eq 0 then begin
        yarr_0 =yarr
        xarr_0  = xarr
        timearr_0  = timearr 
        fluxarr_0  = fluxarr
        fluxerrarr_0  = fluxerrarr
        corrfluxarr_0  = corrfluxarr
     endif
     if a eq 1 then begin
        yarr_1 =yarr
        xarr_1  = xarr
        timearr_1  = timearr 
        fluxarr_1  = fluxarr
        fluxerrarr_1  = fluxerrarr
        corrfluxarr_1  = corrfluxarr
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
  st.Save , '/Users/jkrick/irac_warm/snapshots/yvstime_stare.png'
  st2.Save ,  '/Users/jkrick/irac_warm/snapshots/xvstime_stare.png'
  st3.Save, '/Users/jkrick/irac_warm/snapshots/fluxvstime_stare.png'
  st4.Save, '/Users/jkrick/irac_warm/snapshots/corrfluxvstime_stare.png'
  xy.Save ,  '/Users/jkrick/irac_warm/snapshots/xy_stare.png'


;save these variables
save, /ALL, filename = '/Users/jkrick/irac_warm/snapshots/HD158460_stare_v2.sav'
end
