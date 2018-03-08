pro hd158460_stare
  
  aorname = ['r42506496','r42051584']
    ;read in all the pixel phase correction data
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/pmap_ch2_500x500_0043_maxnorm.fits', pmapdata, pmapheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/psigma_ch2_500x500_0043_111129.fits', psigmadata, psigmaheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/occu_ch2_500x500_0043_111129.fits', occdata, occheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/xgrid_ch2_500x500_0043_111129.fits', xgriddata, xgridheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/ygrid_ch2_500x500_0043_111129.fits', ygriddata, ygridheader

  for a = 0,0 do begin
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
        
                                ;choose 3 pixel aperture, 3-7 pixel background
        abcdflux = f[*,9]       ;put it back in same nomenclature
        fs = fs[*,9]
        x_center = x3
        y_center = y3
        
;correct for pixel phase effect
        corrresult = interp2d(pmapdata,xgriddata,ygriddata,x_center,y_center,/regular)
        pmaperr = interp2d(psigmadata,xgriddata,ygriddata,x_center,y_center,/regular)
        occresult = interp2d(occdata,xgriddata,ygriddata,x_center,y_center,/regular)
                  ;print out some summaries
        ;print, 'x, y, flux, correction', mean(x_center), mean(y_center), mean(abcdflux), mean(corrresult), format = '(A, F10.4, F10.4, F10.4, F10.4)'


        ;start by applying the correction to all 65 subarray images per frame
        corrflux = abcdflux / corrresult
        corrfluxerr = sqrt(fs^2 + pmaperr^2)
     
        ;then nan out those which do not have a high enough occupation
        for ncor = 0, n_elements(abcdflux) - 1 do begin
        
           if occresult(ncor) lt 100. then begin ; 1E6 then begin
              ;print, 'out of range', fitsname(i), ncor
              corrflux(ncor) = abcdflux(ncor) / alog10(-1)
              corrfluxerr(ncor)= abcdflux(ncor)/ alog10(-1)
           endif
           
        endfor
     

        if i eq 0 then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
           corrfluxarr = corrflux[1:*]
           corrfluxerrarr = corrfluxerr[1:*]

        endif else begin
           xarr = [xarr, x_center[1:*]]
           yarr = [yarr, y_center[1:*]]
           fluxarr = [fluxarr, abcdflux[1:*]]
           fluxerrarr = [fluxerrarr, abcdflux[1:*]]
           corrfluxarr = [corrfluxarr, corrflux[1:*]]
           corrfluxerrarr = [corrfluxerrarr, corrfluxerr[1:*]]

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

         ;fill in that structure
     snapshots[a] ={snapob,  ra_ref,dec_ref,xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a)}

     
  endfor                        ;for each AOR
  
 
  save, snapshots, filename='/Users/jkrick/irac_warm/snapshots/stare_corr.sav'
  undefine, snapshots

;save these variables
save, /ALL, filename = '/Users/jkrick/irac_warm/snapshots/HD158460_stare_v3.sav'
end
