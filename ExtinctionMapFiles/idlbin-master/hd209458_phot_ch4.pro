pro hd209458_phot_ch4
 
;hd209458 staring mode 
ra_ref = 330.79488
dec_ref = 18.8843175

  aorname = ['r24649984','r24650240','r24650496']

  for a = 0 , n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/hd209458/'+ string(aorname(a)) + '/ch4/bcd'
     CD, dir                    ; change directories to the correct AOR directory
     command  =  " find . -name '*_bcd.fits' > /Users/jkrick/irac_warm/hd209458/cbcdlist.txt"
     spawn, command
     command2 =  "find . -name '*_bunc.fits' > /Users/jkrick/irac_warm/hd209458/bunclist.txt"
     spawn, command2

     readcol,'/Users/jkrick/irac_warm/hd209458/cbcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/hd209458/bunclist.txt',buncname, format = 'A', /silent

 ;    yarr = fltarr(n_elements(fitsname)*63)
 ;    xarr = fltarr(n_elements(fitsname)*63)
 ;    fluxarr = fltarr(n_elements(fitsname)*63)
 ;    fluxerrarr = fltarr(n_elements(fitsname)*63)
 ;    corrfluxarr = fltarr(n_elements(fitsname)*63)
;     corrfluxerrarr = fltarr(n_elements(fitsname)*63)
     for i =0.D, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        frametime = sxpar(header, 'FRAMTIME')
        bmjd_obs = sxpar(header, 'BMJD_OBS')
        utcs_obs = sxpar(header, 'UTCS_OBS')

        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

        sclkarr = fltarr(64)
        bmjdarr = fltarr(64)
        utcsarr = fltarr(64)
        for nt = 0, 64 - 1 do begin
           sclkarr[nt] = (sclk_obs ) + 0.5*frametime + frametime*nt
           bmjdarr[nt]= bmjd_obs + 0.5*frametime + frametime*nt
           utcsarr[nt]= utcs_obs + 0.5*frametime + frametime*nt
        endfor

        ;do the photometry
        ;get_centroids,fitsname(i), testt, testdt, testx_center, testy_center, testabcdflux, testxs, testys, testfs, testb, /WARM, /APER, APRAD = [3],ra = ra_ref, dec = dec_ref, /silent    

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
  
                    
        if a eq 0 and  i eq 0 then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
           timearr = sclkarr[1:*]
           bmjd = bmjdarr[1:*]
           utcs = utcsarr[1:*]
        endif else begin
           xarr = [xarr, x_center[1:*]]
           yarr = [yarr, y_center[1:*]]
           fluxarr = [fluxarr, abcdflux[1:*]]
           fluxerrarr = [fluxerrarr, fs[1:*]]
          timearr = [timearr, sclkarr[1:*]]
          bmjd = [bmjd, bmjdarr[1:*]]
           utcs = [utcs, utcsarr[1:*]]
       endelse
        
     endfor                     ; for each fits file in the AOR
 ;    xy=m.plot(xarr, yarr, '6r1.',color = 'black', thick = 2,xrange = [14.5,15.5], yrange = [14.5,15.5]);,multi_index = a)
 ;    st = m.plot(timearr, yarr,'6r1', yrange = [14.5, 15.5], xrange = [0, 0.6]) ;, multi_index = n_elements(aorname) + a) ;
     
 ;    print, 'x, y', mean(xarr), mean(yarr)
 ;    print, 'fluxarr', mean(fluxarr)
 ;    print, 'corrfluxarr', mean(corrfluxarr,/nan)
     ;fill in that structure
;     print, 'n xarr, fluxarr, nfits', n_elements(xarr), n_elements(fluxarr), nfits
;     snapshots[a] ={planetob,  ra_ref,dec_ref,xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a)}

  endfor                        ;for each AOR
  
  ;test
 ; print, snapshots[0].ra, snapshots[0].xcen[26]
  
 ; save, snapshots, filename='/Users/jkrick/irac_warm/hd209458/hd209458_corr.sav'
 ; undefine, snapshots
  
;sort by time
a = sort(timearr)
sorttime = timearr[a]
sortxarr = xarr[a]
sortyarr = yarr[a]
sortfluxarr = fluxarr[a]
sortfluxerrarr = fluxerrarr[a]
sortbmjd = bmjd[a]
sortutcs = utcs[a]

save, /all, filename = '/Users/jkrick/irac_warm/hd209458/hd209458_ch4_new.sav'
end
