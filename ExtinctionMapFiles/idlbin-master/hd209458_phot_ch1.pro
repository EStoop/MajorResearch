pro hd209458_phot_ch1
  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]

;hd209458 staring mode 
ra_ref = 330.79488
dec_ref = 18.8843175
startfits = 90.D
  aorname = ['r41628416','r41628672','r41628928','r41629184','r41629440'];,'r41629696','r41629952','r41630208','r41630464']


  ;read in all the pixel phase correction data
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch1_500x500_0043_120110/pmap_ch1_500x500_0043_120110.fits', pmapdata, pmapheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch1_500x500_0043_120110/psigma_ch1_500x500_0043_120110.fits', psigmadata, psigmaheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch1_500x500_0043_120110/occu_ch1_500x500_0043_120110.fits', occdata, occheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch1_500x500_0043_120110/xgrid_ch1_500x500_0043_120110.fits', xgriddata, xgridheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch1_500x500_0043_120110/ygrid_ch1_500x500_0043_120110.fits', ygriddata, ygridheader

  for a = 0 , n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/hd209458/'+ string(aorname(a)) + '/ch1/bcd'
     CD, dir                    ; change directories to the correct AOR directory
     command  =  " find . -name '*_bcd.fits' > /Users/jkrick/irac_warm/hd209458/cbcdlist.txt" ; subarray has no cbcd, cbunc
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
     for i =startfits, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        frametime = sxpar(header, 'FRAMTIME')
        bmjd_obs = sxpar(header, 'BMJD_OBS')
        utcs_obs = sxpar(header, 'UTCS_OBS')

        if i eq startfits then sclk_0 = sxpar(header, 'SCLK_OBS')

        sclkarr = fltarr(64)
        bmjdarr = fltarr(64)
        utcsarr = fltarr(64)
        for nt = 0, 64 - 1 do begin
           sclkarr[nt] = (sclk_obs ) + 0.5*frametime + frametime*nt
           bmjdarr[nt]= bmjd_obs + 0.5*frametime + frametime*nt
           utcsarr[nt]= utcs_obs + 0.5*frametime + frametime*nt
        endfor


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
  
        ;print,'abcdflux', abcdflux
      ;correct for pixel phase effect based on pmaps from Jim
        ;corrflux = pixel_phase_correct_gauss(abcdflux,x_center,y_center,2,  '3_12_20')
        corrresult = interp2d(pmapdata,xgriddata,ygriddata,x_center,y_center,/regular)
        pmaperr = interp2d(psigmadata,xgriddata,ygriddata,x_center,y_center,/regular)
        occresult = interp2d(occdata,xgriddata,ygriddata,x_center,y_center,/regular)
               
        ;print out some summaries
;        print, 'x, y, flux, correction', mean(x_center), mean(y_center), mean(abcdflux), mean(corrresult), format = '(A, F10.4, F10.4, F10.4, F10.4)'
        ;start by applying the correction to all 65 subarray images per frame
        corrflux = abcdflux / corrresult
        corrfluxerr = sqrt(fs^2 + pmaperr^2)
     
        ;then nan out those which do not have a high enough occupation
        for ncor = 0, n_elements(abcdflux) - 1 do begin
        
           if occresult(ncor) lt 40. then begin ; 1E6 then begin
              ;print, 'out of range', fitsname(i), ncor
              corrflux(ncor) = abcdflux(ncor) / alog10(-1)
              corrfluxerr(ncor)= abcdflux(ncor)/ alog10(-1)
           endif
           
        endfor

        if a eq 0 and  i eq startfits then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
           corrfluxarr = corrflux[1:*]
           corrfluxerrarr = corrfluxerr[1:*]
           timearr = sclkarr[1:*]
           bmjd = bmjdarr[1:*]
           utcs = utcsarr[1:*]
        endif else begin
           xarr = [xarr, x_center[1:*]]
           yarr = [yarr, y_center[1:*]]
           fluxarr = [fluxarr, abcdflux[1:*]]
           fluxerrarr = [fluxerrarr, fs[1:*]]
           corrfluxarr = [corrfluxarr, corrflux[1:*]]
           corrfluxerrarr = [corrfluxerrarr, corrfluxerr[1:*]]
           timearr = [timearr, sclkarr[1:*]]
          bmjd = [bmjd, bmjdarr[1:*]]
           utcs = [utcs, utcsarr[1:*]]
        endelse

        ;print, 'meany',i, mean(y_center,/nan)

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
sortcorrfluxarr = corrfluxarr[a]
sortcorrfluxerrarr = corrfluxerrarr[a]
sortbmjd = bmjd[a]
sortutcs = utcs[a]

save, /all, filename = '/Users/jkrick/irac_warm/hd209458/hd209458_ch1_new.sav'
end
