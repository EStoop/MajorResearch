pro pcrs_inc_c_v3
  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]

  ;aorname = ['r44254976','r44255232','r44255488','r44255744','r44256000','r44256256','r44256512','r44256768','r44257024','r44257280'] ;

;actually these are pcrsbias ch2_sub 
;aorname = ['r44407808','r44408064','r44408320','r44408576','r44408832','r44409088','r44409344','r44409600','r44409856','r44410112']

; and these are pmpa data
aorname = ['r44464128','r44463872','r44463616','r44463360','r44463104','r44462848','r44462592','r44462336','r44462080','r44461824','r44461568','r44461056','r44460800','r44460544','r44460032']
  nfits = 242*63.
pmap= replicate({pmapob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' '},n_elements(aorname))
 
  ;read in all the pixel phase correction data
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/pmap_ch2_500x500_0043_111129.fits', pmapdata, pmapheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/psigma_ch2_500x500_0043_111129.fits', psigmadata, psigmaheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/occu_ch2_500x500_0043_111129.fits', occdata, occheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/xgrid_ch2_500x500_0043_111129.fits', xgriddata, xgridheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/ygrid_ch2_500x500_0043_111129.fits', ygriddata, ygridheader

  for a = 0,   n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/pmap/' + string(aorname(a)) 
     CD, dir                    ; change directories to the correct AOR directory
     command  =  " ls ch2/bcd/*_bcd.fits > /Users/jkrick/irac_warm/pmap/cbcdlist.txt"
     spawn, command
     command2 =  "ls ch2/bcd/*bunc.fits > /Users/jkrick/irac_warm/pmap/bunclist.txt"
     spawn, command2

     readcol,'/Users/jkrick/irac_warm/pmap/cbcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/pmap/bunclist.txt',buncname, format = 'A', /silent

     yarr = fltarr(n_elements(fitsname)*63)
     xarr = fltarr(n_elements(fitsname)*63)
     timearr = fltarr(n_elements(fitsname)*63)
     fluxarr = fltarr(n_elements(fitsname)*63)
     fluxerrarr = fltarr(n_elements(fitsname)*63)
     corrfluxarr = fltarr(n_elements(fitsname)*63)
     corrfluxerrarr = fltarr(n_elements(fitsname)*63)
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
        
           if occresult(ncor) lt 30. then begin ; 1E6 then begin
              print, 'out of range', fitsname(i), ncor
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
           fluxerrarr = [fluxerrarr, fs[1:*]]
           corrfluxarr = [corrfluxarr, corrflux[1:*]]
           corrfluxerrarr = [corrfluxerrarr, corrfluxerr[1:*]]
        endelse
        
     endfor                     ; for each fits file in the AOR
     ;fill in that structure
     print, 'ne xarr', n_elements(xarr), n_elements(fluxarr), n_elements(pmap[a].xcen)
    pmap[a] ={pmapob,  ra_ref,dec_ref,xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a)}

  endfor                        ;for each AOR
  
  ;test
  print, pmap[0].ra, pmap[0].xcen[26]
  
  save, pmap, filename='/Users/jkrick/irac_warm/pmap/pmap_corr.sav'
  undefine, pmap
  

end
