pro hd7924_pmap
  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]

;these are the snapshot mode AORs HD158460
;remove those with no pmap coverage
;  aorname = ['r44497152','r44497664','r44497920','r44498176','r44498432','r44498688','r44498944','r44499200','r44499456','r44500224','r44500480','r44500736', 'r44500992', 'r44501248', 'r44501504' ]

;these are the staring mode AORs HD158460
;  aorname = ['r42506496','r42051584']


;koi-69 with pcrspeakup
;aorname = [   'r44448512',  'r44448768']

;hd7924 with peakup
aorname = ['r44605184']

;m = pp_multiplot(multi_layout=[n_elements(aorname),2]);,global_xtitle='Test X axis title',global_ytitle='Test Y axis title')
  ;nfits = 75*63.  ; for snapshots
  ;nfits = 1330*63. ; for staring mode
  nfits = 169. * 63; for koi 69
  nfits = 3400.*63 ; for hd7924

  snapshots = replicate({planetob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', back:fltarr(nfits)},n_elements(aorname))

  ;read in all the pixel phase correction data
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/pmap_ch2_500x500_0043_111129.fits', pmapdata, pmapheader;
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/psigma_ch2_500x500_0043_111129.fits', psigmadata, psigmaheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/occu_ch2_500x500_0043_111129.fits', occdata, occheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/xgrid_ch2_500x500_0043_111129.fits', xgriddata, xgridheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_500x500_0043_111129/ygrid_ch2_500x500_0043_111129.fits', ygriddata, ygridheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/pmap_ch2_0p4s_500x500_0043_111206.fits', pmapdata, pmapheader
 fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/psigma_ch2_0p4s_500x500_0043_111206.fits', psigmadata, psigmaheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/occu_ch2_0p4s_500x500_0043_111206.fits', occdata, occheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/xgrid_ch2_0p4s_500x500_0043_111206.fits', xgriddata, xgridheader
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/ygrid_ch2_0p4s_500x500_0043_111206.fits', ygriddata, ygridheader

  for a = 0 , n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/snapshots/hd7924/' + string(aorname(a)) 
     CD, dir                    ; change directories to the correct AOR directory
     command  =  " ls ch2/bcd/*_bcd.fits > /Users/jkrick/irac_warm/snapshots/hd7924/cbcdlist.txt"
     spawn, command
     command2 =  "ls ch2/bcd/*bunc.fits > /Users/jkrick/irac_warm/snapshots/hd7924/bunclist.txt"
     spawn, command2

     readcol,'/Users/jkrick/irac_warm/snapshots/hd7924/cbcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/snapshots/hd7924/bunclist.txt',buncname, format = 'A', /silent

     yarr = fltarr(n_elements(fitsname)*63)
     xarr = fltarr(n_elements(fitsname)*63)
     timearr = fltarr(n_elements(fitsname)*63)
     fluxarr = fltarr(n_elements(fitsname)*63)
     fluxerrarr = fltarr(n_elements(fitsname)*63)
     backarr = fltarr(n_elements(fitsname)*63)
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
        b = b[*,2]

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
        corrfluxerr = sqrt((fs^2)*(corrresult^2) + (pmaperr^2)*(abcdflux^2))
     
        ;then nan out those which do not have a high enough occupation
        ;for ncor = 0, n_elements(abcdflux) - 1 do begin
        
        ;   if occresult(ncor) lt 7. then begin ; 1E6 then begin
        ;      ;print, 'out of range', fitsname(i), ncor
        ;      corrflux(ncor) = abcdflux(ncor) / alog10(-1)
        ;      corrfluxerr(ncor)= abcdflux(ncor)/ alog10(-1)
        ;   endif
           
        ;endfor

        if i eq 0 then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
           backarr = b[1:*]
           corrfluxarr = corrflux[1:*]
           corrfluxerrarr = corrfluxerr[1:*]
        endif else begin
           xarr = [xarr, x_center[1:*]]
           yarr = [yarr, y_center[1:*]]
           fluxarr = [fluxarr, abcdflux[1:*]]
           fluxerrarr = [fluxerrarr, fs[1:*]]
           backarr = [backarr, b[1:*]]
           corrfluxarr = [corrfluxarr, corrflux[1:*]]
           corrfluxerrarr = [corrfluxerrarr, corrfluxerr[1:*]]
        endelse
        
     endfor                     ; for each fits file in the AOR
 ;    xy=m.plot(xarr, yarr, '6r1.',color = 'black', thick = 2,xrange = [14.5,15.5], yrange = [14.5,15.5]);,multi_index = a)
 ;    st = m.plot(timearr, yarr,'6r1', yrange = [14.5, 15.5], xrange = [0, 0.6]) ;, multi_index = n_elements(aorname) + a) ;
     
 ;    print, 'x, y', mean(xarr), mean(yarr)
 ;    print, 'fluxarr', mean(fluxarr)
 ;    print, 'corrfluxarr', mean(corrfluxarr,/nan)
     ;fill in that structure

     snapshots[a] ={planetob,  ra_ref,dec_ref,xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a), backarr}

  endfor                        ;for each AOR
  
  ;test
  print, snapshots[0].ra, snapshots[0].xcen[26]
  
  save, snapshots, filename='/Users/jkrick/irac_warm/snapshots/hd7924/hd7924_corr.sav'
  undefine, snapshots
  

end
