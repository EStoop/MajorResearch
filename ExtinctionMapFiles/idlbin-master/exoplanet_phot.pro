pro exoplanet_phot
  colorarr = [ 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE', 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE', 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]
;these are the snapshot mode AORs HD158460
;remove those with no pmap coverage
;  aorname = ['r44497152','r44497664','r44497920','r44498176','r44498432','r44498688','r44498944','r44499200','r44499456','r44500224','r44500480','r44500736', 'r44500992', 'r44501248', 'r44501504' ]

;these are the staring mode AORs HD158460
;  aorname = ['r42506496','r42051584']


;koi-69 with pcrspeakup
;aorname = [   'r44448512',  'r44448768']

;hd7924 with peakup
;aorname = ['r44605184']

;pcrs_test
;aorname = ['0048360448', '0046475776', '0048072192', '0048359936','0046482176','0048073728'];,'0046352128'];, '0048073728']
  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_fits/pmap_ch2_0p1s_x4_500x500_0043_121120.fits', pmapdata, pmapheader
  c = contour(pmapdata, /fill, n_levels = 21, rgb_table = 0, xtitle = 'X (pixel)', ytitle = 'Y (pixel)', title = planetname, aspect_ratio = 1, xrange = [-500,1000], yrange = [-500,1500])

;wasp62
;aorname = ['r48685056','r48685312','r48685568','r48685824','r48686080','r48686336','r48686592','r48686848','r48687104','r48687360','r48688640','r48688896','r48689152','r48689408','r48689664','r48689920','r48690176','r48690432','r48690688','r48696064','r48696320','r48696576','r48696832','r48697088','r48697344','r48697600','r48700416','r48700672','r48700928','r48701184','r48701440','r48701696','r48701952','r48702208','r48702464','r48702720','r48706304','r48706560','r48706816','r48707072','r48707328','r48707584','r48707840','r48708096','r48708352','r48708608']

 aorname = ['r48704512','r48693504']

;m = pp_multiplot(multi_layout=[n_elements(aorname),2]);,global_xtitle='Test X axis title',global_ytitle='Test Y axis title')
  ;nfits = 75*63.  ; for snapshots
  ;nfits = 1330*63. ; for staring mode
  ;nfits = 10. * 63; for koi 69
 ; nfits = 10.*63 ; for hd7924

 ; snapshots = replicate({planetob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' '},n_elements(aorname))

  ;read in all the pixel phase correction data
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/pmap_ch2_0p4s_500x500_0043_111206.fits', pmapdata, pmapheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/psigma_ch2_0p4s_500x500_0043_111206.fits', psigmadata, psigmaheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/occu_ch2_0p4s_500x500_0043_111206.fits', occdata, occheader
 ; fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/xgrid_ch2_0p4s_500x500_0043_111206.fits', xgriddata, xgridheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/ygrid_ch2_0p4s_500x500_0043_111206.fits', ygriddata, ygridheader

  for a = 0 , n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/pcrs_planets/wasp62/' + string(aorname(a)) + '/ch2/bcd/'
     CD, dir                    ; change directories to the correct AOR directory
     command  =  " ls *_bcd.fits > /Users/jkrick/irac_warm/pcrs_planets/wasp62/bcdlist.txt"
     spawn, command
     command2 =  "ls *_bunc.fits > /Users/jkrick/irac_warm/pcrs_planets/wasp62/bunclist.txt"
     spawn, command2

     readcol,'/Users/jkrick/irac_warm/pcrs_planets/wasp62/bcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/pcrs_planets/wasp62/bunclist.txt',buncname, format = 'A', /silent

 ;    yarr = fltarr(n_elements(fitsname)*63)
 ;    xarr = fltarr(n_elements(fitsname)*63)
 ;    timearr = fltarr(n_elements(fitsname)*63)
 ;    fluxarr = fltarr(n_elements(fitsname)*63)
 ;    fluxerrarr = fltarr(n_elements(fitsname)*63)
 ;    corrfluxarr = fltarr(n_elements(fitsname)*63)
 ;    corrfluxerrarr = fltarr(n_elements(fitsname)*63)

     xarr = fltarr(10)
     yarr = fltarr(10)
     timearr = fltarr(10)
     fluxarr = fltarr(10)
     fluxerrarr = fltarr(10)

     for i =0.D, 12. do begin;  n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

 ;       for nt = 0, 63 - 1 do begin
;           timearr[i*63+nt] = ((sclk_obs - sclk_0) + 0.4*nt)/60./60.
 ;          timearr[i*63+nt] = ((sclk_obs - sclk_0) + 0.4*nt)/60./60.
;        endfor

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
;        corrresult = interp2d(pmapdata,xgriddata,ygriddata,x_center,y_center,/regular)
;        pmaperr = interp2d(psigmadata,xgriddata,ygriddata,x_center,y_center,/regular)
;        occresult = interp2d(occdata,xgriddata,ygriddata,x_center,y_center,/regular)
               
        ;print out some summaries
;        print, 'x, y, flux, correction', mean(x_center), mean(y_center), mean(abcdflux), mean(corrresult), format = '(A, F10.4, F10.4, F10.4, F10.4)'


 ;;       ;start by applying the correction to all 65 subarray images per frame
 ;       corrflux = abcdflux / corrresult
;        corrfluxerr = sqrt(fs^2 + pmaperr^2)
     
        ;then nan out those which do not have a high enough occupation
 ;       for ncor = 0, n_elements(abcdflux) - 1 do begin
        
;           if occresult(ncor) lt 40. then begin ; 1E6 then begin
;              ;print, 'out of range', fitsname(i), ncor
;              corrflux(ncor) = abcdflux(ncor) / alog10(-1)
;              corrfluxerr(ncor)= abcdflux(ncor)/ alog10(-1)
;           endif
           
;        endfor

        if i eq 0 then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
 ;          corrfluxarr = corrflux[1:*]
 ;          corrfluxerrarr = corrfluxerr[1:*]
        endif else begin
           xarr = [xarr, x_center[1:*]]
           yarr = [yarr, y_center[1:*]]
           fluxarr = [fluxarr, abcdflux[1:*]]
           fluxerrarr = [fluxerrarr, fs[1:*]]
;           corrfluxarr = [corrfluxarr, corrflux[1:*]]
;           corrfluxerrarr = [corrfluxerrarr, corrfluxerr[1:*]]
        endelse
        
     endfor                     ; for each fits file in the AOR
 ;    xy=m.plot(xarr, yarr, '6r1.',color = 'black', thick = 2,xrange = [14.5,15.5], yrange = [14.5,15.5]);,multi_index = a)
 ;    st = m.plot(timearr, yarr,'6r1', yrange = [14.5, 15.5], xrange = [0, 0.6]) ;, multi_index = n_elements(aorname) + a) ;
     
 ;    print, 'x, y', mean(xarr), mean(yarr)
 ;    print, 'fluxarr', mean(fluxarr)
 ;    print, 'corrfluxarr', mean(corrfluxarr,/nan)
     ;fill in that structure

;     snapshots[a] ={planetob,  ra_ref,dec_ref,xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a)}


;plots 
;     print, 'xarr', xarr
;     print, 'yarr', yarr
     xcen500 = 500.* (xarr- 14.5)
     ycen500 = 500.* (yarr - 14.5)
;     print, 'xcen', xcen500
 ;    print, 'ycen', ycen500
     an = plot(xcen500, ycen500, '6r1s', sym_size = 0.2,   sym_filled = 1, color = colorarr[a],/overplot)



  endfor                        ;for each AOR
  
  ;test
 ; print, snapshots[0].ra, snapshots[0].xcen[6]
  
;  save, snapshots, filename='/Users/jkrick/irac_warm/snapshots/hd7924/hd7924_corr.sav'
;  undefine, snapshots
  
;  an.save, dirname+'position.png'

end
