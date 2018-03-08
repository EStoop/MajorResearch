pro phot_55cnc_ch2_42000384
  

;55cnc staring mode ch2
ra_ref = 133.14754
dec_ref = 28.330195

  aorname = ['r42000384'] ;ch2
  nfits =6230.*63             

  AOR55cnc = replicate({ob55cnc, ra:0D,dec:0D,xcen:dblarr(nfits),ycen:dblarr(nfits),flux:dblarr(nfits),fluxerr:dblarr(nfits), corrflux:dblarr(nfits), corrfluxerr:dblarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', bmjdarr:dblarr(nfits), utcsarr:dblarr(nfits), bkgd:dblarr(nfits), bkgderr:dblarr(nfits), xgcntrd:dblarr(nfits), ygcntrd:dblarr(nfits), xcntrd:dblarr(nfits), ycntrd:dblarr(nfits)},n_elements(aorname)) ;, xfwhmarr:dblarr(nfits), yfwhmarr:dblarr(nfits


  for a = 0,   n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/pcrs_planets/55cnc/'+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  =  "find ch2/bcd -name '*_bcd.fits' > /Users/jkrick/irac_warm/pcrs_planets/55cnc/cbcdlist.txt"
     spawn, command
     command2 =  "find  ch2/bcd -name '*bunc.fits' > /Users/jkrick/irac_warm/pcrs_planets/55cnc/bunclist.txt"
     spawn, command2

     readcol,'/Users/jkrick/irac_warm/pcrs_planets/55cnc/cbcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/pcrs_planets/55cnc/bunclist.txt',buncname, format = 'A', /silent

     lc = dblarr(n_elements(fitsname)*64)
     count = 0.D
     for i =0.D, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        frametime = sxpar(header, 'FRAMTIME')
        bmjd_obs = sxpar(header, 'BMJD_OBS')
        utcs_obs = sxpar(header, 'UTCS_OBS')

        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

        sclkarr = dblarr(64)
        bmjdarr = dblarr(64)
        utcsarr = dblarr(64)
        for nt = 0., 64 - 1 do begin
           sclkarr[nt] = (sclk_obs ) + 0.5*frametime + frametime*nt
           bmjdarr[nt]= bmjd_obs + 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt
           utcsarr[nt]= utcs_obs + 0.5*(frametime/60./60./24.)  + (frametime/60./60./24.) *nt
           ;print, 'bmjd stuff', nt, frametime, ' ', bmjd_obs,  ' ', bmjd_obs + 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt,' ', bmjdarr[nt],  format = '(A, F10.2, F10.2,A, F0, A, F0,A,F0)'
        endfor
;        for na = 0, n_elements(bmjdarr) - 1 do print, 'bmjdarr ', bmjdarr(na), format = '(A, F0)'

        ;do the photometry
   
        fits_read, fitsname(i), im, h
        fits_read, buncname(i), unc, hunc
                                ; fits_read, covname, covdata, covheader
        get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                     x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                     x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                     xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                     xfwhm, yfwhm, bb, /WARM
       ; print, 'test return fwhm', xfwhm, yfwhm, x3
      ;choose 3 pixel aperture, 3-7 pixel background
        abcdflux = f[*,9]       ;put it back in same nomenclature
        fs = fs[*,9]
        x_center = x3
        y_center = y3
;        back =  b[*,2]  ;for 3, 3-7
        back =  bb[*,2];  use the box_centroider background as a test
        backerr = bs[*,2]
        ;xfw = xfwhm
        ;yfw = yfwhm
        x_gcntrd = xg
        y_gcntrd = yg
        x_cntrd = xh
        y_cntrd = yh

;--------------------------------
      ;correct for pixel phase effect based on pmaps from Jim
        
        corrflux = iracpc_pmap_corr(abcdflux, x_center, y_center, 2, /threshold_occ, threshold_val = 20)
        corrfluxerr = fs        ;leave out the pmap err for now
     
 ;--------------------------------

        if  i eq 0 then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
           corrfluxarr = corrflux[1:*]
           corrfluxerrarr = corrfluxerr[1:*]
           timearr = sclkarr[1:*]        
           bmjd = bmjdarr[1:*]
           utcs = utcsarr[1:*]
           backarr = back[1:*]
           backerrarr = backerr[1:*]
           xgcntrdarr = x_gcntrd[1:*]
           ygcntrdarr = y_gcntrd[1:*]
           xcntrdarr = x_cntrd[1:*]
           ycntrdarr = y_cntrd[1:*]
  ;         xfwhmarr = xfw
 ;         yfwhmarr = yfw

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
           backarr = [backarr, back[1:*]]
           backerrarr = [backerrarr, backerr[1:*]]
           xgcntrdarr = [xgcntrdarr, x_gcntrd[1:*]]
           ygcntrdarr = [ygcntrdarr,y_gcntrd[1:*]]
           xcntrdarr = [xcntrdarr,x_cntrd[1:*]]
           ycntrdarr = [ycntrdarr,y_cntrd[1:*]]
;           xfwhmarr = [xfwhmarr, xfw[1:*]]
;           yfwhmarr = [yfwhmarr, yfw[1:*]]
        endelse


        ;make a light curve of a single pixel.
;        lc(count) = im[29,5]
;        count = count + 1

     endfor                     ; for each fits file in the AOR
     p1 = plot(xarr, yarr, '1.', xtitle = 'x', ytitle = 'y')
     p2 = plot(bmjd - bmjd(0), yarr, '1.', xtitle = 'time', ytitle = 'y')
     p3 = plot(bmjd- bmjd(0), corrfluxarr, '1.', xtitle = 'time', ytitle = 'corrflux')
     p4 = plot(timearr - timearr(0), xarr, '1,', xtitle = 'time', ytitle = 'X ')
;     p4 = plot(timearr, lc, '1.', xtitle = ' Time', ytitle = 'SB MJy/sr')
     ;fill in that structure
;     print, 'n xarr, fluxarr, nfits', n_elements(xarr), n_elements(fluxarr), nfits
     AOR55cnc[a] ={ob55cnc,  ra_ref,dec_ref,xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a),bmjd,utcs, backarr, backerrarr, xgcntrdarr, ygcntrdarr, xcntrdarr, ycntrdarr};, xfwhmarr, yfwhmarr


;-------------------------------------------------------
;do some background checks
;-------------------------------------------------------

;first check for NANs.
        nan = where(finite(im) lt 1,nancount)
        if nancount gt 0 then print,  'nancount', nancount


        

  endfor                        ;for each AOR
  
  ;test
  print, AOR55cnc[0].ra, AOR55cnc[0].xcen[26]
  
  save, AOR55cnc, filename='/Users/jkrick/irac_warm/pcrs_planets/55cnc/55cnc_phot_ch2_42000384.sav'
;  undefine, AOR55cnc
  
;sort by time
;a = sort(timearr)
;sorttime = timearr[a]
;sortxarr = xarr[a]
;sortyarr = yarr[a]
;sortfluxarr = fluxarr[a]
;sortfluxerrarr = fluxerrarr[a]
;sortcorrfluxarr = corrfluxarr[a]
;sortcorrfluxerrarr = corrfluxerrarr[a]
;sortbmjd = bmjd[a]
;sortutcs = utcs[a]

;save, /all, filename = '/Users/jkrick/irac_warm/pcrs_planets/55cnc/55cnc_phot_ch2.sav'
end


  ;read in all the pixel phase correction data
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/pmap_ch2_0p4s_500x500_0043_111206.fits', pmapdata, pmapheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/psigma_ch2_0p4s_500x500_0043_111206.fits', psigmadata, psigmaheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/occu_ch2_0p4s_500x500_0043_111206.fits', occdata, occheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/xgrid_ch2_0p4s_500x500_0043_111206.fits', xgriddata, xgridheader
;  fits_read, '/Users/jkrick/irac_warm/pmap/pmap_ch2_0p4s_500x500_0043_111206/ygrid_ch2_0p4s_500x500_0043_111206.fits', ygriddata, ygridheader

        ;corrflux = pixel_phase_correct_gauss(abcdflux,x_center,y_center,2,  '3_12_20')
 
;        corrresult = interp2d(pmapdata,xgriddata,ygriddata,x_center,y_center,/regular)
;        pmaperr = interp2d(psigmadata,xgriddata,ygriddata,x_center,y_center,/regular)
;        occresult = interp2d(occdata,xgriddata,ygriddata,x_center,y_center,/regular)
;        corrflux = abcdflux / corrresult
;        corrfluxerr = sqrt(fs^2 + pmaperr^2)

       ;then nan out those which do not have a high enough occupation
;        for ncor = 0, n_elements(abcdflux) - 1 do begin
        
;           if occresult(ncor) lt 40. then begin ; 1E6 then begin
;              ;print, 'out of range', fitsname(i), ncor
;              corrflux(ncor) = abcdflux(ncor) / alog10(-1)
;              corrfluxerr(ncor)= abcdflux(ncor)/ alog10(-1)
;           endif
           
;        endfor
