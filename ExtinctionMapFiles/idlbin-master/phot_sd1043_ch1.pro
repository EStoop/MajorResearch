pro phot_sd1043_ch1
 ;from PID80179

; staring mode ch1
;simbad coords
;ra_ref = 160.89616667
;dec_ref = 12.22083333

;coords from the image
ra_ref = 160.89641
dec_ref = 12.219967

  aorname = ['r45758208', 'r43336448'] ;ch2
  nfits = 3802.     

  AORsd1043 = replicate({obsd1043, ra:0D,dec:0D,xcen:dblarr(nfits),ycen:dblarr(nfits),flux:dblarr(nfits),fluxerr:dblarr(nfits), corrflux:dblarr(nfits), corrfluxerr:dblarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', bmjdarr:dblarr(nfits), utcsarr:dblarr(nfits), bkgd:dblarr(nfits), bkgderr:dblarr(nfits)},n_elements(aorname)) ;, xfwhmarr:dblarr(nfits), yfwhmarr:dblarr(nfits

 

  for a = 0,   n_elements(aorname) - 1 do begin
     nfits = [3802., 3802.]

     xarr = fltarr(nfits(a))
     yarr = fltarr(nfits(a))
     fluxarr = fltarr(nfits(a))
     fluxerrarr = fltarr(nfits(a))
     corrfluxarr = fltarr(nfits(a))
     corrfluxerrarr = fltarr(nfits(a))
     timearr = dblarr(nfits(a))
     bmjd = dblarr(nfits(a))
     utcs = dblarr(nfits(a))
     backarr = fltarr(nfits(a))
     backerrarr =fltarr(nfits(a))
    
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/pcrs_planets/sd1043/'+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  =  "find ch1/bcd -name '*_cbcd.fits' > /Users/jkrick/irac_warm/pcrs_planets/sd1043/cbcdlist.txt"
     spawn, command
     command2 =  "find  ch1/bcd -name '*_cbunc.fits' > /Users/jkrick/irac_warm/pcrs_planets/sd1043/bunclist.txt"
     spawn, command2

     readcol,'/Users/jkrick/irac_warm/pcrs_planets/sd1043/cbcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/pcrs_planets/sd1043/bunclist.txt',buncname, format = 'A', /silent

     for i =0.D, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        frametime = sxpar(header, 'FRAMTIME')
        bmjd_obs = sxpar(header, 'BMJD_OBS')
        utcs_obs = sxpar(header, 'UTCS_OBS')
        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

;        sclkarr = dblarr(64)
;        bmjdarr = dblarr(64)
;        utcsarr = dblarr(64)
;        for nt = 0., 64 - 1 do begin
;           sclkarr[nt] = (sclk_obs ) + 0.5*frametime + frametime*nt
;           bmjdarr[nt]= bmjd_obs + 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt
;           utcsarr[nt]= utcs_obs + 0.5*(frametime/60./60./24.)  + (frametime/60./60./24.) *nt
;           ;print, 'bmjd stuff', nt, frametime, ' ', bmjd_obs,  ' ', bmjd_obs + 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt,' ', bmjdarr[nt],  format = '(A, F10.2, F10.2,A, F0, A, F0,A,F0)'
;        endfor
;        for na = 0, n_elements(bmjdarr) - 1 do print, 'bmjdarr ', bmjdarr(na), format = '(A, F0)'

        ;do the photometry
   
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
        back =  b[*,1] ;for 3, 11-15.5 ;b[*,2]  
        backerr = bs[*,1]
        ;xfw = xfwhm
        ;yfw = yfwhm
;--------------------------------
      ;correct for pixel phase effect based on pmaps from Jim
        
        corrflux = iracpc_pmap_corr(abcdflux, x_center, y_center, 1, /threshold_occ, threshold_val = 20)
        corrfluxerr = fs        ;leave out the pmap err for now
     
 ;--------------------------------

           xarr[i] = x_center
           yarr[i] = y_center
           fluxarr[i] = abcdflux
           fluxerrarr[i] = fs
           corrfluxarr[i] = corrflux
           corrfluxerrarr[i] = corrfluxerr
           timearr[i] = sclk_obs  
           bmjd[i] = bmjd_obs
           utcs[i] = utcs_obs
           backarr[i] = back
           backerrarr[i] = backerr
  ;         xfwhmarr = xfw
 ;         yfwhmarr = yfw

      
     endfor                     ; for each fits file in the AOR
     p1 = plot(xarr, yarr, '1.')
     p2 = plot(bmjd - bmjd(0), yarr, '1.')
     p3 = plot(bmjd- bmjd(0), corrfluxarr, '1.')

     ;fill in that structure
;     print, 'n xarr, fluxarr, nfits', n_elements(xarr), n_elements(fluxarr), nfits
     AORsd1043[a] ={obsd1043,  ra_ref,dec_ref,xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a),bmjd,utcs, backarr, backerrarr};, xfwhmarr, yfwhmarr

  endfor                        ;for each AOR
  
  ;test
  ;print, AORsd1043[0].ra, AOR55cnc[0].xcen[26]
  
  save, AORsd1043, filename='/Users/jkrick/irac_warm/pcrs_planets/sd1043/sd1043_phot_ch1.sav'
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
