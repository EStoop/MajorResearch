pro phot_wasp14_ch2
  

ra_ref = 218.276458
dec_ref = 21.894722

  aorname = ['r45426688',  'r45428224',  'r45428480',  'r45428736',  'r45428992'] ;ch2
  nfits = 360L*63L    
  print, 'nfits',nfits
  AORwasp14 = replicate({wasp14ob, ra:0D,dec:0D,xcen:dblarr(nfits),ycen:dblarr(nfits),flux:dblarr(nfits),fluxerr:dblarr(nfits), corrflux:dblarr(nfits), corrfluxerr:dblarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', bmjdarr:dblarr(nfits), utcsarr:dblarr(nfits), bkgd:dblarr(nfits), bkgderr:dblarr(nfits)},n_elements(aorname))


  for a = 0, n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     dir = '/Users/jkrick/irac_warm/pcrs_planets/wasp14/'+ string(aorname(a) ) 
     CD, dir                    ; change directories to the correct AOR directory
     command  =  "find ch2/bcd -name '*_bcd.fits' > /Users/jkrick/irac_warm/pcrs_planets/wasp14/bcdlist.txt"
     spawn, command
     command2 =  "find  ch2/bcd -name '*bunc.fits' > /Users/jkrick/irac_warm/pcrs_planets/wasp14/bunclist.txt"
     spawn, command2

     readcol,'/Users/jkrick/irac_warm/pcrs_planets/wasp14/bcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/pcrs_planets/wasp14/bunclist.txt',buncname, format = 'A', /silent

     for i =0.D,  n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        frametime = sxpar(header, 'FRAMTIME')
        bmjd_obs = sxpar(header, 'BMJD_OBS')
        utcs_obs = sxpar(header, 'UTCS_OBS')

        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

        sclkarr = dblarr(64)
        bmjdarr = fltarr(64)
        utcsarr = fltarr(64)
        for nt = 0., 64 - 1 do begin
           sclkarr[nt] = (sclk_obs ) + 0.5*frametime + frametime*nt
 ;          bmjdarr[nt]= bmjd_obs + 0.5*frametime + frametime*nt
 ;          utcsarr[nt]= utcs_obs + 0.5*frametime + frametime*nt
           bmjdarr[nt]= bmjd_obs + 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt
           utcsarr[nt]= utcs_obs + 0.5*(frametime/60./60./24.)  + (frametime/60./60./24.) *nt
        endfor
         
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
        back = b[*,2]
        backerr = bs[*,2]

;--------------------------------
      ;correct for pixel phase effect based on pmaps from Jim
        file_suffix = ['500x500_0043_120409.fits','0p1s_x4_500x500_0043_120531.fits']
        ;corrflux = iracpc_pmap_corr(abcdflux, x_center, y_center, 2, /threshold_occ, threshold_val = 20)
        corrflux = iracpc_pmap_corr(abcdflux,x_center,y_center,2,FILE_SUFFIX=file_suffix,/threshold_occ, threshold_val = 20)
        corrfluxerr = fs        ;leave out the pmap err for now
     
 ;--------------------------------

        if i eq 0 then begin
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
         endelse
        
     endfor                     ; for each fits file in the AOR
     p1 = plot(xarr, yarr, '1.')
     p2 = plot(timearr - timearr(0), yarr, '1.')
     p3 = plot(timearr- timearr(0), corrfluxarr, '1.')

     ;fill in that structure
    print, 'n xarr, fluxarr, nfits', n_elements(xarr), n_elements(fluxarr), nfits
    if n_elements(xarr) lt nfits then begin
       xarr2 = dblarr(nfits)
       yarr2 = dblarr(nfits)
       fluxarr2= dblarr(nfits)
       fluxerrarr2= dblarr(nfits)
       corrfluxarr2= dblarr(nfits)
       corrfluxerrarr2= dblarr(nfits)
       timearr2= dblarr(nfits)
       bmjd2 = dblarr(nfits)
       utcs2 = dblarr(nfits)
       back2 = dblarr(nfits)
       backerr2 = dblarr(nfits)

       xarr2[0:n_elements(xarr) - 1] = xarr
       yarr2[0:n_elements(xarr) - 1] = yarr
       fluxarr2[0:n_elements(xarr) - 1] = fluxarr
       fluxerrarr2[0:n_elements(xarr) - 1] = fluxerrarr
       corrfluxarr2[0:n_elements(xarr) - 1] = corrfluxarr
       corrfluxerrarr2[0:n_elements(xarr) - 1] = corrfluxerrarr
       timearr2[0:n_elements(xarr) - 1] = timearr
       bmjd2[0:n_elements(xarr) - 1] = bmjd
       utcs2[0:n_elements(xarr) - 1] = utcs
       back2[0:n_elements(xarr) - 1] = backarr
       backerr2[0:n_elements(xarr) - 1] = backerrarr

       AORwasp14[a] ={wasp14ob,  ra_ref,dec_ref,xarr2,yarr2,fluxarr2,fluxerrarr2, corrfluxarr2, corrfluxerrarr2, sclk_0, timearr2,aorname(a), bmjd2,utcs2,back2, backerr2}
    endif else begin
       AORwasp14[a] ={wasp14ob,  ra_ref,dec_ref,xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a),bmjd, utcs,backarr, backerrarr}
    endelse

  endfor                        ;for each AOR
  
  ;test
  print, AORwasp14[0].ra, AORwasp14[0].xcen[26]
  
  save, AORwasp14, filename='/Users/jkrick/irac_warm/pcrs_planets/wasp14/wasp14_phot_ch2.sav'
  undefine, AORwasp14
  
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

;save, /all, filename = '/Users/jkrick/irac_warm/pcrs_planets/wasp33/wasp33_phot_ch2.sav'
end



       ;then nan out those which do not have a high enough occupation
;        for ncor = 0, n_elements(abcdflux) - 1 do begin
        
;           if occresult(ncor) lt 40. then begin ; 1E6 then begin
;              ;print, 'out of range', fitsname(i), ncor
;              corrflux(ncor) = abcdflux(ncor) / alog10(-1)
;              corrfluxerr(ncor)= abcdflux(ncor)/ alog10(-1)
;           endif
           
;        endfor
