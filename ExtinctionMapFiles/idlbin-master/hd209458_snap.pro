pro hd209458_snap
  colorarr = [  'SILVER'   , 'SKY_BLUE'  , 'SLATE_BLUE' , 'SLATE_GREY' , ' red',  'PLUM' ,  'SPRING_GREEN'  ,  'STEEL_BLUE'   ,  'green', 'TEAL' ,  'THISTLE' ,  'TOMATO'  ,  'TURQUOISE' ,  'VIOLET' , 'YELLOW'  ,  'YELLOW_GREEN', 'PAPAYA_WHIP'  ,' PEACH_PUFF','   PERU',   ' POWDER_BLUE' ]

 aorname = ['r45188864','r45189120','r45189376','r45189632','r45189888','r45190144','r45190400','r45190656','r45190912','r45191168','r45191424','r45191680','r45191936','r45192192','r45192704','r45195264','r45192960','r45193216','r45193472','r45193984','r45193728','r45195520','r45194240','r45194496','r45194752','r45195008','r45196288','r45195776','r45197312','r45196032','r45196544','r45196800','r45197056','r45197568','r45197824','r45198080','r45192448']



;m = pp_multiplot(multi_layout=[n_elements(aorname),2]);,global_xtitle='Test X axis title',global_ytitle='Test Y axis title')

 nfits = 60*63.
 print, 'nfits', nfits
 snapshots = replicate({snapob, ra:0D,dec:0D,xcen:fltarr(nfits),ycen:fltarr(nfits),flux:fltarr(nfits),fluxerr:fltarr(nfits), corrflux:fltarr(nfits), corrfluxerr:fltarr(nfits), sclktime_0:0D, timearr:dblarr(nfits),aor:' ', peakpix:fltarr(nfits), bmjd:dblarr(nfits), utcs:dblarr(nfits)},n_elements(aorname))


  for a =0,  n_elements(aorname) - 1 do begin
     print, 'working on ',aorname(a)
     
     bcddir = '/Users/jkrick/irac_warm/hd209458/' + string(aorname(a)) +'/ch2/bcd/'
     rawdir = '/Users/jkrick/irac_warm/hd209458/' + string(aorname(a)) + '/ch2/raw/'
 
     CD, bcddir                    ; change directories to the correct AOR directory
     command  =  " ls *_bcd.fits > /Users/jkrick/irac_warm/hd209458/bcdlist.txt"
     spawn, command
     command2 =  "ls *_bunc.fits > /Users/jkrick/irac_warm/hd209458/bunclist.txt" ;this really should be the bunc file
     spawn, command2
     CD, rawdir                    ; change directories to the correct AOR directory
     command  =  " ls *dce.fits > /Users/jkrick/irac_warm/hd209458/rawlist.txt"  ;this really should be the raw file
     spawn, command
 
     readcol,'/Users/jkrick/irac_warm/hd209458/rawlist.txt',rawname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/hd209458/bcdlist.txt',fitsname, format = 'A', /silent
     readcol,'/Users/jkrick/irac_warm/hd209458/bunclist.txt',buncname, format = 'A', /silent

     yarr = fltarr(n_elements(fitsname)*63)
     xarr = fltarr(n_elements(fitsname)*63)
     timearr = fltarr(n_elements(fitsname)*63)
     fluxarr = fltarr(n_elements(fitsname)*63)
     fluxerrarr = fltarr(n_elements(fitsname)*63)
     corrfluxarr = fltarr(n_elements(fitsname)*63)
     corrfluxerrarr = fltarr(n_elements(fitsname)*63)
     peakarr = fltarr(n_elements(fitsname)*63)
     print, 'xarr', n_elements(xarr)
     
     for i =0.D, n_elements(fitsname) - 1 do begin ;read each bcd file, find centroid, keep track
        ;print, 'working on ', fitsname(i)         
                                ;now were did it really point to?
        cd, bcddir
        header = headfits(fitsname(i)) ;
        sclk_obs= sxpar(header, 'SCLK_OBS')
        frametime = sxpar(header, 'FRAMTIME')
        bmjd_obs = double(sxpar(header, 'BMJD_OBS'))
        ;print, ' bmjdobs', bmjd_obs, format = '(A, F0)'
        utcs_obs = double(sxpar(header, 'UTCS_OBS'))
        ;print, 'frametime', frametime
        if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

        sclkarr = fltarr(64)
        bmjdarr = dblarr(64)
        utcsarr = dblarr(64)
        for nt = 0, 64 - 1 do begin
           sclkarr[nt] = (sclk_obs ) + 0.5*frametime + frametime*nt
           bmjdarr[nt]= bmjd_obs + (0.5*frametime + frametime*nt) / 60./60./24.   ;in days
           utcsarr[nt]= utcs_obs + (0.5*frametime + frametime*nt)/ 60./60./24. ;in days
           ;print, 'bmjdinside ', bmjd_obs,' ', bmjd_obs + (0.5*frametime + frametime*nt) / 60./60./24.,' ', bmjdarr[nt], format = '(A, F0, A,F0,A, F0)'
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
;--------------------------------
      ;correct for pixel phase effect based on pmaps from Jim
        
        corrflux = iracpc_pmap_corr(abcdflux, x_center, y_center, 2, /threshold_occ, threshold_val = 20)
        corrfluxerr = fs        ;leave out the pmap err for now
     
 ;--------------------------------    
        ;print out some summaries
;        print, 'x, y, flux, correction', mean(x_center), mean(y_center), mean(abcdflux), mean(corrresult), format = '(A, F10.4, F10.4, F10.4, F10.4)'

        ;start by applying the correction to all 65 subarray images per frame
;        corrflux = abcdflux / corrresult
;        corrfluxerr = sqrt(fs^2 + pmaperr^2)
     

        ;----------------------------
        ;add info about the peak pixel flux
        ;print, 'working on raw', rawname(i)
        cd, rawdir
        fits_read, rawname(i), rawdata, rawheader
        mc = reform(rawdata, 32, 32, 64)
        img = mc
; convert to real
        img=img*1.
; fix the problem with unsigned ints
        fix=where((img LT 0),count)
        if (count GT 0) then img[fix]=img[fix]+65536.
; flip the InSb
        img=65535.-img

        peaks = img[15, 16,*]
       ;array of 64 peak pixel values is  img[15,16]
       ; print, 'mean peaks', mean(peaks, /nan)
        ;----------------------------

    ; print, 'finished raw'

       if i eq 0 then begin
           xarr = x_center[1:*]
           yarr = y_center[1:*]
           fluxarr = abcdflux[1:*]
           fluxerrarr = fs[1:*]
           corrfluxarr = corrflux[1:*]
           corrfluxerrarr = corrfluxerr[1:*]
           peakarr = peaks[1:*]
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
           ;print, 'filling peakarr'
           peakarr = [peakarr, peaks[1:*]]
           ;print, 'done filling peakarr'
           timearr = [timearr, sclkarr[1:*]]
           bmjd = double([bmjd, bmjdarr[1:*]])
           utcs = [utcs, utcsarr[1:*]]
        endelse
      ;print, 'finished filling arrays'

     endfor                     ; for each fits file in the AOR
     ;fill in that structure
     help, bmjd
     
     snapshots[a] ={snapob,  ra_ref,dec_ref,xarr,yarr,fluxarr,fluxerrarr, corrfluxarr, corrfluxerrarr, sclk_0, timearr,aorname(a), peakarr, bmjd, utcs}
    ; for bi = 0, n_elements(bmjd) - 1 do  print,'bmjd', snapshots[a].bmjd[bi], format = '(A, F0)'
     ;print, 'finished filling structure'
  endfor                        ;for each AOR
  
  ;test
  ;print, snapshots[0].ra, snapshots[0].xcen[26]
  
  save, snapshots, filename='/Users/jkrick/irac_warm/hd209458/hd209458_snap_082812.sav'
  undefine, snapshots
  

end
