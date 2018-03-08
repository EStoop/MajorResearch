function do_calstar_photometry, ch, dirname, processing
  common track_block
  ft1 = systime(1)

  if processing eq 'SHA' then begin
     fullname = strcompress(dirname+'ch'+string(ch),/remove_all)
     exist_dir = chk_dir(fullname)
  endif
  
;;  if exist_dir eq 1 then begin  ;got a live one
  
     cd, dirname
     pwd
     if ch eq 1 then command ="find ./bcd/00*/ -name 'IRAC.1*bcd_fp.fits' > /Users/jkrick/irac_warm/calstars/allch1bcdlist.txt "
     if ch eq 2 then command ="find ./bcd/00*/ -name 'IRAC.2*bcd_fp.fits' > /Users/jkrick/irac_warm/calstars/allch2bcdlist.txt "
     if processing eq 'SHA' then begin
        if ch eq 1 then command ="find ./ch1/bcd/ -name 'SPITZER_I1*_bcd.fits' > /Users/jkrick/irac_warm/calstars/allch1bcdlist.txt "
        if ch eq 2 then command ="find ./ch2/bcd/ -name 'SPITZER_I2*_bcd.fits' > /Users/jkrick/irac_warm/calstars/allch2bcdlist.txt "
     endif
     
     spawn, command  
     readcol,strcompress('/Users/jkrick/irac_warm/calstars/allch'+string(ch)+'bcdlist.txt',/remove_all), fitsname, format = 'A', /silent
;  print, '      read in bcdlist', systime(1) - ft1
     nfits =  n_elements(fitsname) 
     
     ;;set up storage arrays
     ;;assume all sub array then cut down later
     xarr = fltarr(nfits*64)
     starnamearr = strarr(nfits*64)
     timearr = dblarr(nfits*64)
     yarr = xarr 
     fluxarr = xarr
     fluxerrarr = xarr
     backarr = xarr
     backerrarr = xarr
     corrfluxarr = xarr
     raarr = xarr
     decarr = xarr
     nparr = xarr
     xfwhmarr = xarr
     yfwhmarr = xarr
     raarr = timearr
     decarr = timearr
     bmjdarr = timearr    
     procverarr = starnamearr
     
     startfits = 0L
     stopfits =  nfits - 1
     i = 0L
     for c = startfits, stopfits do begin
                                ;print, '      starting another fits file', systime(1) - ft1
        
        header = headfits(fitsname(c)) ;
        NAXIS= sxpar(header, 'NAXIS')
        AORLABEL= sxpar(header, 'AORLABEL')
        
;       print, i, strmid(AORLABEL, 0, 12), naxis
        if strmid(AORLABEL, 0, 12) eq 'IRAC_calstar'  then begin ;got a good one  
           
;read in the image, get important info from the header
           fits_read,fitsname(c), im, h
           inter = strmid(fitsname(c), 0, 46)
           uncname = strcompress(inter + '*sig_dntoflux.fits',/remove_all)
           if processing eq 'SHA' then begin
              inter = strmid(fitsname(c), 0, 41)
              uncname = strcompress(inter + '*bunc.fits',/remove_all)
           endif
           
;;        print, 'uncname',i,  uncname, fitsname(i)
           fits_read, uncname, unc, hunc, /no_abort ; so it won't crash if the file isn't there but should use the last unc file.
           chnlnum = sxpar(h, 'CHNLNUM')
           ra_ref = sxpar(h, 'RA_REF')
           dec_ref = sxpar(h, 'DEC_REF')
           sos_ver = sxpar(h, 'SOS_VER')
           sclk_obs = sxpar(h, 'SCLK_OBS')
           aintbeg = sxpar(header, 'AINTBEG')
           atimeend = sxpar(header, 'ATIMEEND')
           bmjd_obs = sxpar(header, 'BMJD_OBS')
           proc_ver = sxpar(header, 'CREATOR') ;which processing version
           framtime = sxpar(header, 'FRAMTIME')
           starname = strmid(AORLABEL, 13, 7)
           campaign = sxpar(header, 'CAMPAIGN')
;make sure the target is on the frame
           ADXY, h, ra_ref, dec_ref, xcen, ycen
           if naxis eq 2 then begin
              xmax =250
              ymax = 250
           endif else begin
              xmax = 30
              ymax = 30
           endelse
                                ;if starname eq '1812095' then print, 'working on 1812095'
           if xcen gt 5 and ycen gt 5 and xcen lt xmax and ycen lt ymax then begin
;           print, 'starting photometry ' , i, c, fitsname(c), starname, uncname
                                ;do photometry
;           print, '      about to start get_centroids', systime(1) - ft1
              
              get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                           x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                           x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                           xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                           xfwhm, yfwhm, /WARM
;           print, '      finished get_centroids', systime(1) - ft1
              
                                ;choose the requested pixel aperture =
                                ;3.0 pixels to match
                                ;pixphasecorrectgauss
              abcdflux = f[*,4]      
;           print, 'flux', x3, y3, abcdflux
              fs = fs[*,4]
                                ; 3-7 pixel background
              back = b[*,0]
              backerr = bs[*,0]
              
                                ;make a correction for pixel phase 
              if naxis eq 3 then corrflux = pixel_phase_correct_gauss(abcdflux,x3,y3,chnlnum, '3_3_7',/subarray) $ 
              else corrflux = pixel_phase_correct_gauss(abcdflux,x3,y3,chnlnum, '3_3_7')
              
                                ;corrflux = abcdflux*ppcorr
                                ;apply array dependent correction
              arraycorr = photcordata(x3, y3) 
              corrflux= corrflux * arraycorr
;           print, '      finished both corrections', systime(1) - ft1
              
                                ;save them
              starnamearr[c] = strmid(AORLABEL, 13, 7)
              if naxis eq 3 then begin  
                 xarr[i] = x3
                 yarr[i] = y3
                 fluxarr[i] = abcdflux
                 fluxerrarr[i] = fs
                 corrfluxarr[i] = corrflux
                 backarr[i] = back
                 backerrarr[i] = backerr
                 nparr[i] = np
                 xfwhmarr[i] = xfwhm
                 yfwhmarr[i] = yfwhm
                 raarr[i:i+63] = ra_ref
                 decarr[i:i+63] = dec_ref
                 starnamearr[i:i+63] = starname
                 procverarr[i:i+63] = proc_ver
                 deltatime = (atimeend - aintbeg) / 64.D ; real time for each of the 64 frames
                 nt = dindgen(64)
                 timearr[i] = sclk_obs  + (deltatime*nt)/60./60./24.D ; 0.5*frametime + frametime*nt
                 bmjdarr[i] = bmjd_obs + (deltatime*nt)/60./60./24.D  ; 0.5*(frametime/60./60./24.) + (frametime/60./60./24.)*nt
                 i = i + 64
              endif 
              if naxis eq 2 then begin ; and i eq 0 then begin
                 xarr[i] = x3
                 yarr[i]  =  y3
                 fluxarr[i]  =  abcdflux
                 fluxerrarr[i]  =  fs
                 corrfluxarr[i]  = corrflux
                 timearr[i]  = sclk_obs
                 bmjdarr[i]  = bmjd_obs
                 backarr[i]  =  back
                 backerrarr[i]  = backerr
                 nparr[i]  = np
                 xfwhmarr[i] = xfwhm
                 yfwhmarr[i] = yfwhm
                 raarr[i] = ra_ref
                 decarr[i] = dec_ref
                 starnamearr[i] = starname
                 procverarr[i] = proc_ver
                 i = i + 1
              endif
                                ;print, 'xarr', xarr[0:100]
           endif                ; if the target is on the frame
           
        endif                   ; if  calstar
     endfor                     ; for each fits image
     print, 'final i', i, ' ', nfits
     final = i
     xarr = xarr[0:final-1] 
     yarr = yarr[0:final-1] 
     fluxarr = fluxarr[0:final-1] 
     corrfluxarr = corrfluxarr[0:final-1]
     fluxerrarr = fluxerrarr[0:final-1] 
     backarr = backarr[0:final-1] 
     backerrarr = backerrarr[0:final-1] 
     timearr = timearr[0:final-1] 
     bmjdarr = bmjdarr[0:final-1] 
     nparr = nparr[0:final-1] 
     xfwhmarr = xfwhmarr[0:final-1] 
     yfwhmar = yfwhmarr[0:final-1] 
     starnamearr = starnamearr[0:final-1]
     procverarr = procverarr[0:final-1]
     raarr = raarr[0:final-1] 
     decarr = decarr[0:final-1] 
     
     if i gt 0 then begin
        campaign_int = fix(strmid(campaign, 5))
        if campaignarr eq !NULL then campaignarr = campaign_int else campaignarr = [campaignarr, campaign_int]
     endif
;;  endif  ; does the directory exist?
     return, 0
  end
  
