pro phot_S19_2, chname, normalize = normalize

  ;;phot_S19_2, 'ch1'

  if chname eq 'ch1' then startnum = 0 else startnum = 1
  t1 = systime(1)
  dirname =  '/Users/jkrick/external/S19.2repro/'
  cd, dirname
  
;;for array dependant photometric correction warm
  fits_read, '/Users/jkrick/irac_warm/calstars/arrayloccorr/ch1_photcorr_ap_5.fits', photcor_ch1_s191, photcorhead_ch1
  fits_read, '/Users/jkrick/irac_warm/calstars/arrayloccorr/ch2_photcorr_ap_5.fits', photcor_ch2_s191, photcorhead_ch2
  fits_read, '/Users/jkrick/irac_warm/calstars/arrayloccorr/S19_2/ch1_photcorr_rj_warm.fits', photcor_ch1_s192, photcorhead_ch1
  fits_read, '/Users/jkrick/irac_warm/calstars/arrayloccorr/S19_2/ch2_photcorr_rj_warm.fits', photcor_ch2_s192, photcorhead_ch2

;;for pixel phase correction S19.2
  S192params = {DELTAF_X:[0.045065252D,0.019332239D], DELTAF_Y:[0.051759456D,0.013922875D], X0:[0.13051077D,0.061947404D], Y0:[0.043896764D,0.11991437D], SIGMA_X:[0.17571139D,0.17625694D], SIGMA_Y:[0.15926883D,0.15129765D], F0:[0.95961049D,0.98622255D], THETA:[0.,0.]}


;;set bin_sizes up front
  bsize = 1E-2

  aorname = '0045816576'
  command = 'ls -d 00*/ > aorlist.txt'
;  spawn, command
  readcol, 'aorlist.txt', aorname, format = 'A', /silent
  count19_2 = 0
  medarr = fltarr(n_elements(aorname))
  starname = strarr(n_elements(aorname))

  ;;aperture correction for [3,3-7] from Jim
  aper_corr = [1.1233D,1.1336D]
  counta = 0
  for a = 0, n_elements(aorname) -1 do begin
     undefine, fitsname
     undefine, buncname
     print, 'working on aorname ', aorname(a) , a, n_elements(aorname)
     cd, aorname(a)
     command1 = 'ls IRAC.1*.bcd_fp.fits > ch1bcdnames.txt'
     command2 = 'ls IRAC.2*.bcd_fp.fits > ch2bcdnames.txt'
     command3 =  'ls IRAC.1*sig_dntoflux.fits > ch1buncnames.txt'
     command4 =  'ls IRAC.2*sig_dntoflux.fits > ch2buncnames.txt'
     spawn, command1 & spawn, command2 & spawn, command3 & spawn, command4

     readcol, strcompress(chname + 'bcdnames.txt',/remove_all), fitsname, format = 'A', /silent, count = bcdcount
     readcol, strcompress(chname + 'buncnames.txt',/remove_all), buncname, format = 'A', /silent
;     print, 'n count', bcdcount
     if bcdcount gt 0 then begin   ; have data in that channel
        header = headfits(fitsname(0))  
        ra_ref = sxpar(header, 'RA_REF')
        dec_ref = sxpar(header, 'DEC_REF')
        AORLABEL= sxpar(header, 'AORLABEL')
        chnlnum = sxpar(header, 'CHNLNUM')
        naxis = sxpar(header, 'NAXIS')
        
        if strmid(AORLABEL, 0, 12) eq 'IRAC_calstar' then begin ;and NAXIS lt 3 then begin 
           print, 'inside calstar', strmid(AORLABEL, 13, 7)
           ;;find truth fluxes
           ;;taken from Jim's evernote note "How to Correct for
           ;;Position Dependent Effects"
           CASE strmid(AORLABEL, 13, 7) OF
              'NPM1p64': truth = [63.234,  41.131]
              'KF06T1_': truth = [13.289, 8.000]
              'HD16545': truth = [649.541, 422.219]
              'NPM1p60': truth = [38.328,  24.806]
              '1812095': truth = [8.711,   5.677]
              'KF08T3_': truth = [11.908,   7.262]
              'KF09T1_': truth = [166.614, 102.136]
              'KF06T2_': truth = [10.256,   6.174]
              'NPM1p67': truth = [807.331, 482.708]
              'NPM1p68': truth = [556.244, 327.106]
              'NPM1p66': truth = [141.357,  82.665]
              ELSE: BREAK       ; we don't have a truth flux for other stars
           ENDCASE

           ;;setup arrays to save data
           if NAXIS lt 3 then begin
              fluxarr = fltarr(bcdcount)
           endif else begin
              fluxarr = fltarr(bcdcount*64)
           endelse

           c = 0
           for i = 0, n_elements(fitsname) - 1 do begin
              ;;make sure it is on the frame
              header = headfits(fitsname(i))  
              ADXY, header, ra_ref, dec_ref, xcen, ycen
                 
              if xcen gt 5 and ycen gt 5 and xcen lt 250 and ycen lt 250 then begin
                 count19_2++
                 fits_read, fitsname(i), im, h
                 fits_read, buncname(i), unc, hunc, /no_abort
                 get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                              x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                              x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                              xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                              xfwhm, yfwhm, /WARM,/silent
                 bcdflux = f[*,3]      ; I think this is [3,3-7]
                 ;;make a correction for pixel phase 
;                 print, 'n fluxes', n_elements(f[*,3]), n_elements(x3), n_elements(y3)
                 corrflux = pixel_phase_correct_gauss(f[*,3],x3,y3,chnlnum, '3_3_7',params = S192params)

                 ;;apply array dependent correction
;                 arraycorr = fltarr(n_elements(x3))
;                 for ipc = 0, n_elements(x3) - 1 do begin
;                    print, 'chnlnum', ipc, chnlnum, x3(ipc), y3(ipc)
;                    if chnlnum eq '1' then arraycorr[ipc] = photcor_ch1_s192(x3(ipc), y3(ipc)) else $
;                       arraycorr[ipc] = photcor_ch2_s192(x3(ipc), y3(ipc))
;                 endfor
                 if chnlnum(0) eq '1' then arraycorr = photcor_ch1_s192[x3, y3] else arraycorr = photcor_ch2_s192[x3, y3]
                 corrflux= corrflux * arraycorr
                 ;;save them
                 ;;consider normalizing by truth fluxes
;              xcenarr[c]  = x3 &  ycenarr[c] = y3 & fluxerrarr[c] =
;              fs[0] & backarr[c]= b[0] 
                 fluxarr[c] =  corrflux ;bcdflux;

                 if NAXIS lt 3 then begin
                    c = c + 1
                 endif else begin
                    c = c + 64
                 endelse

              endif             ;if the star is on the frame
           endfor               ; for each fits image
           ;;remove known nans = off the frame
           fluxarr = fluxarr[0:c-1]

           ;;apply aperture correction
           fluxarr = fluxarr*aper_corr(startnum)

           ;;normalize
            if keyword_set(normalize) then begin
               ;;normfactor = mean(fluxarr, /nan)
               normfactor = truth(chnlnum - 1)/ 1E3
;               print, 'before', fluxarr, 'norm', normfactor
               fluxarr = fluxarr / normfactor
;               print, 'after', fluxarr
            endif

           ;;fill an array of all calstars in this channel
           if (allfluxarr eq !NULL) gt 0 then begin
              allfluxarr = fluxarr
           endif else begin
              allfluxarr = [allfluxarr, fluxarr]
           endelse

           ;;plot median normalized value of each individual star
           medarr[counta] = median(fluxarr)
           starname[counta] = strmid(AORLABEL, 13, 7)
           counta = counta + 1

        endif                   ; if aorlabel eq irac_calstar

     endif                      ; if this channel has any data     

     
     cd, dirname                ;reset for the next AOR
  endfor                        ; for each AOR
     
  ;;make a histogram of the delta fluxes.
  good = where(finite(allfluxarr) gt 0)
  allfluxarr = allfluxarr(good)
  print, 'allfluxarr', allfluxarr
  plothist, allfluxarr, xhist, yhist, bin = bsize, /noplot, /noprint
  p1 = barplot(xhist, yhist,  fill_color = 'purple', index = 0, nbars = 2, xtitle = 'Normalized Corrected Flux',$
               ytitle = 'Number', transparency = 50,yrange = [0,30], xrange = [0.8,0.95], title = chname) ;
  start = [median(allfluxarr),0.05,count19_2]
  noise = fltarr(n_elements(xhist))
  noise[*] = 1.0
  gresult= MPFITFUN('mygauss',xhist,yhist, noise, start)  
  xp = findgen(100) / 300. + 0.80
  p1 = plot(xp, gresult(2)/sqrt(2.*!Pi) * exp(-0.5*((xp - gresult(0))/gresult(1))^2.), color = 'purple', thick = 2, overplot = p1)

  ;;find the mean/median/mode and plot that too
  p1 = plot([median(allfluxarr), median(allfluxarr)], [20,25], color = 'purple', thick=2, overplot = p1)

  medarr = medarr[0:counta -1]
  starname = starname[0:counta -1]
  xint = indgen(n_elements(medarr))
  pcompare = plot(xint, medarr, '1s', xtitle = 'star', ytitle = 'Median(Flux/Truth)', title = chname, $
                  xminor = 0, sym_filled = 1, yrange = [0.95,1.05] , xtickname = starname, $
                  xtext_orientation = 45, xrange = [0, n_elements(medarr) - 1])
  bad = where(medarr lt 0.9)
  medarr[bad] = alog10(-1)
  meanval = mean(medarr, /nan)
  pcompare = plot(findgen(20), fltarr(20) + meanval, thick = 2, color = 'gold', overplot = pcompare)
  pcompare = plot(findgen(20), fltarr(20) + 1.0, thick = 2, color = 'grey', overplot = pcompare)

  for ic = 0, counta - 1 do print, starname[ic], medarr[ic]
;;---------------------------------------------------------------------------------------------
;;---------------------------------------------------------------------------------------------
  ;;do photometry on the S19.1 stuff
  ;;taken mostly from derive_fluxconv_6
  print, '--------------------------------------'

  ;;what size aperture are we using for photometry?
;  aprad = 5
  apcorr = [1.049,1.050]
  ;;for array dependant photometric correction
  fits_read, '/Users/jkrick/external/S19.2repro/S19.1SHA/ch1_photcorr_ap_5.fits', photcor_warm_ch1, photcorhead_warm_ch1
  fits_read, '/Users/jkrick/external/S19.2repro/S19.1SHA/ch2_photcorr_ap_5.fits', photcor_warm_ch2, photcorhead_warm_ch2

  dirloc = ['/Users/jkrick/external/S19.2repro/S19.1SHA/']
  cd, dirloc
  command = 'ls -d */ > /Users/jkrick/external/S19.2repro/S19.1SHA/aorlist.txt'
;  spawn, command
  readcol, '/Users/jkrick/external/S19.2repro/S19.1SHA/aorlist.txt', aorname, format = 'A', /silent
  

  for ch = startnum, startnum do begin
;     print, 'Channel', ch+1
     count19_1 = 0
     
     for a = 0, n_elements(aorname) - 1 do begin
        print, 'working on aorname ',aorname(a)
        aname = strcompress(dirloc + aorname(a),/remove_all)
        cd, aname
        pwd
        corrected_fluxarr = fltarr(10000)
        j = 0
        
        if ch eq 0 then begin
           command1  =  "ls ch1/bcd/*_bcd.fits > /Users/jkrick/irac_warm/list.txt"
           command2  =  "ls ch1/bcd/*_bunc.fits > /Users/jkrick/irac_warm/unclist.txt"
           spawn, command1 & spawn, command2
        endif
        if ch eq 1 then begin
           command1  =  "ls ch2/bcd/*_bcd.fits > /Users/jkrick/irac_warm/list.txt"
           command2  =  "ls ch2/bcd/*_bunc.fits > /Users/jkrick/irac_warm/unclist.txt"
           spawn, command1 & spawn, command2
        endif
        readcol,'/Users/jkrick/irac_warm/list.txt', fitsname, format = 'A', /silent, count = bcdcount
        readcol,'/Users/jkrick/irac_warm/unclist.txt', uncname, format = 'A', /silent

        if bcdcount gt 0 then begin ; have data in that channel
           header = headfits(fitsname(0))  
           ra_ref = sxpar(header, 'RA_REF')
           dec_ref = sxpar(header, 'DEC_REF')
           AORLABEL= sxpar(header, 'AORLABEL')
           chnlnum = sxpar(header, 'CHNLNUM')
           naxis = sxpar(header, 'NAXIS')
           
           if strmid(AORLABEL, 0, 12) eq 'IRAC_calstar' and NAXIS lt 3 then begin 
              print, 'inside S19.1 calstar'
              ;;find truth fluxes
              ;;taken from Jim's evernote note "How to Correct for
              ;;Position Dependent Effects"
              CASE strmid(AORLABEL, 13, 7) OF
                 'NPM1p64': truth = [63.234,  41.131]
                 'KF06T1_': truth = [13.289, 8.000]
                 'HD16545': truth = [649.541, 422.219]
                 'NPM1p60': truth = [38.328,  24.806]
                 '1812095': truth = [8.711,   5.677]
                 'KF08T3_': truth = [11.908,   7.262]
                 'KF09T1_': truth = [166.614, 102.136]
                 'KF06T2_': truth = [10.256,   6.174]
                 'NPM1p67': truth = [807.331, 482.708]
                 'NPM1p68': truth = [556.244, 327.106]
                 'NPM1p66': truth = [141.357,  82.665]
                 ELSE: BREAK ; we don't have a truth flux for other stars
              ENDCASE
              print, 'S19.1', strmid(AORLABEL, 13, 7), truth(0)
              
              for i =1, n_elements(fitsname) - 1 do begin
                 fits_read,fitsname(i), bcddata, header
                 fits_read, uncname(i), unc, uheader
                 adxy, header, ra_ref, dec_ref, bcdx, bcdy
              
                 if bcdx gt 10 and bcdy gt 10 and bcdy lt 240 and bcdx lt 240 then begin
                    count19_1++
                    ;;do the photometry
                    get_centroids_for_calstar_jk,bcddata, header, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                                 x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                                 x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                                 xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                                 xfwhm, yfwhm, /WARM,/silent
                    bcdflux = f[3]
;               get_centroids, fitsname(i), t, dt, bcdxcen, bcdycen, bcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = aprad, RA=ra(r), DEC=dec(r),/silent
                    
                    ;;make a correction for pixel phase 
                                ;pixap = strcompress(string(aprad) + '_12_20',/remove_all)
                    corrected_flux = pixel_phase_correct_gauss(f[3],x3,y3,ch+1,'3_3_7')
                    
                    ;;apply array dependent correction
                    photcor_ch1 = photcor_warm_ch1(x3, y3)
                    photcor_ch2 = photcor_warm_ch2(x3, y3)
                    if ch eq 0 then corrected_flux = corrected_flux * photcor_ch1
                    if ch eq 1 then corrected_flux = corrected_flux * photcor_ch2
                    
                    ;;aperture correction for small aperture
                                ;if aprad eq 5 then
;                    corrected_flux = corrected_flux * apcorr(ch)
                    
                                ; fluxes off by a factor of 100 (mJy)     
                    corrected_fluxarr[j] =   corrected_flux ;bcdflux; *1000. ; in mJy
                    
                    j = j + 1
                 endif
                 
              endfor            ;end for each fits image
           ;;remove known nans = off the frame
           corrected_fluxarr = corrected_fluxarr[0:j-1]

           ;;normalize
           if keyword_set(normalize) then begin
              ;normfactor = mean(corrected_fluxarr, /nan)
              normfactor = truth(chnlnum - 1)/1E3
              corrected_fluxarr = corrected_fluxarr / normfactor
           endif


;           print, 'corrected_fluxarr', corrected_fluxarr

           ;;fill an array of all calstars in this channel
           if (allfluxarrold eq !NULL) gt 0 then begin
              allfluxarrold = corrected_fluxarr
           endif else begin
              allfluxarrold = [allfluxarrold, corrected_fluxarr]
           endelse

           endif  ;calstar
        endif  ; has data in that channel

     endfor                     ;end for each aorname

     ;;make a histogram of the delta fluxes.
     good = where(finite(allfluxarrold) gt 0)
     allfluxarrold = allfluxarrold(good)
     print, 'S19.1 allfluxarr', allfluxarrold
     plothist, allfluxarrold, xhist, yhist, bin = bsize, /noplot, /noprint
     p1 = barplot(xhist, yhist,  fill_color = 'orange red', index = 1, nbars = 2, overplot = p1, transparency = 50) ;,
     start = [median(allfluxarrold),0.05,count19_1]
     noise = fltarr(n_elements(xhist))
     noise[*] = 1.0
     gresult= MPFITFUN('mygauss',xhist,yhist, noise, start)  
     p1 = plot(xp, gresult(2)/sqrt(2.*!Pi) * exp(-0.5*((xp - gresult(0))/gresult(1))^2.), color = 'orange red', thick = 2, overplot = p1)

     ;;find the mean/median/mode and plot that too
     p1 = plot([median(allfluxarrold), median(allfluxarrold)], [20, 25], color = 'orange red', thick=2, overplot = p1)
     
     ;;labels
     t = text(0.91, 20, 'S19.1', color = 'orange red', overplot = p1,/data)
     t = text(0.91, 18, 'S19.2', color = 'purple', overplot = p1,/data)
  endfor  ;for each channel
  
  print, 'count 19.1, 19.2', count19_1, ', ' , count19_2
  

 ;;what is the difference between the two fluxes
  ;;assumes these are in the same order, which they should be
  delta = allfluxarr - allfluxarrold
  print, 'delta', delta
  plothist, delta, xhist, yhist, bin = bsize/2., /noplot, /noprint
  p2 = barplot(xhist, yhist, fill_color = 'light blue', xtitle = 'Delta Photometry', ytitle = 'Number', $
               title = chname, xrange = [-0.1, 0.1])
  start = [median(delta),0.01,count19_1]
  noise = fltarr(n_elements(xhist))
  noise[*] = 1.0
  gresult= MPFITFUN('mygauss',xhist,yhist, noise, start)  
  xp = findgen(100) / 300 - 0.1
  p1 = plot(xp, gresult(2)/sqrt(2.*!Pi) * exp(-0.5*((xp - gresult(0))/gresult(1))^2.), thick = 2, overplot = p2)
  text2= text(.03, 15, string(gresult(0)),/data, overplot = p2)

  print, 'time check', systime(1) - t1, n_elements(delta)

end

     

  
; 45816064, 45817088, 45815296, 45816320, 45817344, 45815552,  45817600, 45815808, 45816832, 45821440


