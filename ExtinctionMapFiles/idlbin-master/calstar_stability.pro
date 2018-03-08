pro calstar_stability

 t1 = systime(1)

<<<<<<< HEAD
  dirloc = ['/Volumes/iracdata/flight/IWIC/']
  
;start by getting all of the fits files from ch1, then later sort based on  naxis and exptime and maybe star name, then add ch2
  cd, dirloc
;  command ="find ./IRAC03*/bcd/00*/ -name 'IRAC.1*bcd_fp.fits' >> /Users/jkrick/irac_warm/calstars/allch1bcdlist.txt "
;  spawn, command
;  command2 ="find ./IRAC03*/bcd/00*/ -name 'IRAC.1*sig_dntoflux.fits' >>  /Users/jkrick/irac_warm/calstars/allch1unclist.txt "
;  spawn, command2
  
  for ch = 0, 0 do begin
     
     readcol,strcompress('/Users/jkrick/irac_warm/calstars/allch'+string(ch + 1)+'bcdlist.txt',/remove_all), fitsname, format = 'A', /silent
;     readcol,strcompress('/Users/jkrick/irac_warm/calstars/allch'+ string(ch + 1)+'unclist.txt',/remove_all), uncname, format = 'A', /silent
     
     print, 'nfits', n_elements(fitsname) 
;set up storage arrays
     xcenarr = fltarr(n_elements(fitsname))
     ycenarr = xcenarr
     starnamearr = strarr(n_elements(fitsname))
     timearr = xcenarr
     fluxarr = xcenarr
     fluxerrarr = xcenarr
     backarr = xcenarr
     corrfluxarr = xcenarr
     raarr = xcenarr
     decarr = xcenarr
     
;for array dependant photometric correction warm
     fits_read, '/Users/jkrick/irac_warm/calstars/arrayloccorr/ch1_photcorr_ap_5.fits', photcor_ch1, photcorhead_ch1
     fits_read, '/Users/jkrick/irac_warm/calstars/arrayloccorr/ch2_photcorr_ap_5.fits', photcor_ch2, photcorhead_ch2

;filename saver
     if ch eq 1 then openw, outlun, '/Users/jkrick/irac_warm/calstars/ch1_calstar_bcdlist.txt', /get_lun
     if ch eq 2 then openw, outlun, '/Users/jkrick/irac_warm/calstars/ch2_calstar_bcdlist.txt', /get_lun

     startfits = 0L
     stopfits = n_elements(fitsname) - 1
     c = 0L
     for i= startfits, stopfits do begin
        ;track progress
<<<<<<< HEAD
        if i mod 10000 eq 0 then print, 'fits file ', i, ' ',  systime(1) - t1
;      print, 'working on fitsname', fitsname(i)
        header = headfits(fitsname(i)) ;
        NAXIS= sxpar(header, 'NAXIS')
        
        AORLABEL= sxpar(header, 'AORLABEL')
;       print, i, strmid(AORLABEL, 0, 12), naxis
                                ;cut on full array calstars only (for now)
<<<<<<< HEAD
        if strmid(AORLABEL, 0, 12) eq 'IRAC_calstar' and NAXIS lt 3 then begin            ;got a good one

;add this name to a list of filenames
;that are calstars to use in future
;iterations of this code
           printf, outlun, fitsname(i)
           
           fits_read,fitsname(i), im, h
           inter = strmid(fitsname(i), 0, 59)
           uncname = strcompress(inter + 'sig_dntoflux.fits',/remove_all)
;           print, 'uncname',i,  uncname, fitsname(i)
           fits_read, uncname, unc, hunc, /no_abort  ; so it won't crash if the file isn't there but should use the last unc file.
           
           chnlnum = sxpar(h, 'CHNLNUM')
           ra_ref = sxpar(h, 'RA_REF')
           dec_ref = sxpar(h, 'DEC_REF')
           sos_ver = sxpar(h, 'SOS_VER')

          ;make sure it is on the frame
           ADXY, h, ra_ref, dec_ref, xcen, ycen
           if xcen gt 5 and ycen gt 5 and xcen lt 250 and ycen lt 250 then begin
              timearr[c]  = sxpar(h, 'SCLK_OBS')
              raarr[c] = ra_ref
              decarr[c] = dec_ref
              starnamearr[c] = strmid(AORLABEL, 13, 7)
              
              get_centroids_for_calstar_jk,im, h, unc, ra_ref, dec_ref,  t, dt, hjd, xft, x3, y3, $
                                           x5, y5, x7, y7, xg, yg, xh, yh, f, b, x3s, y3s, x5s, y5s, $
                                           x7s, y7s, fs, bs, xp3, yp3, xp5, yp5, xp7, yp7, xp3s, yp3s, $
                                           xp5s, yp5s, xp7s, yp7s, fp, fps, np, flag, ns, sf, $
                                           xfwhm, yfwhm, /WARM
              
                                ;make a correction for pixel phase 
;              arraycorr = pixel_phase_correct_gauss(f[0],x3,y3,ch+1, '3_3_7')
;              corrflux = f[0]*arraycorr

                                ;apply array dependent correction
;              if chnlnum eq '1' then photcorr = photcor_ch1(x3, y3) else photcorr = photcor_ch2(x3, y3)
;              corrflux= corrflux * photcorr
              
                                ;save them
              xcenarr[c]  = x3 &  ycenarr[c] = y3
              fluxarr[c] = f[0] & fluxerrarr[c] = fs[0] & backarr[c]= b[0] ;& corrfluxarr[c] = corrflux
              c = c + 1
           endif                ; if the target is on the frame
           
        endif                   ; if full array calstar
     endfor                     ; for each fits image
     print, 'final c', c
     xcenarr = xcenarr[0:c-1] 
     ycenarr = ycenarr[0:c-1] 
     fluxarr = fluxarr[0:c-1] 
;     corrfluxarr = corrfluxarr[0:c-1]
     fluxerrarr = fluxerrarr[0:c-1] 
     backarr = backarr[0:c-1] 
     timearr = timearr[0:c-1] 
     starnamearr = starnamearr[0:c-1] 
     raarr = raarr[0:c-1] 
     decarr = decarr[0:c-1] 
                                ;save the variables for plotting seperately
     save,  xcenarr,  ycenarr,  starnamearr,  timearr,  fluxarr,  fluxerrarr,  backarr,  raarr,  decarr , filename =strcompress( '/Users/jkrick/irac_warm/calstars/allch'+string(ch + 1)+'phot.sav',/remove_all) ;corrfluxarr, 
     print, 'time check', systime(1) - t1
     free_lun, outlun
  endfor ; for each channel
 
end
