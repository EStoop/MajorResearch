pro calstar_rawfluxes


  dirloc = ['/Volumes/iracdata/flight/IWIC/IRAC038200/']
  
  cd, dirloc
  command ="find ./bcd/00*/ -name 'IRAC.1*bcd_fp.fits' >> /Users/jkrick/irac_warm/calstars/pc134ch1bcdlist.txt "
;  spawn, command
  command ="find ./bcd//00*/ -name 'IRAC.2*bcd_fp.fits' >> /Users/jkrick/irac_warm/calstars/pc134ch2bcdlist.txt "
;  spawn, command
  command ="find ./raw//00*/ -name 'IRAC.1*mipl.fits' >> /Users/jkrick/irac_warm/calstars/pc134ch1rawlist.txt "
;  spawn, command
  command ="find ./raw//00*/ -name 'IRAC.2*mipl.fits' >> /Users/jkrick/irac_warm/calstars/pc134ch2rawlist.txt "
;  spawn, command

  for ch = 0, 1 do begin
     print, '----------------------------------------'
     print, 'working on ch ', ch + 1
     
     ;fitsname = prepare_fits('/Users/jkrick/irac_warm/calstars/allch1bcdlist.txt')
     if ch eq 0 then begin
        readcol, '/Users/jkrick/irac_warm/calstars/pc134ch1bcdlist.txt', fitsname, format = 'A', /silent
        readcol, '/Users/jkrick/irac_warm/calstars/pc134ch1rawlist.txt', rawname, format = 'A', /silent
     endif else begin
        readcol, '/Users/jkrick/irac_warm/calstars/pc134ch2bcdlist.txt', fitsname, format = 'A', /silent
        readcol, '/Users/jkrick/irac_warm/calstars/pc134ch2rawlist.txt', rawname, format = 'A', /silent
     endelse
     
     print, 'nfits, nraaw', n_elements(fitsname) , n_elements(rawname)
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
     peakpixarr = xcenarr

     startfits = 0L
     stopfits = n_elements(fitsname) - 1
     c = 0L
     for i= startfits, stopfits do begin
        header = headfits(fitsname(i)) ;
        NAXIS= sxpar(header, 'NAXIS')
        
        AORLABEL= sxpar(header, 'AORLABEL')
;       print, i, strmid(AORLABEL, 0, 12), naxis
                                ;cut on full array calstars only (for now)
        if strmid(AORLABEL, 0, 12) eq 'IRAC_calstar' and NAXIS lt 3 then begin            ;got a good one
           
           fits_read,fitsname(i), im, h
           inter = strmid(fitsname(i), 0, 48)
           uncname = strcompress(inter + 'sig_dntoflux.fits',/remove_all)
;           print, 'uncname',i,  uncname, fitsname(i)
           fits_read, uncname, unc, hunc, /no_abort  ; so it won't crash if the file isn't there but should use the last unc file.

           chnlnum = sxpar(h, 'CHNLNUM')
           ra_ref = sxpar(h, 'RA_REF')
           dec_ref = sxpar(h, 'DEC_REF')
           sos_ver = sxpar(h, 'SOS_VER')
           exptime = sxpar(h, 'EXPTIME')

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

;read in the raw data file and get the DN level of the peakpixel      
              fits_read, rawname(i), rawdata, rawheader
;              rawdata = readfits(rawname(i), rawheader)
              barrel = fxpar(rawheader, 'A0741D00')
              fowlnum = fxpar(rawheader, 'A0614D00')
              pedsig = fxpar(rawheader, 'A0742D00')
              ichan = fxpar (rawheader, 'CHNLNUM')
              
                                ;use Bill's code to conver to DN
;              dewrap2, rawdata, ichan, barrel, fowlnum, pedsig, 0, rawdata
     
                                ;or use Jim's code to convert to DN
              rawdata = irac_raw2dn(rawdata,ichan,barrel,fowlnum)
;              rawdata = reform(rawdata, 32, 32, 64)
              peakpixarr[c] = max(rawdata[x3-1:x3+1, y3-1:y3+1])
              print, AORLABEL, peakpixarr[c], exptime
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
     save,  xcenarr,  ycenarr,  starnamearr,  timearr,  fluxarr,  fluxerrarr,  backarr,  raarr,  decarr , peakpixarr,  filename =strcompress( '/Users/jkrick/irac_warm/calstars/pc134ch'+string(ch + 1)+'phot.sav',/remove_all) ;corrfluxarr, 
     print, 'time check', systime(1) - t1
  endfor ; for each channel
 
end
