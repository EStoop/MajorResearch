pro calstar_vs_time

fits_read, '/Users/jkrick/iwic/ch1_photcorr_gauss.fits', photcor_warm_ch1, photcorhead_warm_ch1
fits_read, '/Users/jkrick/iwic/ch2_photcorr_gauss.fits', photcor_warm_ch2, photcorhead_warm_ch2
;fits_read, '/Volumes/irac_drive/ch1_photcorr_warm.fits', photcor_warm_ch1, photcorhead_warm_ch1
;fits_read, '/Volumes/irac_drive/ch2_photcorr_warm.fits', photcor_warm_ch2, photcorhead_warm_ch2

name = ['1812095_cal', 'HD165459', 'NPM1p60','HD184837','HD156896' ]
ra = [273.04000000 ,270.628089 , 261.2178,294.318477,260.231135 ]
dec = [63.49508333, 58.627264, 60.430817, -24.669625,-20.359752]
real_flux_ch1 = [8.68,647.38,38.20, 40.74,50.40]
real_flux_ch2 = [5.66,421.19,24.74, 26.28,32.54]


for ch = 0, 0 do begin
   print, 'Channel', ch+1
   
;pc1: 37652224 ,37651456 ,37651968
;pc2: 37816064,37815552,37816320
;pc3: 37851904 ,37851136 ,37851648
;pc4: 37914112,37913856 ,37913344
   
   dirname = [37652224 ,37651456 ,37651968,37816064,37815552,37816320,37851904 ,37851136 ,37851648, 37914112,37913856 ,37913344,  37997056 ,37993728 ,37993984 ,37994496 ,37996800  ,37994240  ,37994752 ,37995008 ,37995264 ,37995520 ,38052608 ,38048256 ,38048512 ,38048768 ,38049024 ,38049280 ,38049536 ,38049792 ,38051840 ,38052352 ,38048000 ,38050048 ,38047744 ,38050304 ,38091776 ,38065920 ,38069760 ,38071808 ,38089472 ,38088704 ,38087168 ,38086912 ,38090240 ,38071296 ,38071552 ,38068992 ,38087936 ,38091008 ,38087680 ,38066176 ,38072064 ,38070272 ,38070016 ,38067712 ,38088192 ,38089216 ,38087424 ,38070784 ,38070528 ,38071040 ,38067968 ,38089728 ,38088448 ,38091520 ,38069248 ,38089984 ,38080000 ,38067200 ,38139648,38157824,38157568,38157056,38140416,38139904,38140160,38141184,38140672,38141440,38402304,38171136,38171904,38403072,38171648,38402816,38419200,38418432,38418944,38686976,38686720,38686208,38705920,38705664,38705152]

   if ch eq 0 then dirname = strcompress('/Users/jkrick/irac_warm/calstars/r'+string(dirname) +'/ch1/bcd',/remove_all)
   if ch eq 1 then dirname = strcompress('/Users/jkrick/irac_warm/calstars/r'+string(dirname) +'/ch2/bcd',/remove_all)
;   if ch eq 0 then dirname = strcompress('/Volumes/irac_drive/irac_warm/calstars/r'+string(dirname) +'/ch1/bcd',/remove_all)
;   if ch eq 1 then dirname = strcompress('/Volumes/irac_drive/irac_warm/calstars/r'+string(dirname) +'/ch2/bcd',/remove_all)
   
;set up for keeping track of flux, time, and position for each calstar for each aor.
   numoobjects = 3 ;only three primary calstars for now
    calobject = replicate({calob, flux:fltarr(n_elements(dirname)),flux_pixelcorr:fltarr(n_elements(dirname)), flux_arraycorr:fltarr(n_elements(dirname)), time:fltarr(n_elements(dirname)), pos:fltarr(n_elements(dirname)), flux_e:fltarr(n_elements(dirname)), percent_e:fltarr(n_elements(dirname))},numoobjects)
 ;   calobject = replicate({calob, flux:fltarr(n_elements(dirname)*10),flux_pixelcorr:fltarr(n_elements(dirname)*10), flux_arraycorr:fltarr(n_elements(dirname)*10), time:fltarr(n_elements(dirname)*10), pos:fltarr(n_elements(dirname)*10)},numoobjects)

   for r = 0,  n_elements(ra) - 3 do begin ;only use primary calstars    
      print, 'working on star', r
      j = 0
      
      for d = 0, n_elements(dirname) - 1 do begin
         command = 'ls '+dirname(d) +'/*_cbcd.fits > /Users/jkrick/list.txt'
         command2 =  'ls '+dirname(d) +'/*bcd_fp.fits >> /Users/jkrick/list.txt'
         spawn, command
         spawn, command2
         readcol,'/Users/jkrick/list.txt', fitsname, format = 'A', /silent
         
         
         
         for i =0, n_elements(fitsname) - 1 do begin
            fits_read,fitsname(i), bcddata, bcdheader
            clock = sxpar(bcdheader, 'SCLK_OBS')
            
            adxy, bcdheader, ra(r), dec(r), bcdx, bcdy
            if bcdx gt 20 and bcdy gt 20 and bcdy lt 240 and bcdx lt 240 then begin
               
                                ;do the photometry
               get_centroids, fitsname(i), t, dt, bcdxcen, bcdycen, bcdflux, xs, ys, fs, b, $
                              /WARM, /APER, RA=ra(r), DEC=dec(r),/silent
               calobject[r].time[j] = clock
               calobject[r].flux[j] = bcdflux * 1000.
               
               gain = 3.7
               exptime =  sxpar(bcdheader, 'EXPTIME') 
               fluxconv =  sxpar(bcdheader, 'FLUXCONV') 
               fluxconv =  .1198;.1443;
               sbtoe = gain*exptime/fluxconv

               ;flux in electrons
               ;convert to microJy
               eflux = bcdflux *1E6
               ;convert to Mjy/sr
               eflux = eflux / 34.98
               ;convert to e
               eflux = eflux *sbtoe

               ;good guess at photon statistics error
               phot_err = sqrt(eflux)
               calobject[r].flux_e[j] = eflux
               calobject[r].percent_e[j] = sqrt(eflux) / eflux

                                ;make a correction for pixel phase 
               ;corrected_flux = correct_pixel_phase(ch+1,bcdxcen,bcdycen,bcdflux)
               corrected_flux = pixel_phase_correct_gauss(bcdflux,bcdxcen,bcdycen,ch+1)

               calobject[r].flux_pixelcorr[j] = corrected_flux * 1000.

                                ;apply array dependent correction
               photcor_ch1 = photcor_warm_ch1(bcdxcen, bcdycen)
               photcor_ch2 = photcor_warm_ch2(bcdxcen, bcdycen)
               if ch eq 0 then corrected_flux = corrected_flux * photcor_ch1
               if ch eq 1 then corrected_flux = corrected_flux * photcor_ch2
               calobject[r].flux_arraycorr[j] = corrected_flux * 1000.
                              
               ;keep track of the position
               if bcdxcen lt 130 + 10 and bcdxcen gt 130 - 10 and bcdycen lt 130 + 10 and bcdycen gt 130 - 10 then begin
                  calobject[r].pos[j] = 1
                  print, 'pos1'
               endif

               if bcdxcen lt 65 + 10 and bcdxcen gt 65 - 10 and bcdycen lt 190 + 10 and bcdycen gt 190 - 10 then calobject[r].pos[j] = 2
               if bcdxcen lt 190 + 10 and bcdxcen gt 190 - 10 and bcdycen lt 190 + 10 and bcdycen gt 190 - 10 then calobject[r].pos[j] = 3
               if bcdxcen lt 190 + 10 and bcdxcen gt 190 - 10 and bcdycen lt 65 + 10 and bcdycen gt 65 - 10 then calobject[r].pos[j] = 4
               if bcdxcen lt 65 + 10 and bcdxcen gt 65 - 10 and bcdycen lt 65 + 10 and bcdycen gt 65 - 10 then calobject[r].pos[j] = 5
               j = j + 1
            endif   
         endfor                 ;end for each fits image
      endfor                    ; end for each directory = AOR  
      calobject[r].time = calobject[r].time[0:j-1]
      calobject[r].flux = calobject[r].flux[0:j-1]
      calobject[r].flux_e = calobject[r].flux_e[0:j-1]
      calobject[r].percent_e = calobject[r].percent_e[0:j-1]
      calobject[r].flux_pixelcorr = calobject[r].flux_pixelcorr[0:j-1]
      calobject[r].flux_arraycorr = calobject[r].flux_arraycorr[0:j-1]


   endfor                       ; end for each ra, for each star 
endfor                          ;for each channel

save, calobject, filename='/Users/jkrick/idlbin/calobject.sav'
meanflux = mean(calobject[0].flux[0:j-1])
plot, calobject[0].time[0:j-1], calobject[0].flux[0:j-1], psym = 2, xrange = [calobject[0].time[0], calobject[0].time[j-1]], yrange=[meanflux - 0.3*meanflux, meanflux + 0.3*meanflux]

undefine, calobject
end
