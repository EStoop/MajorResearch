pro derive_fluxconv

;for array dependant photometric correction
fits_read, '/Users/jkrick/iwic/photcorr/ch1_photcorr_rj.fits', photcor_cryo_ch1, photcorhead_cryo_ch1
fits_read, '/Users/jkrick/iwic/ch1_photcorr_warm.fits', photcor_warm_ch1, photcorhead_warm_ch1

fits_read, '/Users/jkrick/iwic/photcorr/ch2_photcorr_rj.fits', photcor_cryo_ch2, photcorhead_cryo_ch2
fits_read, '/Users/jkrick/iwic/ch2_photcorr_warm.fits', photcor_warm_ch2, photcorhead_warm_ch2


name = ['1812095_cal', 'HD165459', 'NPM1p60','HD184837','HD156896' ]
ra = [273.04000000 ,270.628089 , 261.2178,294.318477,260.231135 ]
dec = [63.49508333, 58.627264, 60.430817, -24.669625,-20.359752]
real_flux_ch1 = [8.68,647.38,38.20, 40.74,50.40]
real_flux_ch2 = [5.66,421.19,24.74, 26.28,32.54]


for ch = 0, 1 do begin
   print, 'Channel', ch+1
   ratioarr = fltarr(n_elements(ra))
   count = 0
;   dirname = [37993728,37996800,37994240,37993984,37997056,37994496,38052608,38048256,38048512,38048768,38049024,38049280,38049536,38049792,38051840,38052352]

dirname = [37652224 ,37651456 ,37651968 ,37816320 ,37815552 ,37816064 ,37851904 ,37851136 ,37851648 ,37997056 ,37993728 ,37993984 ,37913856 ,37994496 ,37996800 ,37914112 ,37994240 ,37913344 ,37994752 ,37995008 ,37995264 ,37995520 ,38052608 ,38048256 ,38048512 ,38048768 ,38049024 ,38049280 ,38049536 ,38049792 ,38051840 ,38052352 ,38048000 ,38050048 ,38047744 ,38050304 ,38091776 ,38065920 ,38069760 ,38071808 ,38089472 ,38088704 ,38087168 ,38086912 ,38090240 ,38071296 ,38071552 ,38068992 ,38087936 ,38091008 ,38087680 ,38066176 ,38072064 ,38070272 ,38070016 ,38067712 ,38088192 ,38089216 ,38087424 ,38070784 ,38070528 ,38071040 ,38067968 ,38089728 ,38088448 ,38091520 ,38069248 ,38089984 ,38080000 ,38067200 ,38139648,38157824,38157568,38157056,38140416,38139904,38140160,38141184,38140672,38141440]

      if ch eq 0 then dirname = strcompress('/Users/jkrick/irac_warm/calstars/r'+string(dirname) +'/ch1/bcd',/remove_all)
      if ch eq 1 then dirname = strcompress('/Users/jkrick/irac_warm/calstars/r'+string(dirname) +'/ch2/bcd',/remove_all)


      for r = 0,  n_elements(ra) - 1 do begin
         print, 'working on ra', ra(r), dec(r)
         fluxarr = fltarr(10000)
         j = 0
         for d = 0, n_elements(dirname) - 1 do begin
 ;           print, 'working on ' ,dirname(d)
            command = 'ls '+dirname(d) +'/*_cbcd.fits > /Users/jkrick/irac_warm/list.txt'
            command2 =  'ls '+dirname(d) +'/*bcd_fp.fits >> /Users/jkrick/irac_warm/list.txt'
            spawn, command
            spawn, command2
            readcol,'/Users/jkrick/irac_warm/list.txt', fitsname, format = 'A', /silent
            
            for i =1, n_elements(fitsname) - 1 do begin
               fits_read,fitsname(i), bcddata, bcdheader
               adxy, bcdheader, ra(r), dec(r), bcdx, bcdy
               
               if bcdx gt 20 and bcdy gt 20 and bcdy lt 240 and bcdx lt 240 then begin

                                ;do the photometry
                  get_centroids, fitsname(i), t, dt, bcdxcen, bcdycen, bcdflux, xs, ys, fs, b, /WARM, /APER, RA=ra(r), DEC=dec(r),/silent

                             ;make a correction for pixel phase 
                  corrected_flux = correct_pixel_phase(ch+1,bcdxcen,bcdycen,bcdflux)

                      ;apply array dependent correction
                  photcor_ch1 = photcor_warm_ch1(bcdxcen, bcdycen)
                  photcor_ch2 = photcor_warm_ch2(bcdxcen, bcdycen)
;                  print, 'photcor', photcor_ch1, photcor_ch2
                 if ch eq 0 then corrected_flux = corrected_flux * photcor_ch1
                 if ch eq 1 then corrected_flux = corrected_flux * photcor_ch2
                 
                 ; fluxes off by a factor of 100 (mJy)     
                 fluxarr(j) = corrected_flux*1000. ; in mJy

                 j = j + 1
               endif
               
            endfor              ;end for each fits image
         endfor                 ; end for each directory = AOR

         fluxarr = fluxarr[0:j-1]
         
         ;clip any high sigma values out of the mean
        meanclip, fluxarr, meanflux, sigmaflux, clipsig = 3

         if ch eq 0 then begin
            ratioarr[count] =meanflux/ real_flux_ch1(r)
            print,'mean',  meanflux, real_flux_ch1(r),  meanflux/ real_flux_ch1(r)
         endif

         if ch eq 1 then begin
            ratioarr[count] = meanflux / real_flux_ch2(r)
;            print,'mean',  meanflux, real_flux_ch2(r),  meanflux/ real_flux_ch2(r)
         endif

         count = count+1
      endfor                    ; end for each ra, for each star
   
   print, 'ratio', ratioarr, 'mean', mean(ratioarr)
endfor  ;for each channel

end
