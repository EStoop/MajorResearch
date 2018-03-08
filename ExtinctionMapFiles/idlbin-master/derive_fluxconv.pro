pro derive_fluxconv

real_flux_ch1 = [8.68,647.38,38.20];, 9.364];,13.93,10.53]
real_flux_ch2 = [5.66,421.19,24.74];, 9.355];,7.80,5.99]
;for array dependant photometric correction
fits_read, '/Users/jkrick/iwic/photcorr/ch1_photcorr_rj.fits', photcor_cryo_ch1, photcorhead_cryo_ch1
fits_read, '/Users/jkrick/iwic/ch1_photcorr_warm.fits', photcor_warm_ch1, photcorhead_warm_ch1

fits_read, '/Users/jkrick/iwic/photcorr/ch2_photcorr_rj.fits', photcor_cryo_ch2, photcorhead_cryo_ch2
fits_read, '/Users/jkrick/iwic/ch2_photcorr_warm.fits', photcor_warm_ch2, photcorhead_warm_ch2

;for warm pixel phase correction
;restore, '/Users/jkrick/iwic/pixel_phase_img_ch1.sav'   
;restore, '/Users/jkrick/iwic/pixel_phase_img_ch2.sav'   


;pc5 + pc6
;dirname_1 = [37993728,37996800,37994240];,37994752,37995264]
;dirname_1 = [37993728,37996800,37994240,38052608,38048256,38048512,38048768,38049024,38049280,38049536,38049792,38051840,38052352]
;dirname_2 = [37993984,37997056,37994496];,37995008,37995520]
;dirname_2 = [37993984,37997056,37994496,38052608,38048256,38048512,38048768,38049024,38049280,38049536,38049792,38051840,38052352]

name = ['1812095_cal', 'HD165459', 'NPM1p60']
ra = [273.04000000 ,270.628089 , 261.2178];, 260.231135];, 269.49375000,269.65791667]
dec = [63.49508333, 58.627264, 60.430817];,-20.359752];,66.87472222,66.78111111]


for ch = 0, 1 do begin
   print, 'Channel', ch+1
   ratioarr = fltarr(n_elements(ra))
   count = 0
;   for run = 0, 1 do begin
;      if run eq 0 then dirname = dirname_1
;      if run eq 1 then dirname = dirname_2
   dirname = [37993728,37996800,37994240,37993984,37997056,37994496,38052608,38048256,38048512,38048768,38049024,38049280,38049536,38049792,38051840,38052352]

      if ch eq 0 then dirname = strcompress('/Users/jkrick/irac_warm/calstars/r'+string(dirname) +'/ch1/bcd',/remove_all)
      if ch eq 1 then dirname = strcompress('/Users/jkrick/irac_warm/calstars/r'+string(dirname) +'/ch2/bcd',/remove_all)


      for r = 0,  n_elements(ra) - 1 do begin
         print, 'working on ra', ra(r), dec(r)
         fluxarr = fltarr(100)
         j = 0
         for d = 0, n_elements(dirname) - 1 do begin
 ;           print, 'working on ' ,dirname(d)
            command = 'ls '+dirname(d) +'/*_cbcd.fits > /Users/jkrick/irac_warm/ch1list.txt'
            spawn, command
            readcol,'/Users/jkrick/irac_warm/ch1list.txt', fitsname, format = 'A', /silent
            
            for i =1, n_elements(fitsname) - 1 do begin
               fits_read,fitsname(i), bcddata, bcdheader
               adxy, bcdheader, ra(r), dec(r), bcdx, bcdy
               
               if bcdx gt 20 and bcdy gt 20 and bcdy lt 240 and bcdx lt 240 then begin
                  
;                  cntrd, bcddata, bcdx, bcdy, bcdxcen, bcdycen, 5
  
                  ;conver the data from MJy/sr to electrons
;;                  gain = 3.7
;;                  exptime =  sxpar(bcdheader, 'EXPTIME') 
;;                  fluxconv =  sxpar(bcdheader, 'FLUXCONV') 
;;                  sbtoe = gain*exptime/fluxconv
;;                 bcddata = bcddata*sbtoe

                                ;do the photometry
                  get_centroids, fitsname(i), t, dt, bcdxcen, bcdycen, bcdflux, xs, ys, fs, b, /WARM, /APER, RA=ra(r), DEC=dec(r),/silent
                  print, 'x',bcdxcen, 'y',bcdycen
                  ;aper, bcddata, bcdxcen, bcdycen, bcdflux, errap, sky, skyerr, 1, [3], [3,7],/NAN,/flux,/silent
                  ;;aper, bcddata, bcdxcen, bcdycen, bcdflux, errap, sky, skyerr, 1, [10], [12,20],/NAN,/flux,/silent,/exact

                  ;convert back to Mjy/sr
                  ;;bcdflux = bcdflux /sbtoe

                                ;aperture correction
                  ;if ch eq 0 then bcdflux = bcdflux *1.124 ;for ch1
                  ;if ch eq 1 then bcdflux = bcdflux *1.127 ;for ch2

                                ;color correction?
               ;if r le 2 then bcdflux = bcdflux * 1.019                  ;1.001 for ch2   A stars
               ;if r gt 2 then bcdflux = bcdflux * 1.021                  ;1.007 for ch2  K stars
               

                                ;make a correction for pixel phase 
                  corrected_flux = correct_pixel_phase(ch+1,bcdxcen,bcdycen,bcdflux)
;                  xphase = (bcdxcen mod 1) - NINT(bcdxcen mod 1)
;                  yphase = (bcdycen mod 1) - NINT(bcdycen mod 1)
;                  p = sqrt(xphase^2 + yphase^2)
;                  if ch eq 0 then restore, '/Users/jkrick/iwic/pixel_phase_img_ch1_old.sav'   
;                  if ch eq 1 then restore, '/Users/jkrick/iwic/pixel_phase_img_ch2_old.sav'                     
;                  interpolated_relative_flux = interp2d(relative_flux,x_center,y_center,xphase,yphase,/regular)
;                  corrected_flux =bcdflux / interpolated_relative_flux

                      ;apply array dependent correction
                  photcor_ch1 = photcor_warm_ch1(bcdxcen, bcdycen)
                  photcor_ch2 = photcor_warm_ch2(bcdxcen, bcdycen)
;                  print, 'photcor', photcor_ch1, photcor_ch2
                 if ch eq 0 then corrected_flux = corrected_flux * photcor_ch1
                 if ch eq 1 then corrected_flux = corrected_flux * photcor_ch2

                 
                                ;put into microjy
;;                  corrected_flux = corrected_flux * 34.9837
;                  print, 'corrected_flux', corrected_flux / 1000.
;;                  fluxarr(j) = corrected_flux / 1000. ; in mJy
                  j = j + 1
               endif
               
            endfor              ;end for each fits image
         endfor                 ; end for each directory = AOR

         fluxarr = fluxarr[0:j-1]
         
;try throwing out the high value in fluxarr ..... this sends things in
;the wrong direction, aka pulls down the ratio so the discrepancy is
;higher.
;but this is roughly what Patrick does.
         fluxsort = fluxarr[sort(fluxarr)]
         fluxsort = fluxsort[0:n_elements(fluxsort) - 2]
;         print, 'fluxarr', fluxarr
;         print, 'fluxsort', fluxsort

         meanclip, fluxarr, mean, sigma, clipsig = 2

         if ch eq 0 then ratioarr[count] =mean/ real_flux_ch1(r)
         if ch eq 1 then ratioarr[count] = mean / real_flux_ch2(r)
         count = count+1
      endfor                    ; end for each ra, for each star
;   endfor     ;for eachh run
   
   print, 'ratio', ratioarr, 'mean', mean(ratioarr)
endfor  ;for each channel

end
