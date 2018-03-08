pro derive_fluxconv

;for array dependant photometric correction
fits_read, '/Users/jkrick/iwic/photcorr_cryo/ch1_photcorr_rj.fits', photcor_cryo_ch1, photcorhead_cryo_ch1
fits_read, '/Users/jkrick/iwic/ch1_photcorr_gauss.fits', photcor_warm_ch1, photcorhead_warm_ch1

fits_read, '/Users/jkrick/iwic/photcorr_cryo/ch2_photcorr_rj.fits', photcor_cryo_ch2, photcorhead_cryo_ch2
fits_read, '/Users/jkrick/iwic/ch2_photcorr_gauss.fits', photcor_warm_ch2, photcorhead_warm_ch2


name = ['1812095_cal', 'HD165459', 'NPM1p60','HD184837','HD156896' ]
ra = [273.04000000 ,270.628089 , 261.2178];,294.318477,260.231135 ]
dec = [63.49508333, 58.627264, 60.430817];, -24.669625,-20.359752]
real_flux_ch1 = [8.68,647.38,38.20, 40.74,50.40]
real_flux_ch2 = [5.66,421.19,24.74, 26.28,32.54]


for ch = 0, 1 do begin
   print, 'Channel', ch+1
   ratioarr = fltarr(n_elements(ra))
   count = 0


   for r = 0, n_elements(ra) - 1 do begin
      print, 'working on ra', ra(r), dec(r)
      fluxarr = fltarr(10000)
      j = 0
      
;pc5 through pc17
      dirloc = ['/Users/jkrick/iracdata/flight/IWIC/calstars','/Users/jkrick/irac_warm/calstars/r38419200','/Users/jkrick/irac_warm/calstars/r38418432','/Users/jkrick/irac_warm/calstars/r38686976','/Users/jkrick/irac_warm/calstars/r38418944','/Users/jkrick/irac_warm/calstars/r38686208','/Users/jkrick/irac_warm/calstars/r38686720','/Users/jkrick/irac_warm/calstars/r38705920','/Users/jkrick/irac_warm/calstars/r38705664','/Users/jkrick/irac_warm/calstars/r38705152', '/Users/jkrick/irac_warm/calstars/r38818048','/Users/jkrick/irac_warm/calstars/r38817280','/Users/jkrick/irac_warm/calstars/r38817792','/Users/jkrick/irac_warm/calstars/r38925056','/Users/jkrick/irac_warm/calstars/r38924288','/Users/jkrick/irac_warm/calstars/r38924800', '/Users/jkrick/irac_warm/calstars/r39099392','/Users/jkrick/irac_warm/calstars/r39099136','/Users/jkrick/irac_warm/calstars/r39098624']


      for c = 0, n_elements(dirloc) - 1 do begin ;for each directory with calstars in it
         cd, dirloc(c)
         if ch eq 0 then begin
            if c eq 0 then command  =  "find . -name 'IRAC.1*bcd_fp.fits' > /Users/jkrick/irac_warm/list.txt"
            if c gt 0 then command =  "find . -name 'SPITZER_I1*_bcd.fits' > /Users/jkrick/irac_warm/list.txt"
            spawn, command
         endif
         
         if ch eq 1 then begin
            if c eq 0 then command  =  "find . -name 'IRAC.2*bcd_fp.fits' > /Users/jkrick/irac_warm/list.txt"
            if c gt 0 then command =  "find . -name 'SPITZER_I2*_bcd.fits' > /Users/jkrick/irac_warm/list.txt"
            spawn, command
         endif
         
         readcol,'/Users/jkrick/irac_warm/list.txt', fitsname, format = 'A', /silent
         
         for i =1, n_elements(fitsname) - 1 do begin
;            print, 'working on fitsna', fitsname(i)
            fits_read,fitsname(i), bcddata, bcdheader
            adxy, bcdheader, ra(r), dec(r), bcdx, bcdy
               
            if bcdx gt 20 and bcdy gt 20 and bcdy lt 240 and bcdx lt 240 then begin

                                ;do the photometry
               get_centroids, fitsname(i), t, dt, bcdxcen, bcdycen, bcdflux, xs, ys, fs, b, /WARM, /APER, RA=ra(r), DEC=dec(r),/silent

                             ;make a correction for pixel phase 
               ;corrected_flux = correct_pixel_phase(ch+1,bcdxcen,bcdycen,bcdflux)
               corrected_flux = pixel_phase_correct_gauss(bcdflux,bcdxcen,bcdycen,ch+1)

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
            
         endfor                 ;end for each fits image
      endfor    ;end for each directory

         
         fluxarr = fluxarr[0:j-1]
         
         ;clip any high sigma values out of the mean
         meanclip, fluxarr, meanflux, sigmaflux, clipsig = 3

         if ch eq 0 then begin
            ratioarr[count] =meanflux/ real_flux_ch1(r)
            print, 'mean',  meanflux, real_flux_ch1(r),  meanflux/ real_flux_ch1(r)
         endif

         if ch eq 1 then begin
            ratioarr[count] = meanflux / real_flux_ch2(r)
            print,'mean',  meanflux, real_flux_ch2(r),  meanflux/ real_flux_ch2(r)
         endif

         count = count+1
      endfor                    ; end for each ra, for each star
   
   print, 'ratio', ratioarr, 'mean', mean(ratioarr)
endfor  ;for each channel

end
