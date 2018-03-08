pro derive_fluxconv
!p.multi = [0, 2, 3]

;what size aperture are we using for photometry?
aprad = 5
apcorr = [1.049,1.050]
;for array dependant photometric correction
filename1 = strcompress('/Users/jkrick/iwic/ch1_photcorr_ap_' + string(aprad) + '.fits',/remove_all)
fits_read,filename1, photcor_warm_ch1, photcorhead_warm_ch1

filename2 = strcompress('/Users/jkrick/iwic/ch2_photcorr_ap_' + string(aprad) + '.fits',/remove_all)
fits_read, filename2, photcor_warm_ch2, photcorhead_warm_ch2


name = ['1812095_cal', 'HD165459', 'NPM1p60','HD184837','HD156896' ]
ra = [273.04000000 ,270.628089 , 261.2178,294.318477,260.231135 ]
dec = [63.49508333, 58.627264, 60.430817, -24.669625,-20.359752]
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
      
 dirloc = ['/Users/jkrick/iracdata/flight/IWIC/S18.16/calstars']

      for c = 0, n_elements(dirloc) - 1 do begin ;for each directory with calstars in it
         cd, dirloc(c)
         if ch eq 0 then begin
            command1  =  "find . -name 'IRAC.1*bcd_fp.fits' > /Users/jkrick/irac_warm/list.txt"
            spawn, command1
         endif
         
         if ch eq 1 then begin
            command1  =  "find . -name 'IRAC.2*bcd_fp.fits' > /Users/jkrick/irac_warm/list.txt"
            spawn, command1
         endif
         
         readcol,'/Users/jkrick/irac_warm/list.txt', fitsname, format = 'A', /silent
         
         for i =1, n_elements(fitsname) - 1 do begin
;            print, 'working on fitsna', fitsname(i)
            fits_read,fitsname(i), bcddata, bcdheader
            adxy, bcdheader, ra(r), dec(r), bcdx, bcdy
               
            if bcdx gt 20 and bcdy gt 20 and bcdy lt 240 and bcdx lt 240 then begin

                                ;do the photometry
               get_centroids, fitsname(i), t, dt, bcdxcen, bcdycen, bcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = aprad, RA=ra(r), DEC=dec(r),/silent

                             ;make a correction for pixel phase 
                                ;corrected_flux =
                                ;correct_pixel_phase(ch+1,bcdxcen,bcdycen,bcdflux)
               pixap = strcompress(string(aprad) + '_12_20',/remove_all)
               corrected_flux = pixel_phase_correct_gauss(bcdflux,bcdxcen,bcdycen,ch+1,pixap)

                      ;apply array dependent correction
               photcor_ch1 = photcor_warm_ch1(bcdxcen, bcdycen)
               photcor_ch2 = photcor_warm_ch2(bcdxcen, bcdycen)
;                  print, 'photcor', photcor_ch1, photcor_ch2
               if ch eq 0 then corrected_flux = corrected_flux * photcor_ch1
               if ch eq 1 then corrected_flux = corrected_flux * photcor_ch2
               
               ;aperture correction for small aperture
               if aprad eq 5 then corrected_flux = corrected_flux * apcorr(ch)

                 ; fluxes off by a factor of 100 (mJy)     
               fluxarr(j) = corrected_flux*1000. ; in mJy
               
              
               j = j + 1
            endif
            
         endfor                 ;end for each fits image
      endfor    ;end for each directory

         
      fluxarr = fluxarr[0:j-1]
      print, 'fluxarr',fluxarr
      print, 'stddev', stddev(fluxarr)
      ;plothist, fluxarr, bin = 0.1
         ;clip any high sigma values out of the mean
      meanclip, fluxarr, meanflux, sigmaflux, clipsig = 3
      MMM, fluxarr,skymde, sigma, skew
           
;fit a gaussian
      plothist, fluxarr, xhist, yhist, bin = 0.05,/noplot
      start = [median(fluxarr),0.1,15.]
      noise = fltarr(n_elements(xhist))
      noise[*] = 1.0
      gresult= MPFITFUN('mygauss',xhist,yhist, noise, start)  

      if ch eq 0 then begin
         ratioarr[count] =gresult(0)/ real_flux_ch1(r) 
         print, 'mean',  gresult(0), real_flux_ch1(r),  gresult(0)/ real_flux_ch1(r)
      endif

      if ch eq 1 then begin
         ratioarr[count] = gresult(0) / real_flux_ch2(r) 
         print,'mean',  gresult(0), real_flux_ch2(r),  gresult(0)/ real_flux_ch2(r)
      endif


      count = count+1
   endfor                       ; end for each ra, for each star
   
   print, 'ratio', ratioarr, 'mean', mean(ratioarr)
endfor  ;for each channel

end
