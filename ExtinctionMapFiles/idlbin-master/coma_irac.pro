pro coma_irac
;start testing this on ch2 since there it is obvious where the bright star effects the background levels near the plume.

  dirloc = '/Users/jkrick/plumes/coma/spitzer/irac/'
 
 aorname =['r3859456/', 'r3859712/', 'r3859968/'] 


fext = [0.91, 0.94]
for ch = 0, 1 do begin
   meanarr = fltarr(5000)
   sigmaarr = fltarr(5000)
   timearr = fltarr(5000)
   delayarr = fltarr(5000)
   zodyarr = fltarr(5000)
   d = 0

   for a = 0, n_elements(aorname) - 1 do begin
      d = 0
      cd, strcompress(dirloc + aorname(a) + 'ch' + string(ch + 1) + '/bcd',/remove_all)
      spawn, 'pwd'
      if ch eq 0 then command =  "find . -name 'SPITZER_I1*_cbcd.fits' > bcdlist_ch1.txt"
      if ch eq 1 then command =  "find . -name 'SPITZER_I2*_cbcd.fits' > bcdlist_ch2.txt"
      spawn, command
      
      ;read in a list of all cbcds
      readcol,strcompress( 'bcdlist_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A';, /silent
           
      bigim = fltarr(256, 256, n_elements(fitsname))
      c = 0

      for i =1, n_elements(fitsname) - 1 do begin  ;don't use the zero frame
         print, 'working on', fitsname(i)
         
         ;run SExtractor
         cmd = 'sex ' + fitsname(i) + ' -c ../../../default.sex'
         spawn, cmd
         fits_read, fitsname(i), data, header
         mjydata = data
         fits_read, 'segmentation.fits', segdata, segheader

         ;write out image with all objects turned to NaN's
         m = where(segdata gt 0)
         data(m) = alog10(-1)   ;set the objects to NaN's.
         newname = strcompress(strmid(fitsname(i), 0, 37) + 'nan.fits',/remove_all)
         fits_write, newname, data, segheader

        ;keep all images together to make median image for self-dark
         bigim(0,0,c) = data
         c = c + 1
    
        ; do a better background determination and subtraction
         gain  = sxpar(header, 'GAIN')
         flux_conv  = sxpar(header, 'FLUXCONV')
         exptime  = sxpar(header, 'EXPTIME')
         
         bkg = measure_bkg(125,125, data,  gain, exptime, flux_conv)
         
         print, 'bkg', bkg


        ;keep track of mean background levels, and time of
        ;the exposures.  
        timearr(d)  = sxpar(header, 'SCLK_OBS')
        zody  = sxpar(header, 'ZODY_EST')
        zodyarr(d) = zody
        meanarr(d) = mean(mjydata,/nan) - zody ;/ fext(ch)   ;correcting for extended source
        sigmaarr(d) = stddev(data, /nan)
        delayarr(d)  = sxpar(header, 'FRAMEDLY')
        d = d + 1

      endfor   ;for each fits image

      ;generate median image per AOR
      medarr, bigim, meddark

      ;only want to subtract structure in the dark, not the levels.
      ;so subtract the mean from the images so mean = 0
      print, "mean", mean(meddark,/nan)
      meddark = meddark - mean(meddark,/nan)

      fits_write, strcompress( 'med_ch' + string(ch + 1) + '.fits', /remove_all), meddark, imheader

      ;now go back and dark subtract all the cbcd's.
      ;for use in the mosaicing.
      ;also subtract the zody

      for i =1, n_elements(fitsname) - 1 do begin
         fits_read, fitsname(i), data, header
         data = data - meddark

         zody  = sxpar(header, 'ZODY_EST')
         data = data - zody;/fext(ch)

         ;ff effect
         ;NORMALIZE the effect to delay time of 19 = no correction.
         ;so I want to add (19) - (x) to the value of x to get the corrected value
         ;xdelay = sxpar(header, 'FRAMEDLY')
         ;ffcorr = at19 - ( result(0) + result(1)*alog10(xdelay))

         ;if ch eq 0 then data = data + ffcorr
         
         ;now add some overall value so that the levels are above zero.
         data = data + 0.5


         fits_write, strcompress(strmid(fitsname(i), 0, 37) + 'ndark.fits',/remove_all), data, header
      endfor   ;for each fits image


      if ch eq 0 and a eq 0 then q = plot( findgen(n_elements(meanarr)), meanarr/10.,  xtitle = 'Frame Number', ytitle = 'mean levels', xrange = [0,220], yrange = [-0.1, 0.4],/nodata)
      if ch eq 0 and a gt 0 then q = plot( findgen(n_elements(meanarr)), meanarr, /overplot)
  ;    if ch eq 1 and a eq 0 then q1 = plot( findgen(n_elements(meanarr)), meanarr/10., color = 'red', /overplot)
      if ch eq 1 and a gt 0 then q1 = plot( findgen(n_elements(meanarr)), meanarr+0.3, color = 'red', /overplot)

      if ch eq 0 and a gt 0 then q2 = plot(findgen(n_elements(zodyarr)), zodyarr , /overplot)
      if ch eq 1 and a gt 0 then q2 = plot(findgen(n_elements(zodyarr)), zodyarr , color = 'red', /overplot)

   endfor ;for each AOR

   ;save these to work with them later
   ;instead of re-running SExtractor

   delayarr = delayarr[0:d-1]
   meanarr = meanarr[0:d-1]

;   if ch eq 0 then p = plot( delayarr, meanarr, '2*', xtitle = 'delay time (s)', ytitle = 'mean levels(Mjy/sr)')
;   if ch eq 1 then p1 = plot( delayarr, meanarr, '2sr', /overplot)



;   save, /variables, filename =strcompress( '/Users/jkrick/plumes/coma/spitzer/irac/ch1/ndark.sav',/remove_all)

endfor  ;for each channel




end

function measure_bkg, xcen, ycen, data, gain, exptime, flux_conv
 ; print, gain, exptime, flux_conv
;In order to find an accurate background I first randomly choose
;regions in the frame.  I then make a histogram of the values in those
;regions.  I fit the histogram with a gaussian and take the mean of
;the gaussian as the mean background value.

  ;need to convert the data into electrons
  ; print, 'examples', data[100,100], data[150,150]
 data = data * gain*exptime/flux_conv
 ; print, 'examples', data[100,100], data[150,150]
  nrand = 100
  ;choose random locations that are not too close to the edges of the frame
  xbkg_radius = xcen - 1 < 36 < (254 - xcen)
  ybkg_radius = ycen - 1 < 36 < (254 - ycen)
  x = randomu(S,nrand) *(xbkg_radius * 2) + xcen - xbkg_radius
  y = randomu(S,nrand) *(ybkg_radius * 2) + ycen - ybkg_radius
 
  bkgd = fltarr(nrand)
  for ti = 0, n_elements(x) - 1 do begin  ; for each random background region
     pixval = data[x(ti)-2:x(ti)+2,y(ti)-2:y(ti)+2]
     ga = where(finite(pixval) gt 0)
     if n_elements(ga) gt .85*n_elements(pixval) then  begin
        bkgd(ti) =  total(data[x(ti)-2:x(ti)+2,y(ti)-2:y(ti)+2],/nan) ; / n_elements(ga)
 ;       print, 'bkgd(ti)', bkgd(ti)
     endif

  endfor
  ta = where(bkgd ne 0)
  bkgd = bkgd(ta)
  ;print, 'n_ta', n_elements(ta)
  if n_elements(ta) lt 10 then begin  ; if there aren't enough regions to fit a gaussian.
     print, 'arent enough regions to fit a gaussian'
     result = [mean(bkgd), stddev(bkgd), 2.0]
  endif else begin

;do some sigma clipping on bkgd before running it through plothist

     skpix = MEDCLIP(bkgd, med, sigma, clipsig = 3, converge_num = 0.03)
     bkgd = skpix
;fit a gaussian to this distribution
     binsize = 200           
     plothist, bkgd, xhist, yhist, bin = binsize, xrange = [-10000,10000] ,/noplot ;make a histogram of the data
;     print, 'plothist xhist', xhist
;     print, 'plothist yhist', yhist

     start = [mean(bkgd),stddev(bkgd), 2.]
    ; print, 'start', start
     noise = fltarr(n_elements(yhist))
     noise[*] = 1  
                                ;equally weight the values
     result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) ;  fit a gaussian to the histogram sorted data
     bkg = result(0)
;     oplot, xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.)
  endelse

     return, result
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
