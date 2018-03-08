pro virgo_galex_test
;measure accurate photometry and background errors for Galex data of
;ICL plume in the virgo cluster

;make sure this works on a galaxy.  Measure ngc 4407 per gil de paz.

!P.multi = [0,2,1]
greencolor = FSC_COLOR("Green", !D.Table_Size-4)

for ch = 1, 1 do begin  ;for fuv and nuv
   if ch eq 0 then begin
      fits_read, '/Users/jkrick/Virgo/Galex/VirgoCdeg-fd-int.fits', data, header
      fits_read, '/Users/jkrick/Virgo/Galex/VirgoCdeg-fd-rrhr.fits', rrdata, rrheader
      fits_read, '/Users/jkrick/Virgo/Galex/fd_seg.fits', segdata, segheader
      zpt = 18.82
      flux_zpt = 1.4E-15
      explow = 3000
      exphigh = 5000
      binsize=0.000005
      AUV = 7.9*0.031
      lambda = 1516 ;angstroms
      print, '****FUV*****'
   endif

   if ch eq 1 then begin
      fits_read, '/Users/jkrick/Virgo/Galex/VirgoCdeg-nd-int.fits', data, header
      fits_read, '/Users/jkrick/Virgo/Galex/VirgoCdeg-nd-rrhr.fits', rrdata, rrheader
      fits_read, '/Users/jkrick/Virgo/Galex/fd_seg.fits', segdata, segheader ; ***********
      zpt = 20.08
      flux_zpt = 2.06E-16
      explow = 29000
      exphigh = 32000
      binsize=0.0001
      AUV = 8.0*0.031
      lambda = 2267 ;angstroms
      print, '****NUV*****'
   endif

;start with quick rough photometry
 ;need to change the aperture and also
 ;get the sky from the histogram, not
 ;like this.
;   xcen = 1874.
;   ycen = 2146.

;   object = (total(data[xcen-37:xcen+37,ycen-75:ycen+75])) 


;roi photometry
   naxis1 = fxpar(header, 'naxis1')
   naxis2 = fxpar(header, 'naxis2')

;   plotimage, xrange=[1,naxis1], yrange = [1,naxis2],  bytscl(data, min =min(data[0:500]), max = max(data[0:500])) ,/preserve_aspect,  ncolors=8

   px = [1829.4444,1874.7778,1923.8889,1908.7778,1837,1818.1111]
   py = [2213.7778,2210,2138.2222,2036.2222,2021.1111,2126.8889]

   px = [1823.8197,1912.5556,1923.8889,1905,1837,1818.1111]
   py = [2226.9581,2247.7778,2138.2222,2024.8889,2006,2070.2222]

   px=[1823.8257,1886.1119,1897.4453,1874.7786,1837.0068,1818.1095]
   py=[2226.9573,2228.8889,2134.4444,2024.8889,2005.9996,2070.2224]

   ;in ra and dec this time, since that is what Mihos needs
  raroi = [186.889,186.86234,186.8575,186.86722,186.88338,186.89146]
  decroi = [13.228459,13.229257,13.189904,13.144259,13.136393,13.163154]

  adxy, header, raroi, decroi, px2, py2

   photroi = Obj_New('IDLanROI', px, py) 
;   draw_roi,photroi             ;,/line_fill,color='444444'x 

;create a mask out of the roi
   mask = photroi->ComputeMask(Dimensions = [naxis1, naxis2], mask_rule = 1)

;apply the mask to make all data outside the roi = 0
   onlyplume = data*(mask GT 0)
   goodarea = n_elements(where(mask gt 0))
   sq2 = fix(sqrt(goodarea) / 2) ; what is the equivalent square with the same area

;   plotimage, xrange=[1,naxis1], yrange = [1,naxis2],  bytscl(onlyplume, min =min(data[0:500]), max = max(data[0:500])) ,/preserve_aspect,  ncolors=8

   ;do photometry with an elliptical aperture.
   sma = 2.3 ;arcmin
   smb = 1.5 ; arcmin
   pa = 60

   ;start with just a circle.
   ra = 186.63469
   dec = 12.61105
   adxy, header, ra, dec, xcen, ycen
   print, 'ra, dec', ra, dec, xcen, ycen
   cntrd, data, xcen, ycen, x, y
   aper, data, x, y, mags, errap, sky, skyerr, 1, [60,92], [150,240],/NAN
   print, 'mags with zpt = 25', mags
   
   ;yea photometry
;   object = (total(onlyplume)) 


;get the sky background below

;look at multiple regions of that same size around the frame, not on
;the plume, their stddev is the important number.

;mask regions with vastly different exptime
   a = where(rrdata lt explow)
   b = where(rrdata gt exphigh)
   nan = alog10(-1)
   data[a] = nan
   data[b] = nan

;mask regions with object detections
   a = where(segdata gt 0)
   data[a] = nan

;mask the plume region
   data[1821:2146,1866:2259] = nan

;randomly choose centers to look at background regions
   nrand = 400
   x = randomu(S,nrand) * 3200. + 100.
   y = randomu(S,nrand) * 3100. + 200.

   bkgd = fltarr(nrand)
   for i = 0, n_elements(x) - 1 do begin
      pixval = data[x(i)-sq2:x(i)+sq2,y(i)-sq2:y(i)+sq2]
      a = where(finite(pixval) gt 0)
      if n_elements(a) gt .85*n_elements(pixval) then  bkgd(i) =  total(data[x(i)-sq2:x(i)+sq2,y(i)-sq2:y(i)+sq2],/nan) / n_elements(a)
   endfor
   a = where(bkgd ne 0)
;print, 'n good background regions', n_elements(a)
   bkgd = bkgd(a)

   plothist, bkgd, xhist, yhist, bin = binsize,/noplot
;print,'mean, stddev before gauss', mean(bkgd), stddev(bkgd)

;fit a gaussian to this distribution
   start = [mean(bkgd),stddev(bkgd), 80.]
   noise = fltarr(n_elements(yhist))
   noise[*] = 1                                            ;equally weight the values
   result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) ;./quiet    

;now plot the histogram of stddevs of all unmasked pixels
   plothist, bkgd, xjunk, yjunk, bin = binsize, xtitle = 'background levels', ytitle='Number', $
             xrange=[-8*result(1) + result(0),result(0) + 8*result(1)],$
             title = "FUV"      ;, thick=3, xthick=3,ythick=3,charthick=3

   xyouts, result(0) + 3.5*result(1), max(yhist) / 2.  + 0.08*max(yhist), 'mean' + string(result(0)) ;, charthick=3
   xyouts, result(0) + 3.5*result(1), max(yhist) / 2. ,  'sigma' + string(result(1))                 ;, charthick=3
;xyouts, result(0) + 3.5*result(1), max(yhist) / 2. + 2*0.08*max(yhist), 'N pixels used' + string(count), charthick=3

;plot the fitted gaussian and print results to plot

   xarr = findgen(10000)/ 1000000.

   oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
          color = greencolor, thick=3
   ;print, 'gauss fit result', result


   ;coming back to the photometry....
   sky = result(0); total(data[1678.-37.:1678.+37.,2251.-75.:2251.+75.])
   flux =  object / (37.*2.*2.*75.) - sky

;   print, 'flux ', flux, ' counts/s over whole aperture', object, sky
;   print, 'flux', flux/ (37.*2.*2.*75.)  , ' counts/s/pix'
   print, 'flux', flux, ' counts/s/pix'


   ;account for galactic extinction
   ;first take couts/s/pix into mag
   m = -2.5*alog10(flux) +zpt
   m = m - AUV    ;apply extinction
   flux = 10^((m - zpt) / (-2.5));back into counts/s/pix

   arcsec_pixel = 1.5
;   SB = zpt -2.5*alog10(flux/11325.) + 5.*alog10(arcsec_pixel)
   SB = zpt -2.5*alog10(flux) + 5.*alog10(arcsec_pixel)

   print, 'SB ', SB, ' AB mag/sqarcsec'
;   print, 'flux', flux_zpt * flux / (37*2*75*2), ' erg/s/cm2/A'
   flambda = flux_zpt * flux 
   print, 'flux', flambda , ' erg/s/cm2/A'
   
   ;the extra 10% is from the error in the absolute photometry
   errorlambda = flux_zpt*result(1)*1.1
   print, 'error', errorlambda, ' erg/s/cm2/A'
   
   ;convert these into Fnu for eric (instead of Flambda)
   c = 299792458 ;m/s
   c = c / 1E-10 ; angstroms / s
   fnu = flambda * (lambda^2) / c
   print, 'flux ', fnu, ' erg/s/cm2/Hz'
   
   errornu = errorlambda*(lambda^2 )/ c
   print, 'error', errornu, ' erg/s/cm2/Hz'

;try fitting a poisson distribution
;function doesn't work with values below 1
;aka .01 factorial doesn't mean anything.
;so fake it and multiply by a large value, now not per second, but per
;day

;xhist = xhist *2*3600.
;yhist = yhist *2*3600.
;print,'xhist', xhist
;save, filename = '/Users/jkrick/idlbin/xhist', xhist
;plot, xhist, yhist
;start = [18.,1E8]
;result= MPFITFUN('poisson',xhist,yhist, noise, start) ;./quiet    
;oplot, findgen(100),result(1)*exp(-result(0))*( result(0)^(findgen(100))) / (factorial(findgen(100))) , $
;                color = greencolor, thick=3, linestyle = 1




;display the resulting masked image
   if ch eq 0 then fits_write, '/Users/jkrick/Virgo/Galex/fd_masked.fits', fuvdata, fuvheader
   if ch eq 1 then fits_write, '/Users/jkrick/Virgo/Galex/nd_masked.fits', fuvdata, fuvheader
   
endfor

end
