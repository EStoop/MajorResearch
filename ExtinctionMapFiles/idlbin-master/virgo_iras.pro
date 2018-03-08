pro virgo_iras
;measure accurate photometry and background errors for Galex data of
;ICL plume in the virgo cluster
ps_open, filename='/Users/jkrick/Virgo/IRAS/background_high.ps',/portrait,/square,/color

!P.multi = [0,1,2]
!P.thick = 1
!P.charthick = 1
greencolor = FSC_COLOR("Green", !D.Table_Size-4)

for ch = 0, 1 do begin ; for 60 and 100
   if ch eq 0 then begin
      fits_read, '/Users/jkrick/Virgo/IRAS/IRIS_band4/0001_186.44478000_12.66209000_I253B4H0.fits', data, header
      binsize=0.01
      lambda = 100 ;microns
      print, '****100*****'
   endif

   if ch eq 1 then begin
      fits_read, '/Users/jkrick/Virgo/IRAS/IRIS_band3/0001_186.44478000_12.66209000_I253B3H0.fits', data, header
      binsize=0.01
      lambda = 60 ;microns
      print, '****60*****'
   endif


;roi photometry
   naxis1 = fxpar(header, 'naxis1')
   naxis2 = fxpar(header, 'naxis2')
    
;   plotimage, xrange=[1,naxis1], yrange = [1,naxis2],  bytscl(data, min =min(data[0:500]), max = max(data[0:500])) ,/preserve_aspect,  ncolors=8

  
   ;in ra and dec this time, since that is what Mihos needs
   ;this is the bright, central part of the plume
  raroi = [186.889,186.86234,186.8575,186.86722,186.88338,186.89146]
  decroi = [13.228459,13.229257,13.189904,13.144259,13.136393,13.163154]

  ;the full extended plume
   ;raroi = [186.889,186.79603,186.80092,186.84619,186.8559,186.88338,186.89146]
   ;decroi = [13.228459,13.254413,13.193027,13.167863,13.14583,13.136393,13.163154]

   ;but the IRAS images are in B1950, so
   ;precess them back to FK4 coords
  bprecess, raroi, decroi, ra_1950, dec_1950
  
  adxy, header, ra_1950, dec_1950, px2, py2
  
  photroi = Obj_New('IDLanROI', px2, py2) 
;   draw_roi,photroi             ;,/line_fill,color='444444'x 

;create a mask out of the roi
   mask = photroi->ComputeMask(Dimensions = [naxis1, naxis2], mask_rule =2)

;apply the mask to make all data outside the roi = 0
   onlyplume = data*(mask GT 0)
   goodarea = n_elements(where(mask gt 0))
   sq2 = fix(sqrt(goodarea) / 2) ; what is the equivalent square with the same area

   print, 'goodarea', goodarea
   ;plotimage, xrange=[1,naxis1], yrange = [1,naxis2],  bytscl(onlyplume, min =min(data[0:500]), max = max(data[0:500])) ,/preserve_aspect,  ncolors=8

   ;yea photometry
   object = (total(onlyplume)) 


;get the sky background below

;look at multiple regions of that same size around the frame, not on
;the plume, their stddev is the important number.

;mask regions with vastly different exptime
;   a = where(rrdata lt explow)
 ;  b = where(rrdata gt exphigh)
 ;  nan = alog10(-1)
 ;  data[a] = nan
;   data[b] = nan

;mask regions with object detections
;   a = where(segdata gt 0)
;   data[a] = nan

;mask the plume region
;   data[1821:2146,1866:2259] = nan

;randomly choose centers to look at background regions
   nrand = 100
   x = randomu(S,nrand) * 15. + 2.
   y = randomu(S,nrand) * 20. +34.

   
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
   start = [mean(bkgd),stddev(bkgd), 10.]
   noise = fltarr(n_elements(yhist))
   noise[*] = 1                                            ;equally weight the values
   result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) ;./quiet    

;now plot the histogram of stddevs of all unmasked pixels
   plothist, bkgd, xjunk, yjunk, bin = binsize, xtitle = 'Background (Mjy/sr)', ytitle='Number', $
             xrange=[-8*result(1) + result(0),result(0) + 8*result(1)];,$
 ;           thick=3, xthick=3,ythick=3,charthick=3, xticks = 2

   ;xyouts, result(0) + 3.5*result(1), max(yhist) / 2.  + 0.08*max(yhist), 'mean' + string(result(0)) ;, charthick=3
   ;xyouts, result(0) + 3.5*result(1), max(yhist) / 2. ,  'sigma' + string(result(1))                 ;, charthick=3
;xyouts, result(0) + 3.5*result(1), max(yhist) / 2. + 2*0.08*max(yhist), 'N pixels used' + string(count), charthick=3

;plot the fitted gaussian and print results to plot

   xarr = findgen(1000)/ 100.

   oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
           thick=5, linestyle =  2    ;color = greencolor,
   ;print, 'gauss fit result', result


   ;coming back to the photometry....
   sky = result(0); total(data[1678.-37.:1678.+37.,2251.-75.:2251.+75.])
   flux =  (object / goodarea) - sky

;   print, 'flux ', flux, ' counts/s over whole aperture', object, sky
;   print, 'flux', flux/ (37.*2.*2.*75.)  , ' counts/s/pix'
   print, 'flux, sky, sigma', flux, sky, result(1)


   ;account for galactic extinction
   ;first take couts/s/pix into mag
   ;m = -2.5*alog10(flux) +zpt
   ;m = m - AUV    ;apply extinction
   ;flux = 10^((m - zpt) / (-2.5));back into counts/s/pix

   ;arcsec_pixel = 1.5
;   SB = zpt -2.5*alog10(flux/11325.) + 5.*alog10(arcsec_pixel)
   ;SB = zpt -2.5*alog10(flux) + 5.*alog10(arcsec_pixel)

   ;print, 'SB ', SB, ' AB mag/sqarcsec'
;   print, 'flux', flux_zpt * flux / (37*2*75*2), ' erg/s/cm2/A'

   ;now want total flux in counts/s (with
   ;extinction correction already in there
   ;flux = flux*goodarea
   ;print, 'flux', flux, ' counts/s'

   ;f lambda
   ;flambda = flux_zpt * flux 
   ;print, 'flux', flambda , ' erg/s/cm2/A'
   
   ;the extra 10% is from the error in the absolute photometry
   ;error = result(1)*1.1*goodarea
   ;errorlambda = flux_zpt*error
   ;print, 'error', errorlambda, ' erg/s/cm2/A'
   
   ;convert these into Fnu for eric (instead of Flambda)
   ;go into AB mags first (from the galex zpts)
   ;abmag = -2.5*alog10(flux) + zpt
   ;print, 'AB mag', abmag

   ;fnu = magab_to_flux(abmag)
   ;print, 'flux ', fnu, ' erg/s/cm2/Hz'

   ;magerr = -2.5*alog10(flux - error) + zpt - abmag
   ;print, 'magerr', magerr
   ;fnuerr =  fnu - magab_to_flux(abmag + magerr) 
   ;print, 'error', fnuerr, ' erg/s/cm2/Hz'

   ;not working?
 ;c = 299792458 ;m/s
   ;c = c / 1E-10 ; angstroms / s
   ;fnu = flambda * (lambda^2) / c
   ;print, 'flux ', fnu, ' erg/s/cm2/Hz'
   
  ; errornu = errorlambda*(lambda^2 )/ c
  ; print, 'error', errornu, ' erg/s/cm2/Hz'


  

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
   if ch eq 0 then fits_write, '/Users/jkrick/Virgo/IRAS/I100_masked.fits', onlyplume, header
   if ch eq 1 then fits_write, '/Users/jkrick/Virgo/IRAS/I60_masked.fits', onlyplume, header
   
endfor
;ps_close, /noprint,/noid

end
