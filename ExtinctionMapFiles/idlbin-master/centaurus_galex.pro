pro centaurus_galex
;measure accurate photometry and background errors for Galex data of
;ICL plume in the virgo cluster
ps_open, filename='/Users/jkrick/plumes/centaurus/galex/background.ps',/portrait,/square,/color

!P.multi = [0,2,2]
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
nan = alog10(-1)

for ch = 0, 1 do begin  ;for fuv and nuv
   if ch eq 0 then begin
      fits_read, '/Users/jkrick/plumes/centaurus/galex/AIS_432_50432_0001_sv95/AIS_432_sg95-fd-int.fits', data, header
      fits_read, '/Users/jkrick/plumes/centaurus/galex/AIS_432_50432_0001_sv95/AIS_432_sg95-fd-rrhr.fits', rrdata, rrheader
      zpt = 18.82
      flux_zpt = 1.4E-15
      explow = 83
      exphigh = 98
      binsize=0.00005
      AUV = 7.9*0.111 ;galactic extinction in mags
      lambda = 1516 ;angstroms
      print, '****FUV*****'
   endif

   if ch eq 1 then begin
      fits_read, '/Users/jkrick/plumes/centaurus/galex/AIS_432_50432_0001_sv95/AIS_432_sg95-nd-int.fits', data, header
      fits_read, '/Users/jkrick/plumes/centaurus/galex/AIS_432_50432_0001_sv95/AIS_432_sg95-nd-rrhr.fits', rrdata, rrheader
      zpt = 20.08 ;flux should be in photons per second
      flux_zpt = 2.06E-16
      explow = 70
      exphigh = 80
      binsize=0.0001
      AUV = 8.0*0.111  ; ..111 from E(B-V) in schlegel maps
      lambda = 2267 ;angstroms
      print, '****NUV*****'
   endif
;starting off, the units are counts
;   data = data /rrdata
;now the image units are counts/s
;nevermind, this is already done in the image.

;roi photometry
   naxis1 = fxpar(header, 'naxis1')
    naxis2 = fxpar(header, 'naxis2')

   ;in ra and dec 
  raroi1 = [192.56815,192.5636,192.52335,192.44121,192.44233]
  decroi1 = [-41.473335,-41.463567,-41.446162,-41.410885,-41.419382]

  adxy, header, raroi1, decroi1, px1, py1
  
;look at one region for now.
   photroi = Obj_New('IDLanROI', px1, py1) 
;   draw_roi,photroi             ;,/line_fill,color='444444'x 

;create a mask out of the roi
   mask = photroi->ComputeMask(Dimensions = [naxis1, naxis2], mask_rule = 1)

;apply the mask to make all data outside the roi = 0
   onlyplume = data*(mask GT 0)
   goodarea = n_elements(where(mask gt 0))
   sq2 = fix(sqrt(goodarea) / 2) ; what is the equivalent square with the same area

   print, 'goodarea', goodarea
;   plotimage, xrange=[1,naxis1], yrange = [1,naxis2],  bytscl(onlyplume, min =min(data[0:500]), max = max(data[0:500])) ,/preserve_aspect,  ncolors=8

   ;yea photometry
   object = (total(onlyplume)) 

;get the sky background below

;look at multiple regions of that same size around the frame, not on
;the plume, their stddev is the important number.

;mask regions with vastly different exptime
   a = where(rrdata lt explow)
   b = where(rrdata gt exphigh)
   data[a] = nan
   data[b] = nan

;mask regions with object detections
;   a = where(segdata gt 0)
;   data[a] = nan

;generously mask the plume region
   data[1200:1455,1895:2063] = nan

;randomly choose centers to look at background regions
;this needs to be more complicated because it is a round image, actually maybe the exptime requirement will take care of this.
   nrand = 400
   x = randomu(S,nrand) * 3100. + 400.
   y = randomu(S,nrand) * 3000. + 500.

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
   result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) 

;now plot the histogram of stddevs of all unmasked pixels
   plothist, bkgd, xjunk, yjunk, bin = binsize, xtitle = 'Background (counts/s/pix)', ytitle='Number', $
             xrange=[-8*result(1) + result(0),result(0) + 8*result(1)] ,$
            thick=3, xthick=3,ythick=3,charthick=3, xticks = 2

   ;xyouts, result(0) + 3.5*result(1), max(yhist) / 2.  + 0.08*max(yhist), 'mean' + string(result(0)) ;, charthick=3
   ;xyouts, result(0) + 3.5*result(1), max(yhist) / 2. ,  'sigma' + string(result(1))                 ;, charthick=3
;xyouts, result(0) + 3.5*result(1), max(yhist) / 2. + 2*0.08*max(yhist), 'N pixels used' + string(count), charthick=3

;plot the fitted gaussian and print results to plot

   xarr = findgen(10000)/ 1000000.

   oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
           thick=5, linestyle =  2    ;color = greencolor,
   print, 'gauss fit result', result


   ;coming back to the photometry....
   sky = result(0); total(data[1678.-37.:1678.+37.,2251.-75.:2251.+75.])
   flux =  (object / goodarea) - sky

;   print, 'flux ', flux, ' counts/s over whole aperture', object, sky
;   print, 'flux', flux/ (37.*2.*2.*75.)  , ' counts/s/pix'
   print, 'flux', flux, ' counts/s/pix'


   ;account for galactic extinction
   ;first take counts/s/pix into mag
   m = -2.5*alog10(flux) +zpt
   m = m - AUV    ;apply extinction
   flux = 10^((m - zpt) / (-2.5));back into counts/s/pix

   arcsec_pixel = 1.5
;   SB = zpt -2.5*alog10(flux/11325.) + 5.*alog10(arcsec_pixel)
   SB = zpt -2.5*alog10(flux) + 5.*alog10(arcsec_pixel)

   print, 'SB ', SB, ' AB mag/sqarcsec'
;   print, 'flux', flux_zpt * flux / (37*2*75*2), ' erg/s/cm2/A'

   ;now want total flux in counts/s (with
   ;extinction correction already in there
   flux = flux*goodarea
   print, 'flux', flux, ' counts/s'

   ;f lambda
   flambda = flux_zpt * flux 
   print, 'flux', flambda , ' erg/s/cm2/A'
   
   ;the extra 10% is from the error in the absolute photometry
   error = result(1)*1.1*goodarea
   errorlambda = flux_zpt*error
   print, 'error', errorlambda, ' erg/s/cm2/A'
   
   ;convert these into Fnu for eric (instead of Flambda)
   ;go into AB mags first (from the galex zpts)
   abmag = -2.5*alog10(flux) + zpt
   print, 'AB mag', abmag

   fnu = magab_to_flux(abmag)
   print, 'flux ', fnu, ' erg/s/cm2/Hz'

   magerr = -2.5*alog10(flux - error) + zpt - abmag
   print, 'magerr', magerr, flux, error, zpt, abmag
   fnuerr =  fnu - magab_to_flux(abmag + magerr) 
   print, 'error', fnuerr, ' erg/s/cm2/Hz'

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
;   if ch eq 0 then fits_write, '/Users/jkrick/plumes/centaurus/galex/AIS_113_50432_0001_sv95/fd_masked.fits', fuvdata, fuvheader
;   if ch eq 1 then fits_write, '/Users/jkrick/plumes/centaurus/galex/AIS_113_50432_0001_sv95/nd_masked.fits', fuvdata, fuvheader
   
endfor
ps_close, /noprint,/noid

end
