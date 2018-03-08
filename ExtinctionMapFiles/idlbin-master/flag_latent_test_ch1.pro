pro flag_latent_test_ch1


  sq2 = 2
  stopthresh = 4                ;stop correcting when we find 4 in a row below background detection
;this allows for HDR frames to have 2 deep images in a row.
  
  fits_read, '/Users/jkrick/irac_warm/latent/flag_latent/pipeline_test/2849/1/38763008/cbcd_fp_0015.fits', data, header
  currentsclk =  sxpar(header, 'SCLK_OBS') ;keep track of the time stamp
  naxis = sxpar(header, 'NAXIS')           ;figure out if subarray or fullarray
  exptime = sxpar(header, 'EXPTIME')       ;need to know exposure time for unit conversion later
  readmode = sxpar(header, 'READMODE')     ; FULL or SUB array
  
 
            ;how bright was the original source in electrons
  gain = sxpar(header, 'GAIN')
  pixel_scale = abs(sxpar(header, 'PXSCAL2')) ;need the absolute value, sometimes they are negative.
  flux_conv = sxpar(header, 'FLUXCONV')
  known_flux_electron =29955982.0
            
            ;use that to measure the y_intercept coefficient
  y_int = -5614.36 + 880.124*alog10(known_flux_electron)
  print, ' electron, y_int', known_flux_electron, y_int

                                ;measure background
  bkg = measure_bkg(192., 39., data, sq2, gain, exptime, flux_conv)
  print, 'bkg 192, 39', bkg

  ;bkg = measure_bkg(38., 126., data, sq2, gain, exptime, flux_conv)
  ;print, 'bkg 38,126', bkg


end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function measure_bkg, xcen, ycen, data, sq2, gain, exptime, flux_conv

;In order to find an accurate background I first randomly choose
;regions in the frame.  I then make a histogram of the values in those
;regions.  I fit the histogram with a gaussian and take the mean of
;the gaussian as the mean background value.

  ;need to convert the data into electrons
  data = data * gain*exptime/flux_conv

  nrand = 100
  ;choose random locations that are not too close to the edges of the frame
  xbkg_radius = xcen - 1 < 36 < (254 - xcen)
  ybkg_radius = ycen - 1 < 36 < (254 - ycen)
print, 'xcen, ycen', xcen, ycen
  print, 'xbkg_radius', xbkg_radius
  print, 'ybkg_radius', ybkg_radius

  x = randomu(S,nrand) *(xbkg_radius * 2) + xcen - xbkg_radius
  y = randomu(S,nrand) *(ybkg_radius * 2) + ycen - ybkg_radius
 
  print, 'x',x
print, 'y', y

  bkgd = fltarr(nrand)
  for ti = 0, n_elements(x) - 1 do begin  ; for each random background region
     pixval = data[x(ti)-sq2:x(ti)+sq2,y(ti)-sq2:y(ti)+sq2]
     ga = where(finite(pixval) gt 0)
     if n_elements(ga) gt .85*n_elements(pixval) then  begin
        bkgd(ti) =  total(data[x(ti)-sq2:x(ti)+sq2,y(ti)-sq2:y(ti)+sq2],/nan) ; / n_elements(ga)
     endif

  endfor
  ta = where(bkgd ne 0)
  bkgd = bkgd(ta)
  ;print, 'n_ta', n_elements(ta)
  if n_elements(ta) lt 10 then begin  ; if there aren't enough regions to fit a gaussian.
     result = [mean(bkgd), stddev(bkgd), 2.0]
  endif else begin

;do some sigma clipping on bkgd before running it through plothist

     skpix = MEDCLIP(bkgd, med, sigma, clipsig = 3, converge_num = 0.03)
     bkgd = skpix
;fit a gaussian to this distribution
     binsize = 200           
     plothist, bkgd, xhist, yhist, bin = binsize, xrange = [-10000,10000] ;,/noplot ;make a histogram of the data
;     print, 'plothist xhist', xhist
;     print, 'plothist yhist', yhist

     start = [mean(bkgd),stddev(bkgd), 2.]
    ; print, 'start', start
     noise = fltarr(n_elements(yhist))
     noise[*] = 1  
                                ;equally weight the values
     result= MPFITFUN('mygauss',xhist,yhist, noise, start);,/quiet) ; fit a gaussian to the histogram sorted data
     bkg = result(0)
     oplot, xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.)
  endelse

     return, result
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function quick_phot, xcen, ycen, data, sq2, imheader
           ;get rid of NaN's
  data[where(finite(data) lt 1)] = 0
           
          ;convert the data from MJy/sr to electrons
  gain = sxpar(imheader, 'GAIN')
  gain = 3.7                    ;force it to the correct value since currently it is wrong in most headers
  exptime =  sxpar(imheader, 'EXPTIME') 
  fluxconv =  sxpar(imheader, 'FLUXCONV') 
  time = sxpar(imheader, 'SCLK_OBS')
  sbtoe = gain*exptime/fluxconv
  data = data*sbtoe
           
           ;do photometry in a 5x5 box around the object
  flux = total(data[xcen - sq2:xcen +sq2, ycen -sq2:ycen  + sq2])
  return, flux
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mygauss, X, P
return, P(2)/sqrt(2.*!Pi) * exp(-0.5*((X - P(0))/P(1))^2.)

;P(0) = mean
;P(1) = sigma
;P(2) = area
end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function Mjypersrtojansky,  Mjypersr, pixel_scale


; Convert from summed MJy/sr to Jy
scale = pixel_scale^2 ;1.22D * 1.22D
; convert scale from arcsec^2 to sr and scale to Jy
scale = scale * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06

jansky = Mjypersr * scale

return, jansky

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


function mjy_to_electron, mjy, pixel_scale, gain, exptime, flux_conv

jy = mjy/1000.  ;milijanskies to janskies

; Convert from summed  to Jy to Mjy/sr
scale = pixel_scale^2 ;1.22D * 1.22D
; convert scale from arcsec^2 to sr and scale to Jy
scale = scale * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06

Mjypersr = jy / scale

electrons = Mjypersr * gain*exptime/flux_conv


return, electrons

end
