pro coma_acs
;measure accurate photometry and background errors for ACS data of
;ICL plume in the Coma cluster

!P.multi = [0,1,1]
nan = alog10(-1)

for ch = 0, 2 do begin          ;for 435,606,814
   if ch eq 0 then begin
      fits_read, '/Users/jkrick/plumes/coma/hst/acs/acs435_drz.fits', exten_no = 1, data, header
      fits_read, '/Users/jkrick/plumes/coma/hst/acs/acs435_drz.fits', exten_no = 2 ,rrdata, rrheader
      fits_read, '/Users/jkrick/plumes/coma/hst/acs/acs435_segmentation.fits', segdata, segheader
      zpt = 25.68392 ;http://www.stsci.edu/hst/acs/analysis/zeropoints
      photflam = fxpar(header, 'PHOTFLAM')
      A_lambda = .038 ;galactic extinction in mags
      print, '****F435*****'
   endif

   if ch eq 1 then begin
      fits_read, '/Users/jkrick/plumes/coma/hst/acs/acs606_drz.fits', exten_no = 1, data, header
      fits_read, '/Users/jkrick/plumes/coma/hst/acs/acs606_drz.fits', exten_no = 2 ,rrdata, rrheader
      fits_read, '/Users/jkrick/plumes/coma/hst/acs/acs606_segmentation.fits', segdata, segheader
      zpt = 26.50512
      photflam = fxpar(header, 'PHOTFLAM')
      A_lambda = .026 ;galactic extinction in mags
      print, '****F606*****'
   endif
   if ch eq 2 then begin
      fits_read, '/Users/jkrick/plumes/coma/hst/acs/acs814_drz.fits', exten_no = 1, data, header
      fits_read, '/Users/jkrick/plumes/coma/hst/acs/acs814_drz.fits', exten_no = 2 ,rrdata, rrheader
      fits_read, '/Users/jkrick/plumes/coma/hst/acs/acs814_segmentation.fits', segdata, segheader
      zpt = 25.95928
      photflam = fxpar(header, 'PHOTFLAM')
      A_lambda = .017 ;galactic extinction in mags
      print, '****F814*****'
   endif
;starting off, the units are electrons/second 
;convert them to flux units of ergs/s/cm2/Ang
;do this later
  ; data = data *photflam

;mask regions with object detections using segmentation image
; this is not perfect, but a good start
   a = where(segdata gt 0)
   data[a] = nan
  
;roi photometry
   naxis1 = fxpar(header, 'naxis1')
   naxis2 = fxpar(header, 'naxis2')
   
                                ;in ra and dec 
   raroi1 = [194.88847,194.82077,194.83109,194.84896,194.85072,194.89023]
   decroi1 = [28.008323,27.992767,27.989325,27.994105,27.990994,28.000434]
   
   adxy, header, raroi1, decroi1, px1, py1
   
;look at one region for now.
   photroi = Obj_New('IDLanROI', px1, py1) 
;   draw_roi,photroi             ;,/line_fill,color='444444'x 
   
;create a mask out of the roi
   mask = photroi->ComputeMask(Dimensions = [naxis1, naxis2], mask_rule = 1)
   
;apply the mask to make all data outside the roi = 0
   onlyplume = data*(mask GT 0)
   goodarea = n_elements(where(onlyplume gt 0 ))
   sq2 = fix(sqrt(goodarea) / 2) ; what is the equivalent square with the same area
   
   print, 'goodarea', goodarea
   ;plotimage, xrange=[1,naxis1], yrange = [1,naxis2],  bytscl(onlyplume, min =-0.01, max = 0.01) ,/preserve_aspect,  ncolors=8

   ;yea photometry
   object = (total(onlyplume,/nan)) 
   print, 'object', object/goodarea
;-------------------------------
;get the sky background below

;look at multiple regions of that same size around the frame, not on
;the plume, their stddev is the important number.
;mask regions with vastly different exptime
; this is just the edges, and since I am choosing the background by hand, don't need to do this. 

;since there is not much background, will need to specifically choose the background regions.
   rabkg1=[194.84859,194.84551,194.82667,194.85132]
   decbkg1=[28.032145,28.032945,27.999345,28.004148]
   rabkg2=[194.88375,194.85023,194.85204,194.87759]
   decbkg2=[28.014865,28.030705,28.007188,28.012786]
   adxy, header, rabkg1, decbkg1, bx1, by1
   adxy, header, rabkg2, decbkg2, bx2, by2
   
   bkg1roi = Obj_New('IDLanROI', bx1, by1) 
   bkg2roi = Obj_New('IDLanROI', bx2, by2) 
      
   bkg1mask = bkg1roi->ComputeMask(Dimensions = [naxis1, naxis2], mask_rule = 1)
   bkg2mask = bkg2roi->ComputeMask(Dimensions = [naxis1, naxis2], mask_rule = 1)

;apply the mask to make all data outside the roi = 0
   bkg1 = data*(bkg1mask GT 0)
   bkg2 = data*(bkg2mask GT 0)

  ;plotimage, xrange=[1,naxis1], yrange = [1,naxis2],  bytscl(bkg2, min =-0.01, max = 0.01) ,/preserve_aspect,  ncolors=8


   bkg1area = n_elements(where(bkg1mask gt 0))   
   bkg2area = n_elements(where(bkg2mask gt 0))   
   print, 'bkgarea', bkg1area, bkg2area, 'pixels'

   bkg1_perpixel = (total(bkg1,/nan)) / bkg1area
   bkg2_perpixel = (total(bkg2,/nan)) / bkg2area


   print, 'bkg', bkg1_perpixel, bkg2_perpixel,  (bkg1_perpixel+bkg2_perpixel) / 2., 'e-/2/pixel'
   ;print, 'stddev', stddev(bkg1,/nan), stddev(bkg2,/nan), (stddev(bkg1, /nan) + stddev(bkg2,/nan))/2., 'e-/2/pixel'

  error = [bkg1_perpixel, bkg2_perpixel]
   error = stddev(error)
   print, 'error', error ,  'e-/s '

   ;coming back to the photometry....
   sky = (bkg1_perpixel+bkg2_perpixel) / 2.;= result(0); total(data[1678.-37.:1678.+37.,2251.-75.:2251.+75.])
   flux =  (object / goodarea) - sky

   print, 'flux', flux, ' e-/s in one pixel'

   ;account for galactic extinction
   abmag = -2.5*alog10(flux) + zpt
   abmag = abmag - A_lambda    ;apply extinction
   print, 'AB mag', abmag
   ;flux = 10^((m - zpt) / (-2.5));back into flux units

;;------------------------------------------------
;unit roundup
   arcsec_pixel = 0.0499
   SB = zpt -2.5*alog10(flux) + 5.*alog10(arcsec_pixel)

   print, 'SB ', SB, ' AB mag/sqarcsec'
;   print, 'flux', flux_zpt * flux / (37*2*75*2), ' erg/s/cm2/A'

   ;now want total flux in counts/s (with
   ;extinction correction already in there
   flux_erg = flux*goodarea*photflam
   print, 'flux', flux_erg, 'erg/s/cm2/A '
  ; error = error *goodarea*photflam
   ;print, 'error', error, 'erg/s/cm2/A'

   noise_erg = error*goodarea*photflam
   print, 'noise', noise_erg,  'erg/s/cm2/A '
   print, 'SNR', flux_erg/noise_erg

   ;convert these into Fnu for eric (instead of Flambda)
   ;go into AB mags first 
   
   fnu = magab_to_flux(abmag)
   print, 'flux ', fnu, ' erg/s/cm2/Hz'

   magerr = -2.5*alog10(flux - error) + zpt - abmag
   print, 'magerr', magerr, flux, error, zpt, abmag
   fnuerr =  fnu - magab_to_flux(abmag + magerr) 
   print, 'error', fnuerr, ' erg/s/cm2/Hz'

endfor

end
