pro flag_latent

sq2 = 2

;need to generate this file which has all the satstars triggers for
;the entire PAO in it.

;need to read in all imasks's in the entire PAO and make a file
;with the locations of the saturated stars.

;pao 3760
;aorname=['38927616','38927872','38928128','38928384','38926080','38926336','38928640','38926592','38928896','38926848','38929152','38927104','38929408','38927360','38933248','38920192','38925056','38925568','38925312','38925824','38924544','38924800','38923776','38924032','38924288','35418624','35422208','31537152','38763008','38788352','37908480']

dirloc = '/Users/jkrick/irac_warm/latent'
aorname =['38763008']

readcol, '/Users/jkrick/irac_warm/latent/satcorr_phot.tbl', x    ,  y    ,    ra   ,    dec  ,replace,Pratio,Fratio, ferr , hsig , loc_sky, skysig, radius,    raw  ,   model  ,        rectify , skipline = 3


for a = 0,  n_elements(aorname) - 1 do begin
   print, 'aor', aorname(a)
   CD,dirloc+'r' +aorname(a) + '/ch1/bcd'
 

;read in all the individual satcorr_phot files, make them into one
;file with time stamps from sclk

endfor

for a = 0,  n_elements(aorname) - 1 do begin
   print, 'aor', aorname(a)
   CD,dirloc+'r' +aorname(a) + '/ch1/bcd'

  command1 =  ' ls SPITZER_I1*_bcd.fits > bcd.txt'   ;list all images in that directory
   spawn, command1
 
 ;read in the list of images
   readcol, dirloc+'r'+aorname(a)+'/ch1/bcd/bcd.txt', bcdlist, format='A', /silent ;/ch1/bcd
   for f = 0, n_elements(bcdlist) - 1 do begin
      fits_read, bcdlist(f), data, imheader
      print, bcdlist(f)
      sclk =  sxpar(imheader, 'SCLK_OBS')
      lat = where(sclk gt timearr, count)
                                ;lat(count) = the stars which need to
                                ;have their latent tested against the background
      if count gt 0 then begin
           ;measure background
         bkg = measure_bkg(xcen(lat(count)), ycen(lat(count)), data, sq2)
         
         latent_flux = final_yint*exp(-sclk/4.4)
         if latent_flux gt 2*bkg(1) then begin ;or stoplooking lt stopthresh then begin
            print, 'flag the latent'
         endif
      endif
     
   endfor
   
endfor

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function measure_bkg, xcen, ycen, data, sq2
;Maybe this is already done somewhere else?

;In order to find an accurate background I first randomly choose
;regions in the frame.  I then make a histogram of the values in those
;regions.  I fit the histogram with a gaussian and take the mean of
;the gaussian as the mean background value.

  nrand = 100
  ;choose random locations that are not too close to the edges of the frame
  xbkg_radius = xcen - 1 < 36 < (254 - xcen)
  ybkg_radius = ycen - 1 < 36 < (254 - ycen)
  x = randomu(S,nrand) *(xbkg_radius * 2) + xcen - xbkg_radius
  y = randomu(S,nrand) *(ybkg_radius * 2) + ycen - ybkg_radius
 
  bkgd = fltarr(nrand)
  for ti = 0, n_elements(x) - 1 do begin  ; for each random background region
     pixval = data[x(ti)-sq2:x(ti)+sq2,y(ti)-sq2:y(ti)+sq2]
     ga = where(finite(pixval) gt 0)
     if n_elements(ga) gt .85*n_elements(pixval) then  bkgd(ti) =  total(data[x(ti)-sq2:x(ti)+sq2,y(ti)-sq2:y(ti)+sq2],/nan) ; / n_elements(ga)
  endfor
  ta = where(bkgd ne 0)
  bkgd = bkgd(ta)
     
;fit a gaussian to this distribution
  binsize = 200            
  plothist, bkgd, xhist, yhist, bin = binsize,/noplot ;make a histogram of the data
  start = [1000.,200., 2.]
  noise = fltarr(n_elements(yhist))
  noise[*] = 1                                                    ;equally weight the values
  result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet)    ; fit a gaussian to the histogram sorted data
  bkg = result(0)

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


;P(1) = sigma
;P(0) = mean
;P(2) = area
end

;dirloc = '/Users/jkrick/irac_warm/latent/pao3685/'
;aorname = ['38140160','37107968','37115904','37084416','37084928','37120256','37110016','37097984','37093120','32245504','32310528','31980032','32208384','32243712','34771200','34772224','36068864','32030464','32146432','32207360','29439488','29423104','29392384','29381632','29436928','36127744','36182528','35479552','36235008','32634624','36103936','36144896','36123648','36099328','36093952','38163712','38163968','38163200','38163456','38162688','38162944','38164224','38164480','38166272','38166528','38165248','38165504','38164736','38164992','38165760','38166016']          

;open a cal file
;openw, outlunsat, '/Users/jkrick/irac_warm/latent/pao3685/sat.txt', /get_lun

;for a = 0, 1 do begin; n_elements(aorname) - 1 do begin

;   CD,dirloc+'r' +aorname(a) + '/ch1/bcd'
;   command1 =  ' ls SPITZER_I1*imsk.fits >  imask.txt'   ;list all images in that directory
;     spawn, command1
  
  ;read in the list of images
;     readcol, dirloc+'r'+aorname(a)+'/ch1/bcd/imask.txt', imasklist, format='A', /silent ;/ch1/bcd
;     for f = 0, n_elements(imasklist) - 1 do begin
;        fits_read, imasklist(f), data, imheader
        
;need to fix this to get all of them
;        b13 = where(data eq 2^13,c13)
;        b4 = where(data eq 2^4,c4)
;        print, '13', c13
;        if c13 gt 0 then begin
;           for n = 0, n_elements(c13) -1 do begin
;              printf, outlunsat, imasklist(f), '13', b13(n)
;           endfor
;        endif
;        if c4 gt 0 then begin
;           for n = 0, n_elements(c4) -1 do begin
;              printf, outlunsat, imasklist(f), '4', b4(n)
;           endfor
 ;       endif

;     endfor

;endfor

;close, outlunsat
;free_lun, outlunsat


;this isn't working so start with a fake cal file that I can
;use and come back to this.

;timearr = [943394521,943394521,943397881]
;bcdname = ['/Users/jkrick/irac_warm/latent/pao3685/r38163200/ch1/bcd/SPITZER_I1_38163200_0001_0000_1_bcd.fits','/Users/jkrick/irac_warm/latent/pao3685/r38163200/ch1/bcd/SPITZER_I1_38163200_0001_0000_1_bcd.fits', '/Users/jkrick/irac_warm/latent/pao3685/r38162688/ch1/bcd/SPITZER_I1_38162688_0001_0000_1_bcd.fits']
;xcen = [45,201,120]
;ycen = [116,143,201]
;peak_flux = [3.74E8,1.3E8, 7.2E8]
