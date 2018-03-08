pro background_comp
  !P.multi =[0,2,2]
ps_open, filename='/Users/jkrick/irac_warm/background_comp.ps',/portrait,/square,/color

;-----------------------
;example 1
  dirloc = '/Users/jkrick/irac_warm/pc9/IRAC023000/bcd/'
;list of a bunch of AOR's after the AOR with the bright star
;I don't know how to get this information automatically
  aorname = ['0038163456', '0038164480','r32208128/ch1/bcd']
  ;the saturated star is in frame 2
  startframe = 2
  xcen =88.04                ;,162.0,42.,201.0,120.0,59.,218.]
  ycen = 133.69              ;,115.48,116.,143.0,202.,160.,188.]
;--------------------------

sq2 = 2
c = 0
timearr = fltarr(10000)
bkgarr = fltarr(10000)
framarr = fltarr(10000)
zodiarr = fltarr(10000)
sigmarr = fltarr(10000)
predzodiarr = fltarr(10000)
for a =0, n_elements(aorname) - 1 do begin

     CD,dirloc+aorname(a)                       ;+ '/ch1/bcd'
     if a lt 2 then command1 =  ' ls IRAC.1*bcd*.fits >  bcd.txt' ;list all images in that directory
     if a eq 2 then command1 =  ' ls SPITZER_I1*_bcd*.fits >  bcd.txt' ;list all images in that directory
     spawn, command1
  
  ;read in the list of images
     readcol, dirloc+aorname(a)+ '/bcd.txt', bcdlist, format='A', /silent ;/ch1/bcd
                                ;just read in the image with the saturated star on it
     print, 'startframe', bcdlist(startframe+1)


     for f = 0, n_elements(bcdlist) - 1 do begin 
        fits_read, bcdlist(f), data, imheader

  ;set the initial time when the sat star was on the frame
        timearr[c] = sxpar(imheader, 'SCLK_OBS')
        framarr[c] = sxpar(imheader, 'FRAMTIME')
        zodi= sxpar(imheader, 'ZODY_EST')
        flux_conv = sxpar(imheader, 'FLUXCONV')
        gain = sxpar(imheader, 'GAIN')
        rdnoise = sxpar(imheader, 'RONOISE')  ;in electrons already
print, rdnoise
        ;convert to electrons
        zodi = zodi*gain*framarr[c]/flux_conv
        zodiarr[c] =zodi * 25 ;  there are 5x5 summed pixels in the background.
   ;;;;;;;;;
        bkg = measure_bkg( xcen, ycen, data, sq2)
        bkgarr[c] = bkg(0)*gain*framarr[c]/flux_conv
        sigmarr[c] = bkg(1)*gain*framarr[c]/flux_conv
        predzodiarr[c] = sqrt(zodi*25) + rdnoise*25

        c = c + 1
     endfor

  endfor


timearr = timearr[0:c-1]
bkgarr = bkgarr[0:c-1]
framarr = framarr[0:c-1]
zodiarr = zodiarr[0:c-1]
sigmarr = sigmarr[0:c-1]
predzodiarr = predzodiarr[0:c-1]


plot,( timearr - timearr(0)) / 3600, bkgarr,  psym = 2, xtitle = 'time(hrs)', ytitle = 'total background in aperture (electrons)', yrange = [0,3000]
plot,( timearr - timearr(0)) / 3600, sigmarr,  psym = 2, xtitle = 'time(hrs)', ytitle = 'sigma on background distribution (electrons)' , yrange = [0, 1000]
plot,( timearr - timearr(0)) / 3600, zodiarr,  psym = 2, xtitle = 'time(hrs)', ytitle = 'total zodiacal light in aperture (electrons)', yrange = [0,3000]
xyouts, 0, 1300, '12s'
xyouts, 2, 2700, '30s'
xyouts, 7, 300, '2s'
plot,( timearr - timearr(0)) / 3600, sqrt(zodiarr),  psym = 2, xtitle = 'time(hrs)', ytitle = 'sqrt(total zodi in aperture) (electrons)' , yrange = [0, 1000]


print, framarr

ps_close, /noprint,/noid

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
  binsize = 2  
  plothist, bkgd, xhist, yhist, bin = binsize,/noplot ;make a histogram of the data
  start = [0.,10., 2.]
  noise = fltarr(n_elements(yhist))
  noise[*] = 1                                                    ;equally weight the values
  result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) ;/quiet   ; fit a gaussian to the histogram sorted data
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



