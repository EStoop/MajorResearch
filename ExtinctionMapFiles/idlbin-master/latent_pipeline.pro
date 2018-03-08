pro latent_pipeline
  
;given the master_satcorr_phot.tbl from the pipeline, can I recreate the latentflag.tbl
;this is used to check the INT pipeline

sq2 = 2
stopthresh = 4  ;stop correcting when we find 4 in a row below background detection
;this allows for HDR frames to have 2 deep images in a row.

dirloc = '/Users/jkrick/iracdata/S19INT/vers2/'
openw, outlun, '/Users/jkrick/irac_warm/latent/flag_latent/int_aug2011/latent_generated.tbl',/get_lun
;time order of the PAO 3630 from IRAC022700 

aorname = ['r38055168','r38055424','r38055680','r38055936','r38053632','r38053888','r38056192','r38054144','r38056448','r38054400','r38056704','r38054656','r38056960','r38054912','r38057216','r38048256','r38052608','r38053120','r38052864','r38053376','r38052096','r38052352','r38051328','r38051584','r38051840','r32908800','r33220608','r34978304','r32145920','r32459520','r31984896','r32042752','r32329216','r32173056','r32145408','r32309504','r32238848','r32433920','r32176896','r32188160','r32092416','r32410112','r32351232','r32241408','r31974912','r32201472','r32872704','r32840192','r33220352','r34700032','r34745600','r34857728','r32899072','r32335872','r32395520','r32097536','r32145664','r32459264','r31984384','r32042496','r32328960','r32172800','r32144896','r32308992','r35467008','r30481408','r30597632','r31028992','r30951424','r33220096','r37264640','r37256960','r37278720','r37278464','r35373312','r35369216','r37265152','r37265408','r30598912','r37277440','r32187904','r32176640','r32092160','r32433664','r32238592']
;---------------------------------------------------
;need to generate a master file which has all the satstars triggers for
;the entire PAO in it.  Do this by opening each AOR, reading in all
;the satstar triggers, and outputting to the master list.
;----------------------------------------------------

;open file to contain master list of all the saturated stars
mastername = '/Users/jkrick/irac_warm/latent/flag_latent/int_aug2011/master_satcorr_phot.tbl'

 ;read in the formatted list of triggers
fmt = '(F10.2,F10.2,F10.5,F10.5,F10.2,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3,F10.2,F10.2,F10.2,F10.2,D26.8, I10)'
;readfmt, mastername, fmt, x    ,  y    ,    ra   ,    dec  ,replace,Pratio,Fratio, ferr , hsig , loc_sky, skysig, radius,    raw  ,    model  ,  rectify , triggersclk, naxis_trig,skipline = 2

readcol,mastername,x    ,  y    ,    ra   ,    dec  ,replace,Pratio,Fratio, ferr , hsig , loc_sky, skysig, radius,    raw  ,    model  ,  rectify , triggersclk, naxis_trig, format='(F10.2,F10.2,F10.5,F10.5,F10.2,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3,F10.2,F10.2,F10.2,F10.2,D26.8, I10)', skipline = 2

;replace is currently in units of arcseconds, change this to pixels
replace = replace / 1.2

;after the latent has gone below the background we need to stop
;triggering on it.
stopusing = intarr(n_elements(triggersclk))
stopusing(*) = 0

;because idl defines from 0,0 not 1,1
x = x - 1
y = y - 1

for a = 0,  n_elements(aorname) -1 do begin ;for each AOR in the entire PAO...

   aordir = strcompress(dirloc +aorname(a) +'/ch1/bcd/' ,/remove_all)
   CD,aordir;change directories

   print, 'working on AOR ', aordir

   ;command  =  "find . -name '*cbcd.fits' > /Users/jkrick/irac_warm/latent/flag_latent/int_aug2011/cbcdlist.txt"
   command  =  "ls *cbcd.fits > /Users/jkrick/irac_warm/latent/flag_latent/int_aug2011/cbcdlist.txt"
   spawn, command
   readcol,'/Users/jkrick/irac_warm/latent/flag_latent/int_aug2011/cbcdlist.txt',cbcdnames, format = 'A', /silent
 
   for f = 0, n_elements(cbcdnames) - 1 do begin  ; for each image in the directory
      ;print, 'working on ', cbcdnames(f)
      header = headfits(cbcdnames(f))  ;read in the header
      currentsclk =  sxpar(header, 'SCLK_OBS')  ;keep track of the time stamp
      naxis = sxpar(header, 'NAXIS')            ;figure out if subarray or fullarray
      exptime = sxpar(header, 'EXPTIME')  ;need to know exposure time for unit conversion later
      readmode = sxpar(header, 'READMODE') ; FULL or SUB array
      dceid = sxpar(header, 'DCEID')

      ; the stars which need to have their latent tested against the background
      lat = where(currentsclk gt triggersclk and stopusing lt stopthresh, count)
      ;print, 'number of triggers' , count, '-------------'

      for l = 0, count - 1 do begin  ; start this loop if there is at least one trigger
         print, 'x,y', x(lat(l)), y(lat(l))

         ;ignore subarray
         if readmode eq 'SUB' then begin
            ;print, 'subarray'
         endif

         ;if not a subarray...
         if readmode ne 'SUB' then begin

            fits_read, cbcdnames(f), data, header  ;read in the fits image

            ;how bright was the original source in electrons
            gain = sxpar(header, 'GAIN')
            pixel_scale = abs(sxpar(header, 'PXSCAL2'))  ;need the absolute value, sometimes they are negative.
            flux_conv = sxpar(header, 'FLUXCONV')
            known_flux_electron = mjy_to_electron(rectify(lat(l)),  pixel_scale, gain, exptime, flux_conv)
            
            ;use that to measure the y_intercept coefficient
            y_int = -5614.36 + 880.124*alog10(known_flux_electron)
            ;print, 'trigger, mjy, electron, y_int', rectify(lat(l)), known_flux_electron, y_int

                                ;measure background
            bkg = measure_bkg(x(lat(l)), y(lat(l)), data, sq2, gain, exptime, flux_conv)
            ;print, 'bkg', bkgjy 

            ;need time since the star was exposed.
            ltime = (currentsclk - triggersclk(lat(l))) / 3600

            ;this calculates what the predicted flux is of the latent image
            latent_flux = y_int*exp(-ltime/4.4)

            ;print, 'predicted flux', latent_flux, ltime

            if latent_flux gt 2*abs(bkg(1)) then begin  ;if the latent is above the background noise then we need to flag it
                                ;FLAG THE LATENT!
               printf, outlun, aorname(a), dceid ,x(lat(l)), y(lat(l)),  currentsclk ,  triggersclk(lat(l)) ,  10, format = '(A, I10, F10.2, F10.2,D26.8,D26.8, I10)'
                ;printf, outlun, cbcdnames(f), x(lat(l)), y(lat(l)), replace(lat(l)), currentsclk, triggersclk(lat(l)), format = '(A, F10.2, F10.2, F10.2,D26.8,D26.8)'
              ;for x1 =fix(x(lat(l)) - replace(lat(l))) , fix(x(lat(l)) + replace(lat(l))) do begin
               ;   for y1 = fix(y(lat(l)) - replace(lat(l))) , fix(y(lat(l)) + replace(lat(l))) do begin
               ;      print, 'flag ', cbcdnames(f), x1, y1
               ;   endfor
              ; endfor

               ;print, '----------------------------------------'
               stopusing(lat(l)) = 0 ;reset the 'stop using' counter
            endif else begin
               stopusing(lat(l)) = stopusing(lat(l)) + 1
            endelse
            
         endif

      endfor                    ; the number of triggers
      
   endfor  ; end for each image inside the AOR

endfor ; end for each AOR in the PAO

close, outlun
free_lun, outlun


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
  x = randomu(S,nrand) *(xbkg_radius * 2) + xcen - xbkg_radius
  y = randomu(S,nrand) *(ybkg_radius * 2) + ycen - ybkg_radius
 
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
     plothist, bkgd, xhist, yhist, bin = binsize, xrange = [-10000,10000] ,/noplot ;make a histogram of the data
;     print, 'plothist xhist', xhist
;     print, 'plothist yhist', yhist

     start = [mean(bkgd),stddev(bkgd), 2.]
    ; print, 'start', start
     noise = fltarr(n_elements(yhist))
     noise[*] = 1  
                                ;equally weight the values
     result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) ; fit a gaussian to the histogram sorted data
     bkg = result(0)
;     oplot, xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.)
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
