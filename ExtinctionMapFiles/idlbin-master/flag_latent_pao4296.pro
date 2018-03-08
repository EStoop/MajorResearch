pro flag_latent

logdate=string(strmid(systime(0),8,2),format='(i2.2)')+strmid(systime(0),4,3)+strmid(systime(0),22,2)
logfile = '~/irac_warm/latent/flag_latent/idllog_ch1-'+logdate
if file_test(logfile) then print,'Log file found. no logging'
if not file_test(logfile) then print,'Beginning log file...' & journal,logfile

sq2 = 2
stopthresh = 4  ;stop correcting when we find 4 in a row below background detection
;this allows for HDR frames to have 2 deep images in a row.

dirloc = '/Users/jkrick/iracdata/nonstop_pao/'

;open file to contain master list of all the saturated stars
mastername = '/Users/jkrick/irac_warm/latent/flag_latent/pao4296_master_satcorr_phot.tbl'
print, mastername
 ;read in the formatted list of triggers
;fmt = '(F10.2,F10.2,F10.5,F10.5,F10.2,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3,F10.2,F10.2,F10.2,F10.2,D26.8, I10)'
fmt = '(17F0)'

readfmt, mastername, fmt, x    ,  y    ,    ra   ,    dec  ,replace,Pratio,Fratio, ferr , hsig , loc_sky, skysig, radius,    raw  ,   $
         model  ,  rectify , triggersclk, naxis_trig;,/silent

;get rid of repeat triggers (happens in HDR mode)
re_trig = fltarr(n_elements(x))
re_trig[*] = 0
for rec = 1, n_elements(x) - 1 do begin
                                ;if x and y of the trigger are the
                                ;same as the previous trigger within
                                ;0.25 pixels, then we don't
                                ;need to include that trigger
   if x(rec) ge x(rec-1) -0.25 and x(rec) le x(rec-1) +  0.25 and y(rec) ge y(rec-1) -0.25 and y(rec) le y(rec-1) + 0.25 then begin
      re_trig(rec) = 1
      print, 'repeat trigger x(rec),y(rec), x(rec-1), y(rec-1)', x(rec), y(rec), x(rec-1), y(rec-1)
   endif

endfor

ab = where(re_trig lt 1, non)
print, 'number of non-repeat triggers', non

;remake the arrays to disinclude the repeat triggers
x = x(ab) & y = y(ab) & ra = ra(ab) & dec = dec(ab) & replace = replace(ab) & Pratio = Pratio(ab) & Fratio = Fratio(ab) & ferr = ferr(ab) & hsig = hsig(ab) & loc_sky = loc_sky(ab) & skysig = skysig(ab) & radius = radius(ab) & raw = raw(ab) & model = model(ab) & rectify = rectify(ab) & triggersclk = triggersclk(ab) & naxis_trig = naxis_trig(ab) 

;replace is currently in units of arcseconds, change this to pixels
replace = replace / 1.2

;after the latent has gone below the background we need to stop
;triggering on it.
stopusing = intarr(n_elements(triggersclk))
stopusing(*) = 0

;because idl defines from 0,0 not 1,1
x = x - 1
y = y - 1


readcol, '~/irac_warm/latent/flag_latent/cbcdlist.txt', cbcdnames, format = 'A'


;both channels are in this directory! get rid of ch2
cs = 0
ch1_cbcdnames = cbcdnames
for cf = 0, n_elements(cbcdnames) - 1 do begin
   header = headfits(cbcdnames(cf)) ;read in the header
   ch = sxpar(header, 'CHNLNUM')
   if ch lt 2 then begin
      ;print, 'channel', ch, sxpar(header, 'RAWFILE')
      ch1_cbcdnames(cs) = cbcdnames(cf)
      cs = cs + 1
   endif
endfor

ch1_cbcdnames = ch1_cbcdnames[0:cs-1]

;arghhl have to order the images in time order, that is not to be assumed with these file names.
sclkarr = dblarr(n_elements(ch1_cbcdnames))
for mf = 0, n_elements(ch1_cbcdnames) -1  do begin
   header = headfits(ch1_cbcdnames(mf)) ;read in the header
   sclkarr(mf)  = sxpar(header, 'SCLK_OBS')   
endfor

mfa= sort(sclkarr)
sort_cbcd = ch1_cbcdnames(mfa)

;for si = 0, n_elements(sort_cbcd) - 1 do begin
;   header = headfits(sort_cbcd(si)) ;read in the header
;   print, sxpar(header, 'RAWFILE')
;endfor

cd, dirloc
   for f = 0, n_elements(sort_cbcd) - 1 do begin  ; for each image in the directory
      header = headfits(sort_cbcd(f))  ;read in the header
      currentsclk =  sxpar(header, 'SCLK_OBS')  ;keep track of the time stamp
      naxis = sxpar(header, 'NAXIS')            ;figure out if subarray or fullarray
      exptime = sxpar(header, 'EXPTIME')  ;need to know exposure time for unit conversion later
      readmode = sxpar(header, 'READMODE') ; FULL or SUB array
      rawfile = sxpar(header, 'RAWFILE')
      print, 'working on ', sort_cbcd(f), rawfile

      ; the stars which need to have their latent tested against the background
      lat = where(currentsclk gt triggersclk and stopusing lt stopthresh, count)
      print, 'number of triggers' , count, '-------------'

      for l = 0, count - 1 do begin  ; start this loop if there is at least one trigger
         print, 'x,y', x(lat(l)), y(lat(l))


     ;make sure using the correct coords in case of subarray or full array
;         if naxis eq 2 and naxis_trig(lat(l)) eq 2 then goodtogo = 1
;         if naxis eq 3 and naxis_trig(lat(l)) eq 3 then goodtogo = 1
         
;         if naxis eq 3 and naxis_trig(lat(l)) eq 2 and 8 lt x(lat(l)) and x(lat(l)) lt 40 and 216 lt  y(lat(l)) and y(lat(l)) lt 247 then begin
;            x = x - 8
;            y = y - 216
;            goodtogo = 1
;         endif
         
;         if naxis eq 2 and naxis_trig(lat(l)) eq 3 then begin
;            x = x + 8 
;            y = y +216
;            goodtogo = 1
;         endif
         

         ;ignore subarray
         if readmode eq 'SUB' then begin
            print, 'subarray'
         endif

         ;if not a subarray...
         if readmode ne 'SUB' then begin

            fits_read, sort_cbcd(f), data, header  ;read in the fits image

            ;how bright was the original source in electrons
            gain = sxpar(header, 'GAIN')
            pixel_scale = abs(sxpar(header, 'PXSCAL2'))  ;need the absolute value, sometimes they are negative.
            flux_conv = sxpar(header, 'FLUXCONV')
            known_flux_electron = mjy_to_electron(rectify(lat(l)),  pixel_scale, gain, exptime, flux_conv)
            
            ;use that to measure the y_intercept coefficient
            y_int = -5614.36 + 880.124*alog10(known_flux_electron)
            print, 'trigger, mjy, electron, y_int', rectify(lat(l)), known_flux_electron, y_int

                                ;measure background
            bkg = measure_bkg(x(lat(l)), y(lat(l)), data, sq2, gain, exptime, flux_conv)
            print, 'bkg', bkg  ;jy 

            ;need time since the star was exposed.
            ltime = (currentsclk - triggersclk(lat(l))) / 3600

            ;this calculates what the predicted flux is of the latent image
            latent_flux = y_int*exp(-ltime/4.4)

            print, 'predicted flux', latent_flux, ltime

            if latent_flux gt 2*abs(bkg(1)) then begin  ;if the latent is above the background noise then we need to flag it
                                ;FLAG THE LATENT!
               print, sort_cbcd(f), x(lat(l)), y(lat(l)), replace(lat(l)), currentsclk, triggersclk(lat(l)), format = '(A, F10.2, F10.2, F10.2,D26.8,D26.8)'
               ;for x1 =fix(x(lat(l)) - replace(lat(l))) , fix(x(lat(l)) + replace(lat(l))) do begin
               ;   for y1 = fix(y(lat(l)) - replace(lat(l))) , fix(y(lat(l)) + replace(lat(l))) do begin
               ;      print, 'flag ', cbcdnames(f), x1, y1
               ;   endfor
              ; endfor

               print, '----------------------------------------'
               stopusing(lat(l)) = 0 ;reset the 'stop using' counter
            endif else begin
               stopusing(lat(l)) = stopusing(lat(l)) + 1
            endelse
            
         endif

      endfor                    ; the number of triggers
      
   endfor  ; end for each image inside the AOR



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
     result= MPFITFUN('mygauss',xhist,yhist, noise, start) ; fit a gaussian to the histogram sorted data
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
