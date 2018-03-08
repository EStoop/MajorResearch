;+
; NAME:flag_latent
;
; PURPOSE: This code is the idl version that was used to make the flag_latent module in the post-cryo pipeline for ch1 only.  It works across AORs to track and flag all saturated sources causing a latent.  
;
; INPUTS:  Tested on a PAO with aor's as listed inside the code.
;
; OUTPUTS:logfile, master saturated list file which has all saturated sources for the whole PAO, 
;
; MODIFICATION HISTORY:
;Jessica Krick July 2011, initial version
;Jessica Krick Dec 2012, documentation added
;
;-
pro flag_latent
;-----------------------



logdate=string(strmid(systime(0),8,2),format='(i2.2)')+strmid(systime(0),4,3)+strmid(systime(0),22,2)
logfile = '~/irac_warm/latent/flag_latent/idllog_ch1-'+logdate
if file_test(logfile) then print,'Log file found. no logging'
if not file_test(logfile) then print,'Beginning log file...' & journal,logfile

sq2 = 2
stopthresh = 4  ;stop correcting when we find 4 in a row below background detection
;this allows for HDR frames to have 2 deep images in a row.

dirloc = '/Users/jkrick/irac_warm/latent/flag_latent/'

;time order of pc16
;'38927616_ch1','38927872_ch1','38928128_ch1','38928384_ch1',
aorname =['38926080_ch1','38926336_ch1','38928640_ch1','38926592_ch1','38928896_ch1','38926848_ch1','38929152_ch1','38927104_ch1','38929408_ch1','38927360_ch1','38933248_ch1','38920192_ch1','38925056_ch1','38925568_ch1','38925312_ch1','38925824_ch1','38924544_ch1','38924800_ch1','38923776_ch1','38924032_ch1','38924288_ch1','35418624_ch1','35422208_ch1','31537152_ch1','38763008_ch1','38788352_ch1','37908480_ch1'];,'38920448_ch1','37908992_ch1','37909504_ch1','37907200_ch1','38772992_ch1','38778112_ch1','38801920_ch1','31543040_ch1']

;---------------------------------------------------
;need to generate a master file which has all the satstars triggers for
;the entire PAO in it.  Do this by opening each AOR, reading in all
;the satstar triggers, and outputting to the master list.
;----------------------------------------------------

;open file to contain master list of all the saturated stars
mastername = dirloc + 'master_satcorr_phot.tbl'
openw,  masterlun, mastername,/get_lun
                                ;print a header for the log file
printf, masterlun, '|  x    |  y    |    ra   |    dec  |replace|Pratio|Fratio| ferr | hsig | loc sky| skysig| radius|    raw  |   model  |        rectify | sclk | naxis_trig|'
printf, masterlun, '|  -    |  -    |   deg   |    deg  | arcsec|  -   |  -   |  -   |  -   |  MJy/sr| MJy/sr| arcsec|    mJy  |   mJy    |        mJy     | s       | -  |'
ntot = 0
for a = 0,  n_elements(aorname) - 1 do begin
   CD,dirloc +aorname(a)   ; change directories to the correct AOR directory
   print, 'current dir ', dirloc + aorname(a)
   ;use unix to list all the filenames inside each AOR
   command1  =  "find . -name 'satcorr_phot*.tbl' > ./flag_latentlist.txt"
   spawn, command1
   command2  =  "find . -name 'cbcd*.fits' > ./cbcdlist.txt"
   spawn, command2

   ;read in all the individual satcorr_phot files, make them into one file with time stamps from sclk
   readcol,'./flag_latentlist.txt', latentlist, format = 'A', /silent
   readcol,'./cbcdlist.txt',cbcdlist, format = 'A', /silent
   
   for i =0, n_elements(latentlist) - 1 do begin  ;read each file, output it to one common file

      ngood = file_lines(latentlist(i))   ;how many lines are in each file

      if ngood gt 0 then begin  ;if the file is not empty

         readcol,latentlist(i), x    ,  y    ,    ra   ,    dec  ,replace,Pratio,Fratio, ferr , hsig , loc_sky, skysig, radius,    raw  ,   $
                 model  ,  rectify , skipline = 3,/silent

         ;read in each cbcd header, also output the sclk time and the naxis so we
         ;know if it is subarray or full array
         header = headfits(cbcdlist(i))
         sclk = sxpar(header, 'SCLK_OBS')
         naxis_trig = sxpar(header, 'NAXIS')

         for n = 0, ngood - 4 do begin  ;for each trigger in the individual file, output to master file
;            print, 'ngood, ntot',  ngood, '  ', ntot, '  ', latentlist(i)
            ntot = ntot + 1
            
             printf, masterlun, x(n)    ,  y(n)    ,    ra(n)   ,    dec(n)  ,replace(n),Pratio(n),Fratio(n), ferr (n), hsig(n) , $
                    loc_sky(n), skysig(n), radius(n),    raw(n) ,   model(n)  ,  rectify(n) , sclk, naxis_trig,$
                    format = '(F10.2,F10.2,F10.5,F10.5,F10.2,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3,F10.2,F10.2,F10.2,F10.2,D26.8, I10)'

         endfor                 ;end for each trigger
         
      endif                     ; end if the file is not empty
      
   endfor                       ; number of images satcorr_phot lists

endfor                          ;end for each AOR in the PAO

;clean up the open files
close, masterlun
free_lun, masterlun

;-------------------------------
;read in all the triggers, 
;go through each image in each AOR and flag any latent which is above
;the background in that image.
;-------------------------------

;print, 'finished making list, now moving onto reading the triggers'

 ;read in the formatted list of triggers
fmt = '(F10.2,F10.2,F10.5,F10.5,F10.2,F10.3,F10.3,F10.3,F10.3,F10.3,F10.3,F10.2,F10.2,F10.2,F10.2,D26.8, I10)'
readfmt, mastername, fmt, x    ,  y    ,    ra   ,    dec  ,replace,Pratio,Fratio, ferr , hsig , loc_sky, skysig, radius,    raw  ,   $
         model  ,  rectify , triggersclk, naxis_trig,/silent

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

for a = 0, n_elements(aorname) -1 do begin ;for each AOR in the entire PAO...

   CD,dirloc +aorname(a)  ;change directories

   print, 'working on AOR', aorname(a)

   readcol, './cbcdlist.txt', cbcdnames, format = 'A'

   for f = 0, n_elements(cbcdnames) - 1 do begin  ; for each image in the directory
      print, 'working on ', cbcdnames(f)
      header = headfits(cbcdnames(f))  ;read in the header
      currentsclk =  sxpar(header, 'SCLK_OBS')  ;keep track of the time stamp
      naxis = sxpar(header, 'NAXIS')            ;figure out if subarray or fullarray
      exptime = sxpar(header, 'EXPTIME')  ;need to know exposure time for unit conversion later
      readmode = sxpar(header, 'READMODE') ; FULL or SUB array

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

            fits_read, cbcdnames(f), data, header  ;read in the fits image

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
            print, 'bkg', bkgjy 

            ;need time since the star was exposed.
            ltime = (currentsclk - triggersclk(lat(l))) / 3600

            ;this calculates what the predicted flux is of the latent image
            latent_flux = y_int*exp(-ltime/4.4)

            print, 'predicted flux', latent_flux, ltime

            if latent_flux gt 2*abs(bkg(1)) then begin  ;if the latent is above the background noise then we need to flag it
                                ;FLAG THE LATENT!
               print, cbcdnames(f), x(lat(l)), y(lat(l)), replace(lat(l)), currentsclk, triggersclk(lat(l)), format = '(A, F10.2, F10.2, F10.2,D26.8,D26.8)'
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

endfor ; end for each AOR in the PAO


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
case of

endcase
electrons = Mjypersr * gain*exptime/flux_conv


return, electrons

end
