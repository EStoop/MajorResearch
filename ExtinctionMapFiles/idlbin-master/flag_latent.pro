pro flag_latent


;This code determines the duration of latent flagging for non-HDR
;frames during the warm mission.  For channel 1 the flagging is based on an
;empiraclly determined functional form of the latent decay over time.
;Latent images in Ch1 last on the order of hours as they slowly decay
;below the detection threshold.  As such, it is necessary for this
;code to track latents across AORs, and to carefully measure the
;background levels in the images in order to know when the latent
;would no longer be detectable on an individual frame.  The code is
;triggered to start when a saturated star is detected.  The code works
;by first measuring the aperture flux on the saturation corrected image.  With this peak flux measurement, we are then able to
;predict future flux levels of the latent based on an exponential
;function.  To do this, the code systematically reads in every fits
;file in all subsequent AORs, measures the background level, and
;determines if the predicted latent flux is above our defined
;detection threshold of 2 sigma above the background noise.  The flagging is
;stopped when the predicted
;latent flux drops below the detection threshold for three consecutive
;images.  Latent flagging for channel 2 is much simpler.  The latents
;last only a few minutes and therefore we do not need to cross AORs to
;flag them.  There is a simple table of peak flux vs. latent duration.


;want to trigger this code to run if a saturated star has been
;detected.

;--------------------------
;example 1
  dirloc = '/Users/jkrick/irac_warm/pc9/IRAC023000/bcd/'
;list of a bunch of AOR's after the AOR with the bright star
;I don't know how to get this information automatically
  aorname = ['0038163712', '0038163968','0038163200', '0038163456', '0038162688', '0038162944', '0038164224', '0038164480','0038166272','0038166528','0038165248','0038165504','0038164736','0038164992','0038165760','0038166016']
  ;the saturated star is in frame 2
  startframe = 2
  xcen =88.04                ;,162.0,42.,201.0,120.0,59.,218.]
  ycen = 133.69              ;,115.48,116.,143.0,202.,160.,188.]
;--------------------------


;--------------------------
;example 2
  dirloc = '/Users/jkrick/irac_warm/latent/'
;list of a bunch of AOR's after the AOR with the bright star
;I don't know how to get this information automatically
  aorname = ['r38244352', 'r30274048','r30275840','r30276864','r34955008']
  ;the saturated star is in frame 2
  startframe = 51
  xcen = 7.8
  ycen = 81.9
;--------------------------

;--------------------------
;example 3
  dirloc = '/Users/jkrick/irac_warm/latent/'
;list of a bunch of AOR's after the AOR with the bright star
;I don't know how to get this information automatically
  aorname = ['r38217216', 'r38249984','r38268928', 'r30290432', 'r30288640']
  ;the saturated star is in frame 2
  startframe = 173
  xcen = 48.8
  ycen = 6.8
;--------------------------


;--------------------------
;example 1 with different stars
  dirloc = '/Users/jkrick/irac_warm/pc9/IRAC023000/bcd/'
;list of a bunch of AOR's after the AOR with the bright star
;I don't know how to get this information automatically
 aorname = ['0038163200', '0038163456', '0038162688', '0038162944', '0038164224', '0038164480','0038166272','0038166528','0038165248','0038165504','0038164736','0038164992','0038165760','0038166016']
  ;the saturated star is in frame 2
  startframe = 1
  xcen =201.0;,120.0,59.,218.]
  ycen = 143.0;,202.,160.,188.]
;--------------------------


;because idl defines from 0,0 not 1,1
  xcen = xcen - 1
  ycen = ycen - 1

;size of the photometry radius  
  sq2 =2    

;size of the flagging area
;want to consider changing the way this is done
flag_rad = 3

  if ch eq 1 then begin  ;this works only for channel 1
;empirically determined slope from saturated star tests
     slope = 4.4

;how do I know to stop testing for the latent
     stoplooking = 0
     stopthresh = 3             ;want to find 3 in a row below detection
     ;determined empirically from saturated star tests
;  ;before linearity correction
;  known_flux = [5.381E7, 1.2063E8, 3.569E8, 1.0463E9, 2.4542E9]
;  y_int = [570, 1031, 1311, 1662, 2042]
     known_flux = [5.16811E7, 1.2228E8, 3.6333E8, 1.0292E9, 2.5377E9]
     y_int = [570, 830, 1311, 1701, 2027]

;;;;;;;;;;;;;;;
;do photometry on the bright star
;;;;;;;;;;;;;;;
;change to the correct directory where the data are  
     CD,dirloc+aorname(0)                       ;+ '/ch1/bcd'
     command1 =  ' ls IRAC.1*bcd*.fits >  bcd.txt' ;list all images in that directory
     spawn, command1
  
  ;read in the list of images
     readcol, dirloc+aorname(0)+ '/bcd.txt', bcdlist, format='A', /silent ;/ch1/bcd
                                ;just read in the image with the saturated star on it
     print, 'startframe', bcdlist(startframe+1)
     fits_read, bcdlist(startframe+1), data, imheader

  ;do quick photometry on the saturated star
     flux = quick_phot(xcen, ycen, data, sq2,imheader)
     print, 'initial flux', flux
  ;set the initial time when the sat star was on the frame
     time0 = sxpar(imheader, 'SCLK_OBS')

  ;is this subarray or full array
     naxis = sxpar(imheader, 'NAXIS')
     if (naxis eq 2) then begin
        fullarray = 1
        subarray = 0
     endif
     if (naxis eq 3) then begin
        subarray = 1
        fullarray = 0
     endif
     
     if flux gt 5.16E7 then begin ;dimmer objects don't leave enough of a latent to be too concerned
        flag = 1                  ;set the flag!
        stoplooking = 0

;;;;;;;;;;;;;;;
;check remaining frames to really see when the latent goes below the
;noise
;this has to cross AORs
;;;;;;;;;;;;;;;

        for countframe = startframe+2, n_elements(bcdlist) - 1  do begin ;check the rest of the frames in this AOR
           fits_read, bcdlist(countframe), data, imheader
           
  ;do quick photometry on the saturated star
           flux = quick_phot(xcen, ycen, data, sq2,imheader)
           
    ;measure the background in the frame.
     ;;;;;;;;;
           bkg = measure_bkg( xcen, ycen, data, sq2)

 ;is the latent predicted flux greater than 2*noise in the background?
;ie is it a 2sigma detection or not

     ;determine time since the star was on the array
           time = sxpar(imheader, 'SCLK_OBS')
           time = time - time0  ;in units of seconds
           time = time / 3600   ; in unuts of hours
           
                ;generate the output file name from the input filename
           outfile = strcompress(strmid(bcdlist(countframe), 0,31)+ '_mymsk.fits', /remove_all)
                                ;make an image of the correct size with all zeros
           outdata = data
           z = where(outdata ne 0)
           outdata(z) = 0
        

           latent_flux = final_yint*exp(-time/4.4)
           if latent_flux gt 2*bkg(1) or stoplooking lt stopthresh then begin
              ;print, 'flag the pixels which were saturated'
              
           ;set the bits around the saturated star for a latent flag
           ;use the shape that the piepline already has in it
              outdata[xcen - flag_rad :xcen+flag_rad,ycen-flag_rad :ycen+flag_rad] = 1024
              
              flag = 1          ;flag is set
              stoplooking = 0
           endif else begin
              stoplooking  = stoplooking + 1
           endelse
           
           fits_write, outfile, outdata, imheader
           
           ;print, 'time, latent_flux, noise', time, latent_flux,bkg(1), stoplooking
        endfor
        

;;;;;;
;now check the AOR's after this one to see if the latent is
;still aparent.
;;;;;;;;;;;;;;;
        for aorcount =1,  n_elements(aorname) - 1 do begin
           ;print, 'working on aor', aorname(aorcount)
           if stoplooking lt stopthresh then begin ;if the latent is still bright at the end of the previous AOR, then check the next one.
              
;change to the correct directory where the data are  
              CD,dirloc+aorname(aorcount)           ;+ '/ch1/bcd'
              command1 =  ' ls IRAC.1*bcd*.fits >  bcd.txt' ;list all images in that directory
              spawn, command1
              
  ;read in the list of images
              readcol, dirloc+aorname(aorcount)+ '/bcd.txt', bcdlist, format='A', /silent ;'/ch1/bcd'
              for countframe = 1, n_elements(bcdlist) - 1  do begin                       ;ignore the first frame
                 fits_read, bcdlist(countframe), data, imheader
                                ;figure out if it is full array or
        ;subarray and make appropriate coordinate change.
        ; don't continue of the subarray/fullaraay doesn't match with the trigger star.
                 naxis = sxpar(imheader, 'NAXIS')
                 if fullarray eq 1 and naxis eq 2 then goodtogo = 1
                 if fullarray eq 1 and naxis eq 3 and 8 lt xcen and xcen lt 40 and 216 lt  ycen and ycen lt 247 then goodtogo = 1
                 if subarray eq 1 and naxis eq 3 then goodtogo = 1
                 if subarray eq 1 and naxis eq 2 then begin
                    xcen = xcen + 8 
                    ycen = ycen +216
                    goodtogo = 1
                 endif
                 
                ;generate the output file name from the input filename
                 outfile = strcompress(strmid(bcdlist(countframe), 0,31)+ '_mymsk.fits', /remove_all)
                                ;make an image of the correct size with all zeros
                 outdata = data
                 z = where(outdata ne 0)
                 outdata(z) = 0
                 
                 if goodtogo eq 1 and stoplooking lt stopthresh then begin
  ;do quick photometry on the saturated star
                    flux = quick_phot(xcen, ycen, data, sq2,imheader)

    ;measure the background in the frame.
     ;;;;;;;;;
                    bkg = measure_bkg( xcen, ycen, data, sq2)

 ;is the latent predicted flux greater than 2*noise in the background?
;ie is it a 2sigma detection or not

     ;determine time since the star was on the array
                    time = sxpar(imheader, 'SCLK_OBS')
                    time = time - time0 ;in units of seconds
                    time = time / 3600  ; in unuts of hours
                    
                    latent_flux = final_yint*exp(-time/4.4)
                    if latent_flux gt 2*bkg(1)  then begin
                       ;print, 'flag the pixels which were saturated'
                       
           ;set the bits around the saturated star for a latent flag
           ;use the shape that the piepline already has in it
                       outdata[xcen-flag_rad :xcen+flag_rad,ycen-flag_rad :ycen+flag_rad] = 1024

                       flag = 1 ;flag is set
                       stoplooking = 0
                    endif 
                    if latent_flux le 2*bkg(1) and stoplooking lt stopthresh then begin
                       outdata[xcen-flag_rad :xcen+flag_rad,ycen-flag_rad :ycen+flag_rad] = 1024
                       flag = 0
                       stoplooking = stoplooking + 1
                    endif
                    
                    fits_write, outfile, outdata, imheader
                    
                   ; print, 'time, latent_flux, noise', time, latent_flux,bkg(1), stoplooking
                 endif          ;if full array and subarray match
                 
              endfor
              
              
           endif
        endfor
        
     endif                      ;if star is bright enough to leave a significant latent
     
  endif   ;end ch1 

if ch eq 2 then begin

;;;;;;;;;;;;;;;
;do photometry on the bright star
;;;;;;;;;;;;;;;
;change to the correct directory where the data are  
   CD,dirloc+aorname(0)                            ;+ '/ch1/bcd'
   command1 =  ' ls IRAC.2*bcd*.fits >  bcd.txt'   ;list all images in that directory
   spawn, command1
  
  ;read in the list of images
   readcol, dirloc+aorname(0)+ '/bcd.txt', bcdlist, format='A', /silent ;/ch1/bcd
                                ;just read in the image with the saturated star on it
   print, 'startframe', bcdlist(startframe+1)
   fits_read, bcdlist(startframe+1), data, imheader

  ;do quick photometry on the saturated star
   flux = quick_phot(xcen, ycen, data, sq2,imheader)
   print, 'initial flux', flux
  ;set the initial time when the sat star was on the frame
   time0 = sxpar(imheader, 'SCLK_OBS')


   if flux le 1.96E7 then begin
      print, 'flag for 220 seconds'
   endif
   if flux gt 1.96E7 and flux le 1.317E8 then begin
      print, 'flag for 320 seconds'
   endif
   if flux gt 1.317E8 and flux le 3.74E8 then begin
      print, 'flag for 360 seconds'
   endif
   if flux gt 3.74E8  then begin
      print, 'flag for 440 seconds'
   endif

endif   ;end ch2

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



