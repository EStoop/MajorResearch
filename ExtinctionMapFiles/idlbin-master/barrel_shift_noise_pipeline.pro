pro barrel_shift_noise_pipeline, sub, dir_name;, aor_list;, aor_bright, aor_cal
;copied from iwic_210.pro

;example command
;barrel_shift_noise_pipeline, 'sub', '/Users/jkrick/irac_warm/barrel_shift/pipeline_test/'
print,'starting ', systime()

colorname = [ 'black', 'blue',  'teal']

flux_conv = [ 0.1253,0.1469] 

;-------------------------------------------------------------
;-------------------------------------------------------------

;get/make a mask image of the dark field data. = /IRAC/iwic210/dark_mask.fits
fits_read,strcompress(dir_name + '/dark_mask.fits',/remove_all) , objectmaskdata, objectmaskheader

;all full arrays
exptime = ['12','30','f0p4', 'f2', 'f6','h12', 'h0p6_12','h30', 'h1p2','h6','h0p6_6']
;exptime_num = [12, 30, 0.4, 2, 6,12,0.6,30,1.2,6.0,0.6]
xper = [0.12,0.12,2.5,0.35,0.2,0.12, 2.5,0.12, 0.7,0.14, 2.5]
yval = [3000,2500, 2000]
xval = 0.09
;all sub arrays
if sub eq 'sub' then begin
   exptime = ['s0p02', 's0p1','s0p4']
;   exptime_num = [0.02, 0.1, 0.2, 0.4]
   xper = [100,15,8]
   yval = [180, 160,140]
   xval = 1.
endif

shiftname =[ 'shifted_0bit/', 'shifted_1bit/', 'shifted_2bit/']

for exp =0, n_elements(exptime) - 1 do begin

;for full arrays
   if exp eq 0 then xset =  [0,0.15]
   if exp eq 1 then xset =  [0,0.15]
   if exp eq 2 then xset = [0.0,10.0]
   if exp eq 3 then xset = [0.0,1.0]
   if exp eq 4 then xset = [0.0,0.5]
   if exp eq 5 then xset =[0,0.15]
   if exp eq 6 then xset =  [0.0,10.0]
   if exp eq 7 then xset =  [0.0,0.15]
   if exp eq 8 then xset =  [0.0,2.0]
   if exp eq 9 then xset =  [0.0,0.2]
   if exp eq 10 then xset =  [0.0,10.0]
   ;for sub arrays
if sub eq 'sub' then begin
   if exp eq 0 then xset = [40,200]
   if exp eq 1 then xset = [0,20]
   if exp eq 2 then xset = [0,5]
;   if exp eq 3 then xset = [0.0,10.0]
endif
   
   sigma_sigma = fltarr(n_elements(shiftname))

   for ch = 0, 0 do begin       ; 1 do begin       ;for ch1 and ch2

      for shift = 0, n_elements(shiftname) -1  do begin ; 
         print, 'working on ', shiftname(shift)
         print, '-----------'
         print, 'working on noise. ch, exp, shift', ch+ 1, ' ', exptime(exp), ' ', shift

 
         if ch eq 0 then begin
            bcdlist =  strcompress(dir_name +shiftname(shift) + '/bcdlist_ch1_' + exptime(exp)+'s.txt',/remove_all)
         endif else begin
            bcdlist =  strcompress(dir_name +shiftname(shift) + '/bcdlist_ch2_' + exptime(exp)+'s.txt',/remove_all)
         endelse

;XXXif we want to not include the first frame, remove it from this list.
;         readcol, rawlist, rawimage, format="A",/silent
         readcol, bcdlist, bcdimage, format="A",/silent
         print, 'nelements bcdimage', n_elements(bcdimage)
;---------
;setup pointer array now to be filled with data arrays and headers
         bcddataptr  = [Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256))]

         bcdheaderptr = [Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader)]

;need to make these have the correct number of images as read in
;from rawlist
         bcddataptr = bcddataptr[0:n_elements(bcdimage)-1]
         bcdheaderptr = bcdheaderptr[0:n_elements(bcdimage)-1]
;         edataptr = edataptr[0:n_elements(rawimage)-1]
;---------
         medianarr = fltarr(n_elements(bcdimage))
         med_4 = fltarr(4, n_elements(bcdimage))


         for name = 0,  n_elements(bcdimage) -1 do begin
            ;print, 'working on name', name, bcdimage(name)
            fits_read, dir_name+shiftname(shift) +bcdimage(name), bcddata, bcdheader 
            fowlnum = fxpar(bcdheader, 'AFOWLNUM')

            nanpix = where(finite(bcddata) eq 0) ; can't have nan's for sextractor to work.
            bcddata[nanpix] = -1E5
               
; select background pixels based on source mask 
; first change astrometry in header to give the mask image the same x and y of the raw images.
            hastrom, objectmaskdata, objectmaskheader, maskdata, maskheader, bcdheader, missing = 0, errmsg = herr

            naxis1 = fxpar(bcdheader, 'NAXIS1')
            naxis2 = fxpar(bcdheader, 'NAXIS2')
            print, 'naxis', naxis1, naxis2
            ;deal with subarray differently
            if naxis1 lt 100 then begin
               print, 'working on subarray'
;            subarray, rawdata
;               naxis1 = fxpar(bcdheader, 'NAXIS1')
;               naxis2 = fxpar(bcdheader, 'NAXIS2')
               
               ;make a median of all the frames?
 ;              slice = bcddata[*, *, i]
            endif else begin

              ;but not all of the images are covered by the mask.
              ; try masking based on a SExtractor image of the frame itself
               
               ;run SExtractor on the frame XXX
               fits_write, dir_name + '/junk.fits', bcddata, bcdheader
               command = 'sex ' +  dir_name + '/junk.fits'+ ' -c '+ dir_name +'/default.sex'
               spawn, command
               fits_read, 'segmentation.fits', segdata, segheader
               ob = where(segdata gt 0,obcount)

               print, 'just ran sextractor: obcount', obcount

               a = where(maskdata gt 0)  ;this is the dark field mask
               bcddata[a] = alog10(-1) ;masked pixels get NaN's
               bcddata[ob] = alog10(-1)
               bcddata[nanpix] = alog10(-1)
   
;select background pixels based on a good region of the frame
               hextract,bcddata, bcdheader, 13,248,13,248,/silent

            endelse

;measure the median of the background in each image
;median should ignore the NaN values for masked objects
;want to do this per readout 
           
         ;divide into 4 readouts
;            cnvrt_stdim_to_64_256_4, bcddata, bcddata_4
            
;            for read = 0, 3 do begin
;               median_read = median(bcddata_4[*,*,read])
;               bcddata_4[*,*,read] = bcddata_4[*,*,read] - median_read
;               med_4(read,name) = median_read
;            endfor            
         ;convert back to a normal 2D image
;            cnvrt_im_64_256_4_to_stdim, bcddata_4, bcddata

            medianarr(name) = median(bcddata) ;doesn't really need to be a median, mean will be super close because of what I just did.

;XXX ???will need to correct the frames for this median value because the
;noise properties depend on this value

            *bcddataptr[name] = bcddata
            *bcdheaderptr[name] = bcdheader
         endfor  ;for each raw image
;-----------------------------------------------
;measure the stddev within each pixel
         bcddataptr = bcddataptr[0:n_elements(bcdimage)-1]
         bcdheaderptr = bcdheaderptr[0:n_elements(bcdimage)-1]

         count = float(0)
         sigmaarr = fltarr(naxis1*naxis2)
         meanarr = fltarr(naxis1*naxis2)

         if naxis1 gt 100 then begin ; dealing with full array
            for x = 0, 235-1 do begin
               for y = 0,  235 - 1 do begin
                  t = 0
                  arr = fltarr(n_elements(bcddataptr))
                  for i = 0, n_elements(bcddataptr) - 1 do begin
                                ;              print, i, x, y, (*bcddataptr[i])[x,y]
                     if finite((*bcddataptr[i])[x,y]) gt 0 then begin
                        arr(t) = (*bcddataptr[i])[x,y] ;if not masked then include the value
                        t = t + 1
                     endif
                  endfor
                                ;if there are enough measurements for a sigma
                  if t gt n_elements(rawdataptr) - 2 then begin ;XXXXX 2
                     arr = arr[0:t-1]
                     sigma = stddev(arr)
                     
                                ;try sigma with outlier rejection
                     meanclip, arr, clipmean, sigma, clipsig = 4, subs = subarr
;                  if ch eq 1 and sigma gt 0.1 then print,sigma,  'arr', arr
                                ;                print, 'sigma', sigma
                     sigmaarr(count) = sigma
                     meanarr(count) = clipmean
                     count = count + 1
                  endif 
               endfor
            endfor
         endif else begin  ;for sub array data
            for x = 0, naxis1-1 do begin
               for y = 0,  naxis2 - 1 do begin
                  t = 0
                  arr = fltarr(n_elements(bcddataptr)*64)
                  for i = 0, n_elements(bcddataptr) - 1 do begin
                     for j = 1, 64 - 1 do begin  ;ignore the first frame
                                ;if finite((*bcddataptr[i])[x,y,j]) gt 0 then begin
                                ; print,i, j,  (*bcddataptr[i])[x,y,j]
                        arr(t) = (*bcddataptr[i])[x,y,j] ;if not masked then include the value
                        t = t + 1
                                ;   endif
                     endfor     ; for all i
                  endfor        ; for all j
                  
                  if t gt n_elements(rawdataptr) - 2 then begin ;XXXXX 2
                     arr = arr[0:t-1]
                     sigma = stddev(arr)
                     
                                ;try sigma with outlier rejection
                     meanclip, arr, clipmean, sigma, clipsig = 4, subs = subarr
;                  if ch eq 1 and sigma gt 0.1 then print,sigma,  'arr', arr
                                ;                print, 'sigma', sigma
                     sigmaarr(count) = sigma
                     meanarr(count) = clipmean
                     count = count + 1
                  endif 
               endfor           ;for all y
            endfor              ;for all x 
         endelse
          
         sigmaarr = sigmaarr[0:count - 1]
         meanarr = meanarr[0:count - 1]

         save, sigmaarr, filename = '/Users/jkrick/irac_warm/barrel_shift/pipeline_test/sigmaarr.sav'
         save, meanarr, filename = '/Users/jkrick/irac_warm/barrel_shift/pipeline_test/meanarr.sav'

;try using mmm to reject outliers, and use that gaussian for mean etc.
         mmm, sigmaarr, skymod, skysigma, skew
;oplot, findgen(500), (result(2))/sqrt(2.*!Pi) * exp(-0.5*((findgen(500) - (skymod))/(skysigma))^2.), color = redcolor
         print, 'mmm result', skymod, skysigma
         ;if mmm barfs and won't fit the distribution, make it up.
         if skysigma lt 0 then skysigma = stddev(sigmaarr)
;
;plot histogram of stddevs of all unmasked pixels
;don't plot it, just determine the values
;change bin size based on exposure time
         binsize = skymod / 10. ;mean(sigmaarr,/nan) / 10.
         print, 'mean, bin', mean(sigmaarr,/nan), binsize
         goodsig = where(sigmaarr ne 0. )
         print, 'zeros', n_elements(sigmaarr), n_elements(goodsig)
         sigmaarr = sigmaarr[goodsig]
         meanarr = meanarr[goodsig]
         plothist, sigmaarr, xhist, yhist, bin = binsize, /noplot ,/nan
        ; plothist, meanarr, xhist, yhist, bin = binsize, /noplot ,/nan

      ;   for cb = 0, n_elements(xhist) - 1 do print, xhist(cb), yhist(cb)
        ; a = where(xhist lt skymod + 10*skysigma)
         ;print, 'a cutting', n_elements(a), n_elements(xhist)
         ;xhist = xhist(a)
         ;yhist = yhist(a)
         
   
;fit a gaussian to this distribution
;         start = [mean(sigmaarr),stddev(sigmaarr), 100.]
         start = [skymod,skysigma, 5000.]
;         print, 'start', start
;      noise = fltarr(n_elements(sigmaarr))
         noise = fltarr(n_elements(yhist)) 
         noise[*] = 1.0           ;equally weight the values
         low = where(yhist lt 100, lowcount)
         print, 'lowcount', lowcount
         if naxis1 gt 100 then noise[low] = 10.0
;      print, 'n(sigmaarr), n(noise), n(yhist)', n_elements(sigmaarr), n_elements(noise), n_elements(yhist)
         result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet);./quiet    

;make sure sigma value is positive
         result = abs(result)
 
;try calculating mode with Sean's code
 ;        mode_calc, sigmaarr, mode; ,/verbose
 ;        print, 'mode', mode

;try calculating mode with coyote code
         c_mode = coyote_mode(sigmaarr)

;try calculating median
         med = median(sigmaarr)
;         sigma_sigma(pos) = result(1)
         sigma_sigma(shift) = skymod; mode; result(0) ; actually use the mean noise
;now plot the histogram of stddevs of all unmasked pixels
       ;plothist, sigmaarr, xjunk, yjunk, bin = binsize, xtitle = 'Standard deviation', ytitle='Number', $
        ;         xrange=[-8*result(1) + result(0),result(0) + 8*result(1)],$
        ;         title = "channel "+ string(ch+1)+ "    exptime " + string(exptime(exp)) + aorlabel ,$
        ;         thick=3, xthick=3,ythick=3,charthick=3
         if shift eq 0 then  begin
            a0 = plot(xhist, yhist, '3-', color = colorname(shift), xtitle = 'Sigma (Mjy/sr)', ytitle='Number',title = "channel "+ string(ch+1)+ "    exptime " + string(exptime(exp)), xrange =xset)
         endif

         if shift gt 0 then begin
            a0.select
            a0 = plot(xhist, yhist, '3-', color = colorname(shift) ,/overplot)
         endif

;         line_mean = polyline([result(0), result(0)], [3000, 3800],color = 'red',/data)
;         line_t = text(3.0, 2500,'mean',color = 'red',/data )
;         line_mode = polyline([mode, mode], [3000, 3800],color = 'blue',/data)
;         line_t = text(3.0, 2000,'mode',color = 'blue',/data )
;         line_mode = polyline([skymod, skymod], [3000, 3800],color = 'green',/data)
;         line_t = text(3.0, 2200,'mmm_mode',color = 'green',/data )
         
;         xyouts, result(0) + 3.5*result(1), max(yhist) / 2.  + 0.08*max(yhist), 'mean' + string(result(0)), charthick=3
;         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. ,  'sigma' + string(result(1)), charthick=3
;         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. + 2*0.08*max(yhist), 'N pixels used' + string(count), charthick=3
         
 ;;        t2 = text(3.5, 1000. + pos*300. ,  string(result(0)), /current, /data, color = colorname(pos))

;         xarr = findgen(10000)/100.
;         t5 = plot( xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), /overplot, color = colorname(pos))
        
;find the median
         print, 'median', median(sigmaarr)
;         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 1.*0.08*max(yhist), 'median' + string(median(sigmaarr)), charthick=3
   ;      t4 = text(result(0) + 3.5*result(1), max(yhist) / 2. - 1.*0.08*max(yhist), 'median' + string(median(sigmaarr)), /current,/data)
;         line_med = polyline([median(sigmaarr), median(sigmaarr)], [3000, 3800],color = 'purple',/data)
;         line_t = text(3.0, 3000,'median',color = 'purple',/data )

;----------------------------------------------------------
;----------------------------------------------------------
;----------------------------------------------------------
;
      endfor                    ;for each exposure time
      a0 = text( xval, yval(0), 'shift 0'+ strmid(string(mean(sigma_sigma[0])),0,11)+ strmid(string(stddev(sigma_sigma[0])),0,9), color = 'black',/data)
      a0 = text( xval, yval(1), 'shift 1'+ strmid(string(mean(sigma_sigma[1])),0,11)+ strmid(string(stddev(sigma_sigma[1])),0,9), color = 'blue',/data)
      a0 = text( xval, yval(2), 'shift 2'+ strmid(string(mean(sigma_sigma[2])),0,11)+ strmid(string(stddev(sigma_sigma[2])),0,10), color = 'teal',/data)
      
      percentage_change_1 = (mean(sigma_sigma[1]) - mean(sigma_sigma[0])) / mean(sigma_sigma[0])
      percentage_change_2 = (mean(sigma_sigma[2]) - mean(sigma_sigma[0])) / mean(sigma_sigma[0])
;print, 'percentage change 1', percentage_change_1;, percentage_change_2
      a0 = text(xper(exp), yval(1),string(percentage_change_1*100.)+'%', color = 'blue', /data)
      
      a0 = text(xper(exp), yval(2), string(percentage_change_2*100.)+'%',color = 'teal',/data)

      a0.save, strcompress('/Users/jkrick/irac_warm/barrel_shift/pipeline_test/ch'+string(ch+1) + '_' + exptime(exp)+'s.png',/remove_all)

   endfor                       ;for each channel


endfor                          ; for each shift

;print, sigma_sigma
;print, 'mean', mean(sigma_sigma[0:2])
;print, 'sigma', stddev(sigma_sigma[0:2])

;---------------------------------
;make final plots
;frametime = ['s0.02','s0.1','f0.4','s0.4','h0.6_6','h0.6_12','h1.2','f2','h6','f12','h12','f30','h30']
;ch1_pipeline_shift1 = [0.7,1.1,0.1,1.3,0.2,-0.01,0.06,1.5,1.3,0.3,0.04,0.8,-0.2]
;ch1_pipeline_shift2 = [3.2,5.3,2.8,5.6,2.5,2.8,1.1,2,2.5,1.7,2.3,2.1,1.7]
;ch2_pipeline_shift1 = [1.4,1.9,1.3,2.5,-0.5,-0.7,2,2.4,2.1,-1.8,-1.8,1.1,1.1]
;ch2_pipeline_shift2 = [5.4,8.7,5.6,11.6,3.7,3.5,4.7,9.7,11.6,5,4.8,3.2,2.9]

;x = findgen(n_elements(frametime))

;e = plot(x, ch1_pipeline_shift1, '2s', color = orange, yrange = [-10, 15])

print,'finished ', systime()
end

