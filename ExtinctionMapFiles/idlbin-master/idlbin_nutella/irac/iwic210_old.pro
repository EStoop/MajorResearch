pro iwic210

!P.multi = [0,1,2]
ps_open, filename='/Users/jkrick/IRAC/iwic210/noise.ps',/portrait,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)

;get/make a mask image of the dark field data. = /IRAC/iwic210/dark_mask.fits
;	 This is a SExtractor object image
;	 0's for all backgroud pixels, >0's for all object pixels
;	 is not perfect, but is a reasonabe compromise between covering too much area, and getting all object flux masked
;  use the same mask for ch1 and ch2?  definitely use it for all exptimes
fits_read, '/Users/jkrick/IRAC/iwic210/dark_mask.fits', objectmaskdata, objectmaskheader
;---------
exptime = [100,200]

;determine the correct conversion factor from the dark data
;only need to do this once per setpoint
flux_conv = dark_standards()

for ch = 0, 1 do begin   ;for ch1 and ch2

   for exp = 0, n_elements(exptime) - 1 do begin
;Get appropriate raw & bcd data from the central directory.
;    make files with lists of raw images for input; ignore the first frame
;    grab headers from bcd's so that I have ra and dec info
      if ch eq 0 then begin
         rawlist = strcompress('/Users/jkrick/IRAC/iwic210/rawlist_ch1_' + string(exptime(exp))+'s',/remove_all)
         bcdlist =  strcompress('/Users/jkrick/IRAC/iwic210/bcdlist_ch1_'+ string(exptime(exp))+'s',/remove_all)
      endif else begin
         rawlist = strcompress('/Users/jkrick/IRAC/iwic210/rawlist_ch2_'+ string(exptime(exp))+'s',/remove_all)
         bcdlist =  strcompress('/Users/jkrick/IRAC/iwic210/bcdlist_ch2_'+ string(exptime(exp))+'s',/remove_all)
      endelse

;if we want to not include the first frame, remove it from this list.
      readcol, rawlist, rawimage, format="A"
      readcol, bcdlist, bcdimage, format="A"
      
;---------
;setup pointer array now to be filled with data arrays and headers
      rawdataptr  = [Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256))]

      rawheaderptr = [Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader)]

;need to make these have the correct number of images as read in
;from rawlist
      rawdataptr = rawdataptr[0:n_elements(rawimage)-1]
      rawheaderptr = rawheaderptr[0:n_elements(rawimage)-1]
;---------
      medianarr = fltarr(n_elements(rawimage))
      med_4 = fltarr(4, n_elements(rawimage))

;apply the mask to the data
      for name = 0,  n_elements(rawimage) -1 do begin

         fits_read, '/Users/jkrick/IRAC/iwic210/'+rawimage(name), rawdata_int, rawheader
         fits_read, '/Users/jkrick/IRAC/iwic210/' +bcdimage(name), bcddata, bcdheader

;get the flux values converted from integers to DN

         barrel = fxpar(rawheader, 'A0741D00')
         fowlnum = fxpar(rawheader, 'A0614D00')
         pedsig = fxpar(rawheader, 'A0742D00')
         ichan = fxpar (rawheader, 'CHNLNUM')

         dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
   ;pattern is set to 0, I don't know what this is.


;convert from DN/s to MJy/sr

         rawdata = rawdata / exptime(exp) ; first get it into DN/s
         rawdata = rawdata * flux_conv(ch) ;now in Mjy/sr

;rawdata is flipped in y ;flip just to make test output look ok
         hreverse, rawdata, rawheader, 2,/silent

; select background pixels based on source mask 
; first give the mask image the same x and y of the raw images.
         hastrom, objectmaskdata, objectmaskheader, maskdata, maskheader, bcdheader, missing = 100

; then use a where statement to apply the mask to the rawdata
         a = where(maskdata gt 0)
         rawdata[a] = alog10(-1)

;select background pixels based on a good region of the frame
         hextract, rawdata, bcdheader, 13,248,13,248

;measure the median of the background in each image
;median should ignore the NaN values for masked objects
;want to do this per readout 

         ;divide into 4 readouts
         cnvrt_stdim_to_64_256_4, rawdata, rawdata_4

         for read = 0, 3 do begin
            median_read = median(rawdata_4[*,*,read])
            rawdata_4[*,*,read] = rawdata_4[*,*,read] - median_read
            med_4(read,name) = median_read
         endfor


         
         ;mean of all 4 readout medians
;         mean_4 = mean(med_4)

         ;add that mean value back to the entire array
;         rawdata_4 = rawdata_4 + mean_4
                                ;don't add it back in, instead
                                ;can normalize all images to have zero
                                ;background flux which puts them on
                                ;equal footing so there is no ramp in
                                ;background values.

         ;convert back to a normal 2D image
         cnvrt_im_64_256_4_to_stdim, rawdata_4, rawdata

         medianarr(name) = median(rawdata)  ;doesn't really need to be a median, mean will be super close because of what I just did.

;XXX ???will need to correct the frames for this median value because the
;noise properties depend on this value

         *rawdataptr[name] = rawdata
         *rawheaderptr[name] = bcdheader
      endfor


;-----------------------------------------------
;plot the median of the background as a function of image number
;      plot, findgen(n_elements(medianarr)), medianarr, psym = 2, xtitle = 'frame number', ytitle = 'median of background regions',$
;            title = "channel "+string( ch+1)+ "    exptime " + string(exptime(exp)), thick=3, xthick=3,ythick=3,charthick=3

      plot, findgen(n_elements(rawimage)), med_4(0,*),  psym = 1, xtitle = 'frame number', ytitle = 'median of background regions',$
            title = "channel "+string( ch+1)+ "    exptime " + string(exptime(exp)), thick=3, xthick=3,ythick=3,charthick=3,$
            yrange = [mean(med_4(1,*)) - 0.15*(mean(med_4(1,*))), mean(med_4(1,*)) + 0.15*(mean(med_4(1,*)))]

      oplot, findgen(n_elements(rawimage)), med_4(1,*),  psym = 2, thick = 3
      oplot, findgen(n_elements(rawimage)), med_4(2,*),  psym = 4, thick = 3
      oplot, findgen(n_elements(rawimage)), med_4(3,*),  psym = 5, thick = 3

;-----------------------------------------------
;measure the stddev within each pixel
      naxis1 = fxpar(bcdheader, 'NAXIS1')
      naxis2 = fxpar(bcdheader, 'NAXIS2')

      count = float(0)
      sigmaarr = fltarr(naxis1*naxis2)
      for x = 0, naxis1-1 do begin
         for y = 0,  naxis2 - 1 do begin
            t = 0
            arr = fltarr(n_elements(rawdataptr))
            for i = 0, n_elements(rawdataptr) - 1 do begin
               if finite((*rawdataptr[i])[x,y]) gt 0 then begin
                  arr(t) = (*rawdataptr[i])[x,y] ;if not masked then include the value
                  t = t + 1
               endif

            endfor
      ;if there are enough measurements for a sigma
            if t gt 2 then begin
               arr = arr[0:t-1]
               sigma = stddev(arr)
            
        ;try sigma with outlier rejection
 ;        meanclip, arr, clipmean, sigma, clipsig = 4, subs = subarr
            
               sigmaarr(count) = sigma
               count = count + 1
            endif

         endfor

      endfor

      sigmaarr = sigmaarr[0:count - 1]

;plot histogram of stddevs of all unmasked pixels
;don't plot it, just determine the values
      plothist, sigmaarr, xhist, yhist, bin = 0.1, /noplot ; xrange=[-50,200]


;fit a gaussian to this distribution
      start = [mean(sigmaarr),stddev(sigmaarr), 300.]
      noise = fltarr(n_elements(sigmaarr))
      noise[*] = 1              ;equally weight the values

      result= MPFITFUN('gauss',xhist,yhist, noise, start,/quiet)    

;now plot the histogram of stddevs of all unmasked pixels
      plothist, sigmaarr, xjunk, yjunk, bin = 0.1, xtitle = 'Standard deviation', ytitle='Number', $
                xrange=[-8*result(1) + result(0),result(0) + 8*result(1)],$
                title = "channel "+ string(ch+1)+ "    exptime " + string(exptime(exp)),$
                thick=3, xthick=3,ythick=3,charthick=3

      xyouts, result(0) + 3.5*result(1), max(yhist) / 2.,  'mean' + string(result(0)), charthick=3
      xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 0.08*max(yhist),  'sigma' + string(result(1)), charthick=3
      xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 4.*0.08*max(yhist), 'N pixels used' + string(count), charthick=3

;plot the fitted gaussian and print results to plot
      oplot, findgen(1000)/ 10., (result(2))/sqrt(2.*!Pi) * exp(-0.5*((findgen(1000)/10. - (result(0)))/(result(1)))^2.), $
             color = greencolor, thick=3
      print, 'gauss fit result', result


;try using mmm to reject outliers, and use that gaussian for mean etc.
      mmm, sigmaarr, skymod, skysigma, skew
;oplot, findgen(500), (result(2))/sqrt(2.*!Pi) * exp(-0.5*((findgen(500) - (skymod))/(skysigma))^2.), color = redcolor
      print, 'mmm result', skymod, skysigma


;find the median
      print, 'median', median(sigmaarr)
      xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 2.*0.08*max(yhist), 'median' + string(median(sigmaarr)), charthick=3

;measure Q =1/<sigma(i)^2>(i)
      sigsquared = sigmaarr^2
      q = 1 / sigsquared
      q = mean(q)
      q = 1 / sqrt(q)
      print, 'Q ', q
      xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 3.*0.08*max(yhist), 'Q' + string(q), charthick=3

;----------------------------------------------------------
;how many hot pixels are there?
      openw, outlunred, '/Users/jkrick/IRAC/iwic210/highsigma.reg', /get_lun
;printf, outlunred, 'fk5'
      hotcount = fix(0)
      for x = 0, naxis1 - 1 do begin
         for y = 0, naxis2 -2 do begin
            hotpixel = 0
            fin = 0
            for i = 0, n_elements(rawdataptr) - 1 do begin

         ;need result of gaussian fit here so I can't do this in the first loop
               if finite((*rawdataptr[i])[x,y]) gt 0 then begin
                  fin = fin + 1
                  if  (*rawdataptr[i])[x,y] gt 2.0*result(0) + medianarr(i) then begin
                     hotpixel = hotpixel + 1
                  endif
               endif
  
            endfor
      
     ;if that same pixel is hot on all of the images where finite then it must be
                                ;a hot pixel
            if fin gt 3 and hotpixel gt fin - 1 then begin
               hotcount = hotcount + 1
               printf, outlunred, 'circle( ', x, y, ' 4")'
            endif

         endfor
      endfor
      close, outlunred
      free_lun, outlunred

      print, 'number of hot pixels', hotcount
      xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 5.*0.08*max(yhist), 'hot pixels '+ string(hotcount), charthick=3

   endfor                       ;for each exposure time

endfor        ;for each channel

ps_close, /noprint,/noid

end

;-----------------------------------------------
;
;     make a combination of all three positions (18 values per pixel)
;   	  repeat the histogram and measurement values
;

;convert DN to electrons so that all numbers for comparison are in
;electrons
;need to know the gain which changes with Temp & voltage which we are
;changing

;to determine the real gain, use a calibration star taken at the end
;of the AOR.
;do aperture photometry
 
