pro barrel_shift_noise, dir_name;, aor_list;, aor_bright, aor_cal
;copied from iwic_210.pro

;example command
;barrel_shift_noise, '/Users/jkrick/irac_warm/barrel_shift'
print,'starting ', systime()

colorname = [ 'teal', 'teal', 'teal', 'firebrick', 'firebrick', 'firebrick', 'gray', 'gray', 'gray','black']

flux_conv = [ 0.1253,0.1469] 

;-------------------------------------------------------------
;-------------------------------------------------------------

;get/make a mask image of the dark field data. = /IRAC/iwic210/dark_mask.fits
;	 This is a SExtractor object image
;	 0's for all backgroud pixels, >0's for all object pixels
;	 is not perfect, but is a reasonabe compromise between covering too much area, and getting all object flux masked
;  use the same mask for ch1 and ch2?  definitely use it for all exptimes
fits_read,strcompress(dir_name + '/dark_mask.fits',/remove_all) , objectmaskdata, objectmaskheader ;'/Users/jkrick/nutella/IRAC/iwic210/dark_mask.fits'
;---------
exptime = ['0pt4','6','12','0pt6']
exptime_num  = [0.4,6.,12.,0.6]
binsize = [.1, .1, .001, .0001]

noiselimit = [6.3, 12.8] ;e-/s of the zodiacal light

;need to create lists with the appropriate file names by channel,
;raw/bcd, exptime, & position
;filenames, dir_name
; for each position (3 positions)
readcol, dir_name +'/aorlist.txt', aorname, format="A",/silent

sigma_sigma = fltarr(n_elements(aorname))

for pos = 0, n_elements(aorname) -1  do begin; 
print, 'working on aor', aorname(pos)
;filenames appropriate to each position

   for ch =1, 1 do begin; 1 do begin       ;for ch1 and ch2

      for exp = 2, 2 do begin; n_elements(exptime) - 1 do begin
         print, '-----------'
         print, 'working on noise. ch, exp ', ch+ 1, ' ', exptime(exp)


;Get appropriate raw & bcd data from the central directory.
;    make files with lists of raw images for input; ignore the first frame
;    grab headers from bcd's so that I have ra and dec info


;these should also have position names in the file names.
;and these filenames will probably need to be read in from the command
;line.  Or else a file is created by Bill's code in which are
;listed these filenames. or else I cut and past his code into this
;one.
         if pos lt n_elements(aorname) then begin
            if ch eq 0 then begin
               rawlist = strcompress(dir_name + '/raw/'+ aorname(pos)+ '/rawlist_ch1_' + exptime(exp)+'s.txt',/remove_all)
               bcdlist =  strcompress(dir_name + '/bcd/'+ aorname(pos)+ '/bcdlist_ch1_' + exptime(exp)+'s.txt',/remove_all)
            endif else begin
               rawlist = strcompress(dir_name + '/raw/'+ aorname(pos)+ '/rawlist_ch2_' + exptime(exp)+'s.txt',/remove_all)
               bcdlist =  strcompress(dir_name + '/bcd/'+ aorname(pos)+ '/bcdlist_ch2_' + exptime(exp)+'s.txt',/remove_all)
            endelse
         endif else begin
            print, 'working on combined'
            if ch eq 0 then begin
               rawlist = strcompress(dir_name + '/raw/rawlist_ch1_' + exptime(exp)+'s.txt',/remove_all)
               bcdlist =  strcompress(dir_name + '/bcd/bcdlist_ch1_' + exptime(exp)+'s.txt',/remove_all)
            endif else begin
               rawlist = strcompress(dir_name + '/raw/rawlist_ch2_' + exptime(exp)+'s.txt',/remove_all)
               bcdlist =  strcompress(dir_name + '/bcd/bcdlist_ch2_' + exptime(exp)+'s.txt',/remove_all)
            endelse
         endelse

;XXXif we want to not include the first frame, remove it from this list.
         readcol, rawlist, rawimage, format="A",/silent
         readcol, bcdlist, bcdimage, format="A",/silent
         
;         print, 'aorname', aorname(pos)
;        print, 'bcdimage', bcdimage
;---------
;setup pointer array now to be filled with data arrays and headers
         rawdataptr  = [Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256))]
         edataptr  = [Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256))] ;in e-/s

         rawheaderptr = [Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader)]

;need to make these have the correct number of images as read in
;from rawlist
         rawdataptr = rawdataptr[0:n_elements(rawimage)-1]
         rawheaderptr = rawheaderptr[0:n_elements(rawimage)-1]
         edataptr = edataptr[0:n_elements(rawimage)-1]
;---------
         medianarr = fltarr(n_elements(rawimage))
         med_4 = fltarr(4, n_elements(rawimage))

;apply the mask to the data

         for name = 0,  n_elements(rawimage) -1 do begin
 ;           print, 'working on name', name, bcdimage(name)
            if pos lt n_elements(aorname) then begin
               fits_read, dir_name+ '/raw/'+ aorname(pos)+'/'+rawimage(name), rawdata_int, rawheader ;'/Users/jkrick/nutella/IRAC/iwic210/'
               fits_read, dir_name+ '/bcd/'+ aorname(pos)+'/'+bcdimage(name), bcddata, bcdheader ;'/Users/jkrick/nutella/IRAC/iwic210/'
            endif else begin
               fits_read, dir_name+ '/raw/'+rawimage(name), rawdata_int, rawheader ;'/Users/jkrick/nutella/IRAC/iwic210/'
               fits_read, dir_name+ '/bcd/'+bcdimage(name), bcddata, bcdheader ;'/Users/jkrick/nutella/IRAC/iwic210/'
            endelse
               
;get the flux values converted from integers to DN

            barrel = fxpar(rawheader, 'A0741D00')
            fowlnum = fxpar(rawheader, 'A0614D00')
            pedsig = fxpar(rawheader, 'A0742D00')
            ichan = fxpar (rawheader, 'CHNLNUM')
            print, 'details',  ichan, barrel, fowlnum, pedsig
            dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
 ;convert from DN/s to MJy/sr

            rawdata = rawdata / exptime_num(exp) ; first get it into DN/s
            rawdata = rawdata * flux_conv(ch) ;now in Mjy/sr

 ;rawdata is flipped in y ;flip just to make test output look ok
            hreverse, rawdata, rawheader, 2,/silent

; select background pixels based on source mask 
; first change astrometry in header to give the mask image the same x and y of the raw images.
            hastrom, objectmaskdata, objectmaskheader, maskdata, maskheader, bcdheader, missing = 0


                        ;save the image for diagnostics
            newname= '/diag' + strmid(rawimage(name),2, 40)
            fits_write, dir_name+ '/raw/'+ aorname(pos)+newname, rawdata, rawheader


; then use a where statement to apply the mask to the rawdata
;don't apply it to the subarray data; there won't be objects in there anyway.
            if fowlnum eq 2 then begin
;               print, 'working on 0.1s frame, do not apply mask'
                                ;plot up histogram seperately of the
                                ;subarray data for each sub array
;            subarray, rawdata
               naxis1 = fxpar(rawheader, 'NAXIS1')
               naxis2 = fxpar(rawheader, 'NAXIS2')

            endif else begin

              ;but not all of the images are covered by the mask.
              ; try masking based on a SExtractor image of the frame itself
               
               ;run SExtractor on the frame XXX
               nanpix = where(finite(bcddata) eq 0)  ; can't have nan's for sextractor to work.
               bcddata[nanpix] = -1E5
               fits_write, dir_name + '/junk.fits', bcddata, bcdheader
               command = 'sex ' +  dir_name + '/junk.fits'+ ' -c '+ dir_name +'/default.sex'
  ;             print, command
               spawn, command
               fits_read, 'segmentation.fits', segdata, segheader
               ob = where(segdata gt 0,obcount)

  ;             print, 'obcount', obcount

               a = where(maskdata gt 0)  ;this is the dark field mask
               rawdata[a] = alog10(-1) ;masked pixels get NaN's
               rawdata[ob] = alog10(-1)
               rawdata[nanpix] = alog10(-1)
   
;select background pixels based on a good region of the frame
               hextract, rawdata, bcdheader, 13,248,13,248,/silent
               naxis1 = fxpar(bcdheader, 'NAXIS1')
               naxis2 = fxpar(bcdheader, 'NAXIS2')

            endelse


; want to make some data have the units of electrons/second to use in noisy pixel calcuation
            edata = rawdata

            vdet =  fxpar(rawheader, 'A0625E00')
            vrst =  fxpar(rawheader, 'A0624E00')
            bias = vdet - vrst  ; should be a number like 0.6
            gain = fxpar(bcdheader, 'GAIN')
            aorlabel = fxpar(bcdheader, 'AORLABEL')

            edata = (edata/ flux_conv(ch))*gain
            *edataptr[name] = edata


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
            
         ;convert back to a normal 2D image
            cnvrt_im_64_256_4_to_stdim, rawdata_4, rawdata

            medianarr(name) = median(rawdata) ;doesn't really need to be a median, mean will be super close because of what I just did.
            

;XXX ???will need to correct the frames for this median value because the
;noise properties depend on this value

            *rawdataptr[name] = rawdata
            *rawheaderptr[name] = rawheader
         endfor  ;for each raw image
;-----------------------------------------------
;measure the stddev within each pixel
         
         count = float(0)
         sigmaarr = fltarr(naxis1*naxis2)
         meanarr = fltarr(naxis1*naxis2)
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
               if t gt n_elements(rawdataptr) - 2 then begin  ;XXXXX 2
                  arr = arr[0:t-1]
                  sigma = stddev(arr)
            
        ;try sigma with outlier rejection
                  meanclip, arr, clipmean, sigma, clipsig = 4, subs = subarr
       ;           print, 'arr', arr
       ;           print, 'sigma', sigma
                  sigmaarr(count) = sigma
                  meanarr(count) = clipmean
                  count = count + 1
               endif 
            endfor
         endfor

         sigmaarr = sigmaarr[0:count - 1]
         meanarr = meanarr[0:count - 1]
        
         save, sigmaarr, filename = '/Users/jkrick/irac_warm/barrel_shift/sigmaarr.sav'
         save, meanarr, filename = '/Users/jkrick/irac_warm/barrel_shift/meanarr.sav'

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
         binsize = mean(sigmaarr) / 10.
         goodsig = where(sigmaarr ne 0.)
;         print, 'zeros', n_elements(sigmaarr), n_elements(goodsig)
         sigmaarr = sigmaarr[goodsig]
         plothist, sigmaarr, xhist, yhist, bin = binsize, /noplot ,/nan

  ;       print, 'xhist', xhist
  ;       print, 'yhist', yhist

        ; a = where(xhist lt skymod + 10*skysigma)
         ;print, 'a cutting', n_elements(a), n_elements(xhist)
         ;xhist = xhist(a)
         ;yhist = yhist(a)
         
         ;print, 'xhist', xhist
         ;print, 'yhist', yhist

;fit a gaussian to this distribution
;         start = [mean(sigmaarr),stddev(sigmaarr), 100.]
         start = [skymod,skysigma, 5000.]
;         print, 'start', start
;      noise = fltarr(n_elements(sigmaarr))
         noise = fltarr(n_elements(yhist)) 
         noise[*] = 1.0           ;equally weight the values
         low = where(yhist lt 100)
         noise[low] = 10.0
;      print, 'n(sigmaarr), n(noise), n(yhist)', n_elements(sigmaarr), n_elements(noise), n_elements(yhist)
         result= MPFITFUN('mygauss',xhist,yhist, noise, start);./quiet    

;make sure sigma value is positive
         result = abs(result)
 
;try calculating mode with Sean's code
 ;        mode_calc, sigmaarr, mode; ,/verbose
 ;        print, 'mode', mode

;try calculating mdoe with coyote code
         c_mode = coyote_mode(sigmaarr)

;try calculating median
         med = median(sigmaarr)
;         sigma_sigma(pos) = result(1)
         sigma_sigma(pos) = skymod; mode; result(0) ; actually use the mean noise
;now plot the histogram of stddevs of all unmasked pixels
       ;plothist, sigmaarr, xjunk, yjunk, bin = binsize, xtitle = 'Standard deviation', ytitle='Number', $
        ;         xrange=[-8*result(1) + result(0),result(0) + 8*result(1)],$
        ;         title = "channel "+ string(ch+1)+ "    exptime " + string(exptime(exp)) + aorlabel ,$
        ;         thick=3, xthick=3,ythick=3,charthick=3
         if pos eq 0 then  begin
            a0 = plot(xhist, yhist, '3-', color = colorname(pos), xtitle = 'Standard Deviation in Raw Counts', ytitle='Number',title = "channel "+ string(ch+1)+ "    exptime " + string(exptime(exp)), xrange =[0,0.2] )
         endif

         if pos gt 0 then begin
            a0.select
            a1 = plot(xhist, yhist, '3-', color = colorname(pos) ,/overplot)
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

   endfor                       ;for each channel

endfor                          ; for each position

st = text( 0.12, 3000, 'shift 0'+ strmid(string(mean(sigma_sigma[0:2])),0,11)+ strmid(string(stddev(sigma_sigma[0:2])),0,9), color = 'teal',/data)
st2 = text( 0.12, 2500, 'shift 1'+ strmid(string(mean(sigma_sigma[3:5])),0,11)+ strmid(string(stddev(sigma_sigma[3:5])),0,9), color = 'red',/data)
st3 = text( 0.12, 2000, 'shift 2'+ strmid(string(mean(sigma_sigma[6:8])),0,11)+ strmid(string(stddev(sigma_sigma[6:8])),0,10), color = 'gray',/data)

percentage_change_1 = (mean(sigma_sigma[3:5]) - mean(sigma_sigma[0:2])) / mean(sigma_sigma[0:2])
percentage_change_2 = (mean(sigma_sigma[6:8]) - mean(sigma_sigma[0:2])) / mean(sigma_sigma[0:2])
print, 'percentage change 1', percentage_change_1, percentage_change_2
st4 = text(0.12, 3500, 'percentage changes' + string(percentage_change_1*100.)+'%' + string(percentage_change_2*100.)+'%',/data)

print, sigma_sigma
print, 'mean', mean(sigma_sigma[0:2])
print, 'sigma', stddev(sigma_sigma[0:2])

print,'finished ', systime()
end

