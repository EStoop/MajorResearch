pro iwic210, dir_name

;dir_name = '/Users/jkrick/iwic/iwic_recovery4/IRAC017500'
print,'starting ', systime()

!P.multi = [0,1,2]
outputfile = strcompress(dir_name + '/iwic210_100.ps',/remove_all) 
ps_open, filename=outputfile,/portrait,/color

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
graycolor = FSC_COLOR("gray", !D.Table_Size-9)
browncolor = FSC_COLOR("brown", !D.Table_Size-10)

colorname = [ graycolor, redcolor , bluecolor , greencolor , yellowcolor , cyancolor , orangecolor , purplecolor ,  browncolor ]


;-------------------------------------------------------------
;-------------------------------------------------------------
;plot out the noise properties at this setpoint
;-------------------------------------------------------------
;-------------------------------------------------------------

;XXX
filenames, dir_name


;exptime = ['0pt1','0pt4','2','6','12','100','200']
exptime = ['100']; ['12','100']
exptime_num  = [100.]
binsize = [ .001, .0001]
;binsize = [.0001, .0001, .0001, .0001]
biasarr_ch1 = ['B600R3500G3650','B450R3500G3650','B750R3500G3650']
biasarr_ch2 = ['B500R3500G3650','B450R3500G3650','B600R3500G3650']
biasarr_ch1 = [biasarr_ch1, biasarr_ch1, biasarr_ch1]
biasarr_ch2 = [biasarr_ch2, biasarr_ch2, biasarr_ch2]


;exptime = [100,200]
;determine the correct conversion factor from the dark data
;only need to do this once per setpoint

;XXX taking this out for now because stars aren't detected at 50k.
;flux_conv = dark_standards(dir_name)
;print, 'flux_conv', flux_conv
flux_conv = [ 0.1088,0.1388]  ;from headers


;need to create lists with the appropriate file names by channel,
;raw/bcd, exptime, & position
;filenames, dir_name
; for each position (3 positions)
readcol, dir_name + '/aor_list_100', aorname, format="A",/silent

sigmaarr_ch1 = fltarr(9)
meanarr_ch1 = fltarr(9)
temp_ch1 = fltarr(9)
bias_ch1 = strarr(9)
hot_ch1 = fltarr(9)
sigmaarr_ch2 = fltarr(9)
meanarr_ch2 = fltarr(9)
temp_ch2 = fltarr(9)
bias_ch2 = strarr(9)
hot_ch2 = fltarr(9)

;fits_read,strcompress(dir_name + '/dark_mask.fits',/remove_all) , objectmaskdata, objectmaskheader ;'/Users/jkrick/nutella/IRAC/iwic210/dark_mask.fits'

for aor = 0,   n_elements(aorname)  -1 do begin ; 
;print, 'working on aorition', aor
;filenames appropriate to each aorition

   for ch = 0 , 1 do begin      ;for ch1 and ch2

   
      for exp = 0, n_elements(exptime) - 1 do begin
         print, '-----------'
         print, 'working on noise. ch, exp ', ch+ 1, ' ', exptime(exp)
;Get appropriate raw & bcd data from the central directory.
;    make files with lists of raw images for input; ignore the first frame
;    grab headers from bcd's so that I have ra and dec info

;         if aor lt 3 then begin
         if ch eq 0 then begin
            rawlist = strcompress(dir_name + '/raw/'+ aorname(aor)+ '/rawlist_ch1_' + exptime(exp)+'s.txt',/remove_all)
            bcdlist =  strcompress(dir_name + '/bcd/'+ aorname(aor)+ '/bcdlist_ch1_' + exptime(exp)+'s.txt',/remove_all)
         endif else begin
            rawlist = strcompress(dir_name + '/raw/'+ aorname(aor)+ '/rawlist_ch2_' + exptime(exp)+'s.txt',/remove_all)
            bcdlist =  strcompress(dir_name + '/bcd/'+ aorname(aor)+ '/bcdlist_ch2_' + exptime(exp)+'s.txt',/remove_all)
         endelse

              
         readcol, rawlist, rawimage, format="A",/silent

         readcol, bcdlist, bcdimage, format="A",/silent


;---------
;setup pointer array now to be filled with data arrays and headers
         rawdataptr  = [Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)),Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256))]



         rawheaderptr = [Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader),Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader)]

;need to make these have the correct number of images as read in
;from rawlist
         rawdataptr = rawdataptr[0:n_elements(rawimage)-1]
         rawheaderptr = rawheaderptr[0:n_elements(rawimage)-1]

;this doesn't work
;         rawdataptr = REPLICATE(Ptr_New(findgen(256,256)), n_elements(rawimage))
;         rawheaderptr = REPLICATE(Ptr_New(objectmaskheader), n_elements(rawimage))

;---------
         medianarr = fltarr(n_elements(rawimage))
         med_4 = fltarr(4, n_elements(rawimage))

;apply the mask to the data
         for name = 0,  n_elements(rawimage) -1 do begin

            fits_read, dir_name+ '/raw/'+ aorname(aor)+'/'+rawimage(name), rawdata_int, rawheader ;'/Users/jkrick/nutella/IRAC/iwic210/'
            fits_read, dir_name+ '/bcd/'+ aorname(aor)+'/'+bcdimage(name), bcddata, bcdheader                                ;'/Users/jkrick/nutella/IRAC/iwic210/'

;get the flux values converted from integers to DN

            barrel = fxpar(rawheader, 'A0741D00')
            fowlnum = fxpar(rawheader, 'A0614D00')
            pedsig = fxpar(rawheader, 'A0742D00')
            ichan = fxpar (rawheader, 'CHNLNUM')
            waitper = fxpar(rawheader, 'A0615D00')
            AFPAT1B = fxpar(rawheader, 'A0653E00')
            AFPAT2B = fxpar(rawheader, 'A0655E00')
            AFPAT1E = fxpar(rawheader, 'A0662E00')
            AFPAT2E = fxpar(rawheader, 'A0664E00')
            sclk = fxpar(rawheader, 'A0601D00')
            aorid = fxpar(rawheader, 'A0553D00')
                  
            if ch eq 0 then temp_ch1(aor) = AFPAT2B
            if ch eq 1 then temp_ch2(aor) = AFPAT2B

            dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
   ;pattern is set to 0, I don't know what this is.


;convert from DN/s to MJy/sr

            rawdata = rawdata / exptime_num(exp)    ; first get it into DN/s
            rawdata = rawdata * flux_conv(ch)       ;now in Mjy/sr

 ;rawdata is flipped in y ;flip just to make test output look ok
            hreverse, rawdata, rawheader, 2,/silent


; select background pixels based on source mask 
; first change astrometry in header to give the mask image the same x and y of the raw images.
;            hastrom, objectmaskdata, objectmaskheader, maskdata, maskheader, bcdheader, missing = 100


; then use a where statement to apply the mask to the rawdata
;don't apply it to the subarray data; there won't be objects in there anyway.
;            if fowlnum eq 2 then begin
;               print, 'working on 0.1s frame, do not apply mask'
                                ;plot up histogram seperately of the
                                ;subarray data for each sub array
;            subarray, rawdata
            naxis1 = fxpar(rawheader, 'NAXIS1')
            naxis2 = fxpar(rawheader, 'NAXIS2')

;            endif else begin

;               a = where(maskdata gt 0)
;               rawdata[a] = alog10(-1) ;masked pixels get NaN's

;select background pixels based on a good region of the frame
               hextract, rawdata, bcdheader, 13,248,13,248,/silent
               naxis1 = fxpar(bcdheader, 'NAXIS1')
               naxis2 = fxpar(bcdheader, 'NAXIS2')

;            endelse


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
                 
            *rawdataptr[name] = rawdata
            *rawheaderptr[name] = rawheader
                  

         endfor                 ;for each raw image

;-----------------------------------------------
;plot the median of the background as a function of image number
;      plot, findgen(n_elements(medianarr)), medianarr, psym = 2, xtitle = 'frame number', ytitle = 'median of background regions',$
;            title = "channel "+string( ch+1)+ "    exptime " + string(exptime(exp)), thick=3, xthick=3,ythick=3,charthick=3

         plot, findgen(n_elements(rawimage)), med_4(0,*),  psym = 1, xtitle = 'frame number', $
               ytitle = 'median of background regions (MJy/sr)',$
               title = "channel "+string( ch+1)+"  AFPAT2B" + string(AFPAT2B),  thick=3, $
               xthick=3,ythick=3,charthick=3,$
               yrange = [mean(med_4(1,*)) - 0.10*(mean(med_4(1,*))), mean(med_4(1,*)) + 0.10*(mean(med_4(1,*)))]
         
         oplot, findgen(n_elements(rawimage)), med_4(1,*),  psym = 2, thick = 3
         oplot, findgen(n_elements(rawimage)), med_4(2,*),  psym = 4, thick = 3
         oplot, findgen(n_elements(rawimage)), med_4(3,*),  psym = 5, thick = 3

;-----------------------------------------------
;measure the stddev within each pixel
         
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
               if t gt n_elements(rawdataptr) - 2 then begin  ;XXXXX 2
                  arr = arr[0:t-1]
                  sigma = stddev(arr)
            
        ;try sigma with outlier rejection
 ;                 meanclip, arr, clipmean, sigma, clipsig = 4, subs = subarr
            
                  sigmaarr(count) = sigma
                  count = count + 1
               endif 
            endfor
         endfor

         sigmaarr = sigmaarr[0:count - 1]

;try using mmm to reject outliers, and use that gaussian for mean etc.
         mmm, sigmaarr, skymod, skysigma, skew
;oplot, findgen(500), (result(2))/sqrt(2.*!Pi) * exp(-0.5*((findgen(500) - (skymod))/(skysigma))^2.), color = redcolor
         print, 'mmm result', skymod, skysigma


;plot histogram of stddevs of all unmasked pixels
;don't plot it, just determine the values
;change bin size based on exaorure time
         binsize = mean(sigmaarr) / 100.

         plothist, sigmaarr, xhist, yhist, bin = binsize, /noplot 

         a = where(xhist lt skymod + 1*skysigma)
         xhist = xhist(a)
         yhist = yhist(a)


;fit a gaussian to this distribution
;         start = [mean(sigmaarr),stddev(sigmaarr), 100.]
         start = [skymod,skysigma, 80.]
;         print, 'start', start
;      noise = fltarr(n_elements(sigmaarr))
         noise = fltarr(n_elements(yhist))
         noise[*] = 1           ;equally weight the values
;      print, 'n(sigmaarr), n(noise), n(yhist)', n_elements(sigmaarr), n_elements(noise), n_elements(yhist)
         result= MPFITFUN('gauss',xhist,yhist, noise, start,/quiet);./quiet    

;make sure sigma value is positive
         result = abs(result)
;now plot the histogram of stddevs of all unmasked pixels
         plothist, sigmaarr, xjunk, yjunk, bin = binsize, xtitle = 'Standard deviation', ytitle='Number', $
                   xrange=[-8*result(1) + result(0),result(0) + 8*result(1)],$
                   title = "channel "+ string(ch+1)+  "  AFPAT2B" + string(AFPAT2B)  ,$ 
                   thick=3, xthick=3,ythick=3,charthick=3

         xyouts, result(0) + 3.5*result(1), max(yhist) / 2.  + 0.08*max(yhist), 'mean' + string(result(0)), charthick=3
         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. ,  'sigma' + string(result(1)), charthick=3
         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. + 2*0.08*max(yhist), 'N pixels used' + string(count), charthick=3

         xwarm = xjunk
         ywarm = yjunk
         save, xwarm, filename = strcompress('/Users/jkrick/iwic/iwic_recovery3/IRAC017200/xwarm'+string(ichan) +'.sav',/remove_all)
         save, ywarm, filename =strcompress('/Users/jkrick/iwic/iwic_recovery3/IRAC017200/ywarm'+string(ichan) +'.sav',/remove_all)

         if ch eq 0 then meanarr_ch1(aor) = result(0)
         if ch eq 0 then sigmaarr_ch1(aor) = result(1)
         if ch eq 1 then meanarr_ch2(aor) = result(0)
         if ch eq 1 then sigmaarr_ch2(aor) = result(1)

;plot the fitted gaussian and print results to plot

;if exp eq 0 then xarr = findgen(1000)
;if exp eq 1 then xarr = findgen(10000)/ 100.
         if exp eq 0 then xarr = findgen(1000)/ 1000.
         if exp eq 1 then xarr = findgen(1000)/ 1000.



;         oplot, findgen(10000)/ 100., (result(2))/sqrt(2.*!Pi) * exp(-0.5*((findgen(10000)/100. - (result(0)))/(result(1)))^2.), $
         oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
                color = greencolor, thick=3
         print, 'gauss fit result', result

         xwarmfit = xarr
         ywarmfit = (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.)
          save, xwarmfit, filename = strcompress('/Users/jkrick/iwic/iwic_recovery3/IRAC017200/xwarmfit'+string(ichan) +'.sav',/remove_all)
         save, ywarmfit, filename =strcompress('/Users/jkrick/iwic/iwic_recovery3/IRAC017200/ywarmfit'+string(ichan) +'.sav',/remove_all)


;find the median
         print, 'median', median(sigmaarr)
         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 1.*0.08*max(yhist), 'median' + string(median(sigmaarr)), charthick=3

;measure Q =1/<sigma(i)^2>(i)
         sigsquared = sigmaarr^2
         q = 1 / sigsquared
         q = mean(q)
         q = 1 / sqrt(q)
         print, 'Q ', q
         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 2.*0.08*max(yhist), 'Q' + string(q), charthick=3
         
;----------------------------------------------------------
;----------------------------------------------------------
;----------------------------------------------------------
;how many hot pixels are there?
      openw, outlunred, strcompress(dir_name + '/hotpix_' + string(aorname(aor)) + '_'+ string(ichan) + '_'+ string(exptime_num(exp))+'.reg', /remove_all), /get_lun
;      printf, outlunred, 'fk5'
         hotcount = fix(0)
         for x = 1, naxis1 - 1 do begin
            for y = 1, naxis2 -2 do begin
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
      xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 3.*0.08*max(yhist), 'hot pixels '+ string(hotcount), charthick=3
      
      if ch eq 0 then hot_ch1(aor) = hotcount
      if ch eq 1 then hot_ch2(aor) = hotcount
         
   endfor                       ;for each exposure time

   endfor                       ;for each channel

endfor                          ; for each aorition




;now re-order the postscript file
;XXX need to change this once I get all the exposure times.
;newfile = 'psselect -p1,9,17,25,5,13,21,29,2,10,18,26,6,14,22,30,3,11,19,27,7,15,23,31,4,12,20,28,8,16,24,32 ' + dir_name +'/iwic210.ps iwic210_noise.ps'
;spawn, newfile

;for some summary plots


;plot, temp_ch1, meanarr_ch1, psym = 2, xtitle = 'temperature', ytitle = 'mean of stddev', title = 'ch1 100s' $
;      , yrange = [meanarr_ch1(0) - 0.5*meanarr_ch1(0), meanarr_ch1(0) + 0.5*meanarr_ch1(0)], ystyle = 1, $
;      thick = 3, charthick = 3, xthick = 3, ythick = 3
;for j = 0, n_elements(bias_ch1) -1 do xyouts, temp_ch1(j), meanarr_ch1(j), bias_ch1(j)

;plot, temp_ch1, sigmaarr_ch1, psym = 2, xtitle = 'temperature', ytitle = 'sigma of stddev', title = 'ch1 12s' $
;      , yrange = [sigmaarr_ch1(0) - 0.5*sigmaarr_ch1(0), sigmaarr_ch1(0) + 0.5*sigmaarr_ch1(0)], ystyle = 1, $
;      thick = 3, charthick = 3, xthick = 3, ythick = 3

;for j = 0, n_elements(bias_ch1) -1 do xyouts, temp_ch1(j), sigmaarr_ch1(j), bias_ch1(j)

;plot, temp_ch1, hot_ch1, psym = 2,  xtitle = 'temperature', ytitle = 'number of hot pixels ', title = 'ch1 100s' $
;      , yrange = [0,18000], ystyle = 1, $
;      thick = 3, charthick = 3, xthick = 3, ythick = 3
;for j = 0, n_elements(bias_ch1) -1 do xyouts, temp_ch1(j), hot_ch1(j), bias_ch1(j)


;plot, temp_ch2, meanarr_ch2, psym = 2, xtitle = 'temperature', ytitle = 'mean of stddev', title = 'ch2 100s'$
 ;     , yrange = [meanarr_ch2(0) - 0.5*meanarr_ch2(0), meanarr_ch2(0) + 0.5*meanarr_ch2(0)], ystyle = 1, $
 ;     thick = 3, charthick = 3, xthick = 3, ythick = 3

;for j = 0, n_elements(bias_ch2) -1 do xyouts, temp_ch2(j), meanarr_ch2(j), bias_ch2(j)

;plot, temp_ch2, hot_ch2, psym = 2,  xtitle = 'temperature', ytitle = 'number of hot pixels ', title = 'ch2 100s' $
;      , yrange = [hot_ch2(0) - 0.5*hot_ch2(0), hot_ch2(0) + 0.5*hot_ch2(0)], ystyle = 1, $
;      thick = 3, charthick = 3, xthick = 3, ythick = 3
;for j = 0, n_elements(bias_ch2) -1 do xyouts, temp_ch2(j), hot_ch2(j), bias_ch2(j)


;plot, temp_ch2, sigmaarr_ch2, psym = 2, xtitle = 'temperature', ytitle = 'sigma of stddev', title = 'ch2 12s'$
;      , yrange = [sigmaarr_ch2(0) - 0.5*sigmaarr_ch2(0), sigmaarr_ch2(0) + 0.5*sigmaarr_ch2(0)], ystyle = 1, $
;      thick = 3, charthick = 3, xthick = 3, ythick = 3

;for j = 0, n_elements(bias_ch2) -1 do xyouts, temp_ch2(j), sigmaarr_ch2(j), bias_ch2(j)


ps_close, /noprint,/noid

print,'finished ', systime()
end

;write out the image
;            newname = strmid(rawimage(name), 2,36)
;            newname = strcompress(newname + 'wrap.fits', /remove_all)
             ;     print, 'would write to ', dir_name+ '/raw/'+ aorname(aor)+'/'+newname
 ;           fits_write, dir_name+ '/raw/'+ aorname(aor)+'/'+newname, rawdata, rawheader ;'/Users/jkrick/nutella/IRAC/iwic210/'


