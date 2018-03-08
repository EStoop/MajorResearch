pro hot_pix_test, dir_name, aor_list, aor_bright, aor_cal

;dir_name = '/Users/jkrick/iwic/pre_IWIC_2/IRAC016500'
;aor_list = 'aor_t39r35b75_50'

print,'starting ', systime()

!P.multi = [0,1,2]
outputfile = strcompress(dir_name + '/junk_'+strmid(aor_list, 4)+'.ps',/remove_all) 
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

colorname = [ graycolor, redcolor , bluecolor , greencolor , yellowcolor , cyancolor , orangecolor , purplecolor ,  browncolor,graycolor, redcolor , bluecolor , greencolor , yellowcolor , cyancolor , orangecolor , purplecolor ,  browncolor,graycolor, redcolor , bluecolor , greencolor , yellowcolor , cyancolor , orangecolor , purplecolor ,  browncolor ]

;organize all the filenames 
filenames, dir_name, aor_list
;-------------------------------------------------------------
;-------------------------------------------------------------
;look at the bright source test to see if there are latents
;-------------------------------------------------------------
;-------------------------------------------------------------

;run_latent, dir_name, aor_bright

;-------------------------------------------------------------
;-------------------------------------------------------------

;calculate Dn/s/Mjy/sr from the calstar

;flux_conv = run_calstar(dir_name, aor_cal)
flux_conv = [.0873998,.12809]
print, 'flux_conv', flux_conv
;-------------------------------------------------------------
;-------------------------------------------------------------


;plot out the noise properties at this setpoint
;-------------------------------------------------------------
;-------------------------------------------------------------

;get/make a mask image of the dark field data. = /IRAC/iwic210/dark_mask.fits
;	 This is a SExtractor object image
;	 0's for all backgroud pixels, >0's for all object pixels
;	 is not perfect, but is a reasonabe compromise between covering too much area, and getting all object flux masked
;  use the same mask for ch1 and ch2?  definitely use it for all exptimes
fits_read,strcompress(dir_name + '/dark_mask.fits',/remove_all) , objectmaskdata, objectmaskheader ;'/Users/jkrick/nutella/IRAC/iwic210/dark_mask.fits'
;---------
exptime = ['0pt1','0pt4','2','12','100','200']
;exptime = ['0pt1','2','12','100']
exptime_num  = [0.1,0.4,2.,12., 100.,200.]
binsize = [.1, .1, .001, .0001]
;binsize = [.0001, .0001, .0001, .0001]
noiselimit = [6.3, 12.8] ;e-/s of the zodiacal light

;exptime = [100,200]
;determine the correct conversion factor from the dark data
;only need to do this once per setpoint

;XXX taking this out for now because stars aren't detected at 50k.
;flux_conv = dark_standards(dir_name)
;print, 'flux_conv', flux_conv
;flux_conv = [ 0.1088,0.1388]  ;from headers
;flux_conv=[0.0873998   ,  0.128090] ;from calstar

;need to create lists with the appropriate file names by channel,
;raw/bcd, exptime, & position
;filenames, dir_name
; for each position (3 positions)
readcol, dir_name +'/'+ aor_list, aorname, format="A",/silent

hcount = fix(0)
hx= fltarr(1000)
hy= fltarr(1000)
fincount = long(0)
finx = dblarr(256*256.*3)
finy = dblarr(256*256.*3)

for pos = 0,   2  do begin; 
;print, 'working on position', pos
;filenames appropriate to each position

   for ch = 1, 1 do begin       ;for ch1 and ch2

      for exp = 4, 4 do begin; n_elements(exptime) - 1 do begin
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


         if pos lt 3 then begin
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

            if pos lt 3 then begin
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
            
            
            dewrap2, rawdata_int, ichan, barrel, fowlnum, pedsig, 0, rawdata
   ;pattern is set to 0, I don't know what this is.


;convert from DN/s to MJy/sr

            rawdata = rawdata / exptime_num(exp) ; first get it into DN/s
            rawdata = rawdata * flux_conv(ch) ;now in Mjy/sr

 ;rawdata is flipped in y ;flip just to make test output look ok
            hreverse, rawdata, rawheader, 2,/silent


; select background pixels based on source mask 
; first change astrometry in header to give the mask image the same x and y of the raw images.
            hastrom, objectmaskdata, objectmaskheader, maskdata, maskheader, bcdheader, missing = 100


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
               fits_write, strcompress(dir_name +'/test'+strmid(rawimage(name),2),/remove_all), rawdata, bcdheader

               a = where(maskdata gt 0)
               rawdata[a] = alog10(-1) ;masked pixels get NaN's

;select background pixels based on a good region of the frame
 ;              hextract, rawdata, bcdheader, 13,248,13,248,/silent
               naxis1 = fxpar(bcdheader, 'NAXIS1')
               naxis2 = fxpar(bcdheader, 'NAXIS2')

            endelse

; want to make some data have the units of electrons/second to use in noisy pixel calcuation
            edata = rawdata

            vdet =  fxpar(rawheader, 'A0625E00')
            vrst =  fxpar(rawheader, 'A0624E00')
            bias = vdet - vrst  ; should be a number like 0.6
;            print, 'bias', (round(bias*1000.))/1000.
            if (round(bias*100.))/100. eq 0.75 then gain = 3.3
            if (round(bias*100.))/100. eq 0.6 then gain = 3.505
            if (round(bias*100.))/100. eq 0.5  then gain = 3.7
            if (round(bias*100.))/100. eq 0.45 then gain = 3.82
            if (round(bias*100.))/100. eq 0.39 then gain = 3.99
            if (round(bias*100.))/100. eq 0.55 then gain = 3.602;interpolation JK
            if (round(bias*100.))/100. eq 0.4 then gain = 3.94; interpolation JK
            if (round(bias*1000.))/1000. eq 0.426 then gain = 4.0; interpolation JK
            if (round(bias*1000.))/1000. eq 0.526 then gain = 3.64; interpolation JK
            if (round(bias*100.))/100. eq 0.48 then gain = 3.76; interpolation JK

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

            medianarr(name) = median(rawdata) ;doesn't really need to be a median, mean will be super close because of what I just did.

;XXX ???will need to correct the frames for this median value because the
;noise properties depend on this value

            *rawdataptr[name] = rawdata
            *rawheaderptr[name] = bcdheader

            ;output these images for tests
            fits_write, strcompress(dir_name +'/test'+strmid(rawimage(name),2),/remove_all), rawdata, bcdheader
         endfor  ;for each raw image



;-----------------------------------------------
;plot the median of the background as a function of image number
;      plot, findgen(n_elements(medianarr)), medianarr, psym = 2, xtitle = 'frame number', ytitle = 'median of background regions',$
;            title = "channel "+string( ch+1)+ "    exptime " + string(exptime(exp)), thick=3, xthick=3,ythick=3,charthick=3

         plot, findgen(n_elements(rawimage)), med_4(0,*),  psym = 1, xtitle = 'frame number', $
               ytitle = 'median of background regions (MJy/sr)',$
               title = "channel "+string( ch+1)+ "    exptime " + string(exptime(exp)) + "  position " + string(pos), thick=3, $
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
                  meanclip, arr, clipmean, sigma, clipsig = 4, subs = subarr
            
                  sigmaarr(count) = sigma
                  count = count + 1
               endif 
            endfor
         endfor

         sigmaarr = sigmaarr[0:count - 1]
;         save, sigmaarr, filename = '/Users/jkrick/iwic/iwicA/IRAC016300/sigmaarr.sav'

;try using mmm to reject outliers, and use that gaussian for mean etc.
         mmm, sigmaarr, skymod, skysigma, skew
;oplot, findgen(500), (result(2))/sqrt(2.*!Pi) * exp(-0.5*((findgen(500) - (skymod))/(skysigma))^2.), color = redcolor
         print, 'mmm result', skymod, skysigma

;plot histogram of stddevs of all unmasked pixels
;don't plot it, just determine the values
;change bin size based on exposure time
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
                   title = "channel "+ string(ch+1)+ "    exptime " + string(exptime(exp)) + "  position " + string(pos),$
                   thick=3, xthick=3,ythick=3,charthick=3

         xyouts, result(0) + 3.5*result(1), max(yhist) / 2.  + 0.08*max(yhist), 'mean' + string(result(0)), charthick=3
         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. ,  'sigma' + string(result(1)), charthick=3
         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. + 2*0.08*max(yhist), 'N pixels used' + string(count), charthick=3

;plot the fitted gaussian and print results to plot

         if exp eq 0 then xarr = findgen(1000)
         if exp eq 1 then xarr = findgen(10000)/ 100.
         if exp eq 2 then xarr = findgen(10000)/ 100.
         if exp eq 3 then xarr = findgen(1000)/ 1000.
         if exp eq 4 then xarr = findgen(1000)/ 1000.



;         oplot, findgen(10000)/ 100., (result(2))/sqrt(2.*!Pi) * exp(-0.5*((findgen(10000)/100. - (result(0)))/(result(1)))^2.), $
         oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
                color = greencolor, thick=3
         print, 'gauss fit result', result




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

;make smoothed images to work with.
         for i = 0, n_elements(rawdataptr) - 1 do begin
            smoothdata = filter_image(*rawdataptr[i], median=5, /all)
            *rawdataptr[i] = *rawdataptr[i] - smoothdata  ;careful here, overwriting *rawdataptr, can't use that after this runs.
            fits_write, strcompress(dir_name +'/smooth'+strmid(rawimage(i),2),/remove_all), *rawdataptr[i],  *rawheaderptr[i]

            smoothe =  filter_image(*edataptr[i], median=5, /all)
            *edataptr[i] = *edataptr[i] - smoothe  ;careful here, overwriting *rawdataptr, can't use that after this runs.

 ;           print, 'efile', strcompress(dir_name +'/e'+strmid(rawimage(i),2),/remove_all)
;            if pos lt 3 then fits_write, strcompress(dir_name +'/e'+strmid(rawimage(i),2),/remove_all), *edataptr[i],  *rawheaderptr[i]
            

            ;how many pixels are left unmasked on edataptr
            area = where(finite(*edataptr[i]) gt 0)
 ;           print, 'n finite', n_elements(a)
         endfor

;         openw, outlunred, dir_name + '/highsigma_test.reg', /get_lun
         hotcount = fix(0)

         for x = 0, naxis1 - 1 do begin
            for y = 0, naxis2 -2 do begin
               hotpixel = 0
               hpix=0
               fin = 0
               for i = 0, n_elements(rawdataptr) - 1 do begin

         ;need result of gaussian fit here so I can't do this in the first loop
                  if finite((*rawdataptr[i])[x,y]) gt 0 then begin
                     fin = fin + 1
                     if i eq 0 then begin
 ;                       print, fincount
                        finx(fincount) = x;+1
                        finy(fincount) = y;+1
                        fincount = fincount + 1
                     endif

                     if  (*rawdataptr[i])[x,y] gt 10*result(0) + medianarr(i) then begin
                        hotpixel = hotpixel + 1
                     endif

                     if (*edataptr[i])[x,y] gt noiselimit[ch] then begin
                        hpix = hpix + 1
                     endif

                  endif
                  
               endfor
      
     ;if that same pixel is hot on all of the images where finite then it must be
                                ;a hot pixel
               if fin gt 5 and hotpixel gt fin - 1 then begin
                  hotcount = hotcount + 1
;                  printf, outlunred, 'circle( ', x, y, ' 4")'
               endif
               nfin = 5  

               if fin gt nfin and hpix gt fin - 1 then begin
                  hx(hcount) = x ;+1
                  hy(hcount) = y ;+ 1
                  hcount = hcount + 1

;                  printf, outlunred, 'circle( ', x+1, y+1, ' 4")  #color=red'
               endif

            endfor
         endfor
 ;     close, outlunred
 ;     free_lun, outlunred
         
      print, 'number of hot pixels', hotcount
      xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 3.*0.08*max(yhist), 'hot pixels '+ string(hotcount), charthick=3
      print, 'fraction of > '+string(noiselimit(ch)) +'e-/s pixels', hcount / float(n_elements(area)),  float(n_elements(area))

      endfor                    ;for each exposure time

   endfor                       ;for each channel

endfor                          ; for each position


;now get rid of repeats in hx and hy
;first cut off zeros and sort them
hx = hx[0:hcount -1]
hy = hy[0:hcount-1]
hy = hy[sort(hx)]
hx = hx[sort(hx)]

;print, 'hx', hx
;print, 'hy', hy
;now find the uniqe ones
shifts=([-1,1])
indices = where(hx ne shift(hx, 1) or hy ne shift(hy, 1)  , count)
hx = hx(indices)
hy = hy(indices)

openw, outlunred, dir_name + '/hotpix_ch2.reg', /get_lun
for i = 0, n_elements(hx) -1 do printf, outlunred, 'circle( ', hx(i), hy(i), ' 4")  #color=red'
close, outlunred
free_lun, outlunred


;try making an image with finite ones as 1 and others as NANs

data = fltarr(256,256) + alog10(-1)
data(finx, finy) = 1
data(hx, hy) = 0
fits_write, '/Users/jkrick/iwic/ch2_hpix_t29b45r355g361d355.fits', data, *rawheaderptr[0]
print, 'finite', n_elements(finx)
print, 'non first row bad pix', n_elements(hy(where(hy ne 0)))

ps_close, /noprint,/noid


;now re-order the postscript file
;XXX need to change this once I get all the exposure times.
;newfile = 'psselect -p1,9,17,25,5,13,21,29,2,10,18,26,6,14,22,30,3,11,19,27,7,15,23,31,4,12,20,28,8,16,24,32 ' + dir_name +'/iwic210.ps iwic210_noise.ps'
;spawn, newfile

print,'finished ', systime()
end

;cut down the whole list to only include NPM cal stars
;   aorlabelarr = strarr(n_elements(bcdname_ch1))
;   for ic = 0, n_elements(bcdname_ch1) - 1 do begin
;      header = headfits(bcdname_ch1(ic))
;      aorlabelarr(ic) = fxpar(header, 'AORLABEL')
;   endfor

;   a = where(strmatch(aorlabelarr, '*NPM*') eq 1)
;   bcdname_ch1 = bcdname_ch1(a)
   
;   aorlabelarr = strarr(n_elements(bcdname_ch2))
;   for ic = 0, n_elements(bcdname_ch2) - 1 do begin
;      header = headfits(bcdname_ch2(ic))
;      aorlabelarr(ic) = fxpar(header, 'AORLABEL')
;   endfor
   
;   a = where(strmatch(aorlabelarr, '*NPM*') eq 1)
;   bcdname_ch2 = bcdname_ch2(a)
   
   
;   print, n_elements(bcdname_ch1), n_elements(bcdname_ch2)
 
