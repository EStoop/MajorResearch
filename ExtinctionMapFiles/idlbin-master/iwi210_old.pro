pro iwic210, dir_name

;dir_name = '/Users/jkrick/nutella/IRAC/iwic210/r7665664'
print,'starting ', systime()

!P.multi = [0,1,2]
outputfile = strcompress(dir_name + 'iwic210.ps',/remove_all) ;'/Users/jkrick/nutella/IRAC/iwic210/iwic210_noise.ps'
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
;look at the bright source test to see if there are latents
;-------------------------------------------------------------
;-------------------------------------------------------------

;track the locations on the chip on which the bright sources are, and
;see if their aperture fluxes vary over time.

;read in the correct data

;organize input files
command1 = ' find '+dir_name+' -name "*.fits" > ' + dir_name + '/files_test.list'
command2 = 'grep ch1 < '+dir_name+'/files_test.list | grep -v pmask | grep _bcd.fits > '+dir_name+'/ch1_bcd_test.list'
command3 = 'grep ch2 < '+dir_name+'/files_test.list | grep -v pmask | grep _bcd.fits > '+dir_name+'/ch2_bcd_test.list'

a = [command1, command2, command3]
for i = 0, n_elements(a) -1 do spawn, a(i)

;two star positions
ra = [271.015833, 270.93958333]
dec = [66.928333, 66.93416667 ]

gain = [3.3,3.7]


readcol, dir_name+'ch1_bcd_test.list', bcdname_ch1, format="A",/silent
readcol, dir_name+'ch2_bcd_test.list', bcdname_ch2, format="A",/silent
nimages = n_elements(bcdname_ch1)

;two stars to keep track of
xcenarr_1 = (fltarr(nimages))
ycenarr_1 = fltarr(nimages)
xcenarr_2 = (fltarr(nimages))
ycenarr_2 = fltarr(nimages)

;for each channel
for k = 0, 1 do begin
   if k eq 0 then bcdname = bcdname_ch1
   if k eq 1 then bcdname = bcdname_ch2
;find the locations of the stars on all frames
   for i = 0, nimages - 1 do begin
      header = headfits(bcdname(i)) ;read in the header of the image on which we are working
      adxy, header, ra, dec, xcen, ycen
   ;keep track of where the two stars are
      xcenarr_1(i) = xcen(0)
      ycenarr_1(i) = ycen(0)
      xcenarr_2(i) = xcen(1)
      ycenarr_2(i) = ycen(1)

   endfor

   fluxarr_1 = fltarr(nimages,nimages)
   fluxarr_2 = fltarr(nimages,nimages)

;measure the flux at each of those positions
   for i = 0, nimages - 1 do begin
      fits_read, bcdname(i), data, header
   ;star 1
      aper,data, xcenarr_1, ycenarr_1, flux , fluxerr, sky, $
           skyerr, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
;   flux(where(finite(flux) lt 1)) = 1
      fluxarr_1(i,*) = flux
   ;star 2
      aper,data, xcenarr_2, ycenarr_2, flux2 , fluxerr2, sky2, $
           skyerr2, gain(k), [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent
      nanflux = where(finite(flux2) lt 1)

      fluxarr_2(i,*) = flux2
 
   endfor

   ;plot flux as a function of frame number (time)
   ;star 1
   plot, findgen(nimages), fluxarr_1(0,*) / max(fluxarr_1(0,*)), yrange=[1E-4,1], thick = 3, charthick = 3, /ylog,$
         xthick = 3, ythick = 3, title =  strcompress('Ch' + string( k+1) + ' STAR 1'), xtitle = 'Frame Number', ytitle = 'Fractional Flux '
   oplot, findgen(nimages), fluxarr_1(0,*)/ max(fluxarr_1(0,*)), psym = 2, thick = 3
 
   oplot, findgen(nimages), fluxarr_1(1,*) / max(fluxarr_1(1,*)), thick = 3, color = colorname[1]
   oplot, findgen(nimages), fluxarr_1(1,*) / max(fluxarr_1(1,*)), thick = 3, color = colorname[1], psym = 2

  for i =2, nimages - 2 do begin
       if i mod 2 eq 0 then newi = i
      if i  mod 2 eq 1 then newi = i -1
      oplot, findgen(nimages-2), fluxarr_1(i,newi:*) / max(fluxarr_1(i,*)), thick = 3, color = colorname[i]
      oplot, findgen(nimages-2), fluxarr_1(i,newi:*) / max(fluxarr_1(i,*)), thick = 3, color = colorname[i], psym = 2
   endfor

  ;print out the fluxes at the star to see if photometry is stable.
   y = [.013,.02, .03,.05,.08, .13, .2,.3,.5,.8]
   for i = 0, nimages - 2 do xyouts, 6, y(i), max(fluxarr_1(i,*))  

   ;star 2
   plot, findgen(nimages), fluxarr_2(0,*) / max(fluxarr_2(0,*)), yrange= [1E-4,1], thick = 3, charthick = 3, $
         xthick = 3, ythick = 3, title =  strcompress('Ch' + string( k+1) + ' STAR 2'), xtitle = 'frame number', ytitle = 'Fractional Flux ', /ylog
   oplot, findgen(nimages), fluxarr_2(0,*) / max(fluxarr_2(0,*)), psym = 2, thick = 3

   oplot, findgen(nimages), fluxarr_2(1,*) / max(fluxarr_2(1,*)), thick = 3, color = colorname[1]
   oplot, findgen(nimages), fluxarr_2(1,*) / max(fluxarr_2(1,*)), thick = 3, color = colorname[1], psym = 2


   for i =2, nimages - 2 do begin;2
      if i mod 2 eq 0 then newi = i
      if i  mod 2 eq 1 then newi = i -1
      oplot, findgen(nimages-2), fluxarr_2(i,newi:*) / max(fluxarr_2(i,*)), thick = 3, color = colorname[i]
      oplot, findgen(nimages-2), fluxarr_2(i,newi:*) / max(fluxarr_2(i,*)), thick = 3, color = colorname[i], psym = 2
   endfor

   for i = 0, nimages - 2 do xyouts, 6, y(i), max(fluxarr_2(i,*))

endfor
;----------------

;display them in ds9

for j = 0,  1 do begin  ;for each channel
   if j eq 0 then bcdname = bcdname_ch1
   if j eq 1 then bcdname = bcdname_ch2
   command = 'ds9 '
   for i = 0, n_elements(bcdname) -1 do   command = command + ' ' + bcdname(i) + ' -zscale '
   command = command + ' -single'

   print, 'A ds9 image is opening for channel ' + string(j+ 1)
   print, 'Look at the bright sources images for artifacts or latents.'
   print, 'Close the ds9 window to continue.'
   latents = 'no'
   while (latents eq 'no') do begin
;      spawn, command
      read, latents, prompt = 'Are you done looking for artifacts in the ds9 window? (yes or no) '
   endwhile


endfor


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
;exptime = ['0pt1','0pt4','2','6','12','100','200']
exptime = ['0pt1','2','12','100']
exptime_num  = [0.1,2., 12., 100.]
binsize = [.1, .1, .001, .0001]
;binsize = [.0001, .0001, .0001, .0001]

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
readcol, dir_name + '/dir_list_A1', aorname, format="A",/silent


for pos = 0,   n_elements(aorname)   do begin; 
;print, 'working on position', pos
;filenames appropriate to each position

   for ch = 0, 1 do begin       ;for ch1 and ch2

      for exp = 0, n_elements(exptime) - 1 do begin
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

         rawheaderptr = [Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader)]

;need to make these have the correct number of images as read in
;from rawlist
         rawdataptr = rawdataptr[0:n_elements(rawimage)-1]
         rawheaderptr = rawheaderptr[0:n_elements(rawimage)-1]
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

               a = where(maskdata gt 0)
               rawdata[a] = alog10(-1) ;masked pixels get NaN's

;select background pixels based on a good region of the frame
               hextract, rawdata, bcdheader, 13,248,13,248,/silent
               naxis1 = fxpar(bcdheader, 'NAXIS1')
               naxis2 = fxpar(bcdheader, 'NAXIS2')

            endelse

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
;      openw, outlunred, dir_name + 'highsigma.reg', /get_lun
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
                                ;             printf, outlunred, 'circle( ', x, y, ' 4")'
               endif
               
            endfor
         endfor
;      close, outlunred
;      free_lun, outlunred

         print, 'number of hot pixels', hotcount
         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. - 3.*0.08*max(yhist), 'hot pixels '+ string(hotcount), charthick=3
         
      endfor                    ;for each exposure time

   endfor                       ;for each channel

endfor                          ; for each position


ps_close, /noprint,/noid


;now re-order the postscript file
;XXX need to change this once I get all the exposure times.
newfile = 'psselect -p1,9,17,25,5,13,21,29,2,10,18,26,6,14,22,30,3,11,19,27,7,15,23,31,4,12,20,28,8,16,24,32 ' + dir_name +'/iwic210.ps iwic210_noise.ps'
spawn, newfile

print,'finished noise calculations', systime()
;-------------------------------------------------------------
;-------------------------------------------------------------
;calculate Dn/s/Mjy/sr from the calstar
;-------------------------------------------------------------
;-------------------------------------------------------------

;mosaic together cal star data

;XXXXneed to uncomment this line
;need the bcd's here.
junk = run_mopex( '/Users/jkrick/nutella/IRAC/iwic210/NPM1p67.0636/data/r28899328',  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl')

ra = 269.727857
dec = 67.793591

;read in the warm mosaics we just made
fits_read, '/Users/jkrick/nutella/IRAC/iwic210/NPM1p67.0636/data/r28899328/ch1/Combine/mosaic.fits', warmdata_ch1, warmhead_ch1
fits_read, '/Users/jkrick/nutella/IRAC/iwic210/NPM1p67.0636/data/r28899328/ch2/Combine/mosaic.fits', warmdata_ch2, warmhead_ch2

;ch1
adxy, warmhead_ch1, ra, dec, xcen_warm_ch1, ycen_warm_ch1
aper,  warmdata_ch1, xcen_warm_ch1, ycen_warm_ch1, flux_warm_ch1, fluxerr_warm_ch1, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

;ch2
adxy, warmhead_ch2, ra, dec, xcen_warm_ch2, ycen_warm_ch2
aper,  warmdata_ch2, xcen_warm_ch2, ycen_warm_ch2, flux_warm_ch2, fluxerr_warm_ch2, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

;read in the cold mosaics
;same in this case
;l '/Users/jkrick/nutella/IRAC/iwic210/NPM1p67.0636/data/r28899328/ch1/Combine/mosaic.fits'
fits_read,dir_name + 'cal_mosaic_ch1.fits', colddata_ch1, coldhead_ch1
fits_read, dir_name +'cal_mosaic_ch2.fits', colddata_ch2, coldhead_ch2

;ch1
adxy, coldhead_ch1, ra, dec, xcen_cold_ch1, ycen_cold_ch1
aper,  colddata_ch1, xcen_cold_ch1, ycen_cold_ch1, flux_cold_ch1, fluxerr_cold_ch1, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

;ch2
adxy, coldhead_ch2, ra, dec, xcen_cold_ch2, ycen_cold_ch2
aper,  colddata_ch2, xcen_cold_ch2, ycen_cold_ch2, flux_cold_ch2, fluxerr_cold_ch2, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent


;---------------------------------
fluxconv_cold_ch1 = 0.1088  ;from headers
fluxconv_cold_ch2 = 0.1388  ;from headers

fluxconv_warm_ch1  = flux_warm_ch1*fluxconv_cold_ch1 / flux_cold_ch1
fluxconv_warm_ch2  = flux_warm_ch2*fluxconv_cold_ch2 / flux_cold_ch2

print, 'fluxconv_ch1, ch2', fluxconv_warm_ch1, fluxconv_warm_ch2
;-------------------------------------------------------------
;-------------------------------------------------------------


ps_close, /noprint,/noid
print,'finished ', systime()
end


 ;     oplot, findgen(7), fluxarr_2(i, 2:7)/ max(fluxarr_2(i,*)), color = colorname[i], thick = 3
;      oplot, findgen(nimages - i), fluxarr_2(i,nimages-i)/ max(fluxarr_2(i,*)), color = colorname[i], thick = 3
;     oplot, findgen(nimages - i), fluxarr_2(i,nimages-i)/ max(fluxarr_2(i,*)), psym = 2, color = colorname[i], thick = 3
;      oplot, findgen(nimages), fluxarr_2(i,*)/ max(fluxarr_2(i,*)), psym = 2, color = colorname[i], thick = 3

;   oplot, findgen(nimages-2), fluxarr_2(2,2:*) / max(fluxarr_2(2,*)), thick = 3, color = colorname[2]
;   oplot, findgen(nimages-2), fluxarr_2(3,2:*) / max(fluxarr_2(3,*)), thick = 3, color = colorname[3]
;   oplot, findgen(nimages-2), fluxarr_2(4,4:*) / max(fluxarr_2(4,*)), thick = 3, color = colorname[4]
;   oplot, findgen(nimages-2), fluxarr_2(5,4:*) / max(fluxarr_2(5,*)), thick = 3, color = colorname[5]
;   oplot, findgen(nimages-2), fluxarr_2(6,6:*) / max(fluxarr_2(6,*)), thick = 3, color = colorname[6]
;   oplot, findgen(nimages-2), fluxarr_2(7,6:*) / max(fluxarr_2(7,*)), thick = 3, color = colorname[7]
