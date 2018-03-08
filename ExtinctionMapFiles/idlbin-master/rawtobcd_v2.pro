pro rawtobcd_v2

  flux_conv = [ 0.1253,0.1469] 
  dirname = '/Users/jkrick/irac_warm/barrel_shift/pipeline_test/'
  cd,dirname
  command = 'ls shifted_2bit/0039204608/IRAC.1*bcd_fp.fits > ch1_bcdlist.txt'
  spawn, command
  command = 'ls shifted_2bit/0039204608/IRAC.2*bcd_fp.fits > ch2_bcdlist.txt'
  spawn, command
  command =  'ls raw/shifted_2bit/SPITZER_I1_0039204608*_dce.fits > ch1_rawlist.txt'
  spawn, command
  command = 'ls raw/shifted_2bit/SPITZER_I2_0039204608*_dce.fits > ch2_rawlist.txt'
  spawn, command
  
  ch = 1
  readcol, 'ch2_rawlist.txt', rawimage, format="A",/silent
  readcol,'ch2_bcdlist.txt', bcdimage, format="A",/silent
  
  fits_read, '/Users/jkrick/irac_warm/barrel_shift/pipeline_test/r45416960/ch1/raw/SPITZER_I1_0045416960_0003_0000_01_dce.fits', data, header  ;junk raw file just to have the header size
  med_4 = fltarr(4, n_elements(rawimage))
;---------
;setup pointer array now to be filled with data arrays and headers
  rawdataptr  = [Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256))]

  rawheaderptr = [Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header)]

  bcddataptr  = [Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256)), Ptr_New(findgen(256,256))]

  bcdheaderptr = [Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header), Ptr_New(header)]

;need to make these have the correct number of images as read in
;from rawlist
  rawdataptr = rawdataptr[0:n_elements(rawimage)-1]
  rawheaderptr = rawheaderptr[0:n_elements(rawimage)-1]
  bcddataptr = bcddataptr[0:n_elements(bcdimage)-1]
  bcdheaderptr = bcdheaderptr[0:n_elements(bcdimage)-1]
;         edataptr = edataptr[0:n_elements(rawimage)-1]
;---------
  for name = 0,  n_elements(rawimage) -1 do begin
     print, 'working on name', name, rawimage(name)
     fits_read, dirname+rawimage(name), rawdata, rawheader 
     fits_read, dirname+bcdimage(name), bcddata, bcdheader 
     
;get the flux values converted from integers to DN
     
     barrel = fxpar(rawheader, 'A0741D00')
     fowlnum = fxpar(rawheader, 'A0614D00')
     pedsig = fxpar(rawheader, 'A0742D00')
     ichan = fxpar (rawheader, 'CHNLNUM')
     framtime = fxpar(bcdheader, 'FRAMTIME')
     exptime = fxpar(bcdheader, 'EXPTIME')
     naxis1 = fxpar(bcdheader, 'NAXIS1')
     naxis2 = fxpar(bcdheader, 'NAXIS2')

                                ;print, 'details',  ichan, barrel, fowlnum, pedsig
     
                                ;use Bill's code to conver to DN
     dewrap2, rawdata, ichan, barrel, fowlnum, pedsig, 0, rawdata
     
                                ;or use Jim's code to convert to DN
;           rawdata = irac_raw2dn(rawdata,ichan,barrel,fowlnum)
     
                                ;convert from DN/s to MJy/sr
     rawdata = rawdata / exptime   ; first get it into DN/s
     rawdata = rawdata * flux_conv(ch)        ;now in Mjy/sr
     
         ;divide into 4 readouts for raw
     cnvrt_stdim_to_64_256_4, rawdata, rawdata_4          
     for read = 0, 3 do begin
        median_read = median(rawdata_4[*,*,read])
        rawdata_4[*,*,read] = rawdata_4[*,*,read] - median_read
     endfor     
    
         ;convert back to a normal 2D image
     cnvrt_im_64_256_4_to_stdim, rawdata_4, rawdata

        ;divide into 4 readouts for bcd
     cnvrt_stdim_to_64_256_4, bcddata, bcddata_4          
     for read = 0, 3 do begin
        median_read = median(bcddata_4[*,*,read])
        bcddata_4[*,*,read] = bcddata_4[*,*,read] - median_read
     endfor         
         ;convert back to a normal 2D image
     cnvrt_im_64_256_4_to_stdim, bcddata_4, bcddata

     *rawdataptr[name] = rawdata
     *rawheaderptr[name] = rawheader
     *bcddataptr[name] = bcddata
     *bcdheaderptr[name] = bcdheader

                                ;save the image for diagnostics
     newname=dirname + strmid(rawimage(name),0, 39)+'raw2bcd.fits'
     print, 'newname', newname
     fits_write, newname, rawdata, rawheader
  endfor           ;for each image              

;----------------------------------------------
     ;now measure the noise on the  raw2bcd, and bcd next...

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
           if t gt n_elements(rawdataptr) - 2 then begin ;XXXXX 2
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
  
;use mmm to reject outliers, and use that gaussian for mean etc.
     mmm, sigmaarr, skymod, skysigma, skew
     print, 'mmm result', skymod
;plot histogram of stddevs of all unmasked pixels
     binsize = mean(sigmaarr) / 5.
     goodsig = where(sigmaarr ne 0.)
;         print, 'zeros', n_elements(sigmaarr), n_elements(goodsig)
     sigmaarr = sigmaarr[goodsig]
     meanarr = meanarr[goodsig]
     plothist, meanarr, xhist, yhist, bin = binsize, /noplot ,/nan
     a1 = plot(xhist, yhist, '2r', ytitle = 'number', xtitle = 'mean', name = 'raw', title = 'ch2 2s shift2', xrange = [-1,1])

;---------------------------------------------------------------------------------
;---------------------------------------------------------------------------------
;---------------------------------------------------------------------------------
;---------------------------------------------------------------------------------
;bcd

;----------------------------------------------
     ;now measure the noise on the  bcd next...

     count = float(0)
     sigmaarr = fltarr(naxis1*naxis2)
     meanarr = fltarr(naxis1*naxis2)
     for x = 0, naxis1-1 do begin
        for y = 0,  naxis2 - 1 do begin
           t = 0
           arr = fltarr(n_elements(bcddataptr))
           for i = 0, n_elements(bcddataptr) - 1 do begin
              if finite((*bcddataptr[i])[x,y]) gt 0 then begin
                 arr(t) = (*bcddataptr[i])[x,y] ;if not masked then include the value
                 t = t + 1
              endif
           endfor
                                ;if there are enough measurements for a sigma
           if t gt n_elements(bcddataptr) - 2 then begin ;XXXXX 2
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
  
;use mmm to reject outliers, and use that gaussian for mean etc.
     mmm, sigmaarr, skymod, skysigma, skew
     print, 'mmm result', skymod
;plot histogram of stddevs of all unmasked pixels
     binsize = mean(sigmaarr) / 5.
     goodsig = where(sigmaarr ne 0.)
;         print, 'zeros', n_elements(sigmaarr), n_elements(goodsig)
     sigmaarr = sigmaarr[goodsig]
     meanarr = meanarr[goodsig]
     plothist, meanarr, xhist, yhist, bin = binsize, /noplot ,/nan
     a2 = plot(xhist, yhist, '2b', /overplot, name= 'bcd')

     l = legend(target=[a1,a2], position = [0.6, 8000],/data)

end
