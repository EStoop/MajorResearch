pro test_subarray

  bcddataptr  = [Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64)), Ptr_New(findgen(32,32,64))]
  
  bcdheaderptr = [Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader), Ptr_New(objectmaskheader)]
  
  cd, '/Users/jkrick/irac_warm/barrel_shift/pipeline_test/shifted_0bit/'
  bcdlist = 'bcdlist_ch2_s0p02s.txt'
  readcol, bcdlist, bcdimage, format="A",/silent
  
  for name = 0,  n_elements(bcdimage) -1 do begin
     print, 'working on name', name, bcdimage(name)
     fits_read, bcdimage(name), bcddata, bcdheader 

   
     fowlnum = fxpar(bcdheader, 'AFOWLNUM')
     naxis1 = fxpar(bcdheader, 'NAXIS1')
     naxis2 = fxpar(bcdheader, 'NAXIS2')
     print, 'naxis', fowlnum, naxis1, naxis2
     
     nanpix = where(finite(bcddata) eq 0,nnan) ; can't have nan's for sextractor to work.
     print, 'nnan', nnan
     bcddata[nanpix] = -1E5
     
     *bcddataptr[name] = bcddata
     *bcdheaderptr[name] = bcdheader

  endfor

;-----------------------------------------------
;measure the stddev within each pixel
  bcddataptr = bcddataptr[0:n_elements(bcdimage)-1]
  bcdheaderptr = bcdheaderptr[0:n_elements(bcdimage)-1]

  count = float(0)
  sigmaarr = fltarr(naxis1*naxis2)
  meanarr = fltarr(naxis1*naxis2)
  for x = 0, naxis1-1 do begin
     for y = 0,  naxis2 - 1 do begin
        t = 0
        arr = fltarr(n_elements(bcddataptr)*64)
        for i = 0, n_elements(bcddataptr) - 1 do begin
           for j = 0, 64 - 1 do begin
                                ;if finite((*bcddataptr[i])[x,y,j]) gt 0 then begin
             ; print,i, j,  (*bcddataptr[i])[x,y,j]
              arr(t) = (*bcddataptr[i])[x,y,j] ;if not masked then include the value
              t = t + 1
                                ;   endif
           endfor               ; for all i
        endfor                  ; for all j
   
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
     endfor                     ;for all y
  endfor                        ;for all x 
  
  sigmaarr = sigmaarr[0:count - 1]
  meanarr = meanarr[0:count - 1]
  
;try using mmm to reject outliers, and use that gaussian for mean etc.
  mmm, sigmaarr, skymod, skysigma, skew
;oplot, findgen(500), (result(2))/sqrt(2.*!Pi) * exp(-0.5*((findgen(500) - (skymod))/(skysigma))^2.), color = redcolor
  print, 'mmm result', skymod, skysigma
  
  binsize = skymod / 10.        ;mean(sigmaarr,/nan) / 10.
  print, 'mean, bin', mean(sigmaarr,/nan), binsize
  goodsig = where(sigmaarr ne 0. )
  print, 'zeros', n_elements(sigmaarr), n_elements(goodsig)
  sigmaarr = sigmaarr[goodsig]
  plothist, sigmaarr, xhist, yhist, bin = binsize, /noplot ,/nan
  print, 'xhist', xhist
;fit a gaussian to this distribution
  start = [skymod,skysigma, 5000.]
  noise = fltarr(n_elements(yhist)) 
  noise[*] = 1.0                ;equally weight the values
 ; low = where(yhist lt 100, lowcount)
 ; print, 'lowcount', lowcount
 ; noise[low] = 10.0
  result= MPFITFUN('mygauss',xhist,yhist, noise, start) ;./quiet    

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
  a0 = plot(xhist, yhist, '3-',  xtitle = 'Standard Deviation in Raw Counts', ytitle='Number', xrange = [0,200])
  

  end
  
