pro dark_meanvsmode

  ch = 2
  ch1_xrange = [150,450];[ 300, 800]
  ch2_xrange = [150, 250]
  gyrange = [0,8000]
  gxrange = [0,300]
  aordir = '/Users/jkrick/iracdata/flight/IWIC/'
  cd, aordir
  flux_conv = [.1253,.1468]
   
                                ;just work with the last one for now
  if ch eq 1 then filename = './IRAC026200/cal/0041083136/IRAC.1.0041083136.4516290.1.skydark.fits'
  if ch eq 2 then filename = './IRAC026200/cal/0041083136/IRAC.2.0041083136.4516291.1.skydark.fits'
 
;  if ch eq 1 then filename = './IRAC026200/cal/0041082880/IRAC.1.0041082880.4516214.1.skydark.fits'
;  if ch eq 2 then filename = './IRAC026200/cal/0041082880/IRAC.2.0041082880.4516215.1.skydark.fits'
  
;cryo
 ; aordir = '/Users/jkrick/iracdata/flight/calprog/finaldarks/'
 ; cd, aordir

 ; if ch eq 1 then filename = './IRAC009200/cal/0020950528/IRAC.1.0020950528.1826941.9.skydark.fits'
 ; if ch eq 2 then filename = './IRAC009200/cal/0020950528/IRAC.2.0020950528.1826942.6.skydark.fits'
 ; if ch eq 3 then filename = './IRAC009200/cal/0020950528/IRAC.3.0020950528.1826943.6.skydark.fits'
 ; if ch eq 4 then filename = './IRAC009200/cal/0020950528/IRAC.4.0020950528.1826944.6.skydark.fits'

  print, 'working on ', filename
  fits_read, filename, data, header, exten_no = 0
  framtime = sxpar(header, 'FRAMTIME')
  exptime = sxpar(header, 'EXPTIME')
  ch = sxpar(header, 'CHNLNUM')
  gain = sxpar(header, 'GAIN')
  
                                ;need to ignore the uncertainty file
                                ;whichis attached to the data fits
                                ;file, but not a different extension.
  data = data[*,*,0]
  
;get rid of nans
  a = where(finite(data), n_finite)
  data = data[a]
  ;data = data * flux_conv(ch - 1) / framtime
  if ch eq 1 then begin
     binsize = 10.
     gyrange = [0,10000]
     gxrange = ch1_xrange
  endif

  if ch eq 2 then begin
     binsize = 2.
     gyrange = [0, 8000]
     gxrange = ch2_xrange
  endif

  plothist, data, xhist, yhist, bin = binsize, /nan,/noplot ;make a histogram of the data
  a = barplot(xhist, yhist, xtitle = 'Median Pixel Values', ytitle = 'Number', xrange = gxrange, yrange =gyrange)
                                ;use a gaussian to determine the mean background value
  start = [mean(data,/nan), stddev(data,/nan), 4000.0]
  binsize = .1           
  noise = fltarr(n_elements(yhist))
  noise[*] = 1                                                     ;equally weight the                             
  result= MPFITFUN('mygauss',xhist,yhist, noise, start)     ; fit a gaussian to the histogram sorted data
                                ;oplot, xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), linestyle = 2,  color = redcolor     
  
                                ;make a line marking the mean
  b = polyline([result(0), result(0)], [0, 10000], '3r',/data)
  
                                ;calculate mode
  mode_calc, data, mode;,/verbose
  print, 'mode', mode
 ; if ch eq 2 then mode = 199.7
;  if ch eq 1 then mode = 547.622
  c =polyline([mode, mode], [0,10000], '3', /data)

  
  
end
