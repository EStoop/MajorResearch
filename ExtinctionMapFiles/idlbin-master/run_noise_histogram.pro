pro run_noise_histogram

!P.multi = [0,4,4]
greencolor = FSC_COLOR("Green", !D.Table_Size-4)



a='/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159424/ch1/Combine/mosaic.fits'
b= '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159680/ch1/Combine/mosaic.fits'
c= '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159936/ch1/Combine/mosaic.fits'
d= '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037160192/ch1/Combine/mosaic.fits'
ch1_6s = [ a,b,c,d]

ch1_100s = [ '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158400/ch1/Combine/mosaic.fits', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158656/ch1/Combine/mosaic.fits', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158912/ch1/Combine/mosaic.fits','/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159168/ch1/Combine/mosaic.fits']

a='/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159424/ch2/Combine/mosaic.fits'
b= '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159680/ch2/Combine/mosaic.fits'
c= '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159936/ch2/Combine/mosaic.fits'
d='/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037160192/ch2/Combine/mosaic.fits'
ch2_6s = [ a,b,c,d]

ch2_100s = [ '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158400/ch2/Combine/mosaic.fits', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158656/ch2/Combine/mosaic.fits', '/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037158912/ch2/Combine/mosaic.fits','/Users/jkrick/iwic/iwicDAQ/IRAC021800/bcd/0037159168/ch2/Combine/mosaic.fits']


binsize = .01
exptime = ['6','100']

for ch = 0, 1 do begin
   for exp = 0, 1 do begin
      if (ch eq 0 and exp eq 0) then name = ch1_6s 
      if (ch eq 0 and exp eq 1) then name = ch1_100s
      if (ch eq 1 and exp eq 0) then name = ch2_6s 
      if (ch eq 1 and exp eq 1) then name = ch2_100s
      print, 'ch', ch+1, 'exp', exptime(exp)
      for j = 0, n_elements(name) - 1 do begin
         fits_read, name(j), data, header
         data = data[100:600,200:800]
         plothist, data, xhist, yhist,bin = binsize , /noplot 

;fit a gaussian to this distribution
         start = [mean(data),stddev(data), 10000.]
      
         noise = fltarr(n_elements(yhist))
         noise[*] = 1           ;equally weight the values
;      print, 'n(sigmaarr), n(noise), n(yhist)', n_elements(sigmaarr), n_elements(noise), n_elements(yhist)
         result= MPFITFUN('gauss',xhist,yhist, noise, start,/quiet) ;./quiet   
      
      
         plothist, data, xjunk, yjunk, bin = binsize, xtitle = 'pixel values', ytitle='Number', $
                   xrange=[-8*result(1) + result(0),result(0) + 8*result(1)],$
                   title = 'ch' + string(ch + 1) + '   exp' + exptime(exp),  thick = 3
      
         xyouts, result(0) + 3.5*result(1), max(yhist) / 2.  + 0.08*max(yhist), 'mean' + string(result(0))
         xyouts, result(0) + 3.5*result(1), max(yhist) / 2. ,  'sigma' + string(result(1))
      
;plot the fitted gaussian and print results to plot
      
         xarr = findgen(12000)/ 1000.
      
         oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
                color = greencolor, thick=3
         print, 'gauss fit result', result
      endfor

   endfor
endfor




end
