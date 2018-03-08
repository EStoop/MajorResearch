pro background_55cnc
;try to see if I can see background blowups as reported by Demory et al 2011

;55cnc staring mode ch2
  ra_ref = 133.14754
  dec_ref = 28.330195

  cd, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/r39524608/ch2/bcd'
  command  =  "find . -name '*_bcd.fits' > /Users/jkrick/irac_warm/pcrs_planets/55cnc/bcdlist.txt"
;  spawn, command
;  readcol,'/Users/jkrick/irac_warm/pcrs_planets/55cnc/bcdlist.txt',fitsname, format = 'A', /silent
  
;  bigim = fltarr(32, 32, n_elements(fitsname)*64)
;  c = 0L
  
;  for i =0.D,  n_elements(fitsname) - 1 do begin 

;     fits_read, fitsname(i), data, header
     
;     if i eq 0 then sclk_0 = sxpar(header, 'SCLK_OBS')

;     for j = 0, 63  do begin
;        bigim(0,0,c) = data[*,*,j]
;        ;print, 'c, j, data,bigim', c, j, data[15, 15, j], bigim[15,15,c]
;        c = c + 1
;     endfor

;  endfor
  
;save, /all, filename = '/Users/jkrick/irac_warm/pcrs_planets/55cnc/r39524608/bigim.sav'

;------------------------------------------------

restore, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/r39524608/bigim.sav'

meanim = fltarr(32,32)
maxim = fltarr(32,32)
rmsim = fltarr(32,32)

For xpix = 0,31 do begin
  
   for ypix = 0,31 do begin

      x = findgen(c)
      y = bigim[xpix,ypix,*]
                                ;print, n_elements(x), n_elements(y)
                                ; a = plot( x, y, '1.', xtitle = 'frame number', ytitle = 'MJy/sr')
      
; ok and what about the RMS on 10 min timescales
      nbin = 64.*200.
      count = 0
      rmsarr = fltarr(c)
      for si= 0, fix(c/nbin) - 1 do begin
                                ;XXXdo the binning/rmsing...
         rmsarr(count) = stddev(bigim[xpix,ypix,si*nbin:si*nbin + (nbin - 1)])
         count = count + 1
      endfor
      rmsarr = rmsarr[0:count-1]
                                ;print, rmsarr
      
                                ;b = plot(findgen(n_elements(rmsarr)), rmsarr, '1.', xtitle = '10 min', ytitle = 'rms')
;      plothist, rmsarr, xhist, yhist, bin = 0.5,/noplot
      if (xpix gt 12 and xpix lt 18 and ypix gt 12 and ypix lt 18) then plotcolor = 'red' else plotcolor = 'black'
      if (xpix eq 0 and ypix eq 0) then begin
;         p = plot(xhist, yhist, '1-', xtitle = 'RMS', ytitle = 'Number', color = plotcolor)
         p = plot(findgen(n_elements(rmsarr)), rmsarr, '1-', color = plotcolor, xtitle = 'time (10 min bins)', ytitle = 'RMS');, yrange = [0,1000])
      endif else begin
;         p = plot(xhist, yhist, '1-', /overplot, color = plotcolor)
         p = plot(findgen(n_elements(rmsarr)), rmsarr, '1-', color = plotcolor,/overplot)
      endelse
      
;      print, 'x, y, mean rms, max rms', xpix, ypix, mean(rmsarr), max(rmsarr)
      meanim(xpix, ypix) = mean(rmsarr)
      maxim(xpix, ypix) = max(rmsarr)
      rmsim(xpix, ypix) = stddev(bigim[xpix, ypix,*])

   endfor                       ;for each ypix
endfor                          ; for each xpix

im1 = image(meanim, image_dimensions = [32,32], position = [0,0,1,1])
im2 = image(maxim, image_dimensions = [32,32], position = [0,0,1,1])
im3 = image(rmsim, image_dimensions = [32,32], position = [0,0,1,1])

save, /all, '/Users/jkrick/irac_warm/pcrs_planets/55cnc/r39524608/rmsim.sav'
end
