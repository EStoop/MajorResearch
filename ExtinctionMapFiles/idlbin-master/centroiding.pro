pro centroiding
!p.multi = [0,2,1]

;BD+67 1044
ra = 269.727858
dec = 67.793592


;for both channels
for ch = 0,1 do begin
   if ch eq 0 then begin
      dirname= ['/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0037998592', '/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0037998848','/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0037998336','/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038059776','/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038060288','/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0037947648']

      restore, '/Users/jkrick/iwic/pixel_phase_img_ch1_full.sav'   
    
   endif
 if ch eq 1 then begin

    dirname = ['/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038060032', '/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038060544', '/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038060800', '/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038061056', '/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038059520','/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0037947904']

    restore, '/Users/jkrick/iwic/pixel_phase_img_ch2_full.sav'   

 endif

   for d = 0, n_elements(dirname) - 1 do begin
      if d eq 0 and ch eq 0 then command = 'ls '+dirname(d) +'/IRAC.1*.bcd_fp.fits > /Users/jkrick/iwic/ch1list.txt'
      if d gt 0 and ch eq 0 then command =  'ls '+dirname(d) +'/IRAC.1*bcd_fp.fits >> /Users/jkrick/iwic/ch1list.txt'
      if d eq 0 and ch eq 1 then command = 'ls '+dirname(d) +'/IRAC.2*bcd_fp.fits > /Users/jkrick/iwic/ch2list.txt'
      if d gt 0 and ch eq 1 then command =  'ls '+dirname(d) +'/IRAC.2*bcd_fp.fits >> /Users/jkrick/iwic/ch2list.txt'

      spawn, command
   endfor


;read in the list of fits files
   if ch eq 0 then readcol,'/Users/jkrick/iwic/ch1list.txt', fitsname, format = 'A', /silent
   if ch eq 1 then readcol,'/Users/jkrick/iwic/ch2list.txt', fitsname, format = 'A', /silent

;make an object to keep track of numbers

   s_dist = fltarr(n_elements(fitsname))
   cntrd_dist = fltarr(n_elements(fitsname))


   for f = 0, n_elements(fitsname) -1 do begin
      fits_read, fitsname(f), data, header

      ;get rid of NAN's for cntrd and aper!!!
      ;aper can't handle NAN's inside of the aperture!
;      a = where(finite(data) lt 1)
;      data[a] = data[a+1]

      ;apply the correction, see what happens
;      data = data / arraycorr

      adxy, header, ra, dec, x,y
      cntrd, data, x,y,xcen, ycen, 5
      xcen = xcen - fix(xcen)
      ycen = ycen - fix(ycen)
      cntrd_dist(f) = sqrt((xcen - 0.5)^2 + (ycen - 0.5)^2)
;      print, 'cntrd',x,y, xcen, ycen, sqrt((xcen - 0.5)^2 + (ycen - 0.5)^2)
;      aper, data, xcen, ycen, flux, errap, sky, skyerr, 1, [10], [12,20],/NAN,/flux,/silent
 
      get_centroids, fitsname(f), t, dt, sxcen, sycen, flux, xs, ys, fs, b, /WARM, /APER, RA=ra, DEC=dec,/silent
      sxcen = sxcen - fix(sxcen)
      sycen = sycen - fix(sycen)
      s_dist(f) = sqrt((sxcen - 0.5)^2 + (sycen - 0.5)^2)

   ;make a correction for pixel phase 
;      corrected_flux = correct_pixel_phase(ch+1, xcen, ycen, flux)

   endfor
 
   plot, s_dist, cntrd_dist, psym = 3, title = 'ch' + string(ch+1), xtitle = 'distance from pixel center; get_centroids.pro', ytitle = 'distance from pixel center; cntrd.pro'
   oplot, [0,1],[0,1]

endfor





end

