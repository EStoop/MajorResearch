pro array_phot_corr
!p.multi = [0,1,1]

;BD+67 1044
ra = 269.727858
dec = 67.793592

aprad = [2,3,5,8,9,10,11,12]
aprad = [5,8,9,11,12]
aprad = 12
for apcount = 0, n_elements(aprad) - 1 do begin

;for both channels
for ch = 0,1 do begin
   if ch eq 0 then begin
                                ;dirname = ['/Users/jkrick/iwic/iwicDAM/IRAC021500/bcd/0037127168',  '/Users/jkrick/iwic/iwicDAM/IRAC021500/bcd/0037127680']
      dirname= ['/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0037998592', '/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0037998848','/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0037998336','/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038059776','/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038060288','/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0037947648','/Users/jkrick/irac_warm/edge/r38683392/ch1/bcd','/Users/jkrick/irac_warm/edge/r38683904/ch1/bcd']
      
      restore, '/Users/jkrick/iwic/pixel_phase_img_ch1_full.sav'   
;      fits_read, '/Users/jkrick/iwic/ch1_photcorr_warm_gauss_recal.fits', arraycorr, corhead
    
   endif
   if ch eq 1 then begin
;    dirname = ['/Users/jkrick/iwic/iwicDAM/IRAC021500/bcd/0037127424',  '/Users/jkrick/iwic/iwicDAM/IRAC021500/bcd/0037127936', $
;               '/Users/jkrick/iwic/iwicDAM/IRAC021500/bcd/0037128448', '/Users/jkrick/iwic/iwicDAM/IRAC021500/bcd/0037128960',$
;               '/Users/jkrick/iwic/iwicDAM/IRAC021500/bcd/0037126912', '/Users/jkrick/iwic/iwicDAM/IRAC021500/bcd/0037128192',$
;               '/Users/jkrick/iwic/iwicDAM/IRAC021500/bcd/0037128704', '/Users/jkrick/iwic/iwicDAM/IRAC021500/bcd/0037126656']

      dirname = ['/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038060032', '/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038060544', '/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038060800', '/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038061056', '/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0038059520','/Users/jkrick/irac_warm/pc5/IRAC022600/bcd/0037947904','/Users/jkrick/irac_warm/edge/r38683136/ch2/bcd','/Users/jkrick/irac_warm/edge/r38683648/ch2/bcd']
      
      restore, '/Users/jkrick/iwic/pixel_phase_img_ch2_full.sav'   
;    fits_read, '/Users/jkrick/iwic/ch2_photcorr_warm_gauss_recal.fits', arraycorr, corhead

   endif
   
   for d = 0, n_elements(dirname) - 1 do begin
      if d eq 0 and ch eq 0 then command = 'ls '+dirname(d) +'/IRAC.1*.bcd_fp.fits > /Users/jkrick/iwic/ch1list.txt'
      if d gt 0 and ch eq 0 then command =  'ls '+dirname(d) +'/IRAC.1*bcd_fp.fits >> /Users/jkrick/iwic/ch1list.txt'
     ; spawn, command
      if d gt 0 and ch eq 0 then begin 
         command1 =  'ls '+dirname(d) +'/SPITZER_I1*cbcd.fits >> /Users/jkrick/iwic/ch1list.txt'
         spawn, command1
      endif

      if d eq 0 and ch eq 1 then command = 'ls '+dirname(d) +'/IRAC.2*bcd_fp.fits > /Users/jkrick/iwic/ch2list.txt'
      if d gt 0 and ch eq 1 then command =  'ls '+dirname(d) +'/IRAC.2*bcd_fp.fits >> /Users/jkrick/iwic/ch2list.txt'
      spawn, command
      if d gt 0 and ch eq 1 then begin
         command1 =  'ls '+dirname(d) +'/SPITZER_I2*cbcd.fits >> /Users/jkrick/iwic/ch2list.txt'
         spawn, command1
      endif

   endfor
 

;read in the list of fits files
   if ch eq 0 then readcol,'/Users/jkrick/iwic/ch1list.txt', fitsname, format = 'A', /silent
   if ch eq 1 then readcol,'/Users/jkrick/iwic/ch2list.txt', fitsname, format = 'A', /silent

;make an object to keep track of numbers
   bdobject = replicate({bdob, xcen:0D,ycen:0D,aperflux:0D},n_elements(fitsname))
   c = 0
   s_dist = fltarr(n_elements(fitsname))
   cntrd_dist = fltarr(n_elements(fitsname))

   for f = 0, n_elements(fitsname) -1 do begin
      fits_read, fitsname(f), data, header

      ;get rid of NAN's for cntrd and aper!!!
      ;aper can't handle NAN's inside of the aperture!
      a = where(finite(data) lt 1)
      data[a] = data[a+1]

      ;apply the correction, see what happens
;      data = data / arraycorr

;      adxy, header, ra, dec, x,y
;      cntrd, data, x,y,xcen, ycen, 5
;      xcen = xcen - fix(xcen)
;      ycen = ycen - fix(ycen)
 ;     cntrd_dist(f) = sqrt((xcen - 0.5)^2 + (ycen - 0.5)^2)
;      aper, data, xcen, ycen, flux, errap, sky, skyerr, 1, [10], [12,20],/NAN,/flux,/silent
 
      get_centroids, fitsname(f), t, dt, xcen, ycen, flux, xs, ys, fs, b, /WARM, RA=ra, DEC=dec,/aper, APRAD = aprad(apcount),/silent;,  boxwidth=3, backboxwidth=10, boxborder=5; /aper
      ;print, 'x, y', xcen, ycen, format = '(A, F10.2, F10.2)'
      sxcen = xcen - fix(xcen)
      sycen = ycen - fix(ycen)
      s_dist(f) = sqrt((sxcen - 0.5)^2 + (sycen - 0.5)^2)

   ;make a correction for pixel phase 
   ;corrected_flux =correct_pixel_phase(ch+1, xcen, ycen,;flux)   ;pre-linearity correction
      pixap = strcompress(string(aprad(apcount)) + '_12_20',/remove_all)
      corrected_flux = pixel_phase_correct_gauss(flux,xcen,ycen,ch+1,pixap)

   ;need to keep track of the pixel locations and the flux
   ;can't get too near the edge because the aperture
   ;photometry needs an annulus
      if finite(flux) gt 0 then begin
         bdobject[c].xcen = xcen
         bdobject[c].ycen = ycen
         bdobject[c].aperflux = corrected_flux
         c = c + 1

      endif else begin
         print, 'working on ', fitsname(f)
         print, 'not finite', xcen, ycen
         ;print, 'min', min(data[xcen-20:xcen+20, ycen-20:ycen+20])
;         print, 'ap', mean(data[xcen-10:xcen+10, ycen-10:ycen-10])
;         print, 'sky', mean(data[xcen-20:xcen+20, ycen-20:ycen-20])

      endelse

   endfor
   bdobject = bdobject[0:c-1]

;   for j = 0, c - 1 do print, 'xcen', bdobject[j].xcen, 'ycen', bdobject[j].ycen, 'flux', bdobject[j].aperflux


print, 'bdobject.aperflux', bdobject.aperflux
print, 'n_aper',n_elements(bdobject.aperflux)
;-----------------------------------------------------------------------
; binning

   save, filename = '/Users/jkrick/idlbin/test.sav', bdobject

   x1 = bdobject.xcen
   y1 = bdobject.ycen
   z1 = bdobject.aperflux
   z1 = z1 / max(z1)            ; normalize for ease of viewing

;need to bin so that I can average.
; not so easy in two dimensions

   xcenters = [10,22,48,74,100,126,152,178,204,230,246]
   ycenters = [10,22,48,74,100,126,152,178,204,230,246]

   bins = fltarr(121,480)
   count = fltarr(121)
   for i = 0,n_elements(bdobject.xcen) - 1 do begin
      c = 0
      for k = 0, n_elements(xcenters) -1 do begin
         for l = 0, n_elements(ycenters) - 1 do begin
;         print, 'looking', bdobject[i].xcen, bdobject[i].ycen, xcenters(k), ycenters(l)
            if bdobject[i].xcen lt xcenters(k) + 12. and bdobject[i].xcen ge xcenters(k) - 12. and $
               bdobject[i].ycen lt ycenters(l) +12. and bdobject[i].ycen ge ycenters(l) - 12. then begin
;            print, 'match', bdobject[i].xcen, bdobject[i].ycen, xcenters(k), ycenters(l)
               bins(c,count(c)) =  bdobject[i].aperflux ; make a sum of all points here
               count(c) = count(c) + 1
            endif
            c = c + 1
         endfor
      endfor
   endfor
   flux = fltarr(121)
   st = fltarr(121)
   for n = 0, n_elements(flux) - 1 do begin
      new = bins(n,0:count(n)-1)
      flux(n) = mean(new)
      st(n) = stddev(new)
   endfor
   


;store stddev in percentage
   st = st / flux 
   
;try smoothing before fitting
;   flux = smooth(flux,3)

;invert flux
;inverting because that is the previous standard from IRAC
   flux = 1 / flux
  ; flux = flux / max(flux)    ;normalize all at once at the end
      print, 'flux before', flux

 
   y = [ycenters, ycenters, ycenters, ycenters, ycenters, ycenters, ycenters, ycenters, ycenters, ycenters, ycenters] ;9 times 
   x = [xcenters(0), xcenters(0), xcenters(0), xcenters(0), xcenters(0), xcenters(0), xcenters(0), xcenters(0), xcenters(0), xcenters(0), xcenters(0),$
        xcenters(1), xcenters(1), xcenters(1), xcenters(1), xcenters(1), xcenters(1), xcenters(1), xcenters(1), xcenters(1), xcenters(1), xcenters(1), $
        xcenters(2), xcenters(2), xcenters(2), xcenters(2), xcenters(2), xcenters(2), xcenters(2), xcenters(2), xcenters(2), xcenters(2), xcenters(2),$
        xcenters(3), xcenters(3), xcenters(3), xcenters(3), xcenters(3), xcenters(3), xcenters(3), xcenters(3), xcenters(3),xcenters(3), xcenters(3), $
        xcenters(4), xcenters(4), xcenters(4), xcenters(4), xcenters(4), xcenters(4), xcenters(4), xcenters(4), xcenters(4), xcenters(4), xcenters(4),$
        xcenters(5), xcenters(5), xcenters(5), xcenters(5), xcenters(5), xcenters(5), xcenters(5), xcenters(5), xcenters(5), xcenters(5), xcenters(5),$
        xcenters(6), xcenters(6), xcenters(6), xcenters(6), xcenters(6), xcenters(6), xcenters(6), xcenters(6), xcenters(6), xcenters(6), xcenters(6),$
        xcenters(7), xcenters(7), xcenters(7), xcenters(7), xcenters(7), xcenters(7), xcenters(7), xcenters(7), xcenters(7), xcenters(6), xcenters(6),$
        xcenters(8), xcenters(8), xcenters(8), xcenters(8), xcenters(8), xcenters(8), xcenters(8), xcenters(8), xcenters(8), xcenters(8), xcenters(8),$
        xcenters(9), xcenters(9), xcenters(9), xcenters(9), xcenters(9), xcenters(9), xcenters(9), xcenters(9), xcenters(9), xcenters(9), xcenters(9),$
        xcenters(10), xcenters(10), xcenters(10), xcenters(10), xcenters(10), xcenters(10), xcenters(10), xcenters(10), xcenters(10),xcenters(10), xcenters(10)]
   
;test removing bins from the fit
   print, 'n_x', n_elements(x)
   print, 'n_y', n_elements(y)
   print, 'n_f', n_elements(flux)
   sf = stddev(flux)
   mf = median(flux)
   bad = where( flux lt mf - 3*sf or flux gt mf +3*sf,bcount)
   good = where(flux gt mf - 3*sf and flux lt mf + 3*sf,gcount)
   print, 'bad', bad, mf, sf
;keep the bad version of flux, and replace that bad value with
;it's neighbors value
;although we don't need a regular grid for the fit, we do for the display
   flux_plot = flux
   if bcount gt 0 then flux_plot(bad) = flux(bad-1)

   x_fit = x
   y_fit = y
   x = x(good)
   y = y(good)
   flux = flux(good)
   print, 'flux after ', flux
   print, 'flux_plot', flux_plot
;----------------------------------------------------
;now fitting

;   start = [1D, -1E-5, 1E-5, -1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5, 1E-5]
   noise =st
;noise[*] = 1                    ;equally weight the values
   
;   result = mpfit2dfun('polynomial_5', x, y, flux, noise, start)
;   fit =  result(0) + result(1)*x^2 + result(2)*x*y + result(3)*y^2 + $
;          result(4)*x^3 + result(5)*(x^2)*y + result(6) *x*(y^2) + result(7)*y^3 $
;          + result(8)*x^4 + result(9)*(x^3)*y + result(10)*(x^2)*(y^2) + result(11)*x*(y^3) + result(12)*y^4 $
;          + result(13)*x^5 + result(14)*(x^4)*y + result(15)*(x^3)*(y^2) + result(16)*(x^2)*(y^3) + result(17)*x*(y^4) + result(18)*y^5
;----------
   start = [1.0114,-3.5E-6,-6.82E-5,-1.61E-8,1.21E-6, 1.05E-6]
   P = mpfit2dfun('polynomial_hora', x, y, flux, noise, start, status = msg, errmsg = err)
   fit =   P(0) + P(1)*(x_fit-128) + P(2)*(y_fit-128) + P(3)*(x_fit-128)*(y_fit-128) +$
           P(4)*(x_fit-128)^2 + P(5)*(y_fit-128)^2
   print, 'P', P

;-------
;   start = [1.0,.04,50.,30.,120.,120.,0.5]  
;   P = mpfit2dfun('mygauss2d', x, y, flux, noise, start, status = msg, errmsg = err)
;   print, 'status', msg, 'errmsg', err
;   fit =   P(0) + P(1)*exp(-((((x-P(4))*cos(P(6)) - (y - P(5)) *sin(P(6))) / P(2) )^2 + (((x-P(4))*sin(P(6)) + (y-P(5))*cos(P(6))) / P(3) )^2)/2.) 


;use the gaussian fit, so make an array for the full 256x256 fit

;fill the x and y arrays appropriately
y  = findgen(65536)
c = 0l
for i = 0,255 do begin
   for j = 0, 255 do begin
      y(c) = j
      c= c + 1
   endfor
endfor

x  = findgen(65536)
c = 0l
for i = 0,255 do begin
   for j = 0, 255 do begin
      x(c) = i
      c= c + 1
   endfor
endfor

;large_fit =  P(0) + P(1)*exp(-((((x-P(4))*cos(P(6)) - (y - P(5)) *sin(P(6))) / P(2) )^2 + (((x-P(4))*sin(P(6)) + (y-P(5))*cos(P(6))) / P(3) )^2)/2.) 
large_fit =   P(0) + P(1)*(x-128) + P(2)*(y-128) + P(3)*(x-128)*(y-128) +P(4)*(x-128)^2 + P(5)*(y-128)^2

;--------------------------------------------------------
;plotting

   ;set up the array in the right way to be able to use surface
   sflux = fltarr(n_elements(xcenters), n_elements(xcenters))
   print, 'xcen, ycen', n_elements(xcenters), n_elements(ycenters)
   c = 0
   for i = 0, n_elements(xcenters) -1 do begin
      for j = 0, n_elements(ycenters) - 1 do begin
         sflux(i,j) = flux_plot(c)
         c = c + 1
      endfor
   endfor
   help, sflux
   print, 'sflux', sflux
   window,0,retain=2,xsize=500,ysize=500
   surface, sflux,  title='Mesh plot of  data', xtitle='X',   ytitle='Y', ztitle='Z', charsize=2 
   
   window,1,retain=2, title='Surface plot of irregular data',xsize=500,ysize=500
   shade_surf, sflux, title='Surface plot of  data', xtitle='X',   ytitle='Y', ztitle='Z', charsize=2 ;, zrange=[0.8,1.0]
    
   fits_write, strcompress('/Users/jkrick/iwic/ch'+string(ch+1)+'_photcorr_ap_' +string(aprad(apcount)) + '.fits',/remove_all), sflux
;
;----
;want to be able to plot the fit for comparison

   sfit = fltarr(n_elements(xcenters), n_elements(ycenters))
   c = 0
   for i = 0., n_elements(xcenters) -1 do begin
      for j = 0., n_elements(ycenters) - 1 do begin
;      print, i, j, c, fit(c)
         sfit(i,j) = fit(c)
         c = c + 1
      endfor
   endfor
   
   window,2,retain=2,xsize=500,ysize=500
   surface, sfit,title='Mesh plot of fit', xtitle='X',   ytitle='Y', ztitle='Z', charsize=2
   
   window,3,retain=2, title='Surface plot of  data',xsize=500,ysize=500
   shade_surf, sfit, title='Surface plot of fit', xtitle='X',   ytitle='Y', ztitle='Z', charsize=2 ;, zrange=[0.8,1.0]
   
;now make the images the normal 256x256 and write them out
   sfit = fltarr(256,256)
   c = 0L
   for i = 0, 256 -1 do begin
      for j = 0, 256 - 1 do begin
;      print, i, j, c, fit(c)
         sfit(i,j) = large_fit(c)
         c = c + 1
   endfor
   endfor
;-----------------------------------------------

   ;need to normalize,
   central = median(sfit[128-10:128+10,128-10:128+10])
   sfit= sfit  / central ; median(sfit)

   ;make a header with useful information
   mkhdr, shd, sfit

   sxaddpar, shd, 'COMMENT', 'Generated Apr 30, 2010'   
   sxaddpar, shd, 'COMMENT', 'Central and edge data taken during IRAC022600 & 023300'
   sxaddpar, shd, 'COMMENT', 'Uses pixel phase code = pixel_phase_correct_gauss'
   sxaddpar, shd, 'COMMENT', 'Normalized to central 20x20 box'
   
   sxaddpar, shd, 'COMMENT ',' This image corrects an SSC BCD flattened with a zody'                   
   sxaddpar, shd, 'COMMENT','    skyflat to a flatfield appropriate for point sources '                 
   sxaddpar, shd, 'COMMENT','    and a stellar spectrum. Multiply BCD by this image.'                   
   sxaddpar, shd, 'ORIGIN' ,'Spitzer Science Center' ,'Organization generating this FITS file '    
   sxaddpar, shd, 'TELESCOP', 'Spitzer '           ,' SPITZER Space Telescope'                 
   sxaddpar, shd, 'INSTRUME', 'IRAC    '           ,' SPITZER Space Telescope instrument ID'      
   sxaddpar, shd, 'CHNLNUM' ,     ch + 1 ,' 1 digit instrument channel number  '            

   fits_write,  strcompress('/Users/jkrick/iwic/ch'+string(ch+1)+'_photcorr_ap_' +string(aprad(apcount)) + '.fits',/remove_all),sfit,shd
   
;-----------------------------------------------
;compare warm to cold

;   fits_read,strcompress( '/Users/jkrick/iwic/photcorr_cryo/ch'+string(ch+1)+'_photcorr_rj.fits',/remove_all), colddata, coldheader
;   a = where(sfit eq 1.0)
;   print, 'a', a
;   diff = colddata - sfit
   
;   fits_write, strcompress('/Users/jkrick/iwic/ch'+string(ch+1)+'_diff_photcorr_ap_10.fits',/remove_all), diff, coldheader

   ;plot, s_dist, cntrd_dist, psym = 3

endfor

endfor ;apcount




end


;  triangulate, x, y, tr
   
;window,0,retain=2, title='Triangulation of irregular data',xsize=500,ysize=500
;plot, x, y, psym=1, title='Triangulation of irregular data'
;for i=0, n_elements(tr)/3-1 do begin
;    t = [tr(*,i),tr(0,i)]
;    plots, x(t), y(t)
;endfor
   
;   meshed = trigrid(x,y,flux,tr, nx=9, ny=9)
 ;  meshed_shades = byte(255.0*meshed)  
   
;   window,0,retain=2, title='Mesh plot of irregular data',xsize=500,ysize=500
;   surface, meshed, title='Mesh plot of  data', charsize = 3 ;, zrange=[0.8,1.0], ax = 60
   
;window,2,retain=2, title='Surface plot with explicit shades',xsize=500,ysize=500
;shade_surf, meshed, title='Surface plot of irregular data', xtitle='X',   ytitle='Y', ztitle='Z', charsize=3.0, shades=meshed_shades, ax = 50
   
 ;  window,1,retain=2, title='Surface plot of irregular data',xsize=500,ysize=500
;   shade_surf, meshed, title='Surface plot of  data', xtitle='X',   ytitle='Y', ztitle='Z', charsize=3.0 ;, zrange=[0.8,1.0]