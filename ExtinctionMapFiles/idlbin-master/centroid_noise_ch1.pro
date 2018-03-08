pro centroid_noise_ch1
; this code compares the centroiding uncertainties from the different pcrspeakup tests taken after the
;switch to iru 2 in fall 2011

;-----------------------------------------------------
;subarray tests
;-----------------------------------------------------

;second test in Sep

sub_aor = ['0044263936','0044262656','0044262912','0044263168','0044263424','0044264960','0044263680','0044264192','0044264448','0044264704']
yarr = fltarr(n_elements(sub_aor)*64*16)
xarr = fltarr(n_elements(sub_aor)*64*16)
ysarr =fltarr(n_elements(sub_aor)*64*16)
xsarr =fltarr(n_elements(sub_aor)*64*16)
total = 0

for a = 0, n_elements(sub_aor) - 1 do begin
   print, 'working on ', sub_aor(a)
   dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(sub_aor(a)) 
   
   CD, dir                      ; change directories to the correct AOR directory
   command  =  " ls *bcd_fp.fits > /Users/jkrick/irac_warm/pcrsfaint/cbcdlist.txt"
   spawn, command
   
   readcol,'/Users/jkrick/irac_warm/pcrsfaint/cbcdlist.txt',fitsname, format = 'A', /silent
   
   for i =1, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
      
      if i eq 1 then begin
         header = headfits(fitsname(i)) ;
         ra_ref = sxpar(header, 'RA_REF')
         dec_ref = sxpar(header, 'DEC_REF')
         gain = sxpar(header, 'GAIN')
         exptime =  sxpar(header, 'EXPTIME') 
         fluxconv =  sxpar(header, 'FLUXCONV') 
         print, 'ra, dec', ra_ref, dec_ref, gain, exptime, fluxconv
      endif
      
      get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra_ref, dec = dec_ref, /silent     

;figure out flux in electrons
      sbtoe = gain*exptime/fluxconv
               ;convert to microJy
      eflux = abcdflux *1E6
      eunc = fs*1E6
               ;convert to Mjy/sr
      eflux = eflux / 34.98
      eunc = eunc / 34.98
               ;convert to e
      eflux = eflux *sbtoe
      eunc = eunc *sbtoe
      print, 'mean SNR', mean(eflux), mean(eflux)/mean(eunc)

      if i eq 1 and total eq 0 then begin
         xarr = x_center
         yarr = y_center
         xsarr = xs
         ysarr = ys
         total = total + 1
      endif else begin
         xarr = [xarr, x_center]
         yarr = [yarr, y_center]
         xsarr = [xsarr, xs]
         ysarr = [ysarr, ys]
         total = total + 1
      endelse
      
   endfor                       ; for each fits file in the AOR
   
endfor                          ;for each AOR
  ; print, 'xarrs', xsarr
   print, 'nel' ,n_elements(xsarr), total

;plot up a histogram
plothist, xsarr, xsxhist, xsyhist, bin=0.0001, /noprint,/noplot
plothist, ysarr, ysxhist, ysyhist, bin=0.0001, /noprint,/noplot

b = barplot(ysxhist, ysyhist, nbars = 2, index = 0,title = 'ch1', xtitle = 'Y Centroid Uncertainty', ytitle = 'Number',fill_color = 'blue',/ylog, name = 'sub', xrange = [0.004, 0.017])


;-----------------------------------------------------
;now overplot the full array results
;-----------------------------------------------------

full_aor = ['0044258560','0044258816','0044259072','0044258048','0044257792','0044259584','0044259840','0044259328','0044257536','0044258304']
yarr = fltarr(n_elements(full_aor)*319)
xarr = fltarr(n_elements(full_aor)*319)
ysarr =fltarr(n_elements(full_aor)*319)
xsarr =fltarr(n_elements(full_aor)*319)
efluxarr =fltarr(n_elements(full_aor)*319)
euncarr =fltarr(n_elements(full_aor)*319)
total = 0

for a = 0, n_elements(full_aor) - 1 do begin
   print, 'working on ', full_aor(a)
   if a lt 8 then dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(full_aor(a)) 
   if a ge 8 then dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028800/bcd/' + string(full_aor(a)) 

   CD, dir                      ; change directories to the correct AOR directory
   command  =  " ls *bcd_fp.fits > /Users/jkrick/irac_warm/pcrsfaint/cbcdlist.txt"
   spawn, command
   
   readcol,'/Users/jkrick/irac_warm/pcrsfaint/cbcdlist.txt',fitsname, format = 'A', /silent
   
   for i =1, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
 
      if i eq 1 then begin
         header = headfits(fitsname(i)) ;
         ra_ref = sxpar(header, 'RA_REF')
         dec_ref = sxpar(header, 'DEC_REF')
         gain = sxpar(header, 'GAIN')
         exptime =  sxpar(header, 'EXPTIME') 
         fluxconv =  sxpar(header, 'FLUXCONV') 
         print, 'ra, dec', ra_ref, dec_ref, gain, exptime, fluxconv
      endif
      
      get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra_ref, dec = dec_ref, /silent     
;figure out flux in electrons
      sbtoe = gain*exptime/fluxconv
               ;convert to microJy
      eflux = abcdflux *1E6
      eunc = fs*1E6
               ;convert to Mjy/sr
      eflux = eflux / 34.98
      eunc = eunc / 34.98
               ;convert to e
      eflux = eflux *sbtoe
      eunc = eunc *sbtoe
      efluxarr[total] = eflux
      euncarr[total] = eunc
      ;print, 'mean SNR', eflux/ eunc

      if i eq 1 and total eq 0 then begin
         xarr = x_center
         yarr = y_center
         xsarr = xs
         ysarr = ys
         total = total + 1
      endif else begin
         xarr = [xarr, x_center]
         yarr = [yarr, y_center]
         xsarr = [xsarr, xs]
         ysarr = [ysarr, ys]
         total = total + 1
      endelse
      
   endfor                       ; for each fits file in the AOR
   
endfor                          ;for each AOR
  ; print, 'xarrs', xsarr
   print, 'nel' ,n_elements(xsarr), total

;plot up a histogram
plothist, xsarr, xsxhist, xsyhist, bin=0.0001, /noprint,/noplot
plothist, ysarr, ysxhist, ysyhist, bin=0.0001, /noprint,/noplot

b3 = barplot(ysxhist, ysyhist,nbars = 2, index = 0,fill_color = 'coral',/overplot, name = 'full')
print, 'mean SNR ',mean(efluxarr), mean(efluxarr)/mean(euncarr)
;-----------------------------------------------------
;now overplot the full array to sub array results
;-----------------------------------------------------

full_sub_aor = ['0044260096','0044260352','0044261632','0044261376','0044262400','0044261888','0044260864','0044261120','0044260608','0044262144']

yarr = fltarr(n_elements(full_sub_aor)*319)
xarr = fltarr(n_elements(full_sub_aor)*319)
ysarr =fltarr(n_elements(full_sub_aor)*319)
xsarr =fltarr(n_elements(full_sub_aor)*319)
total = 0

for a = 0, n_elements(full_sub_aor) - 1 do begin
   print, 'working on ', full_sub_aor(a)
   if a lt 6 then dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(full_sub_aor(a)) 
   if a ge 6 then dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028800/bcd/' + string(full_sub_aor(a)) 

   CD, dir                      ; change directories to the correct AOR directory
   command  =  " ls *bcd_fp.fits > /Users/jkrick/irac_warm/pcrsfaint/cbcdlist.txt"
   spawn, command
   
   readcol,'/Users/jkrick/irac_warm/pcrsfaint/cbcdlist.txt',fitsname, format = 'A', /silent
   
   for i =1, n_elements(fitsname) - 1 do begin ;read each cbcd file, find centroid, keep track
 
      if i eq 1 then begin
         header = headfits(fitsname(i)) ;
         ra_ref = sxpar(header, 'RA_REF')
         dec_ref = sxpar(header, 'DEC_REF')
         print, 'ra, dec', ra_ref, dec_ref
      endif
      
      get_centroids,fitsname(i), t, dt, x_center, y_center, abcdflux, xs, ys, fs, b, /WARM, /APER, APRAD = [5],ra = ra_ref, dec = dec_ref, /silent     
      if i eq 1 and total eq 0 then begin
         xarr = x_center
         yarr = y_center
         xsarr = xs
         ysarr = ys
         total = total + 1
      endif else begin
         xarr = [xarr, x_center]
         yarr = [yarr, y_center]
         xsarr = [xsarr, xs]
         ysarr = [ysarr, ys]
         total = total + 1
      endelse
      
   endfor                       ; for each fits file in the AOR
   
endfor                          ;for each AOR
  ; print, 'xarrs', xsarr
   print, 'nel' ,n_elements(xsarr), total

;plot up a histogram
plothist, xsarr, xsxhist, xsyhist, bin=0.0001, /noprint,/noplot
plothist, ysarr, ysxhist, ysyhist, bin=0.0001, /noprint,/noplot
for j = 0, n_elements(xshist) - 1 do begin
   print, 'xs', xshist[j], yshist[j]
endfor

b4 = barplot(ysxhist, ysyhist,nbars = 2, index = 1,fill_color = 'orange',/overplot, name = 'full_sub')

;---------------------------------
;and time to make the legend
;---------------------------------
;l = legend(target = [b, b2, b3, b4], position = [.0045, 1000], /data)
t1 = text(0.0125, 1000, 'sub', color = 'blue',/data)
t3 = text(0.009, 600, 'full', color = 'coral',/data)
t4 = text(0.009, 1000, 'full_sub', color = 'orange',/data)

end
