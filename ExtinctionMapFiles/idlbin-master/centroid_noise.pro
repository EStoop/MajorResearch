pro centroid_noise
; this code compares the centroiding uncertainties from the different pcrspeakup tests taken after the
;switch to iru 2 in fall 2011

;-----------------------------------------------------
;now overplot the vw/no_vw tests from wk 509
;-----------------------------------------------------

;second test in Sep

inc_c_aor = ['0044254976','0044255232','0044255488','0044255744','0044256000','0044256256','0044256512','0044256768','0044257024','0044257280']

   yarr = fltarr(n_elements(inc_c_aor)*64*16)
   xarr = fltarr(n_elements(inc_c_aor)*64*16)
   ysarr =fltarr(n_elements(inc_c_aor)*64*16)
   xsarr =fltarr(n_elements(inc_c_aor)*64*16)
   total = 0

for a = 0, n_elements(inc_c_aor) - 1 do begin
   print, 'working on ', inc_c_aor(a)
   dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(inc_c_aor(a)) 
   
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
         print, 'ra, dec', ra_ref, dec_ref
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
      print, 'mean eflux', mean(eflux)

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

b = barplot(ysxhist, ysyhist, nbars = 2, index = 0,title = 'ch2', xtitle = 'Y Centroid Uncertainty', ytitle = 'Number',fill_color = 'blue',/ylog, name = 'inc_c')

;-----------------------------------------------------
;now overplot the vw/no_vw tests from wk 509
;-----------------------------------------------------

vw_aor = ['0044169472','0044173056','0044173312','0044169728','0044169984','0044173568','0044173824','0044170240','0044170496','0044174080','0044174336','0044170752','0044171008','0044174592','0044174848','0044171264','0044171520','0044175104','0044175360','0044171776','0044172032','0044175616','0044175872','0044172288','0044172544','0044176128','0044176384','0044172800']

yarr = fltarr(n_elements(vw_aor)*64*16)
xarr = fltarr(n_elements(vw_aor)*64*16)
ysarr =fltarr(n_elements(vw_aor)*64*16)
xsarr =fltarr(n_elements(vw_aor)*64*16)
total = 0

for a = 0, n_elements(vw_aor) - 1 do begin
   print, 'working on ', vw_aor(a)
   dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(vw_aor(a)) 
   
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

b2 = barplot(ysxhist, ysyhist,nbars = 2, index = 1,fill_color = 'aqua',/overplot, name = 'vw')

;-----------------------------------------------------
;now overplot the full array results
;-----------------------------------------------------

full_aor = ['0044266496','0044265472','0044266240','0044265984','0044265216','0044265728','0044267008','0044267264','0044266752','0044267520']

yarr = fltarr(n_elements(full_aor)*319)
xarr = fltarr(n_elements(full_aor)*319)
ysarr =fltarr(n_elements(full_aor)*319)
xsarr =fltarr(n_elements(full_aor)*319)
efluxarr =fltarr(n_elements(full_aor)*319)
total = 0

for a = 0, n_elements(full_aor) - 1 do begin
   print, 'working on ', full_aor(a)
   if a lt 6 then dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(full_aor(a)) 
   if a ge 6 then dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028800/bcd/' + string(full_aor(a)) 

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
         print, 'ra, dec', ra_ref, dec_ref
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
     ; print, 'mean SNR', eflux/ eunc

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

b3 = barplot(ysxhist, ysyhist,nbars = 2, index = 0,fill_color = 'coral',/overplot, name = 'full')
print, 'mean eflux ', mean(efluxarr)
;-----------------------------------------------------
;now overplot the full array to sub array results
;-----------------------------------------------------

full_sub_aor = ['0044268032','0044267776','0044268800','0044268544','0044268288','0044270080','0044269824','0044269312','0044269056','0044269568']
yarr = fltarr(n_elements(full_sub_aor)*319)
xarr = fltarr(n_elements(full_sub_aor)*319)
ysarr =fltarr(n_elements(full_sub_aor)*319)
xsarr =fltarr(n_elements(full_sub_aor)*319)
total = 0

for a = 0, n_elements(full_sub_aor) - 1 do begin
   print, 'working on ', full_sub_aor(a)
   if a lt 5 then dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028700/bcd/' + string(full_sub_aor(a)) 
   if a ge 5 then dir = '/Users/jkrick/iracdata/flight/IWIC/IRAC028800/bcd/' + string(full_sub_aor(a)) 

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
t1 = text(0.006, 1000, 'inc_c', color = 'blue',/data)
t2 = text(0.006, 2000, 'vw', color = 'aqua',/data)
t3 = text(0.0135, 2000, 'full', color = 'coral',/data)
t4 = text(0.0135, 1000, 'full_sub', color = 'orange',/data)

end
