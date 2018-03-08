pro darkvstime_cryo
ps_open, filename='/Users/jkrick/irac_darks/darksvstime_cryo.ps',/portrait,/square,/color
!P.thick = 3
!P.charthick = 3
!P.multi = [0,1,4]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)

;aordir = '/Users/jkrick/iracdata/flight/IWIC/S18.18.0/darks/'
;dirlist =[ 'IRAC022500/''IRAC022600/','IRAC022700/','IRAC022800/','IRAC022900/','IRAC023000/','IRAC023100/','IRAC023200/','IRAC023300/','IRAC023400/','IRAC023500/','IRAC023600/','IRAC023700/','IRAC023800/','IRAC023900/','IRAC024000/','IRAC024100/','IRAC024200/','IRAC024300/']

aordir = '/Users/jkrick/iracdata/flight/calprog/finaldarks/'

cd, aordir
command = "find . -name '*skydark.fits' > /Users/jkrick/irac_darks/cryodarks.txt"
;spawn, command
print, 'finished making file'

readcol, '/Users/jkrick/irac_darks/cryodarks.txt', fitsname, format = 'A'
ch1_sclkarr = fltarr(n_elements(fitsname))
ch1_meanarr = fltarr(n_elements(fitsname))
ch1_sigmaarr = fltarr(n_elements(fitsname))
ch1_modearr = fltarr(n_elements(fitsname))
ch1_zodiarr = fltarr(n_elements(fitsname))
j = 0

ch2_sclkarr = fltarr(n_elements(fitsname))
ch2_meanarr = fltarr(n_elements(fitsname))
ch2_sigmaarr = fltarr(n_elements(fitsname))
ch2_modearr = fltarr(n_elements(fitsname))
ch2_zodiarr = fltarr(n_elements(fitsname))
k = 0

ch3_sclkarr = fltarr(n_elements(fitsname))
ch3_meanarr = fltarr(n_elements(fitsname))
ch3_sigmaarr = fltarr(n_elements(fitsname))
ch3_modearr = fltarr(n_elements(fitsname))
ch3_zodiarr = fltarr(n_elements(fitsname))
l = 0

ch4_sclkarr = fltarr(n_elements(fitsname))
ch4_meanarr = fltarr(n_elements(fitsname))
ch4_sigmaarr = fltarr(n_elements(fitsname))
ch4_modearr = fltarr(n_elements(fitsname))
ch4_zodiarr = fltarr(n_elements(fitsname))
m = 0

flux_conv = [.1088,.1388, 0.5952, 0.2021]

for a = long(0),  n_elements(fitsname) - 1 do begin
   header = headfits(fitsname(a))
   hdrmode= sxpar(header, 'HDRMODE')
   exptype = sxpar(header, 'EXPTYPE')
   framtime = sxpar(header, 'FRAMTIME')
   ch = sxpar(header, 'CHNLNUM')
   ;caltype = sxpar(header, 'CALTYPE')
   time = sxpar(header, 'SCLK_CAL')
   campaign = sxpar(header, 'CAMPAIGN')
   aorkey = sxpar(header, 'AORKEY')
   ;print, fitsname(a), hdrmode, exptype, framtime

   ;track just the 30s full array darks 
;   if hdrmode eq 0 and framtime eq 30. then begin
   if hdrmode eq 0 and framtime gt 30. and framtime lt 150 and time gt  753985131. and campaign ne 'IRAC011700' then begin 
      if ch eq 4 then print, 'working on ', campaign, aorkey
      fits_read, fitsname(a), data, header, exten_no = 0

      ;need to ignore the uncertainty file
      ;whichis attached to the data fits
      ;file, but not a different extension.
;      z = where(data lt 4, count)
;      if count gt 0 then  data(z) = alog10(-1)
      data = data[*,*,0]

      ;try the zodi
      zodi = sxpar(header, 'ZODY_EST') ; in Mjy/sr

      if ch eq 2 and zodi gt 0.26 then begin
         print, 'high zodi', campaign, zodi
      endif

      ;print, 'data before', data[99], flux_conv(ch-1), framtime
      data = data * flux_conv(ch - 1) / framtime
      ;print, 'data after', data[99], flux_conv(ch-1), framtime
      ;extend = [0.994, 0.937, .772,.737]
      extend = [1.0, 1.0, .772,.737]
 ;     data = data - zodi/extend(ch-1) ;correcting for extended source


  
      ;use a gaussian to determine the background value
      start = [mean(data,/nan), stddev(data,/nan), 4000.0]
      binsize = .1           
      plothist, data, xhist, yhist, bin = binsize, /nan,  yrange = [0, 10000],/noplot;make a histogram of the data
      noise = fltarr(n_elements(yhist))
      noise[*] = 1   ;equally weight the                             
      result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) ; fit a gaussian to the histogram sorted data
                                ;oplot, xhist, result(2)/sqrt(2.*!Pi) * exp(-0.5*((xhist - result(0))/result(1))^2.), linestyle = 2,  color = redcolor
      
     ;use mmm to get the modes, are accurate measures of the peak of the distribution
      mmm, data, mmmmode, mmmsigma, mmmskew 

      if ch eq 1 then begin
         ch1_sclkarr(j) = sxpar(header, 'SCLK_CAL')
         ch1_meanarr(j) = result(0)
         ch1_sigmaarr(j) = result(1)/ sqrt(n_elements(xhist))
         ch1_modearr(j) = mmmmode
         ch1_zodiarr(j) = zodi
         j = j + 1
      endif
      if ch eq 2 then begin
         ch2_sclkarr(k) = sxpar(header, 'SCLK_CAL')
         ch2_meanarr(k) = result(0)
         ch2_sigmaarr(k) = result(1)/ sqrt(n_elements(xhist))
         ch2_modearr(k) = mmmmode
         ch2_zodiarr(k) = zodi
;         print, 'zodi', zodi, campaign
         k = k + 1
      endif
      if ch eq 3 then begin
         ch3_sclkarr(l) = sxpar(header, 'SCLK_CAL')
         ;print, 'ch3', fitsname(a), sxpar(header, 'SCLK_CAL'), result(0), result(1)
         ch3_meanarr(l) = result(0)
         ch3_sigmaarr(l) = result(1)/ sqrt(n_elements(xhist))
         ch3_modearr(l) = mmmmode
         ch3_zodiarr(l) = zodi
         l = l + 1
      endif
      if ch eq 4 then begin
         ch4_sclkarr(m) = sxpar(header, 'SCLK_CAL')
         ch4_meanarr(m) = result(0)
         ch4_sigmaarr(m) = result(1)/ sqrt(n_elements(xhist))
         ch4_modearr(m) = mmmmode
         ch4_zodiarr(m) = zodi
        m = m + 1
      endif
      

   endif

endfor

ch1_sclkarr = ch1_sclkarr[0:j-1]
ch1_meanarr = ch1_meanarr[0:j-1]
ch1_simgaarr = ch1_sigmaarr[0:j-1]
ch1_modearr = ch1_modearr[0:j-1]
ch1_zodiarr = ch1_zodiarr[0:j-1]
ch1_sclksort = ch1_sclkarr[sort(ch1_sclkarr)]
ch1_meansort = ch1_meanarr[sort(ch1_sclkarr)]
ch1_sigmasort = ch1_sigmaarr[sort(ch1_sclkarr)]
ch1_modesort = ch1_modearr[sort(ch1_sclkarr)]
ch1_zodisort = ch1_zodiarr[sort(ch1_sclkarr)]

ch2_sclkarr = ch2_sclkarr[0:k-1]
ch2_meanarr = ch2_meanarr[0:k-1]
ch2_simgaarr = ch2_sigmaarr[0:k-1]
ch2_modearr = ch2_modearr[0:k-1]
ch2_zodiarr = ch2_zodiarr[0:k-1]
ch2_sclksort = ch2_sclkarr[sort(ch2_sclkarr)]
ch2_meansort = ch2_meanarr[sort(ch2_sclkarr)]
ch2_sigmasort = ch2_sigmaarr[sort(ch2_sclkarr)]
ch2_modesort = ch2_modearr[sort(ch2_sclkarr)]
ch2_zodisort = ch2_zodiarr[sort(ch2_sclkarr)]
;print, 'ch2 zodi', ch2_zodisort

ch3_sclkarr = ch3_sclkarr[0:l-1]
ch3_meanarr = ch3_meanarr[0:l-1]
ch3_simgaarr = ch3_sigmaarr[0:l-1]
ch3_modearr = ch3_modearr[0:l-1]
ch3_zodiarr = ch3_zodiarr[0:l-1]
ch3_sclksort = ch3_sclkarr[sort(ch3_sclkarr)]
ch3_meansort = ch3_meanarr[sort(ch3_sclkarr)]
ch3_sigmasort = ch3_sigmaarr[sort(ch3_sclkarr)]
ch3_modesort = ch3_modearr[sort(ch3_sclkarr)]
ch3_zodisort = ch3_zodiarr[sort(ch3_sclkarr)]
;print, 'n_elements', n_elements(ch3_meansort), n_elements(ch3_sclksort)
;print, 'ch3_meanarr', ch3_meanarr

ch4_sclkarr = ch4_sclkarr[0:m-1]
ch4_meanarr = ch4_meanarr[0:m-1]
ch4_simgaarr = ch4_sigmaarr[0:m-1]
ch4_zodiarr = ch4_zodiarr[0:m-1]
ch4_sclksort = ch4_sclkarr[sort(ch4_sclkarr)]
ch4_meansort = ch4_meanarr[sort(ch4_sclkarr)]
ch4_sigmasort = ch4_sigmaarr[sort(ch4_sclkarr)]
ch4_modesort = ch4_modearr[sort(ch4_sclkarr)]
ch4_zodisort = ch4_zodiarr[sort(ch4_sclkarr)]


;ploterror, (ch1_sclksort - ch1_sclksort(0)) / 86400, ch1_meansort / ch1_meansort(1), ch1_sigmasort/ ch1_meansort(1), psym = 2, ytitle = 'Normalized flux', xtitle = 'time in days from pc4 - 22', title = 'Channel 1 F30s'
;ploterror, (ch2_sclksort - ch2_sclksort(0)) / 86400, ch2_meansort / ch2_meansort(1), ch2_sigmasort/ ch2_meansort(1), psym = 4, ytitle = 'Normalized flux', xtitle = 'time in days from pc4 - 22', title = 'Channel 2 F30s'


;plot, (ch1_sclksort - ch1_sclksort(0)) / 86400, ch1_meansort/ ch1_meansort(1) , psym = 2, ytitle = 'Normalized flux', xtitle = 'Time in days starting from IRAC', title = 'Channel 1 F30s', yrange = [0.9,1.5]
;plot, (ch2_sclksort - ch2_sclksort(0)) / 86400, ch2_meansort / ch2_meansort(1),  psym = 4, ytitle = 'Normalized flux', xtitle = 'Time in days starting from pc5', title = 'Channel 2 F30s', yrange = [0.9,1.5]
;plot, (ch3_sclksort - ch3_sclksort(0)) / 86400, ch3_meansort / ch3_meansort(1),  psym = 4, ytitle = 'Normalized flux', xtitle = 'Time in days starting from pc5', title = 'Channel 3 F30s';, yrange = [0.8,1.2]
;plot, (ch4_sclksort - ch4_sclksort(0)) / 86400, ch4_meansort / ch4_meansort(1),  psym = 4, ytitle = 'Normalized flux', xtitle = 'Time in days starting from pc5', title = 'Channel 4 F30s';, yrange = [0.8,1.2]


ps_close, /noprint,/noid

;save, /variables, filename = '/Users/jkrick/irac_darks/cryodarks_100_nozodisub.sav'

end
