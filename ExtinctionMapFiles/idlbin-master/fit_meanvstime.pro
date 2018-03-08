pro fit_meanvstime
;ps_start, filename= '/Users/jkrick/virgo/irac/meanvstime_ch2.ps',/nomatch

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
yellowcolor = FSC_COLOR("yellow", !D.Table_Size-9)
colorarr = [redcolor, orangecolor, yellowcolor, greencolor, cyancolor, bluecolor, purplecolor]

!P.multi = [0,1,1]
!P.Thick = 3
!P.CharThick = 3
!P.Charsize = 1.25 
!X.Thick = 3
!Y.Thick = 3

;restore,  '/Volumes/irac_drive/virgo/irac/ffeffect.sav'
restore,  '/Users/jkrick/virgo/irac/ffeffect.sav'
at19 =  result(0) + result(1)*alog10(19)

  dirloc = '/Users/jkrick/virgo/irac/'
 
 aorname = ['r35320064/','r35320320/','r35320576/','r35320832/','r35321088/','r35321344/','r35321600/','r35321856/','r35322112/','r35322368/','r35322624/','r35322880/','r35323136/','r35323392/','r35323648/','r35323904/','r35324160/','r35324416/','r35324672/','r35324928/','r35325184/','r35325440/','r35325696/','r35325952/','r35326208/']

sortaor = ['r35322368/','r35325952/','r35321856/','r35324928/','r35326208/','r35322112/','r35325440/','r35321344/','r35324416/','r35325696/','r35321600/','r35324672/','r35320832/','r35323904/','r35325184/','r35321088/','r35323648/','r35320576/','r35323136/','r35324160/','r35320320/','r35322880/','r35320064/','r35322624/','r35323392/']

tile1aor = ['r35322368/','r35322112/','r35321600/','r35321088/','r35320320/']

tile2aor = ['r35325952/','r35325440/','r35324672/','r35323648/','r35322880/']

tile3aor = ['r35321856/','r35321344/','r35320832/','r35320576/','r35320064/']

tile4aor = ['r35324928/','r35324416/','r35323904/','r35323136/','r35322624/']

tile5aor = ['r35326208/','r35325696/','r35325184/','r35324160/','r35323392/']


for ch = 0, 0 do begin
   levelarr = fltarr(n_elements(sortaor))
   result1arr = fltarr(n_elements(sortaor))
   result2arr = fltarr(n_elements(sortaor))
   keeptimearr = fltarr(n_elements(sortaor))
   errorarr = fltarr(n_elements(sortaor))

   for a =0 , n_elements(sortaor) - 1 do begin
      meanarr = fltarr(5000)
      sigmaarr = fltarr(5000)
      timearr = fltarr(5000)
      zodiarr = fltarr(5000)
      delayarr = fltarr(5000)
      d = 0
   
      cd, strcompress(dirloc + sortaor(a) + 'ch' + string(ch + 1) + '/bcd',/remove_all)
      spawn, 'pwd'
;      if ch eq 0 then command =  "find . -name 'SPITZER_I1*0_ndark_ff_match.fits' > ndark_ff_ch1.txt"
      if ch eq 0 then command =  "find . -name 'corr_SPITZER_I1*0_ndark_ff.fits' > ndark_ff_ch1.txt"
      if ch eq 1 then command =  "find . -name 'corr_SPITZER_I2*0_ndark.fits' > ndark_ff_ch2.txt"
      spawn, command
      
      ;read in a list of all ndark and ff corrected images
      readcol,strcompress( 'ndark_ff_ch'+ string((ch+1))+'.txt',/remove_all), fitsname, format = 'A';, /silent
           
      for i =3, n_elements(fitsname) - 1 do begin
         ;run SExtractor
         ;spawn, 'which sex'
         cmd = 'sex ' + fitsname(i) + ' -c ../../../default.sex'
         spawn, cmd
         fits_read, fitsname(i), data, header
         fits_read, 'segmentation.fits', segdata, segheader

         ;write out image with all objects turned to NaN's
         m = where(segdata gt 0)
         data(m) = alog10(-1)   ;set the objects to NaN's.
         ;newname = strcompress(strmid(fitsname(i), 0, 32) + 'back.fits',/remove_all)
         ;fits_write, newname, data, segheader

        ;keep track of mean background levels, and time of
        ;the exposures.  
        timearr(d)  = sxpar(header, 'SCLK_OBS')
        m = mean(data,/nan)
        s = robust_sigma(data)
;        if ch eq 0 and m gt 0.2 then begin
;           meanarr(d) = m 
;           sigmaarr(d) = s
;        endif else begin
;           meanarr(d) = meanarr(d-1)
;           sigmaarr(d) = sigmaarr(d-1)
;        endelse
        
        ;ch2
        meanarr(d) = m
        sigmaarr(d) = s

;        meanarr(d) = median(data)
        delayarr(d)  = sxpar(header, 'FRAMEDLY')
        zodiarr(d) = sxpar(header, 'ZODY_EST')

        ;print, 'mean', meanarr(d)
        d = d + 1

      endfor   ;for each fits image
;
      timearr = timearr[0:d-1]
      delayarr=delayarr[0:d-1]
      meanarr = meanarr[0:d-1]
      zodiarr = zodiarr[0:d-1]
      sigmaarr = sigmaarr[0:d-1]

      print, 'zodii', zodiarr

      xt = findgen(n_elements(timearr))
      if a eq 0 then begin
         plot, xt, meanarr,  xtitle = 'Frame Number', ytitle = 'Mean Levels' ,yrange = [0. , 0.2]; yrange = [0.22,0.28];, color = colorarr[a] 
         oplot, xt, zodiarr -0.24 ; +0.12   -0.24
      endif

      if a gt 0 and a lt 7 then begin
         oplot, xt, meanarr;, color = bluecolor; colorarr[a]
         oplot, xt, zodiarr -0.24;+0.12 - 0.24;, color = bluecolor; colorarr[a]
      endif

      if a ge 7 then begin
         oplot, xt, meanarr, color = bluecolor
         oplot, xt, zodiarr -0.24  , color = bluecolor ;+0.12 -0.24
      endif

      noise =sigmaarr
      ;noise[*] = 1              ;equally weight the values
       
;ch2
      if a eq 0 or a eq 5 or a eq 10 or a eq 15 or a eq 20 then begin
         noise[0:10] = 1E6
         noise[60:80] = 1E6
         noise[35:55] = 1E6
      endif
     
      if a eq 1 or a eq 6 or a eq 11 or a eq 16 or a eq 21 then begin
         noise[55:75] = 1E6
         noise[90:110] = 1E6
         noise[10:20] = 1E6
      endif
     
      if a eq 4 or a eq 9 or a eq 14 or a eq 19 or a eq 24 then begin
         noise[15:55] = 1E6
      endif
      
;ch1

;      if a eq 4 or a eq 9 or a eq 14 or a eq 19 or a eq 24 then begin
;         print, 'working on tile 5'
;         noise[8:28] = 1E6
;         noise[39:58] = 1E6
;         noise[80:100] = 1E6
;      endif
      
;      if a eq 1 or a eq 6 or a eq 11 or a eq 16 or a eq 21 then begin  ;tile 2
;         noise[60:90] = 1E6
;     endif

      start = [meanarr(0),0.01,0.01]
 ;     pi = replicate({fixed:0, limited:[0,0], limits:[0.D,0.D]},3)
 ;     pi(0).limited(0) = 1
 ;     pi(0).limits(1) = 0.0
 ;     pi(1).limited(0) = 1
 ;     pi(1).limits(1) = 0.0
 ;     pi(2).limited(0) = 1
 ;     pi(2).limits(1) = 0.0
      result= MPFITFUN('exponential',xt,meanarr, noise, start, perror = osigmaerr);, parinfo = pi) ; fit a function
;      oplot, xt,  result(0) - result(1)*exp(-result(2)*xt), color = bluecolor
      
                               ;keep track of the levels
                               ;later will want to add flux to match levels
      levelarr[a] = result(0)
      result1arr[a] = result(1)
      result2arr[a] = result(2)
      keeptimearr[a] = n_elements(timearr)
      errorarr[a] = osigmaerr(0)

   endfor                       ;for each AOR
   
   
  ; print, 'levelarr', levelarr
   
;try adding flux to get rid of differences
;match to the first one (always the highest I hope, will want to
;check this)

;XXXXthink about how I want to do this
;make it a function or run seperately
;just be careful to only run it once, and also to test that it worked.
 
                                ;save these to work with them later
                                ;instead of re-running SExtractor
;save, /variables, filename = '/Users/jkrick/virgo/irac/fitmean_ch1.sav'

print, 'errorarr  ', errorarr
print, 'errorarr in percent ',errorarr/ levelarr
endfor  ;for each channel

;ps_end

end
