pro plot_firstframe

 ; ps_start, filename = '/Users/jkrick/virgo/irac/ffeffect.ps'
  ;;ps_open, filename= '/Users/jkrick/virgo/irac/ffeffect.ps',/portrait,/square,/color

  ;redcolor = FSC_COLOR("Red", !D.Table_Size-2)
  ;bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
  
  
  restore,  '/Users/jkrick/virgo/irac/ndark_ch1_pc15.sav'
  print, 'ch1 pc15*********************'
  meanarr = meanarr + 0.3
;sort timearr
  c = timearr[sort(timearr)]
  d = meanarr[sort(timearr)]
  
;plot, c, d, psym = 2
  x = findgen(n_elements(c))
;plot, x, d, psym = 2, xrange = [0,3500], xstyle = 1, yrange = [-0.2, 0.1], xtitle = 'frame number', ytitle = 'mean levels(Mjy/sr)'
  
;-----------
  plot, delayarr, meanarr , psym = 6, xtitle = 'Delay Time (s)', ytitle = 'Mean Levels (arbitrary units)', yrange = [0, 0.3], xrange = [5,20], thick = 3, xthick = 3, ythick = 3, charthick = 3
  
;find mean per each bin
  
  h = histogram(delayarr, max = 18, bin = 0.4, reverse_indices=ri, locations= x)
  meanbin = fltarr(n_elements(h))
  sigmabin = fltarr(n_elements(h))

  for j=0,n_elements(meanbin) - 1  do begin
     if ri[j+1] gt ri[j] then begin
     meanbin[j]=median(meanarr[ri[ri[j]:ri[j+1]-1]])
;     sigmabin[j]=stddev(meanarr[ri[ri[j]:ri[j+1]-1]],/nan)
     sigmabin[j] = robust_sigma(meanarr[ri[ri[j]:ri[j+1]-1]])
  endif
  endfor
print, 'meanbin', meanbin
print, 'sigmabin', sigmabin
  oplot, x, meanbin , psym = 2, color = bluecolor, thick = 3
  
;need to fit this ff effect with a function
  xn = x
  start = [1.0,1.0]
  noise = sigmabin
 ; noise[*] = 1                                         ;equally weight the values
  result= MPFITFUN('logarithm',xn,meanbin, noise, start, perror = osigmaerr) ; fit a function
  xnn = findgen(13) + 7
  oplot, xnn, result(0) + result(1)*alog10(xnn), thick = 3;, color = redcolor
  
  print, 'test', 10, result(0) + result(1)*alog10(10)
  print, 'test', 19, result(0) + result(1)*alog10(19)
  print, 'perror', osigmaerr

;NORMALIZE the effect to delay time of 19 = no correction.
;so I want to add (19) - (x) to the value of x to get the corrected value
  
  
  save, result, filename = '/Users/jkrick/virgo/irac/ffeffect_pc15.sav'
;---------------------------------------------------------------------

  
  restore,  '/Users/jkrick/virgo/irac/ndark_ch2_pc15.sav'
  print, 'ch2 pc15*********************'

  meanarr = meanarr + 0.3
;sort timearr
  c = timearr[sort(timearr)]
  d = meanarr[sort(timearr)]
  
;plot, c, d, psym = 2
  x = findgen(n_elements(c))
;  oplot, x, d, psym = 2, color = redcolor ;, yrange = [-0.1, 0.1]
  
;-----------
  oplot, delayarr, meanarr, psym = 2 ,thick = 3;, color = redcolor
  
  h = histogram(delayarr, max = 20, bin = 0.4, reverse_indices=ri, locations= x)
  meanbin = fltarr(n_elements(h))
  sigmabin = fltarr(n_elements(h))

  for j=0,n_elements(meanbin) - 1  do begin
     if ri[j+1] gt ri[j] then begin
     meanbin[j]=median(meanarr[ri[ri[j]:ri[j+1]-1]])
     sigmabin[j]=stddev(meanarr[ri[ri[j]:ri[j+1]-1]])

  endif
  endfor
  
;  oplot, x, meanbin+0.3, psym = 2, color = whitecolor
  
  

;---------------------------------------------------------------------
;--------------------------------------------------------------------


  restore,  '/Users/jkrick/virgo/irac/ndark_ch1_pc16.sav'
  print, 'ch1 pc16*********************'

  meanarr = meanarr + 0.3
;sort timearr
  c = timearr[sort(timearr)]
  d = meanarr[sort(timearr)]
  
;plot, c, d, psym = 2
  x = findgen(n_elements(c))
;plot, x, d, psym = 2, xrange = [0,3500], xstyle = 1, yrange = [-0.2, 0.1], xtitle = 'frame number', ytitle = 'mean levels(Mjy/sr)'
  
;-----------
  oplot, delayarr, meanarr , psym = 5, thick =3;,color = redcolor
  
;find mean per each bin
  
  h = histogram(delayarr, max = 18, bin = 0.4, reverse_indices=ri, locations= x)
  meanbin = fltarr(n_elements(h))
  sigmabin = fltarr(n_elements(h))

  for j=0,n_elements(meanbin) - 1  do begin
     if ri[j+1] gt ri[j] then begin
     meanbin[j]=median(meanarr[ri[ri[j]:ri[j+1]-1]])
     ;sigmabin[j]=stddev(meanarr[ri[ri[j]:ri[j+1]-1]],/nan)
     sigmabin[j] = robust_sigma(meanarr[ri[ri[j]:ri[j+1]-1]])

  endif
  endfor
print, 'meanbin', meanbin
print, 'sigmabin', sigmabin
  
  oplot, x, meanbin , psym = 2, color = bluecolor, thick = 3
  
;need to fit this ff effect with a function
  xn = x
  start = [1.0,1.0]
;  noise = fltarr(n_elements(meanbin))
; noise[*] = 1                                         ;equally weight the values

  noise = sigmabin
  result= MPFITFUN('logarithm',xn,meanbin, noise, start, perror = osigmaerr) ; fit a function
  xnn = findgen(13) + 7
  oplot, xnn, result(0) + result(1)*alog10(xnn), thick = 3;, color = redcolor
  
  print, 'test', 10, result(0) + result(1)*alog10(10)
  print, 'test', 19, result(0) + result(1)*alog10(19)
  print, 'perror', osigmaerr
;NORMALIZE the effect to delay time of 19 = no correction.
;so I want to add (19) - (x) to the value of x to get the corrected value
  
  
  save, result, filename = '/Users/jkrick/virgo/irac/ffeffect_pc16.sav'
;---------------------------------------------------------------------
  
  restore,  '/Users/jkrick/virgo/irac/ndark_ch2_pc16.sav'
  print, 'ch2 pc16*********************'

  meanarr = meanarr + 0.3
;sort timearr
  c = timearr[sort(timearr)]
  d = meanarr[sort(timearr)]
  
;plot, c, d, psym = 2
  x = findgen(n_elements(c))
;  oplot, x, d, psym = 2, color = redcolor ;, yrange = [-0.1, 0.1]
  
;-----------
  oplot, delayarr, meanarr, psym = 1 , thick = 3;, color = redcolor
  
  h = histogram(delayarr, max = 20, bin = 0.4, reverse_indices=ri, locations= x)
  meanbin = fltarr(n_elements(h))
  for j=0,n_elements(meanbin) - 1  do if ri[j+1] gt ri[j] then $
     meanbin[j]=mean(meanarr[ri[ri[j]:ri[j+1]-1]])
  
;  oplot, x, meanbin+0.3, psym = 2, color = whitecolor
 

  ;ps_end
  ;;ps_close, /noprint,/noid

end
