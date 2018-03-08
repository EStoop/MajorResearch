pro plot_checkffcorr
  redcolor = FSC_COLOR("Red");, !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue");, !D.Table_Size-3)

ps_start, filename= '/Users/jkrick/virgo/irac/check_ffcorr.ps',/nomatch
!P.Thick = 3
!P.CharThick = 3
!P.Charsize = 1.25 
!X.Thick = 3
!Y.Thick = 3

!P.multi = [0,1,1]
 redcolor = FSC_COLOR("Red");, !D.Table_Size-2)
  bluecolor = FSC_COLOR("Blue");, !D.Table_Size-3)

  ;restore,  '/Users/jkrick/virgo/irac/ndark_ch1.sav'
  restore,  '/Users/jkrick/virgo/irac/ndark_ff_check_ch1.sav'

a = where(meanarr ne 0)
  meanarr = meanarr[a]
  timearr = timearr[a]
  delayarr = delayarr[a]

;plot, c, d, psym = 2
 ; x = findgen(n_elements(c))
;plot, x, d, psym = 2, xrange = [0,3500], xstyle = 1, yrange = [-0.2, 0.1], xtitle = 'frame number', ytitle = 'mean levels(Mjy/sr)'
  
;-----------
  plot, delayarr, meanarr , psym = 2, xtitle = 'Delay Time (s)', ytitle = 'Mean Levels (arbitrary units)', yrange = [0, 0.3], xrange = [5,20]
;print, 'meanarr', meanarr
;find mean per each bin
  
  h = histogram(delayarr, max = 20, bin = 0.4, reverse_indices=ri, locations= x)
  meanbin = fltarr(n_elements(h))
  for j=0,n_elements(meanbin) - 1  do if ri[j+1] gt ri[j] then $
     meanbin[j]=mean(meanarr[ri[ri[j]:ri[j+1]-1]])
  
  oplot, x, meanbin, psym = 6, color = bluecolor
  xtop = findgen(20)
  ytop = fltarr(20) + meanbin(30)
;oplot, xtop, ytop, color = bluecolor
  
;---------------------------------------------------------------------
 ;restore,  '/Users/jkrick/virgo/irac/ndark_ch2.sav'
  restore,  '/Users/jkrick/virgo/irac/ndark_ff_check_ch1.sav'

a = where(meanarr ne 0)
  meanarr = meanarr[a]
  timearr = timearr[a]
  delayarr = delayarr[a]

;plot, c, d, psym = 2
 ; x = findgen(n_elements(c))
;plot, x, d, psym = 2, xrange = [0,3500], xstyle = 1, yrange = [-0.2, 0.1], xtitle = 'frame number', ytitle = 'mean levels(Mjy/sr)'
  
;-----------
  oplot, delayarr, meanarr, psym = 4
;print, 'meanarr', meanarr
;find mean per each bin
  
  h = histogram(delayarr, max = 20, bin = 0.4, reverse_indices=ri, locations= x)
  meanbin = fltarr(n_elements(h))
  for j=0,n_elements(meanbin) - 1  do if ri[j+1] gt ri[j] then $
     meanbin[j]=mean(meanarr[ri[ri[j]:ri[j+1]-1]])
  
  oplot, x, meanbin, psym = 6, color = bluecolor
  xtop = findgen(20)
  ytop = fltarr(20) + meanbin(30)
;oplot, xtop, ytop, color = bluecolor
  
;---------------------------------------------------------------------
;--------------------------------------------------------------------
    ;------
;plot the old values

;  oldtimearr = oldtimearr[0:n_elements(oldmeanarr)-1]
 
;sort timearr  
;  c = oldtimearr[sort(oldtimearr)]
;  d = oldmeanarr[sort(oldtimearr)]
  
;  x = findgen(n_elements(c))
;  plot, x, d, psym = 2, xrange = [900,1900], xstyle = 1, yrange = [0.2, 0.3], xtitle = 'frame number', ytitle = 'mean levels(Mjy/sr)';xrange=[0,3500]
  
     ;------
;overplot the new values

;timearr = timearr[0:n_elements(meanarr)-1]
;plot, (timearr - timearr(0)) / 3600, meanarr, psym = 2, xtitle = 'time(hrs)', ytitle = 'mean levels(Mjy/sr)', yrange = [0.1, 0.3], xrange = [-500, 100], xstyle = 1
  
;sort timearr
  
;  e1 = timearr[sort(timearr)]
;  f1 = meanarr[sort(timearr)]
  
  ;plot, (c - c[0]) / 3600, d, psym = 2, yrange = [0.2, 0.4], ytitle = 'mean levels', xtitle = 'time';, xrange = [40,70], yrange =[0.1, 0.3]
;  x = findgen(n_elements(e1))
;  oplot, x, f1, psym = 6, color = bluecolor
  

;what are those spikes
;s = where(meanarr gt 0.29)
;print, meanarr(s), delayarr(s), timearr(s)

;---------------------------------------------------------------------
;--------------------------------------------------------------------
;restore, '/Users/jkrick/virgo/irac/zodi.sav'
;oplot, x, sortzodi + 0.1, psym = 7, color = redcolor
;print, 'timearr', n_elements(timearr), e1
;print, 'oldtimearr',n_elements(oldtimearr) ,c
;print, 'meanarr', meanarr

ps_end

end
