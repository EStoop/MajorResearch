pro plot_calstar_vs_time

ps_open, filename='/Users/jkrick/irac_warm/calstars/calstar_vs_time_ch1_test.ps',/portrait,/square,/color


!P.multi = [0,2,3]
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
name = ['1812095_cal', 'HD165459', 'NPM1p60','HD184837','HD156896' ]

restore, '/Users/jkrick/idlbin/calobject.sav'

yr_0 = [6,12]
yr_1 = [600, 650]
yr_2 = [34,44]

for star = 0, 2 do begin

   n = where(calobject[star].flux ne 0 );and calobject[star].flux lt 9.0)
   m = mean(calobject[star].flux(n))
   meanclip, calobject[star].flux(n), m, sigmaflux, clipsig = 1.5

   plothist, calobject[star].flux(n) / m, xhist, yhist, bin = 0.01, xrange = [0.9,1.1], charsize = 1.5, title = name(star) + strcompress(string(m) + 'mJy') + strcompress(string(mean(calobject[star].flux_e(n))) + 'e'), xtitle = 'mean flux', ytitle = 'number', thick = 3
   start = [1.0,0.01, 20.]
   noise = fltarr(n_elements(yhist))
   noise[*] = 1                                                         ;equally weight the values
   result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet)         ;./quiet    

;plot the fitted gaussian and print results to plot
   xarr = findgen(1000)/ 100.      
   oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
          linestyle =  2 , thick = 3      
   xyouts, 1.00, 9, result(1)
   
;now what about the gaussian fits to the corrected flux values
   plothist, calobject[star].flux_pixelcorr(n) / m, xhist, yhist, bin = 0.01, xrange = [0.9,1.1],/noplot
   noise = fltarr(n_elements(yhist))
   noise[*] = 1                                                 ;equally weight the values
   result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet) ;./quiet    
   oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
          linestyle =  2 , color = greencolor, thick  = 3
   xyouts, 1.00, 8, result(1), color = greencolor

   plothist, calobject[star].flux_arraycorr(n) / m, xhist, yhist, bin = 0.01, xrange = [0.9,1.1],/noplot
   noise = fltarr(n_elements(yhist))
   noise[*] = 1                                                 ;equally weight the values
   result= MPFITFUN('mygauss',xhist,yhist, noise, start,/quiet)         ;./quiet    
   oplot, xarr, (result(2))/sqrt(2.*!Pi) * exp(-0.5*((xarr - (result(0)))/(result(1)))^2.), $
          linestyle =  2  , color = redcolor, thick = 3
   xyouts, 1.00, 7, result(1), color = redcolor


;print out the photon noise estimate
   xyouts, 0.9, 9, mean(calobject[star].percent_e(n))


   a = where (calobject[star].pos eq 1,acount)
   b = where (calobject[star].pos eq 2,bcount)
   c = where (calobject[star].pos eq 3,ccount)
   d = where (calobject[star].pos eq 4,dcount)
   e = where (calobject[star].pos eq 5,ecount)
   
   print,'counts', acount, bcount, ccount, dcount, ecount
;   case star of
;      0: yr =yr_0
;      1:yr = yr_1
;      2: yr = yr_2
;   endcase
;   
   plot, calobject[star].time(a), calobject[star].flux(a) / m, psym =2, xrange = [9.32E8, 9.5E8], xstyle = 1, yrange = [0.9,1.1], charsize = 1.5, ytitle = 'normalized flux', xtitle = 'time pc5 - pc14'
   oplot, calobject[star].time(b), calobject[star].flux(b) / m, psym =4
   oplot, calobject[star].time(c), calobject[star].flux(c) / m, psym =5
   oplot, calobject[star].time(d), calobject[star].flux(d) / m, psym =6
   oplot, calobject[star].time(e), calobject[star].flux(e) / m, psym =1
 
   oplot, calobject[star].time(a), calobject[star].flux_pixelcorr(a) / m, psym =2,  color = greencolor
;   plot, calobject[star].time(a), calobject[star].flux_pixelcorr(a) / m, psym =2, xrange = [9.32E8, 9.46E8], xstyle = 1, yrange = [0.9,1.1], charsize = 1.5, color = greencolor, ytitle = 'normalized flux with pixel corrr', xtitle = 'time pc5 - pc11'
   oplot, calobject[star].time(b), calobject[star].flux_pixelcorr(b) / m, psym =4, color = greencolor
   oplot, calobject[star].time(c), calobject[star].flux_pixelcorr(c) / m, psym =5, color = greencolor
   oplot, calobject[star].time(d), calobject[star].flux_pixelcorr(d) / m, psym =6, color = greencolor
   oplot, calobject[star].time(e), calobject[star].flux_pixelcorr(e) / m, psym =1, color = greencolor
 
   oplot, calobject[star].time(a), calobject[star].flux_arraycorr(a) / m, psym =2,  color = redcolor
;   plot, calobject[star].time(a), calobject[star].flux_arraycorr(a) / m, psym =2, xrange = [9.32E8, 9.46E8], xstyle = 1, yrange = [0.9,1.1], charsize = 1.5, color = redcolor, ytitle = 'normalized flux with pixel & array corr', xtitle = 'time pc5 - pc11'
   oplot, calobject[star].time(b), calobject[star].flux_arraycorr(b) / m, psym =4, color = redcolor
   oplot, calobject[star].time(c), calobject[star].flux_arraycorr(c) / m, psym =5, color = redcolor
   oplot, calobject[star].time(d), calobject[star].flux_arraycorr(d) / m, psym =6, color = redcolor
   oplot, calobject[star].time(e), calobject[star].flux_arraycorr(e) / m, psym =1, color = redcolor
 
endfor

ps_close, /noprint,/noid

end
