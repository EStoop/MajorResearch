pro bias_final
ps_open, filename='/Users/jkrick/iwic/start_science2.ps',/portrait,/square,/color
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
blackcolor = FSC_COLOR("black", !D.Table_Size-9)

!P.multi = [0,2,1]
readcol, '/Users/jkrick/iwic/dark_summary/dark_summary_08_18.csv', name,directory,reckey,date,exptime,fowler,waitper,nframes,channel,commandedK,corr_cernoxK,cernoxK,vreset,bias,vdduc,vgg1,new_mean,new_sigma,old_mean,old_sigma,old_noisy,noisypix,flux_conv,latent_strength,latent_decay,format='A,A,I, I, I, I, I,I,I,F10.2,F10.2,F10.2,F10.2,I,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,I, F10.5,F10.2,F10.5,F10.5'

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;temp1
;waiting for corrected diode temps
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
!P.multi = [0,2,2]

ch1_e12b75= where(channel eq 1 and bias eq 750 and exptime eq 12 and new_mean ne 0)
ch2_e12b50= where(channel eq 2 and bias eq 500 and exptime eq 12 and new_mean ne 0)
ch1_e100b75= where(channel eq 1 and bias eq 750 and exptime eq 100 and new_mean ne 0)
ch2_e100b50= where(channel eq 2 and bias eq 500 and exptime eq 100 and new_mean ne 0)

for i = 0, n_elements(ch1_e12b75) - 1 do print, name(ch1_e12b75(i)), corr_cernoxK(ch1_e12b75(i)), new_mean(ch1_e12b75(i)), bias(ch1_e12b75(i)), vgg1(ch1_e12b75(i)), vreset(ch1_e12b75(i)), vdduc(ch1_e12b75(i))
;mean noise

x=cernoxK(ch1_e12b75)
y = new_mean(ch1_e12b75)
plot, x[sort(x)], y[sort(x)] , psym = 2,/ylog, xtitle ='Cernox Temp (K)', $
      ytitle = 'Mean Noise Level (MJy/sr)', title = 'Ch1 Bias=750, Ch2 Bias = 500', xrange=[10,45], yrange=[0.01,1],$
      thick=3, xthick=3, ythick=3, charthick=3
oplot, x[sort(x)], y[sort(x)] 
xyouts, 15.02, .031072, 'C', charthick=3

;plot, corr_cernoxK(ch2_e12b50), new_mean(ch2_e12b50), psym = 2,/ylog, xtitle ='Corrected Cernox Temp (K)', $
;      ytitle = 'Mean Noise Level (MJy/sr)', title = 'Ch2 bias 500', xrange=[10,40], yrange=[0.01,1]
;xyouts, 15.08 ,.04295, 'C'

x = cernoxK(ch2_e12b50)
y = new_mean(ch2_e12b50)
oplot, x[sort(x)], y[sort(x)] , psym = 2, color = bluecolor, thick=3
oplot, x[sort(x)], y[sort(x)] ,  color = bluecolor
xyouts, 15.08 ,.04295, 'C', color = bluecolor, charthick=3

legend, ['Ch1', 'Ch2'], linestyle = [0,0], color = [blackcolor, bluecolor],/left,/top, thick=3, charthick=3

;noisy pix
x = cernoxK(ch1_e100b75)
y = noisypix(ch1_e100b75)
plot, x[sort(x)], y[sort(x)], psym = 2,/ylog, xtitle ='Cernox Temp (K)', $
      ytitle = 'Fraction of Noisy Pixels', title = 'Ch1 Bias=750, Ch2 Bias = 500', xrange=[10,50], yrange=[0.0001,1],$
      thick=3, xthick=3, ythick=3, charthick=3
oplot, x[sort(x)], y[sort(x)]

;oplot, [3,800], [12./(256.*256.), 12./(256.*256.)]
;xyouts, 25, 12./(256.*256.), 'CR in 12s'
;oplot, [3,800], [100./(256.*256.), 100./(256.*256.)]
;xyouts, 25, 100./(256.*256.), 'CR in 100s'
;oplot, [3,800], [200./(256.*256.), 200./(256.*256.)]
;xyouts, 25, 200./(256.*256.), 'CR in 200s'

;plot, corr_cernoxK(ch2_e100b50), noisypix(ch2_e100b50), psym = 2,/ylog, xtitle ='Corrected Cernox Temp (K)', $
;      ytitle = 'Fraction of Noisy Pixels', title = 'Ch2 bias 500',
;      xrange=[10,40], yrange=[0.0001,1]

x =  cernoxK(ch2_e100b50)
y = noisypix(ch2_e100b50)
oplot, x[sort(x)], y[sort(x)],  psym = 2, color = bluecolor, thick=3
oplot, x[sort(x)], y[sort(x)],  color = bluecolor
oplot, [3,800], [12./(256.*256.), 12./(256.*256.)]
xyouts, 12, 12./(256.*256.), 'CR in 12s', charthick=3
oplot, [3,800], [100./(256.*256.), 100./(256.*256.)]
xyouts, 12, 100./(256.*256.), 'CR in 100s', charthick=3
oplot, [3,800], [200./(256.*256.), 200./(256.*256.)]
xyouts, 12, 200./(256.*256.), 'CR in 200s', charthick=3

legend, ['Ch1', 'Ch2'], linestyle = [0,0], color = [blackcolor, bluecolor],/left,/top, thick=3, charthick=3

;latent strength
;a = where(latent_strength(ch2_e12b50) ne 0)
;x = corr_cernoxK(ch2_e12b50(a))
;y = latent_strength(ch2_e12b50(a))
;plot,x[sort(x)], y[sort(x)] , psym = 2, xtitle ='Corrected Cernox Temp (K)', $
;      ytitle = 'Fractional Latent Strength', title = 'Bias 500', xrange=[10,40],$
;      thick=3, xthick=3, ythick=3, charthick=3
;oplot,x[sort(x)], y[sort(x)] , psym = 2, color=bluecolor, thick=3
;oplot,x[sort(x)], y[sort(x)] , color = bluecolor

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;bias
;-----------------------------------------------------------------------------
;---------------------------------------------------------------------------
!P.multi = [0,2,2]

;mean
;---------------------------------------------------------------------------

;ch1 t30
ch1_e2t30 = where(channel eq 1 and corr_cernoxK eq 30.1 and exptime eq 2 and new_mean ne 0  and vgg1 gt -4 and vgg1 lt -3.6)
ch2_e2t30 = where(channel eq 2 and corr_cernoxK eq 30.1 and exptime eq 2 and new_mean ne 0 )
ch1_e12t30 = where(channel eq 1 and corr_cernoxK eq 30.1 and exptime eq 12 and new_mean ne 0  and vgg1 gt -4 and vgg1 lt -3.6)
ch2_e12t30 = where(channel eq 2 and corr_cernoxK eq 30.1 and exptime eq 12 and new_mean ne 0 )
ch1_e100t30 = where(channel eq 1 and corr_cernoxK eq 30.1 and exptime eq 100 and new_mean ne 0  )
ch2_e100t30 = where(channel eq 2 and corr_cernoxK eq 30.1 and exptime eq 100 and new_mean ne 0 )

ch1_e2t30 = where(channel eq 1 and cernoxK eq 29.4 and exptime eq 2 and new_mean ne 0  and vgg1 gt -4 and vgg1 lt -3.6)
ch2_e2t30 = where(channel eq 2 and cernoxK eq 29.9 and exptime eq 2 and new_mean ne 0 )
ch1_e12t30 = where(channel eq 1 and cernoxK eq 29.4 and exptime eq 12 and new_mean ne 0  and vgg1 gt -4 and vgg1 lt -3.6)
ch2_e12t30 = where(channel eq 2 and cernoxK eq 29.9 and exptime eq 12 and new_mean ne 0 )
ch1_e100t30 = where(channel eq 1 and cernoxK eq 29.4 and exptime eq 100 and new_mean ne 0  )
ch2_e100t30 = where(channel eq 2 and cernoxK eq 29.9 and exptime eq 100 and new_mean ne 0 )

for i = 0, n_elements(ch1_e12t30) - 1 do print, 'bias', bias(ch1_e12t30(i)),  new_mean(ch1_e12t30(i)), name(ch1_e12t30(i)), vgg1(ch1_e12t30(i)), vreset(ch1_e12t30(i)), vdduc(ch1_e12t30(i))


x=bias(ch1_e12t30)
y = new_mean(ch1_e12t30)
y = y[sort(x)]
x= x[sort(x)]
plot,  x, y, psym = 2, xrange=[300,700], xtitle ='Bias',thick=3, charthick=3, xthick=3, ythick=3,$
      ytitle = 'Mean Noise Level (MJy/sr)', title = ' 12s', yrange=[0.04,.1]
oplot, x, y

x=bias(ch1_e100t30)
y = new_mean(ch1_e100t30)
oplot,  x[sort(x)], y[sort(x)], psym = 4, thick=3
oplot,  x[sort(x)], y[sort(x)]

x=bias(ch1_e2t30)
y=new_mean(ch1_e2t30)
oplot,  x[sort(x)], y[sort(x)], psym = 6, thick=3
oplot,  x[sort(x)], y[sort(x)]
;--------------------------------------
;ch1 t29

ch1_e100t29 = where(channel eq 1 and corr_cernoxK eq 29.9 and exptime eq 100 and new_mean ne 0  )
ch2_e100t29 = where(channel eq 2 and corr_cernoxK eq 31 and exptime eq 100 and new_mean ne 0 )
ch1_e12t29 = where(channel eq 1 and corr_cernoxK eq 29.9 and exptime eq 12 and new_mean ne 0  )
ch2_e12t29 = where(channel eq 2 and corr_cernoxK eq 31 and exptime eq 12 and new_mean ne 0 )
ch1_e2t29 = where(channel eq 1 and corr_cernoxK eq 29.9 and exptime eq 2 and new_mean ne 0  )
ch2_e2t29 = where(channel eq 2 and corr_cernoxK eq 31 and exptime eq 2 and new_mean ne 0 )


x = bias(ch1_e100t29)
y = new_mean(ch1_e100t29)
oplot,  x[sort(x)], y[sort(x)], psym = 6, thick=3;,color=redcolor
oplot,  x[sort(x)], y[sort(x)];, color = redcolor

x=bias(ch1_e2t29)
y=new_mean(ch1_e2t29)
oplot,  x[sort(x)], y[sort(x)], psym = 6, color = redcolor, thick=3
oplot,  x[sort(x)], y[sort(x)], color = redcolor

x=bias(ch1_e12t29)
y=new_mean(ch1_e12t29)
oplot,  x[sort(x)], y[sort(x)], psym = 6, thick=3;, color = redcolor
oplot,  x[sort(x)], y[sort(x)];, color = redcolor


;legend, ['2s', '12s','100s'], psym=[6,2,4],/right,/top
;legend, ['30K', '29K'], linestyle = [0,0], color = [blackcolor, redcolor],/right,/top
;--------------------------------------
;ch2 t30

x=bias(ch2_e12t30)
y = new_mean(ch2_e12t30)
y = y[sort(x)]
x= x[sort(x)]
oplot,  x, y, psym = 2, color=bluecolor, thick=3;, xrange=[300,800], xtitle ='Bias',thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Mean Noise Level (MJy/sr)', title = 'Ch2 12s', yrange=[0.04,.06]
oplot, x, y, color=bluecolor

x=bias(ch2_e100t30)
y = new_mean(ch2_e100t30)
oplot,  x[sort(x)], y[sort(x)], psym = 4, thick=3
oplot,  x[sort(x)], y[sort(x)]

x=bias(ch2_e2t30)
y=new_mean(ch2_e2t30)
oplot,  x[sort(x)], y[sort(x)], psym = 6, thick=3
oplot,  x[sort(x)], y[sort(x)]


;ch2 t29
x = bias(ch2_e100t29)
y = new_mean(ch2_e100t29)
oplot,  x[sort(x)], y[sort(x)], psym = 4, color = redcolor, thick=3
oplot,  x[sort(x)], y[sort(x)], color = redcolor

x=bias(ch2_e2t29)
y=new_mean(ch2_e2t29)
oplot,  x[sort(x)], y[sort(x)], psym = 6, color = redcolor, thick=3
oplot,  x[sort(x)], y[sort(x)], color = redcolor

x=bias(ch2_e12t29)
y=new_mean(ch2_e12t29)
oplot,  x[sort(x)], y[sort(x)], psym = 6, color = bluecolor, thick=3; redcolor
oplot,  x[sort(x)], y[sort(x)], color = bluecolor; redcolor


;legend, ['2s', '12s','100s'], psym=[6,2,4],/right,/top
legend, ['30K', '29K'], psym = [2,6],/right,/top, thick=3, charthick=3
legend, ['Ch1', 'Ch2'], linestyle = [0,0], color = [blackcolor, bluecolor],/left,/top, thick=3, charthick=3


;noisypix
;---------------------------------------------------------------------------
;ch1
x=bias(ch1_e100t30)
y=noisypix(ch1_e100t30)
plot,x[sort(x)], y[sort(x)], psym = 2,/ylog, xtitle ='Bias', thick=3, charthick=3, xthick=3, ythick=3,$
      ytitle = 'Fraction of Noisy Pixels', xrange=[300,800], yrange=[0.0001,1]
oplot,x[sort(x)], y[sort(x)]

x=bias(ch1_e100t29)
y=noisypix(ch1_e100t29)
oplot, x[sort(x)], y[sort(x)] , psym = 6, thick=3;, color = bluecolor;redcolor
oplot, x[sort(x)], y[sort(x)] ;,  color = bluecolor; redcolor

oplot, [300,800], [12./(256.*256.), 12./(256.*256.)]
xyouts, 600, 12./(256.*256.), 'CR in 12s', charthick=3
oplot, [300,800], [100./(256.*256.), 100./(256.*256.)]
xyouts, 600, 100./(256.*256.), 'CR in 100s', charthick=3
oplot, [300,800], [200./(256.*256.), 200./(256.*256.)]
xyouts, 600, 200./(256.*256.), 'CR in 200s', charthick=3
legend, ['Ch1', 'Ch2'], linestyle = [0,0], color = [blackcolor, bluecolor],/left,/top, thick=3, charthick=3

;ch2
x=bias(ch2_e100t30)
y=noisypix(ch2_e100t30)
oplot,x[sort(x)], y[sort(x)], psym = 2, color=bluecolor, thick=3;,/ylog, xtitle ='Bias', thick=3, charthick=3, xthick=3, ythick=3,$
 ;     ytitle = 'Fraction of Noisy Pixels', title = 'Ch2 ', xrange=[300,800], yrange=[0.0001,1]
oplot,x[sort(x)], y[sort(x)], color=bluecolor

x=bias(ch2_e100t29)
y=noisypix(ch2_e100t29)
oplot, x[sort(x)], y[sort(x)] , psym = 6, color = bluecolor, thick=3
oplot, x[sort(x)], y[sort(x)] ,  color = bluecolor

oplot, [300,800], [12./(256.*256.), 12./(256.*256.)]
xyouts, 600, 12./(256.*256.), 'CR in 12s', charthick=3
oplot, [300,800], [100./(256.*256.), 100./(256.*256.)]
xyouts, 600, 100./(256.*256.), 'CR in 100s', charthick=3
oplot, [300,800], [200./(256.*256.), 200./(256.*256.)]
xyouts, 600, 200./(256.*256.), 'CR in 200s', charthick=3
legend, ['30K', '29K'], psym = [2,6],/right,/top, thick=3, charthick=3

;latent strength
;---------------------------------------------------------------------------

;ch1
;x=bias(ch1_e100t30)
;y=latent_strength(ch1_e100t30)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Bias', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Latent Strength', title = 'Ch1 ', xrange=[300,800]
;oplot,x[sort(x)], y[sort(x)]

;x=bias(ch1_e100t29)
;y=latent_strength(ch1_e100t29)
;oplot, x[sort(x)], y[sort(x)] , psym = 2, color = redcolor
;oplot, x[sort(x)], y[sort(x)] ,  color = redcolor

;legend, ['30K', '29K'], linestyle = [0,0], color = [blackcolor, redcolor],/right,/top

;ch2

;a = where(latent_strength(ch2_e100t30) ne 0)
;x=bias(ch2_e100t30(a))
;y=latent_strength(ch2_e100t30(a))
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Bias', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Fractional Latent Strength',  xrange=[300,800]
;oplot,x[sort(x)], y[sort(x)], psym = 2, color = bluecolor, thick=3
;oplot,x[sort(x)], y[sort(x)], color = bluecolor

;x=bias(ch2_e100t29)
;y=latent_strength(ch2_e100t29)
;oplot, x[sort(x)], y[sort(x)] , psym = 2, color = redcolor
;oplot, x[sort(x)], y[sort(x)] ,  color = redcolor

;legend, ['30K', '29K'], linestyle = [0,0], color = [blackcolor, redcolor],/right,/top


;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;vgg1
;-----------------------------------------------------------------------------
;---------------------------------------------------------------------------
!P.multi = [0,2,2]

ch1_e12t30b45= where(channel eq 1 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 30.1)
ch2_e12t30b45= where(channel eq 2 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 30.1)
ch1_e100t30b45= where(channel eq 1 and bias eq 450 and exptime eq 100 and new_mean ne 0 and corr_cernoxK eq 30.1)
ch2_e100t30b45= where(channel eq 2 and bias eq 450 and exptime eq 100 and new_mean ne 0 and corr_cernoxK eq 30.1)

ch1_e12t29b45= where(channel eq 1 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 29.2)
ch2_e12t29b45= where(channel eq 2 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 29.0)
ch1_e100t29b45= where(channel eq 1 and bias eq 450 and exptime eq 100 and new_mean ne 0 and corr_cernoxK eq 29.2)
ch2_e100t29b45= where(channel eq 2 and bias eq 450 and exptime eq 100 and new_mean ne 0 and corr_cernoxK eq 29.0)


;for i = 0, n_elements(ch1_e12t30b45) - 1 do print, 'vgg1', name(ch1_e12t30b45(i)), bias(ch1_e12t30b45(i)), vgg1(ch1_e12t30;45(i)), vreset(ch1_e12t30b45(i)), vdduc(ch1_e12t30b45(i))
;for i = 0, n_elements(ch2_e12t30b45) - 1 do print, 'vgg1 2', name(ch2_e12t30b45(i)), bias(ch2_e12t30b45(i)), vgg1(ch2_e12t30b45(i)), vreset(ch2_e12t30b45(i)), vdduc(ch2_e12t30b45(i))

;mean
;---------------------------------------------------------------------------
;ch1
x = vgg1(ch1_e12t30b45)
y = new_mean(ch1_e12t30b45)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Vgg1', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Mean Noise Level (MJy/sr)', title = '12s b450', xrange=[-4.2,-3.5], yrange=[.04,.12]
;oplot,x[sort(x)], y[sort(x)]

;x = vgg1(ch1_e100t30b45)
;y = new_mean(ch1_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 6
;oplot,x[sort(x)], y[sort(x)]

x = vgg1(ch1_e12t29b45)
y = new_mean(ch1_e12t29b45)
;oplot,x[sort(x)], y[sort(x)], psym = 6, thick=3;, color = redcolor
;oplot,x[sort(x)], y[sort(x)];,  color = redcolor

;legend, ['30K', '29K'], psym = [2,6],/right,/top, thick=3, charthick=3
;legend, ['Ch1', 'Ch2'], linestyle = [0,0], color = [blackcolor, bluecolor],/left,/top, thick=3, charthick=3

;ch2
x = vgg1(ch2_e12t30b45)
y = new_mean(ch2_e12t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 2, color=bluecolor, thick=3;, xtitle ='Vgg1', thick=3, charthick=3, xthick=3, ythick=3,charsize=1.5,$
;      ytitle = 'Mean Noise Level (MJy/sr)', title = 'Ch2 12s b450', xrange=[-4.2,-3.5], yrange=[.04,.06]
;oplot,x[sort(x)], y[sort(x)], color=bluecolor

;x = vgg1(ch2_e100t30b45)
;y = new_mean(ch2_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 6
;oplot,x[sort(x)], y[sort(x)]

x = vgg1(ch2_e12t29b45)
y = new_mean(ch2_e12t29b45)
;oplot,x[sort(x)], y[sort(x)], psym = 6, color = bluecolor, thick=3
oplot,x[sort(x)], y[sort(x)],  color = bluecolor



;noisypix
;---------------------------------------------------------------------------
;ch1
x = vgg1(ch1_e100t30b45)
y = noisypix(ch1_e100t30b45)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Vgg1', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Fraction of Noisy Pixels', title = 'Bias 450', xrange=[-4.2,-3.5], /ylog, yrange=[.0001,.01]
;oplot,x[sort(x)], y[sort(x)]

x = vgg1(ch1_e100t29b45)
y = noisypix(ch1_e100t29b45)
;oplot,x[sort(x)], y[sort(x)], psym = 6, thick=3;, color = redcolor
;oplot,x[sort(x)], y[sort(x)];,  color = redcolor

;legend, ['30K', '29K'], psym = [2,6],/right,/top, thick=3, charthick=3
;legend, ['Ch1', 'Ch2'], linestyle = [0,0], color = [blackcolor, bluecolor],/left,/top, thick=3, charthick=3

;ch2
x = vgg1(ch2_e100t30b45)
y = noisypix(ch2_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 2, color=bluecolor, thick=3;xtitle ='Vgg1', thick=3, charthick=3, xthick=3, ythick=3,charsize=1.5,$
;      ytitle = 'Fraction of Noisy Pixels', title = 'Ch2  b450', xrange=[-4.2,-3.5], /ylog, yrange=[.0001,.01]
;oplot,x[sort(x)], y[sort(x)], color=bluecolor

x = vgg1(ch2_e100t29b45)
y = noisypix(ch2_e100t29b45)
;oplot,x[sort(x)], y[sort(x)], psym = 6, color = bluecolor, thick=3
;oplot,x[sort(x)], y[sort(x)],  color = bluecolor

;oplot, [-4.2,-3.4], [12./(256.*256.), 12./(256.*256.)]
;xyouts, -4.1, 12./(256.*256.), 'CR in 12s', charthick=3
;oplot, [-4.2,-3.4], [100./(256.*256.), 100./(256.*256.)]
;xyouts, -4.1, 100./(256.*256.), 'CR in 100s', charthick=3
;oplot, [-4.2,-3.4], [200./(256.*256.), 200./(256.*256.)]
;xyouts, -4.1, 200./(256.*256.), 'CR in 200s', charthick=3


;latent_strength
;---------------------------------------------------------------------------
;ch1
;x = vgg1(ch1_e100t30b45)
;y = latent_strength(ch1_e100t30b45)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Vgg1', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Fractional Latent Strength', title = 'Ch1  t30K b450', xrange=[-4.2,-3.5];, /ylog, yrange=[.001,.01]
;oplot,x[sort(x)], y[sort(x)]

;x = vgg1(ch1_e100t29b45)
;y = latent_strength(ch1_e100t29b45)
;oplot,x[sort(x)], y[sort(x)], psym = 2, color = redcolor
;oplot,x[sort(x)], y[sort(x)],  color = redcolor

;legend, ['30K', '29K'], linestyle = [0,0], color = [blackcolor, redcolor],/left,/top

;ch2
x = vgg1(ch2_e100t30b45)
y = latent_strength(ch2_e100t30b45)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Vgg1', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Fractional Latent Strength', title = 'Bias 450', xrange=[-4.2,-3.5];, /ylog, yrange=[.001,.01]
;oplot,x[sort(x)], y[sort(x)], color=bluecolor
;oplot,x[sort(x)], y[sort(x)], color=bluecolor, psym=2, thick=3

;x = vgg1(ch2_e100t29b45)
;y = latent_strength(ch2_e100t29b45)
;oplot,x[sort(x)], y[sort(x)], psym = 2, color = redcolor
;oplot,x[sort(x)], y[sort(x)],  color = redcolor

;legend, ['30K', '29K'], linestyle = [0,0], color = [blackcolor, redcolor],/left,/top

;legend, ['Ch1', 'Ch2'], linestyle = [0,0], color = [blackcolor, bluecolor],/left,/top

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;vreset 
;-----------------------------------------------------------------------------
;---------------------------------------------------------------------------
!P.multi = [0,2,2]

ch1_different1= where(channel eq 1 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 30.1 and vgg1 eq -3.65 and vreset eq -3.6 and vdduc eq -3.5)
ch1_different2= where(channel eq 1 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 30.1 and vgg1 eq -3.65 and vreset eq -3.4 and vdduc eq -3.5)
ch1_same1= where(channel eq 1 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 30.1 and vgg1 eq -3.65 and vreset eq -3.5 and vdduc eq - 3.5)
;ch1_same2= where(channel eq 1 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 30.1 and vgg1 eq -3.68 and vreset eq -3.55 and vdduc eq - 3.55)
ch2_different1= where(channel eq 2 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 30.1 and vgg1 eq -3.6 and vreset eq -3.4 and vdduc eq -3.5)
ch2_same1= where(channel eq 2 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 30.1 and vgg1 eq -3.6 and vreset eq -3.5 and vdduc eq -3.5)

ch1_e12t30b45= where(channel eq 1 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 30.1  and vgg1 eq -3.65 and vdduc eq -3.5)
ch2_e12t30b45= where(channel eq 2 and bias eq 450 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 30.1  and vgg1 eq -3.6 and vdduc eq -3.5)
ch1_e100t30b45= where(channel eq 1 and bias eq 450 and exptime eq 100 and new_mean ne 0 and corr_cernoxK eq 30.1  and vgg1 eq -3.65 and vdduc eq -3.5)
ch2_e100t30b45= where(channel eq 2 and bias eq 450 and exptime eq 100 and new_mean ne 0 and corr_cernoxK eq 30.1  and vgg1 eq -3.6 and vdduc eq -3.5)

;for i = 0, n_elements(ch1_e12t30b45) - 1 do print, 'r1', name(ch1_e12t30b45(i)), bias(ch1_e12t30b45(i)), vgg1(ch1_e12t30b45(i)), vreset(ch1_e12t30b45(i)), vdduc(ch1_e12t30b45(i))
;for i = 0, n_elements(ch2_e12t30b45) - 1 do print, 'r 2', name(ch2_e12t30b45(i)), bias(ch2_e12t30b45(i)), vgg1(ch2_e12t30b45(i)), vreset(ch2_e12t30b45(i)), vdduc(ch2_e12t30b45(i))

;mean
;---------------------------------------------------------------------------
;ch1

x = vreset(ch1_e12t30b45)
y = new_mean(ch1_e12t30b45)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Vreset', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Mean Noise Level (MJy/sr)', title = 'Ch1 12s b450 g365 d35', xrange=[-3.7, -3.3], yrange=[0.03,0.1]
;oplot,x[sort(x)], y[sort(x)]

;x = vreset(ch1_e100t30b45)
;y = new_mean(ch1_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 6
;oplot,x[sort(x)], y[sort(x)]

;legend, ['12s', '100s'], psym = [2,6],/left,/center

;ch2
x = vreset(ch2_e12t30b45)
y = new_mean(ch2_e12t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 2, color=bluecolor, thick=3;, xtitle ='Vreset', thick=3, charthick=3, xthick=3, ythick=3,charsize=1.5,$
;;      ytitle = 'Mean Noise Level (MJy/sr)', title = 'Ch2 12s b450 g36 d35', xrange=[-3.7, -3.3]
;oplot,x[sort(x)], y[sort(x)], color=bluecolor

;x = vreset(ch2_e100t30b45)
;y = new_mean(ch2_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 6
;oplot,x[sort(x)], y[sort(x)]


;legend, ['12s', '100s'], psym = [2,6],/left,/center
;legend, ['Ch1', 'Ch2'], linestyle = [0,0], color = [blackcolor, bluecolor],/left,/top, thick=3, charthick=3



;x =[0,0,1,1]
;y1 = [new_mean(ch1_different1), new_mean(ch1_different2), new_mean(ch1_same1)]
;plot,x, y1, psym = 2, xtitle ='different    same', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Mean Noise Level (MJy/sr)', title = 'Ch1 12s b450', xrange= [-0.1, 1.1], xstyle = 1
;oplot,x,y1

;ch2
;x =[0,1,1]
;y1 = [new_mean(ch2_different1), new_mean(ch2_same1)]
;plot,x, y1, psym = 2, xtitle ='different    same', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Mean Noise Level (MJy/sr)', title = 'Ch1 12s b450', xrange= [-0.1, 1.1], xstyle = 1
;oplot,x,y1


;noisypix
;---------------------------------------------------------------------------
;ch1
x = vreset(ch1_e100t30b45)
y = noisypix(ch1_e100t30b45)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Vreset', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Fraction of Noisy Pixels', title = 'b450 g36 d35', /ylog, yrange=[0.001,0.01], xrange=[-3.7, -3.3]
;oplot,x[sort(x)], y[sort(x)]


;ch2
x = vreset(ch2_e100t30b45)
y = noisypix(ch2_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 2, color=bluecolor, thick=3; xtitle ='Vreset', thick=3, charthick=3, xthick=3, ythick=3,charsize=1.5,$
;;      ytitle = 'Fraction of Noisy Pixels', title = 'Ch2 b450 g36 d35',/ylog, yrange=[0.001,0.01], xrange=[-3.7, -3.3]
;oplot,x[sort(x)], y[sort(x)], color=bluecolor

;oplot, [-3.7, -3.3], [12./(256.*256.), 12./(256.*256.)]
;;xyouts, -3.65, 12./(256.*256.), 'CR in 12s'
;oplot, [-3.7, -3.3], [100./(256.*256.), 100./(256.*256.)]
;xyouts, -3.65, 100./(256.*256.), 'CR in 100s', charthick=3
;oplot, [-3.7, -3.3], [200./(256.*256.), 200./(256.*256.)]
;xyouts, -3.65, 200./(256.*256.), 'CR in 200s', charthick=3


;legend, ['Ch1', 'Ch2'], linestyle = [0,0], color = [blackcolor, bluecolor],/left,/top, thick=3, charthick=3

;latent_strength
;---------------------------------------------------------------------------
;ch1
x = vreset(ch1_e100t30b45)
y = latent_strength(ch1_e100t30b45)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Vreset', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Fractional Latent Strength', title = 'b450 g36 d35', xrange=[-3.7, -3.3];, /ylog, yrange=[.001,.01]
;oplot,x[sort(x)], y[sort(x)], color=bluecolor
;oplot,x[sort(x)], y[sort(x)], color=bluecolor, psym=2, thick=3


;ch2
;x = vreset(ch2_e100t30b45)
;y = latent_strength(ch2_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 2, color=bluecolor; xtitle ='Vreset', thick=3, charthick=3, xthick=3, ythick=3,charsize=1.5,$
;      ytitle = 'Fractional Latent Strength', title = 'Ch2 b450 g36 d35', xrange=[-3.7, -3.3];, /ylog, yrange=[.001,.01]
;oplot,x[sort(x)], y[sort(x)], color=bluecolor

;legend, ['Ch1', 'Ch2'], linestyle = [0,0], color = [blackcolor, bluecolor],/right,/bottom


;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;vdduc
;-----------------------------------------------------------------------------
;--------------------------------------------------------------------------
!P.multi = [0,2,2]


ch1_e12t30b45= where(channel eq 1 and bias eq 750 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 32.1  and vgg1 eq -3.75 and vreset eq -3.6)
ch2_e12t30b45= where(channel eq 2 and bias eq 600 and exptime eq 12 and new_mean ne 0 and corr_cernoxK eq 32.1  and vgg1 eq -3.7 and vreset eq -3.6)
ch1_e100t30b45= where(channel eq 1 and bias eq 750 and exptime eq 100 and new_mean ne 0 and corr_cernoxK eq 32.1  and vgg1 eq -3.75 and vreset eq -3.6)
ch2_e100t30b45= where(channel eq 2 and bias eq 600 and exptime eq 100 and new_mean ne 0 and corr_cernoxK eq 32.1  and vgg1 eq -3.7 and vreset eq -3.6)

for i = 0, n_elements(ch1_e12t30b45) - 1 do print, 'd1', corr_cernoxK(ch1_e12t30b45(i)), name(ch1_e12t30b45(i)), bias(ch1_e12t30b45(i)), vgg1(ch1_e12t30b45(i)), vreset(ch1_e12t30b45(i)), vdduc(ch1_e12t30b45(i))
for i = 0, n_elements(ch2_e12t30b45) - 1 do print, 'd 2', name(ch2_e12t30b45(i)), bias(ch2_e12t30b45(i)), vgg1(ch2_e12t30b45(i)), vreset(ch2_e12t30b45(i)), vdduc(ch2_e12t30b45(i))

;mean
;---------------------------------------------------------------------------
;ch1

x = vdduc(ch1_e12t30b45)
print, 'ch1', x
y = new_mean(ch1_e12t30b45)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Vdduc', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Mean Noise Level (MJy/sr)', title = '12s b750 g375 r36', xrange=[-3.6, -3.0]
;oplot,x[sort(x)], y[sort(x)]

;x = vdduc(ch1_e100t30b45)
;y = new_mean(ch1_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 6
;oplot,x[sort(x)], y[sort(x)]

;legend, ['12s', '100s'], psym = [2,6],/left,/center

;ch2
x = vdduc(ch2_e12t30b45)
print, 'ch2', x
y = new_mean(ch2_e12t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 2,color=bluecolor, thick=3; xtitle ='Vdduc', thick=3, charthick=3, xthick=3, ythick=3,charsize=1.5,$
;;      ytitle = 'Mean Noise Level (MJy/sr)', title = 'Ch2 12s b600 g37 r36', xrange=[-3.6, -3.0]
;oplot,x[sort(x)], y[sort(x)], color=bluecolor

;x = vdduc(ch2_e100t30b45)
;y = new_mean(ch2_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 6
;oplot,x[sort(x)], y[sort(x)]


;legend, ['12s', '100s'], psym = [2,6],/left,/center



;noisypix
;---------------------------------------------------------------------------
;ch1
x = vdduc(ch1_e100t30b45)
y = noisypix(ch1_e100t30b45)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Vdduc', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Fraction of Noisy Pixels', title = 'b750 g375 r36', /ylog, xrange=[-3.6, -3.0], yrange=[0.001,1]
;oplot,x[sort(x)], y[sort(x)]


;ch2
x = vdduc(ch2_e100t30b45)
y = noisypix(ch2_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 2, color=bluecolor , thick=3;xtitle ='Vdduc', thick=3, charthick=3, xthick=3, ythick=3,charsize=1.5,$
;;      ytitle = 'Fraction of Noisy Pixels', title = ' b600 g37 r36',/ylog, xrange=[-3.6, -3.0]
;oplot,x[sort(x)], y[sort(x)], color=bluecolor


;oplot, [-3.6, -3.0], [12./(256.*256.), 12./(256.*256.)]
;xyouts, -3.55, 12./(256.*256.), 'CR in 12s', charthick=3
;oplot,[-3.6, -3.0], [100./(256.*256.), 100./(256.*256.)]
;xyouts, -3.55, 100./(256.*256.), 'CR in 100s', charthick=3
;oplot, [-3.6, -3.0], [200./(256.*256.), 200./(256.*256.)]
;xyouts, -3.55, 200./(256.*256.), 'CR in 200s', charthick=3


;latent_strength
;---------------------------------------------------------------------------
;ch1
x = vdduc(ch1_e100t30b45)
y = latent_strength(ch1_e100t30b45)
;plot,x[sort(x)], y[sort(x)], psym = 2, xtitle ='Vdduc', thick=3, charthick=3, xthick=3, ythick=3,$
;      ytitle = 'Fractional Latent Strength', title = ' b750 g375 r36', xrange=[-3.6, -3.0];, /ylog, yrange=[.001,.01]
;oplot,x[sort(x)], y[sort(x)], color= bluecolor
;oplot,x[sort(x)], y[sort(x)], color= bluecolor, psym=2, thick=3



;ch2
;x = vdduc(ch2_e100t30b45)
;y = latent_strength(ch2_e100t30b45)
;oplot,x[sort(x)], y[sort(x)], psym = 2, color=bluecolor; xtitle ='Vdduc', thick=3, charthick=3, xthick=3, ythick=3,charsize=1.5,$
;      ytitle = 'Fractional Latent Strength', title = 'Ch2 12s b600 g37 r36', xrange=[-3.6, -3.0];, /ylog, yrange=[.001,.01]
;oplot,x[sort(x)], y[sort(x)], color=bluecolor


;------

ps_close, /noprint,/noid

end
