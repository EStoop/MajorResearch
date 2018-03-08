pro bias_final


;this code reads in the excel file of all IWIC210 bias and temperature
;tests.  Then using 'where' commands it selects any sample and plots the
;interesting quantities.

;Use 'new_mean' and 'new_sigma' instead of the 'old'
;versions which have incorrect flux_conv's

ps_open, filename='/Users/jkrick/iwic/test.ps',/portrait,/square,/color
;ps_start, filename='/Users/jkrick/iwic/start_science.ps'

redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
greencolor = FSC_COLOR("Green", !D.Table_Size-4)
yellowcolor = FSC_COLOR("Yellow", !D.Table_Size-5)
cyancolor = FSC_COLOR("cyan", !D.Table_Size-6)
orangecolor = FSC_COLOR("orange", !D.Table_Size-7)
purplecolor = FSC_COLOR("purple", !D.Table_Size-8)
blackcolor = FSC_COLOR("black", !D.Table_Size-9)

!P.multi = [0,2,1]
readcol, '/Users/jkrick/iwic/dark_summary/dark_summary_10_05.csv', name,directory,reckey,date,exptime,fowler,waitper,nframes,channel,commandedK,corr_cernoxK,cernoxK,vreset,bias,vdduc,vgg1,new_mean,new_sigma,old_mean,old_sigma,old_noisy,noisypix,flux_conv,latent_strength,latent_decay,format='A,A,I, I, I, I, I,I,I,F10.2,F10.2,F10.2,F10.2,I,F10.2,F10.2,F10.2,F10.2,F10.2,F10.2,I, F10.5,F10.2,F10.5,F10.5'


a = where(bias eq 450 and exptime eq 100 and new_mean ne 0 and channel eq 1 and vreset eq -3.6 )
print, 'n_elements(a)', n_elements(a)
for i = 0, n_elements(a) - 1 do begin
   print,  format = '(A,F10.1, I, I, F10.3, F10.2, F10.2, F10.5)', name(a(i)), cernoxK(a(i)), channel(a(i)), bias(a(i)), vgg1(a(i)), vreset(a(i)), vdduc(a(i)), new_mean(a(i))
endfor


a = where(bias eq 450 and exptime eq 100 and new_mean ne 0 and channel eq 2 and vreset eq -3.55 )
print, 'n_elements(a)', n_elements(a)
for i = 0, n_elements(a) - 1 do begin
   print,  format = '(A,F10.1, I, I, F10.3, F10.2, F10.2, F10.5)', name(a(i)), cernoxK(a(i)), channel(a(i)), bias(a(i)), vgg1(a(i)), vreset(a(i)), vdduc(a(i)), new_mean(a(i))
endfor

;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
;temp
;-----------------------------------------------------------------------------
;-----------------------------------------------------------------------------
!P.multi = [0,2,2]

ch1_e12b50= where(channel eq 1 and bias eq 500 and exptime eq 12 and new_mean ne 0)
ch2_e12b50= where(channel eq 2 and bias eq 500 and exptime eq 12 and new_mean ne 0)
ch1_e100b50= where(channel eq 1 and bias eq 500 and exptime eq 100 and new_mean ne 0)
ch2_e100b50= where(channel eq 2 and bias eq 500 and exptime eq 100 and new_mean ne 0)

;mean noise

x=cernoxK(ch1_e12b50)
y = new_mean(ch1_e12b50)

print, x, y

plot, x[sort(x)], y[sort(x)] , psym = 2,xtitle ='Cernox Temp (K)', $
      ytitle = 'Mean Noise Level (MJy/sr)', title = 'Bias = 500', xrange=[25,40], yrange=[0.02,.1],$
      thick=3, xthick=3, ythick=3, charthick=3
oplot, x[sort(x)], y[sort(x)] 

oplot, findgen(60), fltarr(60)+0.031072, linestyle = 2, thick = 3


x = cernoxK(ch2_e12b50)
y = new_mean(ch2_e12b50)
oplot, x[sort(x)], y[sort(x)] , psym = 2, color = bluecolor, thick=3
oplot, x[sort(x)], y[sort(x)] ,  color = bluecolor
oplot, findgen(60), fltarr(60)+0.04295, linestyle = 2, thick = 3, color = bluecolor

legend, ['Ch1', 'Ch2', 'cryo 1', 'cryo 2'], linestyle = [0,0,2,2], color = [blackcolor, bluecolor, blackcolor, bluecolor],/left,/top, thick=3, charthick=3


ps_close, /noprint,/noid
;ps_end, /png
end
