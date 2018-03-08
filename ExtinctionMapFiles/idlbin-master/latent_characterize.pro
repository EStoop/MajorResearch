pro latent_characterize
;the question is to measure the background and as at what time the
;latent flux is equal to the background.

;ps_start, filename= '/Users/jkrick/iwic/latent_characterize.ps'
ps_open, filename='/Users/jkrick/iwic/latent_characterize.ps',/portrait,/square,/color

!P.multi = [0,0,1]

!P.thick = 1
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
bluecolor = FSC_COLOR("Blue", !D.Table_Size-3)
blackcolor = FSC_COLOR("Black", !D.Table_Size-4)
;----------------------------------------------------
;----------------------------------------------------
;----------------------------------------------------
;----------------------------------------------------
;----------------------------------------------------
;ch1
;----------------------------------------------------
;----------------------------------------------------
;----------------------------------------------------
;----------------------------------------------------

gain = 3.82
flux_conv = .1277;,.1461] ;from Mark's most recent post on MODA

;----------------------------------------------------
;latent test pc2 hd166780 k = 3.9 100s
;----------------------------------------------------
restore, '/Users/jkrick/IWIC/latent.sav'

time = ch1_x
flux = ch1_y

;only really want the peak and then down, not the three measurements
;at peak
time = time[2:*]
flux = flux[2:*]

;need to reset time
time = time - time(0)

;convert flux which is in Mjy/sr into electrons
exptime = 100.
flux = flux * gain
flux = flux * exptime
flux = flux / flux_conv

plot, 1/flux, time, psym = 2, /xlog,yrange = [-0.1,2],$
      xtitle = '1/flux (electrons)', ytitle = 'time (hrs)', title = 'Ch1'
  
start = [.05,-1E-4]
noise = fltarr(n_elements(flux))
noise[*] = 1                                                ;equally weight the values
result= MPFITFUN('exponential',1/flux,time, noise, start)   ;./quiet    
  
xarr = findgen(10000)/ 1000000.
  
oplot, xarr, result(0)*exp(-(xarr)/(result(1)))

;so at what time does the latent equal background
print, 'time at flux = 0.1', result(0)*exp(-(1/0.1)/(result(1)))


;----------------------------------------------------
;latent test vs78 nsv6890 k=1.7 12s
;----------------------------------------------------

restore,  '/Users/jkrick/iwic/latent_vs78.sav'
time = time[1:*]
flux  = flux[1:*]
time = time - time(0)

;convert flux which is in Mjy/sr into electrons
exptime = 12.
flux = flux * gain
flux = flux * exptime
flux = flux / flux_conv

oplot, 1/flux, time, psym = 6, color = redcolor
start = [.05,-1E-4]
noise = fltarr(n_elements(flux))
noise[*] = 1                                                ;equally weight the values
result= MPFITFUN('exponential',1/flux,time, noise, start)   ;./quiet    
    
oplot, xarr, result(0)*exp(-(xarr)/(result(1))), linestyle = 3, color = redcolor

;so at what time does the latent equal background
print, 'time at flux = 0.1', result(0)*exp(-(1/0.1)/(result(1)))

;----------------------------------------------------
;latent test nsv6114 k = 2.4 30s
;----------------------------------------------------
restore, '/Users/jkrick/iwic/latent_nsv6114.sav'

;convert flux which is in Mjy/sr into electrons
exptime = 30.
flux = flux * gain
flux = flux * exptime
flux = flux / flux_conv

oplot, 1/flux, time, psym = 4, color = bluecolor
start = [.05,-1E-4]
noise = fltarr(n_elements(flux))
noise[*] = 1                                                ;equally weight the values
result= MPFITFUN('exponential',1/flux,time, noise, start)   ;./quiet    
    
oplot, xarr, result(0)*exp(-(xarr)/(result(1))), linestyle = 2, color = bluecolor

;so at what time does the latent equal background
print, 'time at flux = 0.1', result(0)*exp(-(1/0.1)/(result(1)))

legend, ['K=3.9, 100s', 'K=2.4, 30sHDR','K=1.7, 12s'], psym = [2,4,6], color = [blackcolor, bluecolor, redcolor],/left,/top, thick=3, charthick=3

;----------------------------------------------------
;----------------------------------------------------
;----------------------------------------------------
;----------------------------------------------------
;----------------------------------------------------
;ch2
;----------------------------------------------------
;----------------------------------------------------
;----------------------------------------------------
;----------------------------------------------------

gain = 3.82
flux_conv = .1461;from Mark's most recent post on MODA

;----------------------------------------------------
;latent test pc2 hd166780 k = 3.9 100s
;----------------------------------------------------
restore, '/Users/jkrick/IWIC/latent.sav'

time = ch2_x
flux = ch2_y

;only really want the peak and then down, not the three measurements
;at peak
time = time[2:*]
flux = flux[2:*]


;need to reset time
time = time - time(0)

;latent goes negative, then messes up the fits
flux = flux(where(time lt 1.))
time = time(where(time lt 1.))
;convert flux which is in Mjy/sr into electrons
exptime = 100.
flux = flux * gain
flux = flux * exptime
flux = flux / flux_conv

for i = 0, n_elements(time) - 1 do print, time(i), flux(i)
plot, 1/flux, time, psym = 2, /xlog,yrange = [-0.1,2],xrange=[1E-7,1E-2],$
      xtitle = '1/flux (electrons)', ytitle = 'time (hrs)', title = 'Ch2'
  
start = [.01,-4E-4]
noise = fltarr(n_elements(flux))
noise[*] = 1                                                ;equally weight the values
result= MPFITFUN('exponential',1/flux,time, noise, start)   ;./quiet    
  
xarr = findgen(10000)/ 1000000.
  
oplot, xarr, result(0)*exp(-(xarr)/(result(1)))

;so at what time does the latent equal background
print, 'time at flux = 0.1', result(0)*exp(-(1/0.1)/(result(1)))



;----------------------------------------------------
;latent test nsv6114 k = 2.4 30s
;----------------------------------------------------
restore, '/Users/jkrick/iwic/latent_nsv6114_ch2.sav'

;convert flux which is in Mjy/sr into electrons
exptime = 30.
flux = flux * gain
flux = flux * exptime
flux = flux / flux_conv

oplot, 1/flux, time, psym = 4, color = bluecolor
start = [.05,-1E-4]
noise = fltarr(n_elements(flux))
noise[*] = 1                                                ;equally weight the values
result= MPFITFUN('exponential',1/flux,time, noise, start)   ;./quiet    
    
oplot, xarr, result(0)*exp(-(xarr)/(result(1))), linestyle = 2, color = bluecolor

;so at what time does the latent equal background
print, 'time at flux = 0.1', result(0)*exp(-(1/0.1)/(result(1)))

legend, ['K=3.9, 100s', 'K=2.4, 30sHDR'], psym = [2,4], color = [blackcolor, bluecolor],/left,/top, thick=3, charthick=3


;ps_end, /png
ps_close, /noprint,/noid

end
