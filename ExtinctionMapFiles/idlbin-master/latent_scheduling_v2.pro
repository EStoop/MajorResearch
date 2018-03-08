pro latent_scheduling
;redcolor = FSC_COLOR("Red", !D.Table_Size-2)
;ps_open, filename='/Users/jkrick/irac_warm/scheduling/latent_scheduling.ps',/portrait,/square,/color

;Ch1 only

;mag = findgen(20) / 5.
;exptime = [0.02, 0.1,0.4,2,6,12,30,100]
pixel_scale = 1.22
gain = 3.7
flux_conv = .1253

ltime = 0.5 ; in hrs.
;flux_electron = findgen(10000) 

latent_thresh = 400.
;for m = 0, n_elements(flux_electron) - 1 do begin
mag = findgen(25) / 5. - 2.0
;exptime = [0.02, 0.1,0.4,2,6,12,30,100]
pixel_scale = 1.22
gain = 3.7
flux_conv = .1253
latent_flux = fltarr(n_elements(mag))
time = fltarr(n_elements(mag))
flux_electron = fltarr(n_elements(mag))
exptime = [0.02, 0.1,0.4,2,6,12,30,100]
back = [21., 60.,70., 120., 210., 310., 320.,575.] ;DN  
back_e = gain* back

for e = 0, n_elements(exptime) - 1 do begin
   for m = 0, n_elements(mag) - 1 do begin
      flux = magab_to_flux(mag(m)) ;in erg/s/cm2/Hz
      flux = flux / 1E-23          ; in jy
      flux = flux *1000            ; in mjy
      flux_electron(m) = mjy_to_electron(flux, pixel_scale, gain, exptime(e), flux_conv)
      
      y_int =  -5614.36 +880.124*alog10(flux_electron(m)) ;
   ;   latent_flux(m) = y_int*exp(-ltime/4.4)
 ;    print, mag(m), exptime(e), y_int, flux_electron(m)
      
      ;go after time required to get below some flux threshold = the background
      time(m) = -4.4*alog(latent_thresh/y_int)
     ; print, mag(m), y_nt, time(m), flux_electron(m)
   endfor                       ; end for each flux level
   if e eq 0 then begin
      a = plot(mag, time, '1b', xtitle = '3.6 Magnitude', ytitle = 'latent decay time (hrs)', title = 'Ch1', yrange = [0, 11], xticklen = 1)
 ;     t = text(4.0, time(n_elements(time) - 1),string(exptime(e)),/data, color = 'red') 
   endif else begin
      a = plot(mag, time, '1b', /overplot)
      t = text(2.5, time(n_elements(time) - 1),string(exptime(e)),/data, color = 'blue') 
   endelse

endfor                          ; for each exptime
p = polyline([-2,3], [8, 8],color = 'gold',/data, thick = 2)
p = polyline([-2,3], [9, 9],color = 'orange',/data, thick = 2)
p = polyline([-2,3], [10, 10],color = 'red',/data, thick = 2)

end
