pro latent_scheduling
redcolor = FSC_COLOR("Red", !D.Table_Size-2)
ps_open, filename='/Users/jkrick/irac_warm/scheduling/latent_scheduling.ps',/portrait,/square,/color

;Ch1 only

mag = findgen(20) / 5.
exptime = [0.02, 0.1,0.4,2,6,12,30,100]
pixel_scale = 1.22
gain = 3.7
flux_conv = .1253
;start by asking what the flux is at a certain time, then switch to what time to get to a certain flux.
ltime = 24 ; in hrs.
plot,mag, exptime, psym = 3, xtitle = 'Star Magnitude [3.6]',ytitle = 'Exposure Time',/nodata,/ylog, yrange = [0.01, 200]

for m = 0, n_elements(mag) - 1 do begin

   for e = 0, n_elements(exptime) - 1 do begin
      ;need to go from magnitudes to expected electrons;
      flux = magab_to_flux(mag(m)) ;in erg/s/cm2/Hz
      flux = flux / 1E-23          ; in jy
      flux = flux *1000            ; in mjy
      flux_electron = mjy_to_electron(flux, pixel_scale, gain, exptime(e), flux_conv)
      y_int = -5614.36 + 880.124*alog10(flux_electron)
      latent_flux = y_int*exp(-ltime/4.4)
      mjy = electron_to_mjy(latent_flux, pixel_scale, gain, exptime(e), flux_conv)
      print, 'mag', mag(m), 'exptime', exptime(e), 'latent', mjy
      xyouts, mag[m], exptime[e],strmid(string(mjy),7,4),/data, charsize = 0.5, color = redcolor
      ;a = string(latent_flux)
      ;print, 'test', strmid(a, 7, 4)
      ;keep track of these variables

   endfor                       ; end for each exposure time

endfor ; end for each magnitude
ps_close, /noprint,/noid

end
