;function to calcluate noise pixel
function noisepix, im, xcen, ycen, ronoise, gain, exptime, fluxconv,naxis, apradnp, skyradin, skyradout
;  apradnp = 3;4
;  skyradin = 3;10
;  skyradout = 6;12
  print, 'starting noisepix ';, ronoise, mean(xcen), mean(ycen)
  convfac = gain*exptime/fluxconv

  if naxis gt 2 then begin  ; for subarray
     np = fltarr(64,/NOZERO)
     for npj = 0, 63 do begin
        indim = im[*,*,npj]
        indim = indim*convfac
        ; aper requires a 2d array
        aper, indim, xcen[npj], ycen[npj], topflux, topfluxerr, xb, xbs, 1.0, apradnp,[skyradin,skyradout],/flux, /silent, /nan, readnoise = ronoise;,/exact;, setskyval = 0
        aper, indim^2, xcen[npj], ycen[npj], bottomflux, bottomfluxerr, xb, xbs, 1.0,apradnp,[skyradin,skyradout],/flux, /silent, /nan, readnoise = ronoise;,/exact;, setskyval = 0
        
        np[npj] = topflux^2 / bottomflux
     endfor
  endif

  if naxis lt 3 then begin ; for full array
     indim = im*convfac
        
        aper, indim, xcen, ycen, topflux, topfluxerr, xb, xbs, 1.0, apradnp,[skyradin,skyradout],/flux,/exact, /silent, /nan, readnoise = ronoise, setskyval = 0
        aper, indim^2, xcen, ycen, bottomflux, bottomfluxerr, xb, xbs, 1.0,apradnp,[skyradin,skyradout],/flux,/exact, /silent, /nan, readnoise = ronoise, setskyval = 0
        
        np = topflux^2 / bottomflux
  endif

     return, np               
end

