function mjy_to_electron, mjy, pixel_scale, gain, exptime, flux_conv

jy = mjy/1000.  ;milijanskies to janskies

; Convert from summed  to Jy to Mjy/sr
scale = pixel_scale^2 ;1.22D * 1.22D
; convert scale from arcsec^2 to sr and scale to Jy
scale = scale * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06

Mjypersr = jy / scale

electrons = Mjypersr * gain*exptime/flux_conv


return, electrons

end
