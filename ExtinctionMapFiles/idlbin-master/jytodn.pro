function jytodn, jy, flux_conv, exptime, pixel_scale

;flux_conv = .1198, .1443 Mjy/sr  /  Dn/s
;pixel_scale = 1.22


scale = pixel_scale^2 ;1.22D * 1.22D
; convert scale from arcsec^2 to sr and scale to Jy
scale = scale * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06

mjypersr = jy / scale


dn= mjypersr * exptime /flux_conv

return, dn

end
