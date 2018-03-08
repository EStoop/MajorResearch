function dntojansky, dn, flux_conv, exptime, pixel_scale
;flux_conv = .1198, .1447 Mjy/sr  /  Dn/s
;pixel_scale = 1.22

mjypersr = dn*flux_conv/exptime

; Convert from summed MJy/sr to Jy
scale = pixel_scale^2 ;1.22D * 1.22D
; convert scale from arcsec^2 to sr and scale to Jy
scale = scale * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06

jansky = mjypersr * scale

return, jansky

end
