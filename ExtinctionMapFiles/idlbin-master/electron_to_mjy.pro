function electron_to_mjy, electrons, pixel_scale, gain, exptime, flux_conv

sbtoe = gain*exptime/flux_conv
Mjypersr = electrons / sbtoe

scale = pixel_scale^2 ;1.22D * 1.22D
; convert scale from arcsec^2 to sr and scale to Jy
scale = scale * !DPI * !DPI / (3600.D * 3600.D * 180.D * 180.D) * 1.0D+06

jy = Mjypersr * scale

return,jy*1000. ; milijanskies
end
