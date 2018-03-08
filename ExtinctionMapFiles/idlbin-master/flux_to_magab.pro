function flux_to_magab, flux
;flux in units of ergs/s/cm^2/Hz
ab = -2.5*alog10(flux) - 48.6
;print, 'flux ergs/s/cm2/Hz', flux, 'AB mag ', ab

return, ab

end
