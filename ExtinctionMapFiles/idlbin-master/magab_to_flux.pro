function magab_to_flux, magab
;flux in units of ergs/s/cm^2/Hz

flux = 10^((magab + 48.6)/(-2.5))
;print, 'flux in ergs/s/cm2/Hz ', flux, '.  AB mag ', magab
return, flux

end
