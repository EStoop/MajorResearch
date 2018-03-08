pro virgo_optical

;from NED
Av = .102
Ab = .132

;from craig rudick
;for the smaller region
;mv = 15.87
;merrv = .1
;mb = 16.56
;merrb = .09

;for the larger region
mv = 16.51
merrv = .32
mb = 17.09
merrb = .25

;apply galactic extinction
mv = mv - Av
mb = mb - Ab

;if the mags are not AB, then first conver to AB
;http://www.astro.utoronto.ca/~patton/astro/mags.html#conversions
mv = mv -0.044
mb = mb - .163

;convert AB mags  to flux
fluxv = magab_to_flux(mv)
fluxb = magab_to_flux(mb)

;and the errors also convert to flux
fluxerrv =  fluxv - magab_to_flux(mv + merrv) 
fluxerrb = fluxb - magab_to_flux(mb + merrb) 

;final numbers
print, 'flux V in erg/s/cm2/Hz ', fluxv , ' +- ', fluxerrv
print, 'flux B in erg/s/cm2/Hz ', fluxb , ' +- ', fluxerrb



end
