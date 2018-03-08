;*****************************************************************************
;Given the observed wavelength in microns, the flux in Jansky and the redshift,
;return the rest wavelength in microns and the luminosity in nuLnu in
;solar luminosities.
;*****************************************************************************

function nulnu, wave_um, flux_Jy, z

    solar_luminosity = 3.862e33 ;erg/sec/Lsun
    c_cm = 3e10

    parsec_to_cm = 3.086e18 ;cm per parsec
    cm_to_um = 1e4

    c_um = c_cm * cm_to_um

    nu = c_um / wave_um      ;Hz
    Fnu_cgs = flux_Jy * 1e-23 ;erg s^-1 cm^-2 Hz^-1
    nuFnu = double(nu * Fnu_cgs) ;erg/s/cm^2
    
    distance_pc = lum_dist(z)
    Inf_factor = double((parsec_to_cm)^2. / solar_luminosity) ;cm^2 s pc^-2 erg^-1 Lsun
    nuLnu = 4.0 * !PI * (distance_pc)^2. * nuFnu * Inf_factor ;solar luminosities

;    return, {rest_wave:wave_um / (1.+ z), nulnu:nuLnu}
    return, nuLnu
end
