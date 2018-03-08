function Mjysr_to_magsquarc, mjysrpix, pixscale
;convert one SB to another

jy_pix = Mjypersrtojansky(mjysrpix, pixscale)
erg_pix = jy_pix * 1E-23
mag_pix = flux_to_magab(erg_pix)
mag_squarc = mag_pix / ((1.2)^2)
return, mag_squarc

end
