function run_calstar, dir_name, aor_cal
;dir_name = '/Users/jkrick/iwic/IWICA/IRAC016300'

print, 'dir_name in calstar', dir_name + '/'+aor_cal

readcol, dir_name + '/'+aor_cal, calaor, format="A",/silent
print, 'calaor', calaor
;junk = run_mopex_long_calstar(dir_name +'/bcd/'+calaor,  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl', 0.4)

ra = 269.727857
dec = 67.793591

;read in the warm mosaics we just made
fits_read, dir_name + '/ch1/Combine/mosaic.fits', warmdata_ch1, warmhead_ch1
fits_read,dir_name + '/ch2/Combine/mosaic.fits', warmdata_ch2, warmhead_ch2

;ch1
adxy, warmhead_ch1, ra, dec, xcen_warm_ch1, ycen_warm_ch1
aper,  warmdata_ch1, xcen_warm_ch1, ycen_warm_ch1, flux_warm_ch1, fluxerr_warm_ch1, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

;ch2
adxy, warmhead_ch2, ra, dec, xcen_warm_ch2, ycen_warm_ch2
aper,  warmdata_ch2, xcen_warm_ch2, ycen_warm_ch2, flux_warm_ch2, fluxerr_warm_ch2, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

;read in the cold mosaics
;same in this case
;l '/Users/jkrick/nutella/IRAC/iwic210/NPM1p67.0636/data/r28899328/ch1/Combine/mosaic.fits'
fits_read,dir_name + '/cal_mosaic_ch1.fits', colddata_ch1, coldhead_ch1
fits_read, dir_name +'/cal_mosaic_ch2.fits', colddata_ch2, coldhead_ch2

;ch1
adxy, coldhead_ch1, ra, dec, xcen_cold_ch1, ycen_cold_ch1
aper,  colddata_ch1, xcen_cold_ch1, ycen_cold_ch1, flux_cold_ch1, fluxerr_cold_ch1, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

;ch2
adxy, coldhead_ch2, ra, dec, xcen_cold_ch2, ycen_cold_ch2
aper,  colddata_ch2, xcen_cold_ch2, ycen_cold_ch2, flux_cold_ch2, fluxerr_cold_ch2, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent


;---------------------------------
fluxconv_cold_ch1 = 0.1088  ;from headers
fluxconv_cold_ch2 = 0.1388  ;from headers

fluxconv_warm_ch1  = flux_cold_ch1*fluxconv_cold_ch1 / flux_warm_ch1
fluxconv_warm_ch2  = flux_cold_ch2*fluxconv_cold_ch2 / flux_warm_ch2


print, 'fluxconv_ch1, ch2 calstar', fluxconv_warm_ch1, fluxconv_warm_ch2
return, [fluxconv_warm_ch1, fluxconv_warm_ch2]

end
