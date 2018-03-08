function dark_standards, dir_name, aor_list

print, 'starting dark_standards'
;would like to combine 6 BCDs from cold mission dark fields of 100s
;for standards, but will use the mosaiced dark field images from
;Jason.  These will have different noise properties, but I
;don't then have to worry about field rotation.
; is this ok?


;find a few bright, but unsaturated stars, do photometry
;better to use dark field catalog then to try to find them from warm
;photometry

readcol, dir_name + '/stars.txt', star_ra, star_dec, format="A"

;-----------------------------------------------------
;do photometry on warm data

;create warm mosaic
;XXX change this to be the correct new directory
readcol, dir_name + '/'+aor_list, aorname, format="A",/silent


;just look at the first position
print, 'run mopex on', strcompress(dir_name +'/bcd/'+aorname(0),/remove_all)
junk = run_mopex_iwic210(strcompress(dir_name +'/bcd/'+aorname(0),/remove_all),  '/Applications/mopex/mopex-script-env.csh', 'mosaic_ch1.nl', 'mosaic_ch2.nl', 100)

;read in new warm mosaics
fits_read, dir_name +'/ch1/Combine/mosaic.fits', warmdata_ch1,warmhead_ch1
fits_read, dir_name +  '/ch2/Combine/mosaic.fits', warmdata_ch2,warmhead_ch2

;ch1
;keep the stars inside the warm image
adxy, warmhead_ch1, star_ra, star_dec, xcen_warm_ch1, ycen_warm_ch1
good= where(xcen_warm_ch1 gt 20 and ycen_warm_ch1 gt 20 and xcen_warm_ch1 lt 470 and ycen_warm_ch1 lt 470)
xcen_warm_ch1 = xcen_warm_ch1(good)
ycen_warm_ch1 = ycen_warm_ch1(good)
;convert back to ra and dec for using on the cold mission
xyad, warmhead_ch1, xcen_warm_ch1 , ycen_warm_ch1, ra_warm_ch1, dec_warm_ch1



;do photometry
aper,  warmdata_ch1, xcen_warm_ch1, ycen_warm_ch1, flux_warm_ch1, fluxerr_warm_ch1, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

;------
;ch2
adxy, warmhead_ch2, star_ra, star_dec, xcen_warm_ch2, ycen_warm_ch2
;keep the stars inside the warm image
good = where(xcen_warm_ch2 gt 20 and ycen_warm_ch2 gt 20 and xcen_warm_ch2 lt 470 and ycen_warm_ch2 lt 470)
xcen_warm_ch2 = xcen_warm_ch2(good)
ycen_warm_ch2 = ycen_warm_ch2(good)
;convert back to ra and dec for using on the cold mission
xyad, warmhead_ch2, xcen_warm_ch2 , ycen_warm_ch2, ra_warm_ch2, dec_warm_ch2

;do photometry
aper,  warmdata_ch2, xcen_warm_ch2, ycen_warm_ch2, flux_warm_ch2, fluxerr_warm_ch2, sky, $
       skyerr, 3.7, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent



;--------------------------------------------------------

;do photometry on the cold data at the positions of the stars used for
;warm photometry

;read in cold data
fits_read, dir_name+'/mosaic_ch1.fits', colddata_ch1, coldhead_ch1
fits_read, dir_name+ '/mosaic_ch2.fits', colddata_ch2, coldhead_ch2

;ch1
adxy, coldhead_ch1 , ra_warm_ch1, dec_warm_ch1, xcen_cold_ch1, ycen_cold_ch1
aper,  colddata_ch1, xcen_cold_ch1, ycen_cold_ch1, flux_cold_ch1, fluxerr_cold_ch1, sky, $
       skyerr, 3.3, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent

;ch2
adxy, coldhead_ch2, ra_warm_ch2, dec_warm_ch2, xcen_cold_ch2, ycen_cold_ch2
aper,  colddata_ch2, xcen_cold_ch2, ycen_cold_ch2, flux_cold_ch2, fluxerr_cold_ch2, sky, $
       skyerr, 3.7, [5], [15,25], [-100,1000], /nan,/exact,/flux,/silent


;--------------------------------------------------------
fluxconv_cold_ch1 = 0.1088  ;from headers
fluxconv_cold_ch2 = 0.1388  ;from headers

fluxconv_warm_ch1  = flux_cold_ch1*fluxconv_cold_ch1 / flux_warm_ch1
fluxconv_warm_ch2  = flux_cold_ch2*fluxconv_cold_ch2 / flux_warm_ch2


;throw out highest value
;at least one bad value in there, and it won't matter if it was
;a good value
fluxconv_warm_ch1 = fluxconv_warm_ch1(sort(fluxconv_warm_ch1))
fluxconv_warm_ch1 = fluxconv_warm_ch1[0:n_elements(fluxconv_warm_ch1) - 2]
;print, 'median ch1 conversion', median(fluxconv_warm_ch1,/even)

fluxconv_warm_ch2 = fluxconv_warm_ch2(sort(fluxconv_warm_ch2))
fluxconv_warm_ch2 = fluxconv_warm_ch2[0:n_elements(fluxconv_warm_ch2) - 2]
;print, 'median ch2 conversion', median(fluxconv_warm_ch2,/even)


flux_conv = [median(fluxconv_warm_ch1), median(fluxconv_warm_ch2)]
print, 'flux_conv dark field stars', flux_conv, n_elements(fluxconv_warm_ch2)
return, flux_conv
end

