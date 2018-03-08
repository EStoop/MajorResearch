pro virgo_globulars_aper
     
     filename =  '/Users/jkrick/Virgo/IRAC/ch1/ch1_Combine-mosaic/mosaic_bkgd.fits'
     fits_read, filename, data, header
        
;read in the locations of the globulars
     ;readcol, '/Users/jkrick/Virgo/IRAC/globulars/ngvs_ra_dec_g24.txt',  ra_start, dec_start ;, format = '(G10.4, G10.4)'
     ftab_ext, '/Users/jkrick/Virgo/IRAC/globulars/NGVS.ptcat.ugiz.i24.5.v0.2.pps99.iracregion.fits', [6,7], ra_start, dec_start
     adxy, header, ra_start, dec_start, x_start, y_start
     
     print, 'n_elements(x_start)', n_elements(x_start)

;set up the objects
     gclus = replicate({gob, gnum:0I,gra:0D, gdec:0D,ch1_ra:0D, ch1_dec:0D, x_start:0D, y_start:0D, ch1_xcen:0D, ch1_ycen:0D, ch1_flux:0D,ch1_fluxerr:0D, ch1_mag:0D,ch1_magerr:0D, ch1_sky:0D, ch1_skyerr:0D, ch2_flux:0D,ch2_fluxerr:0D, ch2_mag:0D,ch2_magerr:0D, ch2_sky:0D, ch2_skyerr:0D},n_elements(ra_start))

     gclus ={gob, findgen(n_elements(x_start)), ra_start, dec_start, 0.,0.,0.,0.,0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.,0.}


;first mark those globulars off the field of view
     a = where(gclus.x_start gt 0 and gclus.y_start gt 0 and finite(data(gclus.x_start, gclus.y_start)) eq 1)
     gclus[a].ch1_flux = -99
     gclus[a].ch2_flux = -99
     ;x_start = x_start(a)
     ;y_start = y_start(a)
     ;ra_start = ra_start(a)
     ;dec_start = dec_start(a)
     

;find good centroids for the NGVS positions in the IRAC dataset
     cntrd, data, gclus.x_start, gclus.y_start, xcen, ycen, 1.0 ;,/debug
     
;get rid of the ones which return NANs which really are off the field
;just 5 of these
     b = where(finite(xcen) eq 1 and xcen gt 0 and ycen gt 0)
     xcen = xcen(b)
     ycen = ycen(b)
     x_start = x_start(b)
     y_start = y_start(b)
     ra_start = ra_start(b)
     dec_start = dec_start(b)

;how many globular candidates are left?
     print, 'there are ', n_elements(x_start), 'globulars in the FOV'

   
;consider removing the objects where the centroids have moved by more than an arcsec or so.
;this means that the object listed was not able to be centroided on, and so we will get photometry for a neighbor.
;which is not relevant information
;would obviously like to minimize the number of these by altering the cntrd command.

;--------------------------------------------
;run aperture photometry on the NGVS positions
     phpadu = 3.7               ; XXXdon't know what else to set this at?
     aper, data, x_start, y_start, flux_start, ferrap_start, sky_start, skyerr_start, phpadu, [3.16], [12,20],/nan, /exact,/flux,/silent
     

;run some aperture photometry on the new centroids
     phpadu = 3.7               ; XXXdon't know what else to set this at?
     aper, data, xcen, ycen, flux, ferrap, sky, skyerr, phpadu, [3.16], [12,20],/nan, /exact,/flux,/silent

;aperture correction from SWIRE
     apcorr_ch1 = 1.6 ;XXX find the real number
     flux = flux * apcorr_ch1
     flux_start = flux_start * apcorr_ch1
;--------------------------------------------

;how far is the centroided position from the starting position?
     dist = sqrt((xcen - x_start)^2 + (ycen-y_start)^2)
     dist = dist *0.6           ; put in arcseconds, not in pixels
     
;want a cumulative distribution for this; then make a cut by percentage.
;need to sort on distance
     sortdist= dist(sort(dist))
     N1 = n_elements(sortdist)
     f1 = (findgen(N1) + 1.)/ N1
     c = plot( sortdist, f1, '6b2', xtitle='Distance from Original Center (arcsec)', ytitle='Cumulative Fraction', xrange=[0,2], yrange=[0,1])
     
     line_x = [0, 2.0]
     line_y = [0.9, 0.9]
     line5 = polyline(line_x, line_y, thick = 2, color = !color.black,/data)
     
;test plot distance vs. brightness.  Probably the faint ones that are at larger distances
  ; although maybe not since the new source could be any brightness
     t = plot( dist,flux_start - flux, '6r1+', xtitle = 'distance from original position (arcsec)',ytitle = 'flux _start - flux(Mjy/sr)' )
     line_x = [0.9, 0.9]
     line_y = [-1.5, 1.5]
     line5 = polyline(line_x, line_y, thick = 2, color = !color.black,/data)
     
;take out those sources which have distances greater than 0.9arcsec, or 10% of the total sources
     d = where(dist le 0.9, dcount)
     xcen = xcen(d)
     ycen = ycen(d)
     x_start = x_start(d)
     y_start = y_start(d)
     ra_start = ra_start(d)
     dec_start = dec_start(d)
     flux = flux(d)
     ferrap= ferrap(d)
     sky = sky(d)
     skyerr = skyerr(d)

     print, 'there are ', dcount, 'globulars after cutting on distance'
 ;make sure I did this correctly
 ;plothist,  (sqrt((xcen - x_start)^2 + (ycen-y_start)^2))*0.6, xarr, yarr, bin = 0.1
;--------------------------------------------
     
;convert new centroids into ra and dec
     xyad, header, xcen, ycen, ra, dec
     
;make a region file to visually check the centroiding
     openw, outlun, '/Users/jkrick/Virgo/IRAC/globulars/irac_ra_dec_g24.reg', /get_lun
     for r = 0,  n_elements(ra) - 1 do begin
        if finite(ra(r)) lt 1 then ra(r) = 0.0
        if finite(dec(r)) lt 1 then dec(r) = 0.0
        printf, outlun, format ='(A, g10.9, A, g10.9, A)', 'J2000; circle  ',ra_start(r), 'd  ', dec_start(r),    'd 5.0" # color=green'
        printf, outlun,format ='(A, g10.9, A, g10.9, A)', 'J2000; circle  ', ra(r), 'd  ', dec(r),    'd 4.5" # color=red'
     endfor
     close, outlun
     free_lun, outlun
;--------------------------------------------
     
;convert into useful units  
     jy = Mjypersrtojansky(flux,0.6) ; convert from Mjy/sr to jy
     erg = jy * 1E-23                ;conver from jy to erg/s/cm2/Hz
     magab = flux_to_magab(erg)      ;convert from erg/s/cm2/Hz to MagAB
     
;XXXstill need to convert flux errors into mag errors
     
     
;--------------------------------------------
     
  
;need to do the other channel.
   ;use the ch1 positions to measure
   ;simple aperture photometry on ch2

     filename =  '/Users/jkrick/Virgo/IRAC/ch2/ch2_Combine-mosaic/mosaic_bkgd.fits'
     fits_read, filename, data_ch2, header_ch2
     phpadu = 3.7               ; XXXdon't know what else to set this at?
     aper, data_ch2, xcen, ycen, flux_ch2, ferrap_ch2, sky_ch2, skyerr_ch2, phpadu, [3.16], [12,20],/nan, /exact,/flux,/silent
;aperture correction from SWIRE
     apcorr_ch2 = 1.6 ;XXX find the real number
     flux_ch2 = flux_ch2 * apcorr_ch2

;convert into useful units  
     jy_ch2 = Mjypersrtojansky(flux_ch2,0.6) ; convert from Mjy/sr to jy
     erg_ch2 = jy_ch2 * 1E-23                ;conver from jy to erg/s/cm2/Hz
     magab_ch2 = flux_to_magab(erg_ch2)      ;convert from erg/s/cm2/Hz to MagAB
;XXXstill need to convert flux errors into mag errors

;--------------------------------------------


;what is the best way to store all of this?
;need some sort of useful output file
     for i = 0, n_elements(xcen) - 1 do begin
        print, format = '(F10.2, F10.2, F10.5, F10.5, F10.2, F10.4, F10.4, F10.2, F10.2, F10.2)', xcen(i), ycen(i), ra(i), dec(i), magab(i), flux(i), ferrap(i), magab_ch2(i), flux_ch2(i), ferrap_ch2(i)
     endfor


end
