pro virgo_globulars
set_plot, 'ps'
device, filename = '/Users/jkrick/virgo/irac/globulars/flaghist.ps'

readcol,'/Users/jkrick/virgo/irac/globulars/ch1_test.cat',NUMBER,X_WORLD ,Y_WORLD ,X_IMAGE,Y_IMAGE,FLUX_AUTO,MAG_AUTO,MAGERR_AUTO,FLUX_APER,MAG_APER,MAGERR_APER,FLUX_BEST,MAG_BEST,MAGERR_BEST,FWHM_IMAGE,ISOAREA_IMAGE,FLUX_MAX,ELLIPTICITY,CLASS_STAR,FLAGS,THETA_IMAGE,A_IMAGE,B_IMAGE,junk, format="A"

plothist, class_star, xhist, yhist, bin=0.10, xrange = [-1., 1.], xtitle = 'class_star', ytitle = 'number'

plot, isoarea_image, mag_aper, psym = 3, xrange = [1, 3E3], yrange = [10, 32],/xlog, xtitle = 'isoarea_image', ytitle = 'mag_aper'
device, /close

;make a region file to compare with Peng's input region file.
openw, outlun, '/Users/jkrick/Virgo/IRAC/globulars/irac_ra_dec_g24.reg', /get_lun
for r=0, n_elements(X_WORLD) -1 do  printf, outlun, 'J2000; circle  ', X_WORLD(r), '  ', Y_WORLD(r),    ' 6.0" # color=red'

close, outlun
free_lun, outlun

end
