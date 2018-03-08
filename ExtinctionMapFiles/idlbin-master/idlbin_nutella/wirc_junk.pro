pro wircjunk

 readcol,"/Users/jkrick/nutella/palomar/wirc/j/jband_wirc_nohead.txt", number,x_image, y_image,mag_iso,magerr_iso,mag_aper,magerr_aper,mag_auto,magerr_auto,mag_best, magerr_best, kron_rad,junk, junk,isoareaimage, ra,dec,junk,junk,junk, junk,fwhm_image, junk,ell, junk , format="A"

;plot, mag_best, mag_auto, psym = 2, xrange=[10,20], yrange=[10,20]


plothist, mag_best / mag_auto, bin = 0.001

end
