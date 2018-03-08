pro virgo_irac_smooth

  fits_read, '/Users/jkrick/Desktop/mosaic.fits', data, header
  x1  = sxpar(header, 'NAXIS1')
  x2  = sxpar(header, 'NAXIS2')

  hrebin, data, header, newdata, newheader, fix(x1/12), fix(x2/12)

  fits_write, '/Volumes/irac_drive/Virgo/IRAC/ch2/ch2_Combine-mosaic/mosiac_12.fits', data, header

end
