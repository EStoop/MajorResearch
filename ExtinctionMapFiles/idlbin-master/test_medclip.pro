pro test_medclip

  fits_read, '/Users/jkrick/irac_warm/latent/r30288640/ch1/bcd/SPITZER_I1_30288640_0054_0000_1_cbcd.fits', data, header
  
  bkgd = data[90:130,90:130]
  print, 'before medclip', min(bkgd), max(bkgd)
  wsm = MEDCLIP(bkgd, med, sigma, clipsig = 3)
  print, 'back in main code'
  print, 'after medclip', n_elements(wsm), min(bkgd[wsm]), max(bkgd[wsm]), med, sigma
  
  print, bkgd[wsm]
  
end
